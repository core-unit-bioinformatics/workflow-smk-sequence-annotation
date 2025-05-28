#!/usr/bin/env python3

import argparse as argp
import enum
import pathlib as pl
import re

import pandas as pd
import pyranges as pr
import xopen

___prog__ = "dump_labeled_ref_bed.py"
__version__ = "0.1"
__description__ = """
This script takes a normalized PAF alignment file as single input.

It reads a BED file with labeled regions and performs a liftover-like
operation of the BED regions onto the query sequence.

By default, it ignores alignments that are smaller than half the
size of the smallest labeled region.
"""


CIGAR_OPS_REGEXP = "[0-9]+(\\=|X|D|I|M)"


class CIGARstep(enum.Enum):
    TARGET = 0
    QUERY = 1
    BOTH = 2
    IDENT = 3


class CIGARwalker:
    __slots__ = (
        "cigar", "cigar_ops", "move_ops", "orientation",
        "tstart", "tend", "qstart", "qend",
        "target_regions", "query_regions",
        "row_idx",
        "last_start_cigar", "last_start_target", "last_start_query"
    )

    def __init__(self, tstart, tend, qstart, qend, orientation, cigar, row_idx):

        self.cigar_ops = re.compile(CIGAR_OPS_REGEXP)
        self.cigar = cigar
        self.move_ops = {
            "=": CIGARstep.IDENT,
            "M": CIGARstep.BOTH,
            "X": CIGARstep.BOTH,
            "D": CIGARstep.TARGET,
            "I": CIGARstep.QUERY
        }
        self.orientation = orientation
        assert self.orientation in [1,-1]
        self.tstart = tstart
        self.tend = tend
        if self.orientation < 0:
            # reverse orientation, i.e. query is reversed.
            # simplify CIGAR walking by flipping start/end,
            # that is, we are always walking 'left to right'
            self.qstart = qend * -1
            self.qend = qstart * -1
        else:
            self.qstart = qstart
            self.qend = qend
        assert self.qstart < self.qend
        self.row_idx = row_idx

        # saving state for quicker iterations
        self.last_start_cigar = 0
        self.last_start_target = self.tstart
        self.last_start_query = self.qstart

        return None

    def walk_cigar(self):
        """Iterate through CIGAR operations,
        starting from self.last_start_cigar

        Yields a 3-tuple of CIGAR[pos, step, operation]
        where 'pos' is just used to save the state/position
        in the CIGAR string (ensure single iteration/O(n)
        runtime), 'step' is the step size and 'operation'
        is one of the CIGARstep enum types indicating which
        coordinates should be advanced by 'step' bp
        (target or query).
        """

        iter_cigar_ops = self.cigar_ops.finditer(
            self.cigar[self.last_start_cigar:]
        )

        for cigar_op in iter_cigar_ops:
            cigar_pos, _ = cigar_op.span()
            tmp = cigar_op.group(0)
            move_op = self.move_ops[tmp[-1]]
            step = int(tmp[:-1])
            yield cigar_pos, step, move_op
        return

    def exhaustive_walk(self):
        """This is a simplified version of the find_range()
        function that is called in cases where the labeled region
        is larger than the alignment portion. It iterates through
        the entrie CIGAR string and returns the necessary statistics
        """
        iter_t = self.last_start_target
        iter_q = self.last_start_query
        aln_block_length = 0
        aln_block_ident = 0
        max_ident_block = 0
        for cigar_pos, step, move in self.walk_cigar():
            # The alignment block length is computed by
            # simply summing up all CIGAR steps.
            aln_block_length += step
            if move in [CIGARstep.IDENT, CIGARstep.BOTH]:
                iter_t += step
                iter_q += step
            elif move == CIGARstep.TARGET:
                iter_t += step
            else:
                iter_q += step
            if move == CIGARstep.IDENT:
                # The matching block length can only be computed
                # if the CIGAR string uses the =/X operations.
                # Potential TODO: also allow M operation here?
                aln_block_ident += step
                max_ident_block = max(max_ident_block, step)

        if self.orientation < 0:
            # reverse orientation, i.e. query is reversed.
            # simplify CIGAR walking by flipping start/end,
            # that is, we are always walking 'left to right'
            lifted_start = self.qend * -1
            lifted_end = self.qstart * -1
        else:
            lifted_start = self.qstart
            lifted_end = self.qend

        assert lifted_start < lifted_end, f"{lifted_start} - {lifted_end}"
        assert lifted_start >= 0, f"{lifted_start} - {lifted_end}"
        assert lifted_end > 0, f"{lifted_start} - {lifted_end}"

        return lifted_start, lifted_end, aln_block_ident, aln_block_length, max_ident_block

    def find_range(self, find_start, find_end):
        """This function finds the range between
        start and end in target coordinates and
        returns the equivalent range in query
        coordinates (~ performs a liftover).

        This functions corrects for the scenario
        that (i) a small region is lifted over that
        is fully contained in a long alignments;
        (ii) the respective region has been deleted
        in the query.
        """

        iter_t = self.last_start_target
        iter_q = self.last_start_query
        lifted_start = None
        lifted_end = None
        aln_block_length = 0
        aln_block_ident = 0
        max_ident_block = 0
        deleted_region = False
        in_query_region = False

        stepping_in = False
        stepping_out = False
        stepping_over = False
        offset_in = 0
        offset_out = 0
        for cigar_pos, step, move in self.walk_cigar():
            # lookahead - are we stepping over the find range?
            if iter_t + step > find_end and iter_t < find_start:
                # stepping over the find range in one step;
                # this is a special case that is handled at the end
                # of the function
                stepping_over = True
            # lookahead - are we stepping into the find range?
            elif iter_t + step >= find_start and iter_t <= find_start:
                stepping_in = True
                offset_in = find_start - iter_t
                assert offset_in >= 0

                # the following preserves the starting state
                # for subsequent calls to walk_cigar()
                self.last_start_cigar = cigar_pos
                self.last_start_target = iter_t
                self.last_start_query = iter_q

            # lookahead - are we stepping out of the find range?
            elif iter_t + step >= find_end:
                stepping_out = True
                offset_out = iter_t + step - find_end
                assert offset_out >= 0
            else:
                stepping_in = False
                stepping_out = False
                stepping_over = False

            # make the step
            if move in [CIGARstep.IDENT, CIGARstep.BOTH]:
                iter_t += step
                iter_q += step
            elif move == CIGARstep.TARGET:
                iter_t += step
            else:
                iter_q += step

            if self.orientation < 0:
                assert iter_q < 0 or stepping_out or stepping_over, f"{iter_q} / {stepping_out} / {stepping_over}"

            if stepping_in:
                # note for lifted start: step size has already been
                # added to iter_q above, just correct for the offset
                lifted_start = iter_q - (offset_in * self.orientation)
                aln_block_length += (step - offset_in)
                if move in [CIGARstep.IDENT, CIGARstep.BOTH]:
                    aln_block_ident += (step - offset_in)
                    if move == CIGARstep.IDENT:
                        max_ident_block = max(max_ident_block, step - offset_in)

            if lifted_start is not None:
                in_query_region = True

            if in_query_region and not (stepping_in or stepping_out):
                aln_block_length += step
                if move == CIGARstep.IDENT:
                    aln_block_ident += step
                    max_ident_block = max(max_ident_block, step)

            elif in_query_region and stepping_out:
                aln_block_length += (step - offset_out)
                if move == CIGARstep.IDENT:
                    aln_block_ident += (step - offset_out)
                    max_ident_block = max(max_ident_block, step - offset_out)
                lifted_end = iter_q
                # the following accounts for overshooting in case where
                # alignments start exactly at the sequence/contig start
                if self.orientation > 0 or (self.orientation < 0 and lifted_end > 0):
                    lifted_end -= offset_out
                else:
                    lifted_end += offset_out
                break

            elif stepping_over:
                break

        if lifted_start is None:
            # Handle special cases as described in the docstring.
            # Now, the following must hold, otherwise
            # lifted_start should have been set:
            assert (iter_t - step) < find_start and stepping_over
            if move in [CIGARstep.IDENT, CIGARstep.BOTH]:
                # A step in both coordinate spaces suggests
                # that the find range represents a small region,
                # i.e. case (i) in the docstring

                # unroll the step
                iter_t -= step
                iter_q -= step
                # set state-preserving values to
                # what they were before
                self.last_start_cigar = cigar_pos
                self.last_start_target = iter_t
                self.last_start_query = iter_q

                # now manually define the coordinate lift
                # by adjusting the iter_q value
                find_range_size = find_end - find_start
                assert find_range_size > 0
                # not that iter_t has been reset above
                offset_start = find_start - iter_t
                assert offset_start > 0, offset_start
                lifted_start = iter_q + (offset_start)
                lifted_end = lifted_start + (find_range_size)
                if not in_query_region:
                    # this would imply no alignment block / identity block
                    # length info has been collected.
                    # For (very) small regions, this is an overoptimistic
                    # educated guess, but damage should be limited.
                    aln_block_length = find_range_size
                    aln_block_ident = find_range_size
                    max_ident_block = find_range_size
            elif move == CIGARstep.TARGET:
                # A step only in the target coordinate space
                # suggests that the region has been deleted
                # in the query and can thus not be lifted,
                # i.e. case (ii) in the docstring
                iter_t -= step
                self.last_start_cigar = cigar_pos
                self.last_start_target = iter_t
                self.last_start_query = iter_q
                # special values indicating no liftover
                lifted_start = -1
                lifted_end = -1
                deleted_region = True
            else:
                # This cannot have happened because find range
                # works in target coordinate space and stepping
                # over boundaries in that space is the only way
                # of triggering an end of the iteration
                raise RuntimeError(f"Unexpected - query overstepping with CIGAR OP: {step} / {move}")

        assert lifted_start is not None, f"{iter_t} - {find_start}:{find_end} / {stepping_in}"
        assert lifted_end is not None, f"{iter_t} - {find_start}:{find_end} / {stepping_out} / {stepping_over}"

        if self.orientation < 0 and not deleted_region:
            # the following is a workaround for situations where (reverse) alignments
            # start at the beginning of the query sequence, which (rarely) leads to lifted
            # coordinates that run out of the negatives and are thus already positive
            if lifted_start < 0 and lifted_end < 0:
                lifted_start, lifted_end = lifted_end * -1, lifted_start * -1
            else:
                tmp_end = lifted_end
                if lifted_start < 0:
                    lifted_end = lifted_start * -1
                else:
                    lifted_end = lifted_start
                if tmp_end < 0:
                    lifted_start = tmp_end * -1
                else:
                    lifted_start = tmp_end
                assert lifted_start < lifted_end
        if deleted_region:
            pass
        else:
            assert lifted_start >= 0, f"{iter_t} - {iter_q}: {lifted_start} -- {lifted_end} // {find_start} -- {find_end}"
            assert lifted_end > 0, f"{iter_t} - {iter_q} - {find_end}"
            assert lifted_start < lifted_end, f"{lifted_start} - {lifted_end}"
        return lifted_start, lifted_end, aln_block_ident, aln_block_length, max_ident_block


def parse_command_line():

    parser = argp.ArgumentParser(prog=___prog__, description=__description__)

    parser.add_argument(
        "--input-paf", "--norm-paf", "-paf", "-in",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="paf_alignments",
        required=True
    )

    parser.add_argument(
        "--input-bed", "--region-labels", "-labels",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="region_labels",
        required=True
    )

    parser.add_argument(
        "--output-bed", "-bed", "-out",
        type=lambda fp: pl.Path(fp).resolve(),
        dest="output",
        required=True
    )

    parser.add_argument(
        "--version", "-v",
        action="version"
    )

    args = parser.parse_args()

    return args


def read_region_label_file(file_path):

    regions = pd.read_csv(file_path, sep="\t", header=0)
    # need at least: seq - start - end - name
    assert len(regions.columns) > 3
    regions = regions[regions.columns[:4]]
    regions.columns = ["chrom", "start", "end", "name"]
    regions["size"] = regions["end"] - regions["start"]
    return regions


def check_sequence_compatibility(alignments, labeled_regions):

    n_seq_aln = alignments["target_name"].nunique()
    n_seq_regions = labeled_regions[labeled_regions.columns[0]].nunique()

    if n_seq_aln == n_seq_regions == 1:
        # lenient behavior: for a single name, this script can deal with
        # different naming patterns
        alignments["target_name"] = alignments["target_name"].replace(
            {alignments.target_name.iloc[0]: labeled_regions.chrom.iloc[0]},
            inplace=False
        )
    else:
        # if more than one sequence name, the names must match exactly
        seq_aln_names = alignments["target_name"].unique()
        seq_regions_names = labeled_regions[labeled_regions.columns[0]].unique()

        if set(seq_aln_names).union(set(seq_regions_names)) != set(seq_aln_names):
            raise ValueError(
                "Sequences in alignments and labeled regions do not match."
            )
    return


def read_paf_alignment_file(file_path):
    """
    Reads a normalized PAF alignment file and returns
    a DataFrame with the alignments. Hence, the PAF file
    must have a proper header.
    """
    alignments = pd.read_csv(file_path, sep="\t", header=0)
    alignments.sort_values(by=["target_name", "target_start"], inplace=True)
    alignments.reset_index(drop=True, inplace=True)
    return alignments


def compute_overlaps(alignments, labeled_regions):

    region_labels = pr.from_dict(
        {
            "Chromosome": labeled_regions.chrom,
            "Start": labeled_regions.start,
            "End": labeled_regions.end,
            "Name": labeled_regions.name,
            "region_idx": labeled_regions.index.values
        }
    )

    aligned_blocks = pr.from_dict(
        {
            "Chromosome": alignments.target_name,
            "Start": alignments.target_start,
            "End": alignments.target_end,
            "align_idx": alignments.index.values
        }
    )

    broken_alignments = aligned_blocks.intersect(region_labels, how=None)
    # NB here: the 'how' parameter is set to None
    # to ensure that all intersections are returned,
    # i.e. overlapping regions in the source BED
    # file / region annotation may have various overlaps.
    # Think of it as reporting all left-with-right intersections
    # plus all right-with-left intersections.
    # We need to select an intersect only if it does not extend beyond
    # the original source region size.
    # This is done in a scope outside of this function.

    # bring in the the region label information with an
    # additional suffix to avoid name clashes
    broken_alignments = broken_alignments.join(region_labels, suffix="_region_label")

    # subsequent operations require a DataFrame
    # also NB: by construction, the dataframe is sorted by
    # the coordinates of the region labels / target sequence
    broken_alignments = broken_alignments.df
    return broken_alignments


def select_best_fit_overlaps(align_label_overlaps):

    row_select = []
    priority_matches = set()
    for row in align_label_overlaps.itertuples():
        equal_start = row.Start == row.Start_region_label
        smaller_equal_end = row.End <= row.End_region_label
        larger_start = row.Start > row.Start_region_label
        lookup_key = (row.align_idx, row.region_idx)
        if equal_start & smaller_equal_end:
            # perfect match, potantially cut off because alignment ends
            row_select.append(row.Index)
            priority_matches.add(lookup_key)
            continue
        elif larger_start & smaller_equal_end:
            # could be a secondary containment, select only
            # if no priority match exists for this combination
            if lookup_key in priority_matches:
                continue
            row_select.append(row.Index)
            continue
        else:
            # this has the potential of being a dangerous "skip"
            continue

    align_label_overlaps = align_label_overlaps.loc[row_select, :]

    return align_label_overlaps


def perform_liftover(align_label_overlaps, alignments, labeled_regions, min_alignment_threshold):

    last_alignment = None
    lifted_regions = []
    for row in align_label_overlaps.itertuples():
        alignment = alignments.loc[row.align_idx]
        if last_alignment is None or last_alignment != row.align_idx:
            last_alignment = row.align_idx
            cw = CIGARwalker(
                alignment.target_start,
                alignment.target_end,
                alignment.query_start,
                alignment.query_end,
                alignment.align_orient,
                alignment.cg_cigar,
                row.align_idx
            )
        # NB: this needs to check the alignment target coordinates and not the
        # start/end coordinates in the labeled region intersections, which may
        # only cover a portion of the respective alignment
        if alignment.target_start >= row.Start_region_label and alignment.target_end <= row.End_region_label:
            # special case: the alignment is smaller than the labeled region
            # and thus we can simply walk through the CIGAR string
            align_size = alignment.target_end - alignment.target_start
            if align_size < min_alignment_threshold:
                continue
            lifted_block = cw.exhaustive_walk()
        else:
            find_start = max(row.Start, row.Start_region_label)
            find_end = min(row.End, row.End_region_label)
            range_size = find_end - find_start
            if range_size < min_alignment_threshold:
                continue
            try:
                lifted_block = cw.find_range(find_start, find_end)
            except AssertionError as e:
                print(row)
                print(alignment)
                raise

        region_label = labeled_regions.at[row.region_idx, "name"]
        lifted_seq = alignment.query_name
        lifted_start, lifted_end = lifted_block[:2]
        if lifted_start == -1 and lifted_end == -1:
            # region has been deleted in the query
            continue
        aln_matching, aln_block, max_ident = lifted_block[2:]
        # for pct_matching, we need the region size
        region_size = row.End_region_label - row.Start_region_label
        pct_matching = round(aln_matching / region_size * 100, 2)

        lifted_regions.append(
            (
                lifted_seq, lifted_start, lifted_end,
                region_label, 0, alignment.align_orient,
                region_label, aln_matching, aln_block, max_ident,
                pct_matching
            )
        )

    lifted_regions = pd.DataFrame(
        lifted_regions,
        columns=[
            "chrom", "start", "end", "name", "score", "strand",
            "label", "aln_matching", "aln_block", "longest_id_block",
            "pct_matching"
        ]
    )
    lifted_regions["strand"] = lifted_regions["strand"].replace(
        {1: "+", -1: "-"}
    ).astype(str)

    lifted_regions["rank_matching"] = lifted_regions["pct_matching"].rank(
        method="average", ascending=True, pct=True
    )
    lifted_regions["rank_block"] = lifted_regions["longest_id_block"].rank(
        method="average", ascending=True, pct=True
    )

    # this will be used to replace the score column in the BED output
    lifted_regions["mean_rank"] = ((
        lifted_regions["rank_matching"] + lifted_regions["rank_block"]
    ) / 2 * 1000).round(0).astype(int)

    lifted_regions.sort_values(by=["chrom", "start"], inplace=True)

    return lifted_regions


def cluster_lifted_regions(lifted_regions):

    intervals = pr.from_dict(
        {
            "Chromosome": lifted_regions.chrom,
            "Start": lifted_regions.start,
            "End": lifted_regions.end,
            "Name": lifted_regions.name,
            "Strand": lifted_regions.strand,
            "pd_idx": lifted_regions.index.values
        }
    )
    intervals = intervals.cluster(strand="same", by="Name")
    cluster_ids = pd.Series(
        intervals.Cluster.values,
        index=intervals.pd_idx.values,
        name="cluster_id"
    )
    _n_rows = lifted_regions.shape[0]
    # little sanity check
    lifted_regions = lifted_regions.merge(cluster_ids, left_index=True, right_index=True)
    assert lifted_regions.shape[0] == _n_rows

    return lifted_regions


def derive_bed_rows(lifted_regions):

    bed_rows = []
    for cluster_id, regions in lifted_regions.groupby("cluster_id"):
        assert regions.chrom.nunique() == 1
        assert regions.name.nunique() == 1
        assert regions.strand.nunique() == 1

        # collect values...
        seq = regions.chrom.iloc[0]
        start = regions.start.min()
        end = regions.end.max()
        assert start < end, f"{start} - {end}: {regions}"
        name = regions.name.iloc[0]
        score = regions.mean_rank.max()
        assert 0 <= score <= 1000, f"{regions}"
        max_pct_matching = regions.pct_matching.max()
        max_block_id = regions.longest_id_block.max()
        new_label = f"{name}[IDPCT:{max_pct_matching}|IDBLK:{max_block_id}]"
        bed_rows.append(
            (
                seq, start, end, new_label, score,
                name, max_pct_matching, max_block_id, regions.shape[0]
            )
        )
    bed_rows = pd.DataFrame(
        bed_rows,
        columns=[
            "chrom", "start", "end", "name", "score",
            "label", "max_pct_matching", "longest_id_block",
            "merged_alignments"
        ]
    )
    bed_rows.sort_values(by=["chrom", "start"], inplace=True)

    return bed_rows


def main():

    args = parse_command_line()

    labeled_regions = read_region_label_file(args.region_labels)

    min_alignment_threshold = labeled_regions["size"].min() // 2

    alignments = read_paf_alignment_file(args.paf_alignments)

    check_sequence_compatibility(alignments, labeled_regions)

    align_label_overlaps = compute_overlaps(alignments, labeled_regions)

    align_label_overlaps = select_best_fit_overlaps(align_label_overlaps)

    lifted_regions = perform_liftover(align_label_overlaps, alignments, labeled_regions, min_alignment_threshold)

    lifted_regions = cluster_lifted_regions(lifted_regions)

    bed_rows = derive_bed_rows(lifted_regions)

    with xopen.xopen(args.output, "wt") as out_bed:
        out_bed.write("#")
        bed_rows.to_csv(
            out_bed,
            sep="\t",
            header=True,
            index=False,
        )

    return 0


if __name__ == "__main__":
    main()
