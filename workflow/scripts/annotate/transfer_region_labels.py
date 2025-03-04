#!/usr/bin/env python3

import argparse as argp
import collections as col
import enum
import itertools as itt
import pathlib as pl
import re

import pandas as pd
import pandas.errors as pderr
import numpy as np
import numpy.ma as msk

class ShrinkingError(Exception):
    pass

class ClosedRegionError(Exception):
    pass

LabelRegion = col.namedtuple(
    "LabelRegion",
    ["start", "end", "name", "score", "orient", "size", "row_idx", "label_idx"]
)

LabelInterval = col.namedtuple(
    "LabelInterval",
    ["start", "end", "label", "label_idx"]
)

class CIGARstep(enum.Enum):
    TARGET = 0
    QUERY = 1
    BOTH = 2

class CIGARwalker:
    __slots__ = (
        "cigar", "cigar_ops", "move_ops", "orientation",
        "tstart", "tend", "qstart", "qend",
        "target_regions", "query_regions",
        "row_idx"
    )

    def __init__(self, tstart, tend, qstart, qend, orientation, cigar, row_idx):

        self.cigar_ops = re.compile("[0-9]+(\=|X|D|I|M)")
        self.cigar = cigar
        self.move_ops = {
            "=": CIGARstep.BOTH,
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
            self.qstart = qend * -1
            self.qend = qstart * -1
        else:
            self.qstart = qstart
            self.qend = qend
        assert self.qstart < self.qend
        self.target_regions = []
        self.query_regions = []
        self.row_idx = row_idx

        return None

    def walk_cigar(self):

        for cigar_op in self.cigar_ops.finditer(self.cigar):
            tmp = cigar_op.group(0)
            move_op = self.move_ops[tmp[-1]]
            step = int(tmp[:-1])
            yield step, move_op
        return

    def build_region(self, start, end, label, label_idx):

        if start < 0 or end < 0:
            assert self.orientation < 0
            start = abs(start)
            end = abs(end)
            if start < end:
                pass
            else:
                start, end = end, start

        lr = LabelRegion(
            start, end, label, 1000,
            self.orientation, end - start,
            self.row_idx, label_idx
        )
        assert lr.start < lr.end, f"{lr} / {self.row_idx}"
        return lr

    def label_intervals(self, labeled_intervals):

        iter_t = self.tstart
        last_t = self.tstart
        iter_q = self.qstart
        last_q = self.qstart

        current_interval = labeled_intervals[0]
        next_interval = 1
        assert iter_t < current_interval.end

        count_ops = col.Counter()
        for (step, move_op) in self.walk_cigar():
            count_ops[move_op] += step
            if iter_t > current_interval.end:
                # if the alignment balance is way off,
                # skip the region either in the target
                # or in the query
                target_over_query_ratio = self.check_balanced_alignment(count_ops)
                # just shorthand...
                balance = target_over_query_ratio
                while 1:
                    if iter_t < current_interval.end:
                        count_ops = col.Counter()
                        break
                    step_back = iter_t - current_interval.end
                    assert step_back > 0
                    self.target_regions.append(
                        self.build_region(
                            last_t, current_interval.end,
                            current_interval.label,
                            current_interval.label_idx
                        )
                    )
                    last_t = current_interval.end

                    # tricky: cigar ops commonly advance
                    # both target and query, but this is
                    # still not necessarily accurate
                    if balance > 1.5:
                        # do not create this region in the query;
                        # many more 'move' ops for the target than
                        # for the query suggest INS/DUP in target
                        last_q = iter_q
                    else:
                        query_end = iter_q - step_back
                        assert last_q < query_end
                        self.query_regions.append(
                            self.build_region(
                                last_q, query_end,
                                current_interval.label,
                                current_interval.label_idx
                            )
                        )
                        last_q = query_end
                    current_interval = labeled_intervals[next_interval]
                    next_interval += 1
            if move_op == CIGARstep.TARGET:
                iter_t += step
            elif move_op == CIGARstep.QUERY:
                iter_q += step
            else:
                iter_t += step
                iter_q += step

        # add last region for both
        self.target_regions.append(
            self.build_region(
                last_t, iter_t, current_interval.label,
                current_interval.label_idx
            )
        )
        self.query_regions.append(
            self.build_region(
                last_q, iter_q, current_interval.label,
                current_interval.label_idx
            )
        )

        assert iter_t == self.tend, iter_t
        assert abs(iter_q) == abs(self.qend), iter_q

        return

    def get_labeled_regions(self, which):

        if which == "query":
            return sorted(self.query_regions)
        elif which == "target":
            return sorted(self.target_regions)
        else:
            raise ValueError(which)

    def check_balanced_alignment(self, count_ops):

        target = count_ops[CIGARstep.TARGET] + count_ops[CIGARstep.BOTH]
        query = count_ops[CIGARstep.QUERY] + count_ops[CIGARstep.BOTH]
        return target/query


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        "--input", "--norm-paf", "-i",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="norm_paf",
        required=True
    )

    parser.add_argument(
        "--target-labels", "--region-labels",
        "-label-bed", "-l",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="region_labels",
        required=True
    )

    parser.add_argument(
        "--fail-on-empty", "-e",
        action="store_true",
        default=False,
        dest="fail_on_empty",
        help="Fail on empty input instead of generating empty output."
    )

    parser.add_argument(
        "--min-shrinking-fraction", "-f",
        type=float,
        default=0.25,
        dest="min_shrinking_fraction",
        help=(
            "If an aligned region is shrunk, discard it if "
            "the new size is less than this fraction of the "
            "original size. Default: 0.25"
        )
    )

    parser.add_argument(
        "--output", "--out-region-annotation",
        "--bed-out", "-o",
        type=lambda fp: pl.Path(fp).resolve(),
        dest="out_regions",
        required=True,
        help="Output file for annotated regions (6+1 BED-like format)."
    )

    parser.add_argument(
        "--close-gap-size", "-c",
        type=int,
        default=0,
        required=False
    )

    parser.add_argument(
        "--simplified-output", "-s",
        type=lambda fp: pl.Path(fp).resolve(),
        default=None,
        required=False,
        help="Output file for simplified region annotation w/ closed gaps."
    )

    parser.add_argument(
        "--dump-debug-output", "-d",
        "--debug-out", "--debug",
        action="store_true",
        default=False,
        dest="debug_out"
    )

    args = parser.parse_args()

    return args


def read_normalized_paf_file(file_path):

    try:
        table = pd.read_csv(file_path, sep="\t", header=0)
        # reindex starting at one so tracking the region coverage
        # can use value 0 as "not covered"
        # Throughout this script, it is vital to not change the index!
        table.set_index(
            np.arange(1, table.shape[0]+1, dtype=int),
            inplace=True, verify_integrity=True
        )
    except pderr.EmptyDataError:
        table = None

    return table


def check_open_region(target_region, region_cover, min_shrink_fraction):

    select_uncovered = region_cover[target_region.start:target_region.end] == 0
    num_uncovered = select_uncovered.sum()
    if num_uncovered > 0:
        if num_uncovered >= target_region.size:
            new_target_region = [target_region]
        else:
            # try to shrink target region
            # to fit between to other covered regions
            fraction = num_uncovered / target_region.size
            if fraction < min_shrink_fraction:
                raise ShrinkingError
            else:
                region = np.arange(target_region.start, target_region.end, dtype=int)
                assert region[0] == target_region.start
                assert region[-1] == target_region.end - 1
                new_start = region[select_uncovered].min()
                assert new_start >= target_region.start
                # new_end + 1: correct for half-open indexing in numpy;
                # right-exclusive indexing; converting from position to
                # index
                new_end = region[select_uncovered].max() + 1
                assert new_end <= target_region.end

                if not (region_cover[new_start:new_end] == 0).all():
                    # inverting mask => mask everything that is already covered
                    split_alignment = msk.MaskedArray(
                        data=region,
                        mask=~select_uncovered
                    )
                    new_target_region = []
                    old_size = target_region.size
                    for split_region in msk.clump_unmasked(split_alignment):
                        new_start = region[split_region].min()
                        assert new_start >= target_region.start
                        new_end = region[split_region].max() + 1
                        assert new_end <= target_region.end
                        assert (region_cover[new_start:new_end] == 0).all()
                        new_size = new_end - new_start
                        region_score = round(new_size / old_size * 1000)
                        assert region_score <= 1000
                        # ["start", "end", "name", "score", "orient", "size", "row_idx", "label_idx"]
                        new_target_region.append(
                            LabelRegion(
                                new_start, new_end, target_region.name,
                                region_score, target_region.orient,
                                new_size, target_region.row_idx, target_region.label_idx
                            )
                        )
                else:
                    old_size = target_region.size
                    new_size = new_end - new_start
                    assert new_size > 0
                    region_score = round(new_size / old_size * 1000)
                    assert region_score <= 1000
                    new_target_region = [
                        LabelRegion(
                            new_start, new_end, target_region.name,
                            region_score, target_region.orient,
                            new_size, target_region.row_idx, target_region.label_idx
                        )
                    ]
    else:
        raise ClosedRegionError

    return new_target_region


def make_process_record(process_log, indices, anchor_idx, anchor_label, other_label):

    for i in indices:
        assert i not in process_log, f"Collision: {indices} / anchor {anchor_idx}"
        if i == anchor_idx:
            process_log[i] = anchor_label
        else:
            process_log[i] = other_label
    return


def _dump_process_log(process_log, subset=None, reorder=False):

    keys = sorted(process_log.keys()) if reorder else process_log.keys()

    for row_idx in keys:
        if subset is not None:
            if row_idx not in subset:
                continue
        print(row_idx, " --- ", process_log[row_idx])
    return


def get_overlaps_by_region(regions, seqname, interval):
    overlapping = regions.loc[
        regions[regions.columns[0]] == seqname, :
    ].index.overlaps(interval)
    ovl_labels = []
    if overlapping.any():
        ovl_labels = sorted(
            LabelInterval(region.start, region.end, region.name, region.enum_id)
            for region in regions.loc[overlapping, :].itertuples()
        )
    return ovl_labels


def get_overlaps(regions, seqname, interval):
    overlapping = regions.index.overlaps(interval)
    ovl_labels = []
    if overlapping.any():
        ovl_labels = sorted(
            LabelInterval(region.start, region.end, region.name, region.enum_id)
            for region in regions.loc[overlapping, :].itertuples()
        )
    return ovl_labels


def find_region_cover(paf, target_seq, labels, generic_label_lookup, min_shrink_fraction, region_scores, process_log):
    """This function processes all alignments for a target sequence
    ordered by priority. The priority is computed as follows:
    - add new binary column: mapq_nonzero (split alignments into MAPQ == 0 and MAPQ > 0)
    - sort by
    --- mapq_nonzero (desc)
    --- mapq (desc)
    --- align_matching (desc)
    --- align_total (asc)
    --- target_start (asc)
    ---> assign priority from high to low following sort order
    """
    target_subset = paf.loc[
        paf["query_name"] == target_seq, :
    ].copy()
    target_subset["mapq_nonzero"] = 0
    target_subset.loc[target_subset["mapq"] > 0, "mapq_nonzero"] = 1
    target_subset.sort_values(
        ["mapq_nonzero", "mapq", "align_matching", "align_total", "query_start"],
        ascending=[False, False, False, True, True],
        inplace=True
    )
    target_subset["priority"] = np.flip(np.arange(0, target_subset.shape[0], dtype=int))

    # track what has been covered already
    region_cover = np.zeros(target_subset["query_length"].iloc[0], dtype=int)

    min_shrink_score = round(min_shrink_fraction * 1000)

    if generic_label_lookup:
        overlaps = get_overlaps
    else:
        overlaps = get_overlaps_by_region

    labeled_regions = []

    for aln_row in target_subset.itertuples():
        # this now iterates from high to low priority
        aln_iv = pd.Interval(aln_row.target_start, aln_row.target_end, closed="left")
        ovl_labels = overlaps(labels, aln_row.target_name, aln_iv)
        if ovl_labels:
            cw = CIGARwalker(
                aln_row.target_start,
                aln_row.target_end,
                aln_row.query_start,
                aln_row.query_end,
                aln_row.align_orient,
                aln_row.cg_cigar,
                aln_row.Index
            )
            try:
                cw.label_intervals(ovl_labels)
            except AssertionError:
                print(ovl_labels)
                raise

            for labeled_region in cw.get_labeled_regions("query"):

                try:
                    accept_regions = check_open_region(
                        labeled_region, region_cover, min_shrink_fraction)
                except ShrinkingError:
                    pass
                except ClosedRegionError:
                    pass
                else:
                    for accept_region in accept_regions:
                        region_cover[
                            accept_region.start:accept_region.end
                        ] = accept_region.label_idx
                        labeled_regions.append(accept_region)

    return region_cover, labeled_regions


def merge_regions(target_seq, regions):

    start = regions[0].start
    end = regions[-1].end
    cum_score = sum(region.score for region in regions)
    avg_score = int(round(cum_score/len(regions), 0))
    label = regions[0].name
    orient = "+" if regions[0].orient > 0 else "-"
    row_idx = regions[0].row_idx
    return target_seq, start, end, label, avg_score, orient, row_idx


def produce_region_annotation(target_seq, region_cover, labeled_regions):

    masked_regions = msk.masked_not_equal(region_cover, 0)
    # all _covered_ regions are now masked

    sequence_regions = []
    for slice in msk.clump_unmasked(masked_regions):
        sequence_regions.append(
            (target_seq, slice.start, slice.stop, "uncertain", 0, "+")
        )

    labeled_regions = col.deque(sorted(labeled_regions))
    labeled_regions.append(None)

    active = []
    while 1:
        region = labeled_regions.popleft()
        if region is None:
            break
        if not active:
            active.append(region)
            continue
        if region.label_idx == active[-1].label_idx:
            if region.start == active[-1].end:
                active.append(region)
                continue
        sequence_regions.append(
            merge_regions(target_seq, active)
        )
        active = [region]

    if active:
        sequence_regions.append(
            merge_regions(target_seq, active)
        )

    return sequence_regions


def read_region_labels(bed_like):

    labels = pd.read_csv(bed_like, sep="\t", header=0)
    sort_columns = [labels.columns[0], labels.columns[1]]
    labels.sort_values(sort_columns, inplace=True)
    labels["enum_id"] = np.arange(1, labels.shape[0]+1, dtype=int)
    iv_idx = pd.IntervalIndex.from_tuples(
        [(row.start, row.end) for row in labels.itertuples()],
        closed="left"
    )
    labels.set_index(iv_idx, inplace=True)
    return labels


def sanity_check_sequence_names(alignments, labels, args):

    num_label_seq = labels[labels.columns[0]].nunique()
    num_aln_seq = alignments["target_name"].nunique()
    if num_label_seq == num_aln_seq == 1:
        generic_region = True
    else:
        label_seq_names = set(labels[labels.columns[0]].unique())
        aln_seq_names = set(alignments["target_name"].unique())
        if len(label_seq_names.intersection(aln_seq_names)) == 0:
            err_msg = (
                "ERROR:\n"
                f"Region label file (BED) {args.region_labels} "
                f"and alignment file (PAF) {args.norm_paf} "
                "do not share reference/target "
                "sequence names - files seem incompatible."
            )
            raise RuntimeError(err_msg)
        else:
            generic_region = False
    return generic_region


def main():

    process_log = col.OrderedDict()
    region_scores = dict()

    args = parse_command_line()

    paf = read_normalized_paf_file(args.norm_paf)

    if paf is None and args.fail_on_empty:
        raise RuntimeError(f"Empty input PAF: {args.norm_paf}")
    elif paf is None:
        out_regions = pd.DataFrame(
            [],
            columns=["#chrom", "start", "end", "name", "score", "strand", "anchor_row"]
        )
    else:
        labels = read_region_labels(args.region_labels)
        generic_label_lookup = sanity_check_sequence_names(paf, labels, args)

        out_regions = []
        for query_seq in paf["query_name"].unique():
            if query_seq != "chrY_HG01433-J1_DACBCF94":
                continue

            region_cover, labeled_regions = find_region_cover(
                paf, query_seq, labels, generic_label_lookup,
                args.min_shrinking_fraction,
                region_scores, process_log
            )
            sequence_regions = produce_region_annotation(query_seq, region_cover, labeled_regions)
            out_regions.extend(sequence_regions)

        out_regions = pd.DataFrame.from_records(
            out_regions,
            columns=["#chrom", "start", "end", "name", "score", "strand", "anchor_row"]
        )
        out_regions["anchor_row"] = out_regions["anchor_row"].fillna(0, inplace=False)
        out_regions["anchor_row"] = out_regions["anchor_row"].astype(int)
        out_regions.sort_values(["#chrom", "start", "end"], inplace=True)
        print(out_regions.head(50))
        raise

    args.out_regions.parent.mkdir(exist_ok=True, parents=True)
    out_regions.to_csv(args.out_regions, sep="\t", header=True, index=False)

    if False:
        simplify_region_annotation(out_regions, 10000)

    if args.debug_out and not paf is None:
        paf.insert(0, "row_idx_process_log", [process_log[i] for i in paf.index])
        if args.out_regions.suffix == ".gz":
            debug_out = args.out_regions.stem
        else:
            debug_out = args.out_regions
        debug_out = debug_out.with_suffix(".process-debug-log.tsv.gz")
        paf.index.name = "paf_row_idx"
        paf.to_csv(debug_out, sep="\t", header=True, index=True)

    return 0


if __name__ == "__main__":
    main()
