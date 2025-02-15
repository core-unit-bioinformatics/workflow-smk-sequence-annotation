#!/usr/bin/env python3

import argparse as argp
import collections as col
import itertools as itt
import pathlib as pl

import pandas as pd
import numpy as np
import numpy.ma as msk

class ShrinkingError(Exception):
    pass

class ClosedRegionError(Exception):
    pass

TargetRegion = col.namedtuple(
    "TargetRegion",
    ["start", "end", "orient", "size", "score", "anchor_idx"]
)


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        "--input", "--norm-paf", "-i",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="norm_paf",
        required=True
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
        "--bed-out",
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

    table = pd.read_csv(file_path, sep="\t", header=0)
    # reindex starting at one so tracking the region coverage
    # can use value 0 as "not covered"
    # Throughout this script, it is vital to not change the index!
    table.set_index(
        np.arange(1, table.shape[0]+1, dtype=int),
        inplace=True, verify_integrity=True
    )

    return table


def select_anchor_subset(subset, active_rows):

    anchor_region_idx = subset.loc[active_rows, "align_matching"].idxmax()
    anchor_name, anchor_size = subset.loc[
        anchor_region_idx, ["query_name", "query_length"]
    ].values

    select_name = subset["query_name"] == anchor_name
    select_size = subset["query_length"] == anchor_size

    selector = select_name & select_size

    anchor_subset = subset.loc[selector, :]

    return anchor_subset, anchor_region_idx


def determine_target_alignment_overlap(anchor_subset, anchor_idx):
    """The overlap proceeds in a greedy fashion, which is mostly
    relevant in highly repetitive regions where alignments tend
    to be fragmented/scattered.
    """

    # anchor orientation always determines global orientation
    orient = anchor_subset.loc[anchor_idx, "align_orient"]

    start = anchor_subset.loc[anchor_idx, "target_start"]
    end = anchor_subset.loc[anchor_idx, "target_end"]

    selected_regions = set([anchor_idx])

    while 1:

        select_start = anchor_subset["target_start"] < end
        select_end = anchor_subset["target_end"] > start

        select_overlapping = select_start & select_end

        merge_regions = anchor_subset.loc[select_overlapping, :]
        if merge_regions.shape[0] > len(selected_regions):
            start = merge_regions["target_start"].min()
            end = merge_regions["target_end"].max()
            selected_regions = selected_regions.union(
                set(merge_regions.index.values)
            )
        elif merge_regions.shape[0] == len(selected_regions):
            # reached steady state
            tstart = merge_regions["target_start"].min()
            tend = merge_regions["target_end"].max()
            processed_rows = set(merge_regions.index.values)
            break
        else:
            raise

    return tstart, tend, orient, processed_rows


def process_anchor_subset(anchor_subset, anchor_idx):

    if anchor_subset.shape[0] > 1:
        tstart, tend, orient, processed_rows = determine_target_alignment_overlap(
            anchor_subset, anchor_idx
        )
    else:
        tstart = anchor_subset["target_start"].iloc[0]
        tend = anchor_subset["target_end"].iloc[0]
        orient = anchor_subset["align_orient"].iloc[0]
        # turn into set b/c of call above to
        # determine_target_alignment_overlap
        processed_rows = set([anchor_idx])

    region_size = tend - tstart
    assert region_size > 0, f"{tstart} / {tend}"
    target_region = TargetRegion(tstart, tend, orient, region_size, 1000, anchor_idx)

    return target_region, processed_rows


def check_open_region(target_region, region_cover, min_shrink_fraction):

    select_uncovered = region_cover[target_region.start:target_region.end] == 0
    num_uncovered = select_uncovered.sum()
    if num_uncovered > 0:
        if num_uncovered >= target_region.size:
            new_target_region = target_region
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
                new_end = region[select_uncovered].max()
                assert new_end <= target_region.end

                assert (region_cover[new_start:new_end] == 0).all(), f"Fragmented alignment: {target_region}"

                old_size = target_region.size
                new_size = new_end - new_start
                assert new_size > 0
                region_score = round(new_size / old_size * 1000)
                assert region_score <= 1000
                new_target_region = TargetRegion(
                    new_start, new_end, target_region.orient,
                    new_size, region_score, target_region.anchor_idx
                )
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


def find_region_cover(paf, target_seq, min_shrink_fraction, region_scores, process_log):

    # The subsequent code makes use of idxmax
    # to select anchor regions. idxmax returns
    # the _first_ occurrence of the maximal value
    # in a column. Hence, we sort here by
    # >>> align_matching, mapq, align_total, target_start
    # to enable breaking ties in that order.
    # If alignments tie nevertheless, the result
    # is defined by Panda's sorting algorithm
    target_subset = paf.loc[
        paf["target_name"] == target_seq, :
    ].copy()
    target_subset.sort_values(
        ["align_matching", "mapq", "align_total", "target_start"],
        ascending=[False, False, True, True],
        inplace=True
    )

    # track what has been covered already
    region_cover = np.zeros(target_subset["target_length"].iloc[0], dtype=int)

    active_rows = target_subset.index.values
    processed_rows = set()

    min_shrink_score = round(min_shrink_fraction * 1000)

    while 1:

        anchor_subset, anchor_idx = select_anchor_subset(
            target_subset, active_rows
        )
        target_region, proc_rows = process_anchor_subset(
            anchor_subset, anchor_idx
        )
        try:
            target_region = check_open_region(
                target_region, region_cover, min_shrink_fraction
            )
        except ShrinkingError:
            make_process_record(
                process_log, proc_rows, anchor_idx,
                f"SKIP-ANCHOR-SHRINK-{anchor_idx}", f"SKIP-MERGE-SHRINK-{anchor_idx}"
            )
            for p in proc_rows:
                region_scores[p] = min_shrink_score
        except ClosedRegionError:
            make_process_record(
                process_log, proc_rows, anchor_idx,
                f"SKIP-ANCHOR-CLOSED-{anchor_idx}", f"SKIP-MERGE-CLOSED-{anchor_idx}"
            )
            for p in proc_rows:
                region_scores[p] = 0
        else:
            # multiply by orientation -> revs are < 0
            region_cover[
                target_region.start:target_region.end
            ] = target_region.anchor_idx * target_region.orient
            region_scores[anchor_idx] = target_region.score
            make_process_record(
                process_log, proc_rows, anchor_idx,
                f"SET-ANCHOR-{anchor_idx}", f"SET-MERGE-{anchor_idx}"
            )
        processed_rows = processed_rows.union(proc_rows)
        active_rows = sorted(
            set(active_rows) - processed_rows
        )
        if len(active_rows) == 0:
            break

    return region_cover


def extract_region_label(idx, paf):
    label = paf.loc[abs(idx), "query_name"]
    return label


def produce_region_annotation(target_seq, region_cover, region_scores, paf):

    masked_regions = msk.masked_not_equal(region_cover, 0)
    # all _covered_ regions are now masked

    sequence_regions = []
    for slice in msk.clump_unmasked(masked_regions):
        sequence_regions.append(
            (target_seq, slice.start, slice.stop, "uncertain", 0, "+")
        )

    masked_regions = msk.masked_equal(region_cover, 0)
    # all _uncovered_ regions are now masked
    for slice in msk.clump_unmasked(masked_regions):
        sub = masked_regions[slice.start:slice.stop]
        idx = sub.data[0]
        if (sub.data == idx).all():
            region_label = extract_region_label(idx, paf)
            orient = "+" if idx > 0 else "-"
            sequence_regions.append(
                (
                    target_seq, slice.start, slice.stop,
                    region_label, region_scores[abs(idx)],
                    orient, abs(idx)
                )
            )
        else:
            uniq_values, uniq_starts = np.unique(sub.data, return_index=True)
            offset = slice.start
            for (start1, idx1), (start2, idx2) in itt.pairwise(sorted(zip(uniq_starts, uniq_values))):
                region_start = offset + start1
                region_end = offset + start2
                assert idx1 != idx2
                assert (masked_regions[region_start:region_end] == idx1).all()
                region_label = extract_region_label(idx1, paf)
                orient = "+" if idx1 > 0 else "-"
                sequence_regions.append(
                    (
                        target_seq, region_start, region_end,
                        region_label, region_scores[abs(idx1)],
                        orient, abs(idx1)
                    )
                )

    return sequence_regions


def adapt_end_entries(regions, adapt):
    if adapt == "start":
        regions.loc[1, "start"] = 0
        regions = regions.loc[regions.index[1:], :].copy()
    elif adapt == "end":
        last_index = regions.index[-1]
        next_to_last = regions.index[-2]
        regions.loc[next_to_last, "end"] = regions.loc[last_index, "end"]
        regions = regions.loc[regions.index[:-1], :].copy()
    else:
        raise
    return regions


def simplify_region_annotation(regions, gap_size_threshold):


    for seq in regions["#chrom"].unique():
        simplified_regions = []
        # NB: the regions dataframe has already been
        # sorted in the calling scope
        sub = regions.loc[regions["#chrom"] == seq, :].reset_index(drop=True, inplace=False)
        sub["size"] = sub["end"] - sub["start"]
        if sub["name"].iloc[0] == "uncertain" and sub["size"].iloc[0] < gap_size_threshold:
            sub = adapt_end_entries(sub, "start")
        if sub["name"].iloc[-1] == "uncertain" and sub["size"].iloc[-1] < gap_size_threshold:
            sub = adapt_end_entries(sub, "end")

        for name, name_regions in sub.groupby("name"):
            if name == "uncertain":
                continue
            if name_regions.shape[0] == 1:
                # drop region size and accept as is
                region_spec = tuple(name_regions.iloc[0].values[:-1])
                simplified_regions.append(region_spec)
            else:
                print(name_regions)
                raise



def main():

    process_log = col.OrderedDict()
    region_scores = dict()

    args = parse_command_line()

    paf = read_normalized_paf_file(args.norm_paf)

    out_regions = []
    for target_seq in paf["target_name"].unique():

        region_cover = find_region_cover(
            paf, target_seq, args.min_shrinking_fraction,
            region_scores, process_log
        )
        sequence_regions = produce_region_annotation(target_seq, region_cover, region_scores, paf)
        out_regions.extend(sequence_regions)

    out_regions = pd.DataFrame.from_records(
        out_regions,
        columns=["#chrom", "start", "end", "name", "score", "strand", "anchor_row"]
    )
    out_regions["anchor_row"] = out_regions["anchor_row"].fillna(0, inplace=False)
    out_regions["anchor_row"] = out_regions["anchor_row"].astype(int)
    out_regions.sort_values(["#chrom", "start", "end"], inplace=True)

    args.out_regions.parent.mkdir(exist_ok=True, parents=True)
    out_regions.to_csv(args.out_regions, sep="\t", header=True, index=False)

    if True:
        simplify_region_annotation(out_regions, 10000)

    if args.debug_out:
        pass

    return 0


if __name__ == "__main__":
    main()
