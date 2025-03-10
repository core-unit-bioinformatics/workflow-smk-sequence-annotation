#!/usr/bin/env python3

import argparse as argp
import collections as col
import itertools as itt
import pathlib as pl

import pandas as pd
import xopen

IVCLOSED = "left"

Region = col.namedtuple("Region", ["chrom", "start", "end", "name"])


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--annotations", "--tables",
        "-at", "-a", "-t",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        nargs="+",
        required=True,
        dest="annotation_tables"
    )

    parser.add_argument(
        "--annotation-labels", "--labels",
        "-l",
        type=str,
        nargs="+",
        required=True,
        dest="annotation_labels"
    )

    parser.add_argument(
        "--intersect", "-i",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        required=True,
        dest="intersect"
    )

    parser.add_argument(
        "--drop-lowconf",
        type=str,
        default=["uncertain"],
        nargs="*",
        dest="drop_lowconf"
    )

    parser.add_argument(
        "--output", "--out-bed",
        "-o",
        type=lambda fp: pl.Path(fp).resolve(strict=False),
        dest="output",
        required=True
    )

    args = parser.parse_args()

    return args


def check_header(filepath):

    with xopen.xopen(filepath) as table:
        first_line = table.readline().strip().split()
        try:
            _ = int(first_line[1])  # potential start
            skip_rows = 0
        except ValueError:
            skip_rows = 1
        except IndexError:
            # this should imply: empty input file
            skip_rows = 0
    return skip_rows


def load_annotation(filepath, label, drop_lowconf):

    skip_rows = check_header(filepath)

    df = pd.read_csv(
        filepath, sep="\t", header=None, skiprows=skip_rows,
        usecols=[0,1,2,3,4,5], names=["chrom", "start", "end", "name", "score", "strand"]
    )

    if df.empty:
        return dict()

    if len(drop_lowconf) > 0:
        label_bad = lambda label: any(lowconf in label for lowconf in drop_lowconf)
        df = df.loc[
            ~df["name"].apply(label_bad), :
        ].copy()

    #df["name"] += "|" + df["strand"]

    partitions = dict()
    for seq, sub in df.groupby("chrom"):
        ividx = pd.IntervalIndex.from_arrays(sub["start"], sub["end"], closed=IVCLOSED)
        s = pd.Series(sub["name"].values, index=ividx)
        partitions[(seq, label)] = s
    return partitions


def prepare_annotation_files(annot_files, file_labels, drop_lowconf):

    label_is_file = pl.Path(file_labels[0]).is_file()

    if len(annot_files) == 1 and len(file_labels) == 1 and label_is_file:
        # assume both listings
        with open(annot_files[0]) as listing:
            load_files = listing.read().strip().split()
        with open(file_labels[0]) as listing:
            assign_labels = listing.read().strip().split()
    else:
        load_files = annot_files
        assign_labels = file_labels
    assert len(load_files) == len(assign_labels)

    label_lut = dict()
    for annot_file, annot_label in zip(load_files, assign_labels):
        fp = pl.Path(annot_file).resolve(strict=True)
        ann = load_annotation(fp, annot_label, drop_lowconf)
        label_lut.update(ann)

    return label_lut


def check_obo_errors(regions):

    regions["size"] = regions["end"] - regions["start"]
    size_one = regions.index[regions["size"] == 1]

    drop_rows = []
    new_rows = []
    for obo_idx in size_one:
        if obo_idx < 1:
            continue
        if obo_idx == regions.shape[0] - 1:
            break
        before_idx = obo_idx - 1
        after_idx = obo_idx + 1
        if regions.at[before_idx, "chrom"] != regions.at[after_idx, "chrom"]:
            continue
        if regions.at[before_idx, "name"] == regions.at[after_idx, "name"]:
            chrom = regions.at[obo_idx, "chrom"]
            start = regions.at[before_idx, "start"]
            end = regions.at[after_idx, "end"]
            name = regions.at[before_idx, "name"]
            size = end - start
            new_region = (chrom, start, end, name, size)

            new_rows.append(new_region)

            drop_rows.append(before_idx)
            drop_rows.append(obo_idx)
            drop_rows.append(after_idx)

    if drop_rows:
        regions.drop(drop_rows, axis=0, inplace=True)
    if new_rows:
        new_rows = pd.DataFrame.from_records(
            new_rows, columns=["chrom", "start", "end", "name", "size"]
        )
        regions = pd.concat([regions, new_rows], axis=0, ignore_index=False)
        regions.sort_values(["chrom", "start"], inplace=True)

    regions.drop("size", axis=1, inplace=True)
    regions.reset_index(drop=True, inplace=True)
    return regions


def main():

    args = parse_command_line()

    label_lut = prepare_annotation_files(
        args.annotation_tables, args.annotation_labels,
        args.drop_lowconf
    )
    intersect = pd.read_csv(args.intersect, sep="\t", header=0)

    out_regions = []

    for seq, subset in intersect.groupby("chrom"):
        for row in subset.itertuples():
            iv = pd.Interval(row.start, row.end, closed=IVCLOSED)
            label_sources = row.list.split(",")
            collected_labels = set()
            for label_source in label_sources:
                annotation = label_lut[(seq, label_source)]
                select = annotation.index.overlaps(iv)
                if not select.any() and len(collected_labels) == 0:
                    # may happen if low-confidence labels
                    # were dropped when loading annotation tables
                    collected_labels.add("UNCERTAIN")
                elif select.any():
                    labels = set(annotation[select].values)
                    collected_labels = collected_labels.union(labels)
                else:
                    # no overlap (e.g., low-confidence label),
                    # but other labels were already collected
                    # beforehand
                    pass
            final_label = sorted(collected_labels)
            if not final_label:
                final_label = "UNCERTAIN"
            elif len(final_label) > 1:
                final_label = ",".join(filter(lambda l: l != "UNCERTAIN", final_label))
            else:
                final_label = final_label[0]
            region = Region(seq, row.start, row.end, final_label)
            out_regions.append(region)

    out_regions = sorted(out_regions)

    # sanity checking
    for a, b in itt.pairwise(out_regions):
        if a.chrom != b.chrom:
            continue
        assert a.end == b.start

    out_regions = pd.DataFrame.from_records(
        out_regions, columns=["chrom", "start", "end", "name"]
    )

    out_regions = check_obo_errors(out_regions)
    out_regions.rename({"chrom": "#chrom"}, axis=1, inplace=True)

    args.output.parent.mkdir(exist_ok=True, parents=True)
    out_regions.to_csv(args.output, sep="\t", header=True, index=False)

    return 0



if __name__ == "__main__":
    main()
