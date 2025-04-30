#!/usr/bin/env python3

import argparse as argp
import pathlib as pl
import re

import pandas as pd
import pyranges as pr

___prog__ = "dump_region_db_bed.py"
__version__ = "0.1"
__description__ = """
This script takes a normalized PAF alignment file as single input.
That PAF file is expected to hold alignment information for
a region database aligned against a sequence of interest (the target).
In this setting, the region database is assumed to represent a number
of labeled sequences as individual entries in a FASTA file (the query).
This script computes basic statistics on the alignment quality:
(i) pct. matching positions in the alignment relative to the query length
and (ii) and the largest block of identity matches if a CIGAR string
is present in the PAF file. Note that, if the CIGAR string was not
computed with the =/X operation symbols (for identity or mismatch),
but only M (match), then this value is taken as proxy. The score column in
the output BED file is computed as the percentile rank scaled to 0-1000
of the average of (i) and (ii) or just (i) if no CIGAR string is present.
The usecase for this script is to label the target sequence with the
sequence labels from the region database (the query).

Alignments/regions overlapping in the query will be merged in a name/label-
and strand-specific way. The statistics (i) and (ii) are then merged
for the resulting BED entry by taking the maximal value of all merged
alignment entries.

This script does otherwise no filtering to avoid discarding "low-quality"
mappings in highly repetitive regions.
"""


def parse_command_line():

    parser = argp.ArgumentParser(prog=___prog__, description=__description__)

    parser.add_argument(
        "--input-paf", "--norm-paf", "-paf", "-in",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="paf_alignments",
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


def get_max_identity_block_size(cigar_string):

    max_block = 0
    for block in re.finditer("[0-9]+\=", cigar_string):
        max_block = max(max_block, int(block.group(0)[:-1]))

    if max_block == 0:
        # issue warning? CIGAR was not produced with =/X operations
        # so take 'matching' as proxy
        for block in re.finditer("[0-9]+M", cigar_string):
            max_block = max(max_block, int(block.group(0)[:-1]))
    return max_block


def add_cluster_ids(paf):

    # NB below: the 'Name' column is built
    # from the query names, i.e., from the names
    # of the FASTA sequences in the region
    # database.
    target_iv = pr.from_dict(
        {
            "Chromosome": paf.target_name,
            "Start": paf.target_start,
            "End": paf.target_end,
            "Strand": paf.strand,
            "Name": paf.query_name,
            "pd_idx": paf.index.values
        }
    )
    # NB last entry: we carry the index of the original
    # Pandas dataframe with us to enable merging back the
    # cluster IDs into the PAF dataframe.

    # merge by identical strand and name
    target_iv = target_iv.cluster(strand="same", by="Name")

    # Series just holds Cluster IDs as computed
    # by PyRanges
    cluster_ids = pd.Series(
        target_iv.Cluster.values,
        index=target_iv.pd_idx.values,
        name="cluster_id"
    )

    # this: just a sanity check given the dev stage of PyRanges
    _n_rows = paf.shape[0]

    paf = paf.merge(cluster_ids, left_index=True, right_index=True)

    # ensure nothing was lost
    assert paf.shape[0] == _n_rows

    return paf


def build_bed_regions(paf):

    bed_regions = []
    for cluster_id, alignments in paf.groupby("cluster_id"):
        # must be given because PyRanges was clustering by name
        assert alignments.query_name.nunique() == 1
        seq = alignments.target_name.iloc[0]
        name = alignments.query_name.iloc[0]
        strand = alignments.strand.iloc[0]
        start = alignments.target_start.min()
        end = alignments.target_end.max()
        assert start < end
        # this here: cf. description; avoiding the overhead
        # of calculating that exactly, so in the resulting BED
        # region, the interpretation of these values is that
        # ... over all regions merged
        # ... there is at least one query/target alignment [block]
        # with that score / pct. identity / id. block length
        score = alignments.mean_rank.max()
        max_pct = alignments.pct_matching.max()
        max_block = alignments.max_id_block.max()
        # this overloaded label is in preparation of merging these
        # BED entries with regions/annotations from other sources,
        # i.e. to retain some information for later filtering
        new_label = f"{name}[IDPCT:{max_pct}|IDBLK:{max_block}]"
        bed_regions.append(
            (seq, start, end, new_label, score, strand, name, max_pct, max_block, alignments.shape[0])
        )

    bed_df = pd.DataFrame.from_records(
        bed_regions,
        columns=[
            "chrom", "start", "end", "name", "score", "strand",
            "label", "max_pct_align_match", "longest_id_block",
            "merged_alignments"
        ],
    )

    bed_df.sort_values(["chrom", "start", "end"], inplace=True)
    # ensure proper header row when dumped to file
    bed_df.rename({"chrom": "#chrom"}, axis=1, inplace=True)

    return bed_df


def main():

    args = parse_command_line()

    paf = pd.read_csv(args.paf_alignments, sep="\t")

    # these two columns must exist
    paf["pct_matching"] = round(paf.align_matching / paf.query_length * 100, 2)
    paf["rank_matching"] = paf.pct_matching.rank(method="average", ascending=True, pct=True)

    # the cigar string column may exist
    if "cg_cigar" in paf.columns:
        paf["max_id_block"] = paf.cg_cigar.apply(get_max_identity_block_size)
        paf["rank_block"] = paf.max_id_block.rank(method="average", ascending=True, pct=True)
        paf["mean_rank"] = ((paf.rank_matching + paf.rank_block) / 2 * 1000).round(0).astype(int)
    else:
        paf["mean_rank"] = (paf.rank_matching * 1000).round(0).astype(int)

    # this replacement is done only for PyRanges, which only accepts
    # string values for the strand information
    paf["strand"] = paf.align_orient.replace({1: "+", -1: "-", 0: "."}).astype(str)

    paf = add_cluster_ids(paf)
    bed_regions = build_bed_regions(paf)

    args.output.parent.mkdir(exist_ok=True, parents=True)
    bed_regions.to_csv(args.output, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
