#!/usr/bin/env python3

import argparse as argp
import dataclasses as dcl
import pathlib as pl
import re
import sys

import pandas as pd
pd.set_option('future.no_silent_downcasting', True)
import xopen


@dcl.dataclass
class HmmerTableRow:
    """Field names defined according to
    http://eddylab.org/software/hmmer/Userguide.pdf
    (v3.4 / pages 69f)
    """
    target_name: str
    target_accession: str
    query_name: str
    query_accession: str
    hmm_hit_start: int
    hmm_hit_end: int
    target_hit_start: int
    target_hit_end: int
    target_env_start: int
    target_env_end: int
    target_length: int
    strand: str
    evalue: str  # avoid ugly floating point conversions
    bit_score: float
    bias: float
    description_of_target: str


class HmmerTable:
    __slots__ = (
        "header", "metadata", "datatypes", "na_fields",
        "data", "first_data_line",
        "run_complete", "verbose", "ignore_incomplete"
    )

    def __init__(self, quiet=True, ignore_incomplete=False):
        self.run_complete = False
        self.metadata = []
        self.data = []
        self.first_data_line = None
        self.verbose = not quiet
        self.ignore_incomplete = ignore_incomplete
        return None

    def add_row(self, line_num, table_row):

        assert isinstance(table_row, str)
        if table_row.startswith("# target name"):
            self.add_table_header(table_row)
        elif re.match("^#(\\s+)?\\-+", table_row) is not None:
            # separating row, do nothing
            return
        elif table_row == "# [ok]":
            self.run_complete = True
        elif table_row == "#":
            # blank/separating line
            return
        elif table_row.startswith("# "):
            self.metadata.append(table_row)
        else:
            if self.first_data_line is None:
                self.first_data_line = line_num
            self.add_data_row(line_num, table_row)
        return

    def add_table_header(self, header_row):

        singletons = "(accession|hmmfrom|alifrom|envfrom|strand|E\\-value|score|bias)"
        composed = "([a-z]+\\s{1}[a-z]+(\\s{1}target)?)"
        header_expr = re.compile(f"({singletons}|{composed})")

        _expected_columns = 16
        identified_columns = 0

        column_names = [field.name for field in dcl.fields(HmmerTableRow)]

        column_types = [field.type for field in dcl.fields(HmmerTableRow)]

        na_fields = [
            idx for idx, field in enumerate(dcl.fields(HmmerTableRow))
            if field.type is str and field.name not in ["target_name", "query_name", "strand", "evalue"]
        ]

        recognized_header_fields = set()
        for column_pos, mobj in enumerate(header_expr.finditer(header_row), start=0):
            raw_column_name = mobj.group(0)
            try:
                clean_column_name = column_names[column_pos]
                assert clean_column_name not in recognized_header_fields
            except IndexError:
                raise ValueError(f"Do not recognize table header field: {raw_column_name}")
            except AssertionError:
                raise ValueError(f"Duplicate matching of header field: {raw_column_name}")
            if self.verbose:
                sys.stderr.write(
                    f"\nMatching input column {column_pos} with header '{raw_column_name}' "
                    f"to new name '{clean_column_name}'\n"
                )
            identified_columns += 1
            recognized_header_fields.add(clean_column_name)

        if identified_columns != _expected_columns:
            raise ValueError(
                f"Number of fields mismatch: {_expected_columns} vs {identified_columns} "
                f" --- {sorted(recognized_header_fields)}"
            )

        self.header = column_names
        self.datatypes = column_types
        self.na_fields = na_fields
        return

    def add_data_row(self, line_num, data_row):
        """This function adds a 'data' row with the following transformations:
        1) replace '-' by 'n/a' in fields where '-' does not indicate a complement match

        NB: this function does _not_ change the coordinates.
        By experimentation, the coordinate intervals reported by HMMER
        are 1-based, half-open. In order to turn these coordinates into
        0-based, half-open (= BED-compatible intervals, Python slicing etc.),
        the >start< coordinate must be decremented by 1 (start -= 1)
        """
        fields = data_row.split()
        fields = [
            field.replace("-", "n/a") if column_pos in self.na_fields else field
            for column_pos, field in enumerate(fields)
        ]
        fields = [
            column_type(value) for column_type, value in zip(self.datatypes, fields)
        ]
        if len(fields) != len(self.header):
            raise ValueError(f"Data row {line_num} does not match header: {self.header} vs {fields} fields")

        self.data.append(HmmerTableRow(*tuple(fields)))
        return

    def check_complete(self):
        if self.run_complete:
            if self.verbose:
                sys.stderr.write(f"\nHMMER table complete with {len(self.data)} records\n")
        elif self.ignore_incomplete:
            pass
        else:
            raise RuntimeError(
                f"HMMER table looks incomplete/corrupt (# metadata records: {len(self.metadata)})"
            )
        return

    def get_metadata(self):

        if self.run_complete:
            md_block = "\n".join(self.metadata) + "\n"
        else:
            md_block = "# empty-file-no-metadata"
        return md_block

    def to_dataframe(self, keep_hmmer_format=False, score_threshold=None):

        df = pd.DataFrame.from_records(
            [dcl.astuple(table_row) for table_row in self.data],
            columns=self.header
        )
        if not keep_hmmer_format:
            df["strand"] = df["strand"].replace(
                {"+": 1, "-": -1, ".": 0},
                inplace=False
            ).astype(int)
            # switch coordinates for negative strand
            # from HMMER manual:
            # "strand: The strand on which the hit was found (“-" when alifrom>ali to)."
            negative_strand = df["strand"] < 0
            df.loc[negative_strand, ["target_hit_start", "target_hit_end"]] = df.loc[
                negative_strand, ["target_hit_end", "target_hit_start"]
            ].values
            df.loc[negative_strand, ["target_env_start", "target_env_end"]] = df.loc[
                negative_strand, ["target_env_end", "target_env_start"]
            ].values
            assert (df["target_hit_start"] < df["target_hit_end"]).all()
            assert (df["target_env_start"] < df["target_env_end"]).all()
        assert not pd.isnull(df).any(axis=0).any()

        if score_threshold is not None:
            df["high_quality_hit"] = 0
            select_hiq = df["bit_score"] > score_threshold
            df.loc[select_hiq, "high_quality_hit"] = 1

        return df

    def to_bedlike(self, score_threshold=None):
        """The coordinate transformation below assumes the following
        base enumeration/coordinate reporting by HMMER:

        0.|.1.2.3.4.5.|.6.7.8 === Convert 0-based: (2-1):6
        1.|.2.3.4.5.6.|.7.8.9 === HMMER forward: 2:6
        A.|.C.C.G.T.T.|.G.C.A === Forward strand / match CCGTT
        ------------------------------------------------------
        T.G.C.|.A.A.C.G.G.|.T === Reverse strand / match AACGG
        9.8.7.|.6.5.4.3.2.|.1 === HMMER reverse: 6:2
        8.7.6.|.5.4.3.2.1.|.0 === Convert 0-based: (2-1):6

        """

        df = self.to_dataframe(keep_hmmer_format=False, score_threshold=score_threshold)

        columns_to_keep = [
            "target_name", "target_hit_start", "target_hit_end",
            "query_name", "bit_score", "strand", "evalue", "high_quality_hit"
        ]

        columns_to_keep = [col for col in columns_to_keep if col in df.columns]
        bedlike = df[columns_to_keep].copy()

        if not df.empty:

            assert (df["target_hit_start"] < df["target_hit_end"]).all()

            bedlike["target_hit_start"] -= 1
            bedlike["strand"] = bedlike["strand"].replace(
                {-1: "-", 1: "+", 0: "."},
                inplace=False
            )
            bedlike.sort_values(
                ["target_name", "target_hit_start", "target_hit_end"],
                inplace=True
            )
            bedlike["bit_score"] = bedlike["bit_score"].round(0).astype(int)

        return bedlike


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--hmmer-table", "-ht",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="hmmer_table",
        required=True,
        help="Path to HMMER output 'table'."
    )

    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        default=False,
        dest="quiet",
        help="Be quiet about ongoing operations."
    )

    parser.add_argument(
        "--add-metadata", "-md",
        action="store_true",
        default=False,
        dest="add_metadata"
    )

    parser.add_argument(
        "--ignore-incomplete",
        action="store_true",
        default=False,
        dest="ignore_incomplete"
    )

    parser.add_argument(
        "--keep-hmmer-format",
        action="store_true",
        default=False,
        dest="keep_hmmer_format"
    )

    parser.add_argument(
        "--add-score-threshold", "-score-t",
        type=int,
        dest="score_threshold",
        default=None,
        help=(
            "If set, add a new binary column 'high_quality_hit' and "
            "label hit as 1 if its bit score is above the threshold "
            "and 0 otherwise. Default: None"
        )
    )

    parser.add_argument(
        "--output-table", "-ot", "-o",
        type=lambda fp: pl.Path(fp).resolve(strict=False),
        dest="output_table",
        required=True
    )

    parser.add_argument(
        "--output-bedlike", "-bed",
        type=lambda fp: pl.Path(fp).resolve(strict=False),
        dest="output_bedlike"
    )

    args = parser.parse_args()

    return args


def main():

    args = parse_command_line()

    with xopen.xopen(args.hmmer_table, "r") as table:
        hmmer_table = HmmerTable(args.quiet, args.ignore_incomplete)
        for ln, line in enumerate(table, start=1):
            hmmer_table.add_row(ln, line.strip())
        hmmer_table.check_complete()

    args.output_table.parent.mkdir(exist_ok=True, parents=True)
    with xopen.xopen(args.output_table, "w") as table:
        if args.add_metadata:
            table.write(hmmer_table.get_metadata())
        df = hmmer_table.to_dataframe(args.keep_hmmer_format, args.score_threshold)
        df.to_csv(table, sep="\t", header=True, index=False)

    if args.output_bedlike is not None:
        bedlike = hmmer_table.to_bedlike(args.score_threshold)
        args.output_bedlike.parent.mkdir(exist_ok=True, parents=True)
        with xopen.xopen(args.output_bedlike, "w") as bedfile:
            _ = bedfile.write("#")
            bedlike.to_csv(bedfile, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
