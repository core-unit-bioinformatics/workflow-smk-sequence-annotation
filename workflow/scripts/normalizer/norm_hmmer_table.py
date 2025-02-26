#!/usr/bin/env python3

import argparse as argp
import dataclasses as dcl
import pathlib as pl
import re
import sys

import pandas as pd
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
        elif re.match("^#\\s+\\-+", table_row) is not None:
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
            if field.type is str and field.name not in ["target_name", "strand", "evalue"]
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

    def to_dataframe(self):

        df = pd.DataFrame.from_records(
            [dcl.astuple(table_row) for table_row in self.data],
            columns=self.header
        )
        return df


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
        "--output-table", "-ot", "-o",
        type=lambda fp: pl.Path(fp).resolve(strict=False),
        dest="output_table"
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
        df = hmmer_table.to_dataframe()
        if df is not None:
            df.to_csv(table, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
