#!/usr/bin/env python3

import argparse as argp
import dataclasses as dcl
import pathlib as pl
import re

import pandas as pd
pd.set_option('future.no_silent_downcasting', True)
import xopen


@dcl.dataclass
class RepeatMaskerTableRow:
    """This follow the RepeatMasker docs from:
    https://www.repeatmasker.org/webrepeatmaskerhelp.html

    NB: the above docs are outdated (?) for the columns
    'repeat_match_start', 'repeat_match_end' and 'repeat_remaining_bp'
    because they are reporting bp prior to the beginning
    of the repeat motif match, and not bp remaining after the end.

    This script has been tested with output from RepeatMasker v4.1.7
    """
    sw_score: int
    pct_sub: float
    pct_del: float
    pct_ins: float
    query_name: str
    query_match_start: int
    query_match_end: int
    query_remaining_bp: int
    match_orientation: int
    repeat_motif_name: str
    repeat_class_family: str
    repeat_match_start: int
    repeat_match_end: int
    repeat_remaining_bp: int
    repeat_match_id: int
    higher_scoring_overlap: int


class RepeatMaskerTable:
    __slots__ = "header", "data", "datatypes", "column_positions"

    def __init__(self):
        self.header = None
        self.datatypes = None
        self.data = []
        return None

    def add_row(self, line_num, table_row):

        columns = table_row.strip().split()
        try:
            _ = int(columns[0])
        except ValueError:
            if line_num > 2:
                # probably wrong
                raise ValueError(f"Invalid line {line_num}: {table_row}")
            self.add_table_header(columns)
        else:
            self.add_data_row(columns)

    def add_table_header(self, columns):

        if self.header is not None:
            return

        check_keywords = re.compile("(SW|score|perc|matching|position|begin|end)")
        # sanity check: does this look like a header line?
        if not any([check_keywords.search(column) is not None for column in columns]):
            raise ValueError(f"Expected header line, but this looks different: {columns}")

        column_names = [field.name for field in dcl.fields(RepeatMaskerTableRow)]
        column_types = [field.type for field in dcl.fields(RepeatMaskerTableRow)]

        self.header = column_names
        self.datatypes = column_types
        self.column_positions = dict(
            (col_name, pos) for pos, col_name in enumerate(column_names, start=0)
        )

        return

    def add_data_row(self, fields):
        """This function adds a 'data' row with the following transformations:
        1) if non-existent, add field 16 with default '.' (asterisk or empty in default output)
        2) strip parenthesis from all coordinate fields (where applicable)
        3) standardize the orientation/complement match column to 1/forward/+ or -1/reverse/- [C]
        4) standardize the 'higher scoring overlap' column to 0/no or 1/yes

        NB: this function does _not_ chnage the coordinates.
        By experimentation, the coordinate intervals reported by RepeatMasker
        are 1-based, half-open. In order to turn these coordinates into
        0-based, half-open (= BED-compatible intervals, Python slicing etc.),
        the >start< coordinate must be decremented by 1 (start -= 1).
        """

        assert self.header is not None

        num_fields = len(fields)
        if not 15 <= num_fields <= 16:
            raise ValueError(f"Expected 15 or 16 fields/cells in table row: {fields}")
        if num_fields == 15:
            fields.append(".")

        # manually normalize a couple of fields
        _col_idx = self.column_positions
        # strip potential parenthesis
        par_fields = [
            "query_remaining_bp", "repeat_match_start",
            "repeat_match_end", "repeat_remaining_bp"
        ]
        for par_col in par_fields:
            col_pos = _col_idx[par_col]
            fields[col_pos] = fields[col_pos].strip("()")

        # norm orientation column
        fields[_col_idx["match_orientation"]] = fields[
            _col_idx["match_orientation"]
        ].replace("+", "1").replace("C", "-1")

        # norm higher scoring overlap column
        fields[_col_idx["higher_scoring_overlap"]] = fields[
            _col_idx["higher_scoring_overlap"]
        ].replace(".", "0").replace("*", "1")

        fields = [
            column_type(value) for column_type, value in zip(self.datatypes, fields)
        ]

        self.data.append(RepeatMaskerTableRow(*tuple(fields)))

        return

    def to_dataframe(self):

        df = pd.DataFrame.from_records(
            [dcl.astuple(row) for row in self.data],
            columns=self.header
        )
        assert (df["query_match_start"] < df["query_match_end"]).all()
        return df

    def to_bedlike(self, dataframe):
        """The coordinate transformation below assumes the following
        base enumeration/coordinate reporting by RepeatMasker:

        0.|.1.2.3.4.5.|.6.7.8 === Convert 0-based: (2-1):6
        1.|.2.3.4.5.6.|.7.8.9 === RepMask forward: 2:6
        A.|.C.C.G.T.T.|.G.C.A === Forward strand / match CCGTT
        ------------------------------------------------------
        T.G.C.|.A.A.C.G.G.|.T === Reverse strand / match AACGG
        9.8.7.|.6.5.4.3.2.|.1 === RepMask reverse: 6:2
        8.7.6.|.5.4.3.2.1.|.0 === Convert 0-based: (2-1):6

        """
        columns_to_keep = [
            "query_name", "query_match_start", "query_match_end",
            "repeat_motif_name", "sw_score", "match_orientation",
            "repeat_class_family", "higher_scoring_overlap",
            "pct_sub", "pct_del", "pct_ins"
        ]

        bedlike = dataframe[columns_to_keep].copy()
        bedlike["match_orientation"] = bedlike["match_orientation"].replace(
            {1: "+", -1: "-", 0: "."},
            inplace=False
        )
        bedlike.sort_values(["query_name", "query_match_start", "query_match_end"], inplace=True)
        return bedlike


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--repeatmasker-table", "-rt",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="repeatmasker_table",
        required=True
    )

    parser.add_argument(
        "--output-table", "-o",
        type=lambda fp: pl.Path(fp).resolve(strict=False),
        dest="output_table",
        required=True
    )

    parser.add_argument(
        "--output-bedlike", "-bed",
        type=lambda fp: pl.Path(fp).resolve(strict=False),
        default=None,
        dest="output_bedlike"
    )

    args = parser.parse_args()

    return args


def main():

    args = parse_command_line()

    with xopen.xopen(args.repeatmasker_table) as table:
        rmtable = RepeatMaskerTable()
        for ln, line in enumerate(table, start=1):
            if not line.strip():
                continue
            rmtable.add_row(ln, line.strip())
        df = rmtable.to_dataframe()

    args.output_table.parent.mkdir(exist_ok=True, parents=True)
    df.to_csv(
        args.output_table,
        sep="\t",
        header=True,
        index=False
    )

    if args.output_bedlike is not None:
        bedlike = rmtable.to_bedlike(df)
        args.output_bedlike.parent.mkdir(exist_ok=True, parents=True)
        with xopen.xopen(args.output_bedlike, "w") as bedfile:
            _ = bedfile.write("#")
            _ = bedlike.to_csv(bedfile, sep="\t", header=True, index=False)

    return 0


if __name__ == "__main__":
    main()
