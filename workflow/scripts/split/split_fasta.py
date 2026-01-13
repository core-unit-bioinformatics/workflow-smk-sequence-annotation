#!/usr/bin/env python3

import argparse as argp
import pathlib as pl
import re

import dnaio


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        "--input", "-i",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="input",
        required=True
    )

    parser.add_argument(
        "--sample", "-s",
        type=str,
        dest="sample",
        required=True
    )

    parser.add_argument(
        "--output", "-o",
        type=lambda fp: pl.Path(fp).resolve(strict=False),
        dest="output",
        required=True
    )

    args = parser.parse_args()

    return args


def main():

    args = parse_command_line()

    args.output.mkdir(exist_ok=True, parents=True)

    sample_name = args.sample.strip("_")
    sample_prefix = sample_name + "_"

    accept_seq_name = re.compile("^[a-z0-9]{1,49}$", flags=re.IGNORECASE)

    sample_table = []
    with dnaio.open(args.input) as fasta:
        for record in fasta:
            seqname = record.name
            if accept_seq_name.match(seqname) is None:
                raise ValueError(
                    f"Malformed sequence name for sample {sample_name}: "
                    f"{seqname} (checked: {accept_seq_name.pattern})"
                )
            split_output_name = sample_prefix + seqname
            # cast to str for sample table output
            seqlength = str(len(record.sequence))

            out_file = args.output.joinpath(
                f"{split_output_name}.fasta"
            ).resolve()
            with dnaio.FastaWriter(out_file) as dump:
                dump.write(record.name, record.sequence)
            sample_table.append((split_output_name, seqlength, str(out_file), str(args.input)))

    sample_table_name = f"{sample_name}.splits.tsv"
    sample_table_file = args.output.joinpath(sample_table_name)
    with open(sample_table_file, "w") as table:
        header = ["sample", "seq_length", "input_path", "source_file"]
        table.write("\t".join(header) + "\n")
        for row in sample_table:
            table.write("\t".join(row) + "\n")

    return 0


if __name__ == "__main__":
    main()
