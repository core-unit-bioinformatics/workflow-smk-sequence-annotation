#!/usr/bin/env python3

import argparse as argp
import pathlib as pl
import subprocess as sp

import filetype as ft
import xopen


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "-i", "--input", "--input-file",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="input_file",
        required=True
    )

    parser.add_argument(
        "-o", "--output", "--output-file",
        type=lambda fp: pl.Path(fp).resolve(strict=True),
        dest="output_file",
        required=True
    )

    parser.add_argument(
        "--no-name-check",
        action="store_true",
        default=False,
        dest="no_name_check"
    )

    args = parser.parse_args()
    return args


def check_sequence_name(seq_name):

    name_ok = True

    # this is a RepeatMasker constraint
    constraint_length = len(seq_name) < 50
    name_ok &= constraint_length

    # potentially more to add here...

    return name_ok


def determine_input_type(file_path):

    kind = ft.guess(file_path)
    if kind is None:
        # assume text/uncompressed
        is_compressed = False
    else:
        is_compressed = True

    # check for FASTA vs FASTQ
    with xopen.xopen(file_path) as seqfile:
        for line in seqfile:
            if not line.strip():
                continue
            header_start = line.strip()[0]
            break

    if header_start == ">":
        seq_format = "fasta"
    elif header_start == "@":
        seq_format = "fastq"
    else:
        raise RuntimeError(
            f"Unknown sequence file format {file_path}\n"
            f"header start char is: {header_start}"
        )

    return seq_format, is_compressed


def exec_sys_call(call):

    if isinstance(call, list):
        cmd = " ".join(list(map(str, call)))
    elif isinstance(call, str):
        cmd = call
    else:
        raise TypeError(f"Cannot handle input type for system call: {call}")

    _ = sp.check_call(cmd, shell=True)

    return


def convert_file(input_file, output_file):

    output_file.parent.mkdir(exist_ok=True, parents=True)
    cmd = ["seqtk", "-A", "-C", "-S", input_file, ">", output_file]
    exec_sys_call(cmd)
    return


def link_file(input_file, output_file):

    output_file.parent.mkdir(exist_ok=True, parents=True)
    cmd = ["ln", input_file, output_file]
    exec_sys_call(cmd)
    return


def create_fasta_index(file_path):

    assert file_path.is_file().resolve(strict=True)
    cmd = ["samtools", "faidx", file_path]
    exec_sys_call(cmd)
    return


def main():

    args = parse_command_line()

    seq_format, is_compressed = determine_input_type(args.input_file)

    if args.no_name_check:
        if is_compressed or seq_format == "fastq":
            convert_file(args.input_file, args.output_file)
        else:
            # plain fasta
            assert not is_compressed
            assert seq_format == "fasta"
            link_file(args.input_file, args.output_file)
    else:
        raise NotImplementedError(
            "Name checking of input sequence files not yet supported"
        )

    # by default, create FASTA index file
    create_fasta_index(args.output_file)

    return 0


if __name__ == "__main__":
    main()
