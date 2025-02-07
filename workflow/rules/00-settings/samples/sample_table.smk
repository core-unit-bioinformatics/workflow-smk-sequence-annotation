import pathlib

import pandas


SAMPLES = None
SAMPLE_INPUT = None
PATH_IDS = None

CONSTRAINT_ALL_SAMPLES = None


def process_sample_sheet():

    SAMPLE_SHEET_FILE = pathlib.Path(config["samples"]).resolve(strict=True)

    SAMPLE_SHEET = pandas.read_csv(
        SAMPLE_SHEET_FILE,
        sep="\t",
        header=0,
        comment="#"
    )
    sample_sheet_header = SAMPLE_SHEET.columns

    sample_column = _find_sample_column(sample_sheet_header)
    input_column = _find_file_input_column(sample_sheet_header)

    header_renamer = {
        sample_column: "sample",
        input_column: "input_path"
    }

    for other_column in sample_sheet_header:
        if other_column in header_renamer:
            continue
        header_renamer[other_column] = other_column.lower()

    SAMPLE_SHEET.rename(header_renamer, axis=1, inplace=True)
    assert SAMPLE_SHEET["sample"].nunique() == SAMPLE_SHEET.shape[0]

    sample_input = dict()
    path_ids = []
    for row in SAMPLE_SHEET.itertuples():
        sample = row.sample
        input_info = dict()

        input_path = row.input_path
        input_files, input_hashes, input_path_ids = _collect_input_files(
            input_path
        )
        input_info["all_files"] = input_files
        input_info["all_hashes"] = input_hashes
        input_info["all_path_ids"] = input_path_ids
        input_info["by_hash"] = dict(
            (fhash, fpath) for fhash, fpath in zip(input_hashes, input_files)
        )
        input_info["by_id"] = dict(
            (pathid, fpath) for pathid, fpath in zip(input_path_ids, input_files)
        )
        sample_input[sample] = input_info
        path_ids.extend(input_path_ids)

    assert len(sample_input) == SAMPLE_SHEET.shape[0]

    global SAMPLES
    SAMPLES = sorted(SAMPLE_SHEET["sample"].values)

    global SAMPLE_INPUT
    SAMPLE_INPUT = sample_input

    global PATH_IDS
    PATH_IDS = path_ids

    global CONSTRAINT_ALL_SAMPLES
    CONSTRAINT_ALL_SAMPLES = _build_constraint(SAMPLES)

    return


process_sample_sheet()
