import collections

"""
All functions in this module are candidates
for inclusion in the template / commons modules.
"""

def check_data_identifier(identifier):
    """This utility function may be useful in other
    modules of the workflow --- do not start function
    name with underscore.
    """
    if re.match("[A-Za-z0-9]+$", identifier) is None:
        raise ValueError(f"Invalid data identifier: {identifier}")
    return


def match_sample_path_id(*wildcards):

    assert wildcards[0][0][0] == "sample"
    assert wildcards[1][0][0] == "path_id"

    sample_wildcards = set(t[1] for t in wildcards[0])
    path_wildcards = set(t[1] for t in wildcards[1])

    other_wildcards = collections.defaultdict(list)
    if len(wildcards) > 2:
        for wildcard_list in wildcards[2:]:
            [
                other_wildcards[name].append(value)
                for name, value in wildcard_list
            ]

    wildcard_combinations = []
    for sample in sample_wildcards:
        all_paths_for_sample = SAMPLE_INPUT[sample]["all_path_ids"]
        for path in all_paths_for_sample:
            if path not in path_wildcards:
                continue
            this_combination = {
                "sample": sample,
                "path_id": path
            }
            if len(other_wildcards) > 0:
                for other_wildcard, other_values in other_wildcards.items():
                    for value in other_values:
                        tmp = dict(this_combination)
                        tmp[other_wildcard] = value
                        wildcard_combinations.append(tmp)
            else:
                wildcard_combinations.append(this_combination)

    return wildcard_combinations


def _find_sample_column(sample_sheet_header):
    """TODO this must be solved in a generic manner
    in the commons modules.
    """
    sample_column = None
    for column in sample_sheet_header:
        if column.lower() == "sample":
            # preferred name
            sample_column = column
        elif column.lower() in ["sample_name", "sample_id"]:
            if sample_column is not None:
                if VERBOSE:
                    err_msg = (
                        "Warning: ambiguous column header in sample sheet\n"
                        f"Selected >sample< column header: {sample_column}\n"
                        f"Current column header: {column}\n"
                    )
                    logerr(err_msg)
                continue
            else:
                sample_column = column
        else:
            pass
    if sample_column is None:
        raise RuntimeError(
            f"Cannot identify sample column in sample sheet header: {sample_sheet_header}"
        )

    return sample_column


def _find_file_input_column(sample_sheet_header):
    """TODO this must be solved in a generic manner
    in the commons modules iff the sample sheet
    has a generic format.

    TODO abstraction: merge w/ function above
    """
    # TODO enum type
    accepted_names = [
        "input", "path", "file",
        "file_path", "filepath",
        "folder", "directory"
    ]
    accepted_names += [
        ("input_" + n).strip("_") for n in accepted_names
        if "input" not in n
    ]
    accepted_names += [n + "s" for n in accepted_names]

    input_column = None
    for column in sample_sheet_header:
        if column.lower() == "input_path":
            # preferred name
            input_column = column
        elif column.lower() in accepted_names:
            if input_column is not None:
                if VERBOSE:
                    err_msg = (
                        "Warning: ambiguous column header in sample sheet\n"
                        f"Selected >input_path< column header: {input_column}\n"
                        f"Current column header: {column}\n"
                    )
                    logerr(err_msg)
                    continue
            input_column = column
        else:
            pass
    if input_column is None:
        raise RuntimeError(
            f"Cannot identify sample column in sample sheet header: {sample_sheet_header}"
        )

    return input_column



def _normalize_sample_sex(sample_sex):
    """NB: even though it is _often_ the case that
    digits indicate sex/karyotype with the smaller
    one identifying males, there is no general rule for that.
    Hence, this function can only normalize character
    abbreviations.
    """

    known_norm = {
        "f": "female",
        "m": "male",
        "u": "any",
        "d": "any",
        "x": "any"
    }

    norm_sex = known_norm.get(sample_sex.lower(), sample_sex.lower())

    if not norm_sex in ["male", "female", "any"]:
        raise ValueError(f"Cannot normalize sample sex: {sample_sex}")

    return norm_sex


def _read_input_files_from_fofn(fofn_path):
    """Read input file listing from
    file of file names.
    """

    input_files = []
    with open(fofn_path, "r") as listing:
        for line in listing:
            if not line.strip():
                continue

            file_path = pathlib.Path(line.strip())
            if file_path.suffix not in INPUT_FILE_EXTENSIONS:
                if VERBOSE:
                    logerr(f"Skipping file w/ unknown extension: {file_path}")
                continue
            try:
                file_path = file_path.resolve(strict=True)
            except FileNotFoundError:
                try:
                    file_path = DATA_ROOT.joinpath(line.strip()).resolve(strict=True)
                except FileNotFoundError:
                    err_msg = "\nERROR\n"
                    err_msg += f"Cannot find file: {line.strip}\n"
                    err_msg += f"Data root is set to: {DATA_ROOT}\n"
                    logerr(err_msg)
                    raise
            input_files.append(file_path)

    if not input_files:
        logerr(f"No file paths found in {fofn_path} w/ extensions: {INPUT_FILE_EXTENSIONS}")
        raise FileNotFoundError(f"No input files loaded from FOFN: {fofn_path}")

    return sorted(input_files)


def _subset_path(full_path):
    """This helper exists to reduce
    the absolute path to a file
    to just the file name and its
    parent. Computing the hash over this
    shortened path means that even an active
    workflow run could be copied to another
    infrastructure where the root path / data root
    is different as long as relative paths
    are kept consistent.
    Premature optimization?
    """
    folder_name = full_path.parent.name
    file_name = full_path.name
    subset_path = f"{folder_name}/{file_name}"
    # if it so happens that the file resides
    # in a root-level location, strip off
    # leading slash
    return subset_path.strip("/")


def _collect_input_files(path_spec):
    """
    Generic function to collect input files
    from one or more paths on the file system.
    """
    input_files = []
    input_hashes = []
    input_path_ids = []
    for sub_input in path_spec.split(","):
        input_path = pathlib.Path(sub_input).resolve()
        if input_path.is_file() and input_path.name.endswith(".fofn"):
            fofn_files = _read_input_files_from_fofn(input_path)
            fofn_hashes = [
                hashlib.sha256(
                    _subset_path(fp).encode("utf-8")
                ).hexdigest() for fp in fofn_files
            ]
            input_files.extend(fofn_files)
            input_hashes.extend(fofn_hashes)
            input_path_ids.extend(
                [fh[:PATH_ID_LENGTH] for fh in fofn_hashes]
            )
        elif input_path.is_file():
            input_hash = hashlib.sha256(
                _subset_path(input_path).encode("utf-8")
            ).hexdigest()
            input_files.append(input_path)
            input_hashes.append(input_hash)
            input_path_ids.append(input_hash[:PATH_ID_LENGTH])
        elif input_path.is_dir():
            collected_files = _glob_collect_files(input_path)
            collected_hashes = [
                hashlib.sha256(
                    _subset_path(fp).encode("utf-8")
                ).hexdigest() for fp in collected_files
            ]
            input_files.extend(collected_files)
            input_hashes.extend(collected_hashes)
            input_path_ids.extend(
                [ch[:PATH_ID_LENGTH] for ch in collected_hashes]
            )
        else:
            raise ValueError(f"Cannot handle input: {sub_input}")

    num_files = len(input_files)
    num_hashes = len(set(input_hashes))
    if num_hashes != num_files:
        # That would indeed imply an actual hash collision (unlikely),
        # or the user used the same path more than once (likely)
        # for this row (sample) in the sample sheet. We consider
        # this an user error and abort.
        err_msg = (
            f"Path specification: {path_spec}\n"
            f"Number of collected files: {num_files}\n"
            f"Number of path hashes: {num_hashes}\n"
            "ERROR: must be 1-to-1\n"
            "Have you specified the same path twice?\n"
        )
        logerr(err_msg)
        raise RuntimeError("Hash collision or duplicate paths")



    return input_files, input_hashes, input_path_ids


def _build_constraint(values):
    escaped_values = sorted(map(re.escape, map(str, values)))
    constraint = "(" + "|".join(escaped_values) + ")"
    return constraint


def _glob_collect_files(folder):

    all_files = set()
    for pattern in INPUT_FILE_EXTENSIONS:
        pattern_files = set(folder.glob(f"**/*.{pattern}"))
        all_files = all_files.union(pattern_files)
    all_files = [f for f in sorted(all_files) if f.is_file()]
    if len(all_files) < 1:
        logerr(f"No files in path {folder} detected w/ extensions: {INPUT_FILE_EXTENSIONS}")
        raise FileNotFoundError(f"No input files found underneath {folder}")
    return all_files
