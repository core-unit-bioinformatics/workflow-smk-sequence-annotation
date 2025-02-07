
# the following settings are unlikely to be
# modified by the user and have thus hard-coded
# defaults

# common file extensions / file suffix
# for sequence files
_DEFAULT_FILE_EXTENSIONS = [
    ".fasta", ".fasta.gz", ".fa", ".fa.gz",
    ".fna", ".fna.gz",
    ".fastq", ".fastq.gz", ".fq", ".fq.gz"
]

# input files w/ an unknown file extension
# will be skipped
INPUT_FILE_EXTENSIONS = config.get("input_file_extensions", _DEFAULT_FILE_EXTENSIONS)
assert isinstance(INPUT_FILE_EXTENSIONS, list)
# the following: Python's pathlib.Path().suffix returns
# the file suffix w/ a leading dot - correct potential
# user input errors
INPUT_FILE_EXTENSIONS = ["." + fext.strip(".") for fext in INPUT_FILE_EXTENSIONS]


# the path id (= hash over file path)
# is shortened to this to increase
# readability. If hash collisions occur,
# the workflow run will be aborted (must be...)
_DEFAULT_PATH_ID_LENGTH = 8
PATH_ID_LENGTH = config.get("path_id_length", _DEFAULT_PATH_ID_LENGTH)
assert isinstance(PATH_ID_LENGTH, int)


# the old-school tools used in this workflow
# may stumble over long or otherwise 'malformed'
# sequence headers (certain characters etc.)
# the following option enables the user to skip
# the lengthy name checking. Input files will then
# just be unzipped if necessary.
SKIP_SEQUENCE_HEADER_NAME_CHECK = config.get("skip_sequence_header_name_check", True)
assert isinstance(SKIP_SEQUENCE_HEADER_NAME_CHECK, bool)


