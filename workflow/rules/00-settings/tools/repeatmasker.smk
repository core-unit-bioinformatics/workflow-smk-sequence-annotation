import pathlib

RUN_REPEATMASKER = config.get("run_repeatmasker", False)

REPEATMASKER_OFFLINE_SETUP = config.get("repeatmasker_offline_setup", False)

REPEATMASKER_DEFAULT_LIBRARY_FILE = config.get("repeatmasker_default_library_file", None)
if REPEATMASKER_DEFAULT_LIBRARY_FILE is None:
    if REPEATMASKER_OFFLINE_SETUP:
        err_msg = "Offline setup for RepeatMasker requires a default library file"
        logerr(err_msg)
        raise RuntimeError(err_msg)
else:
    REPEATMASKER_DEFAULT_LIBRARY_FILE = pathlib.Path(
        REPEATMASKER_DEFAULT_LIBRARY_FILE
    ).resolve(strict=True)

# The following: as of RepeatMasker 4.1.7+
# the Dfam/FamDB has been split into separate partitions.
# Making use of the Dfam/FamDB partitions is not required,
# but if so, the root partition (fist file, index 0)
# always has to be used.
# See here / installations instructions
# https://repeatmasker.org/RepeatMasker/
# "The root ("dfam##_full.0.h5") partition is required if you plan to use Dfam,
# however any combination of additional partitions may also be downloaded and configured."
REPEATMASKER_DFAM_ROOT_PART_FILE = config.get("repeatmasker_dfam_root_part_file", None)
if REPEATMASKER_DFAM_ROOT_PART_FILE is not None:
    REPEATMASKER_DFAM_ROOT_PART_FILE = pathlib.Path(
        REPEATMASKER_DFAM_ROOT_PART_FILE
    ).resolve(strict=True)

REPEATMASKER_DFAM_ADD_PART_FILES = config.get("repeatmasker_dfam_add_part_files", None)
assert REPEATMASKER_DFAM_ADD_PART_FILES is None, "Additional DFAM partitions not supported by this workflow"

REPEATMASKER_SPECIES = config.get("repeatmasker_species", "human")
