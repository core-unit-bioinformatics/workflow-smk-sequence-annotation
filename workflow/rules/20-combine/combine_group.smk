import re
import collisions
import itertools as itt

"""
This is one of two dedicated modules to combine all
labels (module the respective tool not being run)
in a single BED files. The combination can either
be done for all label annotation files (the other module),
or in a grouped manner by looking at the prefix
of the label identifiers (this module) as
specified in the config file.
"""

# first determine all annotation label prefixes
_ANNOTATION_LABEL_PREFIXES = set()

for name in itt.chain(MINIMAP_REGION_DB_NAMES, MINIMAP_LABELED_REFERENCE_NAMES):
    prefix = re.match(r"^[0-9a-zA-Z]+", name)
    if prefix is not None:
        prefix = prefix.group(0)
    if prefix is None or len(prefix) < 3:
        raise ValueError(
            "Label combination grouped by prefix requires valid (3+ alphanumeric characters) "
            f"prefix at the beginning of the label name/identifier: {prefix}"
        )
    _ANNOTATION_LABEL_PREFIXES.add(prefix)

_ANNOTATION_LABEL_PREFIXES = sorted(_ANNOTATION_LABEL_PREFIXES)


_MULTI_ANNOTATION_INPUT_GRP = collections.defaultdict(list)
_MULTI_ANNOTATION_LABELS_GRP = collections.defaultdict(list)

if RUN_HMMER and USE_HMMER_IN_LABEL_COMBINATION:
    # for HMMER, this only uses the high-quality hits
    # NB: HMMER motifs/labels are not grouped by prefix
    for prefix in _ANNOTATION_LABEL_PREFIXES:
        _MULTI_ANNOTATION_INPUT_GRP[prefix].extend(
            sorted(
                expand(
                    rules.compress_subset_hmmer_output.output.bed,
                    motif=HMMER_MOTIF_NAMES,
                    allow_missing=True
                )
            )
        )
        _MULTI_ANNOTATION_LABELS_GRP[prefix].extend(sorted(HMMER_MOTIF_NAMES))


_RUN_AND_USE_MINIMAP_REGIONDB = RUN_MINIMAP_REGIONDB and USE_REGIONDB_IN_LABEL_COMBINATION
_RUN_AND_USE_MINIMAP_LABELREF = RUN_MINIMAP_LABELREF and USE_LABELREF_IN_LABEL_COMBINATION


if _RUN_AND_USE_MINIMAP_REGIONDB and _RUN_AND_USE_MINIMAP_LABELREF:

    for prefix in _ANNOTATION_LABEL_PREFIXES:

        if any(label in MINIMAP_REGION_DB_NAMES for label in MINIMAP_LABELED_REFERENCE_NAMES):
            _disjoin_regiondb = [
                label + ".regiondb" for label in sorted(MINIMAP_REGION_DB_NAMES)
                if label.startswith(prefix)
            ]
            _disjoin_labelref = [
                label + ".labelref" for label in sorted(MINIMAP_LABELED_REFERENCE_NAMES)
                if label.startswith(prefix)
            ]
        else:
            _disjoin_regiondb = sorted(
                label for label in MINIMAP_REGION_DB_NAMES if label.startswith(prefix)
            )
            _disjoin_labelref = sorted(
                label for label in MINIMAP_LABELED_REFERENCE_NAMES if label.startswith(prefix)
            )
        _MULTI_ANNOTATION_INPUT_GRP[prefix].extend(
            sorted(
                expand(
                    rules.dump_region_db_bed.output.bed,
                    region_db=MINIMAP_REGION_DB_NAMES,
                    allow_missing=True
                )
            )
        )
        _MULTI_ANNOTATION_LABELS_GRP[prefix].extend(_disjoin_regiondb)

        _MULTI_ANNOTATION_INPUT_GRP[prefix].extend(
            sorted(
                expand(
                    rules.dump_labeled_ref_bed.output.bed,
                    labelref=MINIMAP_LABELED_REFERENCE_NAMES,
                    allow_missing=True
                )
            )
        )
        _MULTI_ANNOTATION_LABELS_GRP[prefix].extend(_disjoin_labelref)

elif _RUN_AND_USE_MINIMAP_REGIONDB:

    for prefix in _ANNOTATION_LABEL_PREFIXES:
        _selected_region_dbs = [
            label for label in MINIMAP_REGION_DB_NAMES if label.startswith(prefix)
        ]
        _MULTI_ANNOTATION_INPUT_GRP[prefix].extend(
            sorted(
                expand(
                    rules.dump_region_db_bed.output.bed,
                    region_db=_selected_region_dbs,
                    allow_missing=True
                )
            )
        )
        _MULTI_ANNOTATION_LABELS_GRP[prefix].extend(sorted(_selected_region_dbs))

elif _RUN_AND_USE_MINIMAP_LABELREF:

    for prefix in _ANNOTATION_LABEL_PREFIXES:
        _selected_label_refs = [
            label for label in MINIMAP_LABELED_REFERENCE_NAMES if label.startswith(prefix)
        ]
        _MULTI_ANNOTATION_INPUT_GRP[prefix].extend(
            sorted(
                expand(
                    rules.dump_labeled_ref_bed.output.bed,
                    labelref=_selected_label_refs,
                    allow_missing=True
                )
            )
        )
        _MULTI_ANNOTATION_LABELS_GRP[prefix].extend(sorted(_selected_label_refs))

else:
    pass


if _MULTI_ANNOTATION_LABELS_GRP:
    for prefix in _ANNOTATION_LABEL_PREFIXES:
        assert len(set(_MULTI_ANNOTATION_LABELS_GRP[prefix])) == len(_MULTI_ANNOTATION_LABELS_GRP[prefix])
