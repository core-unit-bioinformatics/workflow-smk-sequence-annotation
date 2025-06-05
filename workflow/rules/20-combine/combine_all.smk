"""
This is one of two dedicated modules to combine all
labels (module the respective tool not being run)
in a single BED files. The combination can either
be done for all label annotation files (this module),
or in a grouped manner by looking at the prefix
of the label identifiers (the other module) as
specified in the config file.
"""


_MULTI_ANNOTATION_INPUT_ALL = []
_MULTI_ANNOTATION_LABELS_ALL = []

if RUN_HMMER:
    # for HMMER, this only uses the high-quality hits
    _MULTI_ANNOTATION_INPUT_ALL.extend(
        sorted(
            expand(
                rules.compress_subset_hmmer_output.output.bed,
                motif=HMMER_MOTIF_NAMES,
                allow_missing=True
            )
        )
    )
    _MULTI_ANNOTATION_LABELS_ALL.extend(sorted(HMMER_MOTIF_NAMES))


if RUN_MINIMAP_REGIONDB and RUN_MINIMAP_LABELREF:

    if any(label in MINIMAP_REGION_DB_NAMES for label in MINIMAP_LABELED_REFERENCE_NAMES):
        _disjoin_regiondb = [label + ".regiondb" for label in sorted(MINIMAP_REGION_DB_NAMES)]
        _disjoin_labelref = [label + ".labelref" for label in sorted(MINIMAP_LABELED_REFERENCE_NAMES)]
    else:
        _disjoin_regiondb = sorted(MINIMAP_REGION_DB_NAMES)
        _disjoin_labelref = sorted(MINIMAP_LABELED_REFERENCE_NAMES)
    _MULTI_ANNOTATION_INPUT_ALL.extend(
        sorted(
            expand(
                rules.dump_region_db_bed.output.bed,
                region_db=MINIMAP_REGION_DB_NAMES,
                allow_missing=True
            )
        )
    )
    _MULTI_ANNOTATION_LABELS_ALL.extend(_disjoin_regiondb)

    _MULTI_ANNOTATION_INPUT_ALL.extend(
        sorted(
            expand(
                rules.dump_labeled_ref_bed.output.bed,
                labelref=MINIMAP_LABELED_REFERENCE_NAMES,
                allow_missing=True
            )
        )
    )
    _MULTI_ANNOTATION_LABELS_ALL.extend(_disjoin_labelref)

elif RUN_MINIMAP_REGIONDB:
    _MULTI_ANNOTATION_INPUT_ALL.extend(
        sorted(
            expand(
                rules.dump_region_db_bed.output.bed,
                region_db=MINIMAP_REGION_DB_NAMES,
                allow_missing=True
            )
        )
    )
    _MULTI_ANNOTATION_LABELS_ALL.extend(sorted(MINIMAP_REGION_DB_NAMES))

elif RUN_MINIMAP_LABELREF:
    _MULTI_ANNOTATION_INPUT_ALL.extend(
        sorted(
            expand(
                rules.dump_labeled_ref_bed.output.bed,
                labelref=MINIMAP_LABELED_REFERENCE_NAMES,
                allow_missing=True
            )
        )
    )
    _MULTI_ANNOTATION_LABELS_ALL.extend(sorted(MINIMAP_LABELED_REFERENCE_NAMES))

else:
    pass


if _MULTI_ANNOTATION_LABELS_ALL:
    assert len(set(_MULTI_ANNOTATION_LABELS_ALL)) == len(_MULTI_ANNOTATION_LABELS_ALL)
