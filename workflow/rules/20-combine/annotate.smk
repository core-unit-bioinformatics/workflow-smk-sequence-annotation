
_MULTI_ANNOTATION_INPUT = []
_MULTI_ANNOTATION_LABELS = []

if RUN_HMMER:
    # for HMMER, this only uses the high-quality hits
    _MULTI_ANNOTATION_INPUT.extend(
        sorted(
            expand(
                rules.compress_subset_hmmer_output.output.bed,
                motif=HMMER_MOTIF_NAMES,
                allow_missing=True
            )
        )
    )
    _MULTI_ANNOTATION_LABELS.extend(sorted(HMMER_MOTIF_NAMES))


if RUN_MINIMAP_REGIONDB:
    _MULTI_ANNOTATION_INPUT.extend(
        sorted(
            expand(
                rules.create_annotation_region_db.output.bed,
                region_db=MINIMAP_REGION_DB_NAMES,
                allow_missing=True
            )
        )
    )
    _MULTI_ANNOTATION_LABELS.extend(
        [
            label + ".regiondb" for label in
            sorted(MINIMAP_REGION_DB_NAMES)
        ]
    )


if RUN_MINIMAP_LABELREF:
    _MULTI_ANNOTATION_INPUT.extend(
        sorted(
            expand(
                rules.create_annotation_labeled_ref.output.bed,
                labelref=MINIMAP_LABELED_REFERENCE_NAMES,
                allow_missing=True
            )
        )
    )
    _MULTI_ANNOTATION_LABELS.extend(
        [
            label + ".labelref" for label in
            sorted(MINIMAP_LABELED_REFERENCE_NAMES)
        ]
    )


if _MULTI_ANNOTATION_LABELS:
    assert len(set(_MULTI_ANNOTATION_LABELS)) == len(_MULTI_ANNOTATION_LABELS)


rule bedtools_annotation_multi_intersect:
    input:
        bed_files = _MULTI_ANNOTATION_INPUT
    output:
        table = DIR_PROC.joinpath(
            "20-combine", "annotate", "multiinter",
            "{sample}.{path_id}.annot-isect.tsv.gz"
        )
    conda:
        DIR_ENVS.joinpath("biotools", "interval_tools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    params:
        header=" ".join(_MULTI_ANNOTATION_LABELS)
    shell:
        "bedtools multiinter -header -names {params.header} -i {input.bed_files} | gzip > {output.table}"


rule run_all_combine_annotations:
    input:
        tables = expand(
            rules.bedtools_annotation_multi_intersect.output.table,
            match_sample_path_id,
            sample=SAMPLES,
            path_id=PATH_IDS
        )
