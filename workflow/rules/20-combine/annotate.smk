
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


if RUN_MINIMAP_REGIONDB and RUN_MINIMAP_LABELREF:

    if any(label in MINIMAP_REGION_DB_NAMES for label in MINIMAP_LABELED_REFERENCE_NAMES):
        _disjoin_regiondb = [label + ".regiondb" for label in sorted(MINIMAP_REGION_DB_NAMES)]
        _disjoin_labelref = [label + ".labelref" for label in sorted(MINIMAP_LABELED_REFERENCE_NAMES)]
    else:
        _disjoin_regiondb = sorted(MINIMAP_REGION_DB_NAMES)
        _disjoin_labelref = sorted(MINIMAP_LABELED_REFERENCE_NAMES)
    _MULTI_ANNOTATION_INPUT.extend(
        sorted(
            expand(
                rules.create_annotation_region_db.output.bed,
                region_db=MINIMAP_REGION_DB_NAMES,
                allow_missing=True
            )
        )
    )
    _MULTI_ANNOTATION_LABELS.extend(_disjoin_regiondb)

    _MULTI_ANNOTATION_INPUT.extend(
        sorted(
            expand(
                rules.create_annotation_labeled_ref.output.bed,
                labelref=MINIMAP_LABELED_REFERENCE_NAMES,
                allow_missing=True
            )
        )
    )
    _MULTI_ANNOTATION_LABELS.extend(_disjoin_labelref)

elif RUN_MINIMAP_REGIONDB:
    _MULTI_ANNOTATION_INPUT.extend(
        sorted(
            expand(
                rules.create_annotation_region_db.output.bed,
                region_db=MINIMAP_REGION_DB_NAMES,
                allow_missing=True
            )
        )
    )
    _MULTI_ANNOTATION_LABELS.extend(sorted(MINIMAP_REGION_DB_NAMES))

elif RUN_MINIMAP_LABELREF:
    _MULTI_ANNOTATION_INPUT.extend(
        sorted(
            expand(
                rules.create_annotation_labeled_ref.output.bed,
                labelref=MINIMAP_LABELED_REFERENCE_NAMES,
                allow_missing=True
            )
        )
    )
    _MULTI_ANNOTATION_LABELS.extend(sorted(MINIMAP_LABELED_REFERENCE_NAMES))

else:
    pass


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


localrules: dump_annotation_label_listings
rule dump_annotation_label_listings:
    input:
        bed_files = _MULTI_ANNOTATION_INPUT
    output:
        lst_files = DIR_PROC.joinpath(
            "20-combine", "annotate", "multiinter",
            "{sample}.{path_id}.annot-files.lst"
        ),
        lst_labels = DIR_PROC.joinpath(
            "20-combine", "annotate", "multiinter",
            "{sample}.{path_id}.annot-labels.lst"
        ),
    params:
        labels = _MULTI_ANNOTATION_LABELS
    run:
        n_files = 0
        with open(output.lst_files) as listing:
            for filepath in input.bed_files:
                _ = listing.write(f"{filepath}\n")
                n_files += 1
        n_labels = 0
        with open(output.lst_labels) as listing:
            for label in params.labels:
                _ = listing.write(f"{label}\n")
                n_labels += 1
        assert n_files == n_labels
    # END OF RUN BLOCK


rule relabel_multi_annotation_table:
    input:
        lst_files = rules.dump_annotation_label_listings.output.lst_files,
        lst_labels = rules.dump_annotation_label_listings.output.lst_labels,
        isect = rules.bedtools_annotation_multi_intersect.output.table
    output:
        tmp_bed = temp(
            DIR_PROC.joinpath(
                "20-combine", "annotate", "multiinter", "tmp",
                "{sample}.{path_id}.relabeled.bed"
            )
        ),
        bed = DIR_RES.joinpath(
            "annotations", "combined",
            "{sample}.{path_id}.relabeled.bed.gz"
        )
    conda:
        DIR_ENVS.joinpath("scripts", "pyseq.yaml")
    params:
        script=find_script("add_annot_labels")
    shell:
        "{params.script} -i {input.isect} -a {input.lst_files} -l {input.lst_labels} --out-bed {output.tmp_bed}"
            " && "
        "bgzip -c {output.tmp_bed} > {output.bed}"
            " && "
        "tabix -p bed {output.bed}"


rule run_all_combine_annotations:
    input:
        tables = expand(
            rules.bedtools_annotation_multi_intersect.output.table,
            match_sample_path_id,
            sample=SAMPLES,
            path_id=PATH_IDS
        ),
        combined = expand(
            rules.relabel_multi_annotation_table.output.bed,
            match_sample_path_id,
            sample=SAMPLES,
            path_id=PATH_IDS
        )
