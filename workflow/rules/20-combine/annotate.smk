
### combination path 1:
# multi-intersect all annotation files
# the resulting table is currently 0/1 encoding
# for overlaps; still needs-postprocessing to add labels


rule bedtools_annotation_multi_intersect:
    input:
        bed_files = select_combination_input
    output:
        table = DIR_PROC.joinpath(
            "20-combine", "annotate", "multiinter",
            "{sample}.{path_id}.cmb-{group_prefix}.annot-isect.tsv.gz"
        )
    conda:
        DIR_ENVS.joinpath("biotools", "interval_tools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    params:
        header=lambda wildcards: " ".join(select_combination_labels(wildcards))
    shell:
        "bedtools multiinter -header -names {params.header} -i {input.bed_files} | gzip > {output.table}"


localrules: dump_annotation_label_listings
rule dump_annotation_label_listings:
    input:
        bed_files = select_combination_input
    output:
        lst_files = DIR_PROC.joinpath(
            "20-combine", "annotate", "multiinter",
            "{sample}.{path_id}.cmb-{group_prefix}.annot-files.lst"
        ),
        lst_labels = DIR_PROC.joinpath(
            "20-combine", "annotate", "multiinter",
            "{sample}.{path_id}.cmb-{group_prefix}.annot-labels.lst"
        ),
    params:
        labels = lambda wildcards: " ".join(select_combination_labels(wildcards))
    run:
        n_files = 0
        with open(output.lst_files, "w") as listing:
            for filepath in input.bed_files:
                _ = listing.write(f"{filepath}\n")
                n_files += 1
        n_labels = 0
        with open(output.lst_labels, "w") as listing:
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
                "{sample}.{path_id}.cmb-{group_prefix}.relabeled.bed"
            )
        ),
        bed = DIR_RES.joinpath(
            "annotations", "combined",
            "{sample}.{path_id}.cmb-{group_prefix}.relabeled.bed.gz"
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


### combination path 2: concat all labels
# this is only useful for manually checking label
# precision/correctness at individual loci (e.g. in IGV)


rule concat_multi_annotation_labels:
    input:
        bed_files = select_combination_input
    output:
        tmp_bed = temp(
            DIR_PROC.joinpath(
                "20-combine", "annotate", "multiinter", "tmp",
                "{sample}.{path_id}.cmb-{group_prefix}.concat-annot.bed"
            )
        ),
        bed = DIR_RES.joinpath(
            "annotations", "combined",
            "{sample}.{path_id}.cmb-{group_prefix}.concat.bed.gz"
        )
    conda:
        DIR_ENVS.joinpath("biotools", "interval_tools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt
    shell:
        "zcat {input.bed_files} | cut -f 1-4 | grep -v '^#' |sort -V -k1,1 -k2,3n > {output.tmp_bed}"
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
            path_id=PATH_IDS,
            group_prefix=GROUP_PREFIX_WILDCARDS
        ),
        combined = expand(
            rules.relabel_multi_annotation_table.output.bed,
            match_sample_path_id,
            sample=SAMPLES,
            path_id=PATH_IDS,
            group_prefix=GROUP_PREFIX_WILDCARDS
        ),
        concat = expand(
            rules.concat_multi_annotation_labels.output.bed,
            match_sample_path_id,
            sample=SAMPLES,
            path_id=PATH_IDS,
            group_prefix=GROUP_PREFIX_WILDCARDS
        )
