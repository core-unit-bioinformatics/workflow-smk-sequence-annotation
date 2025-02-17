
rule minimap_align_region_db:
    input:
        fasta = rules.check_input_sequences.output.norm_file,
        region_db = lambda wildcards: DIR_GLOBAL_REF.joinpath(
            MINIMAP_REGION_DB[wildcards.region_db]
        )
    output:
        paf = DIR_PROC.joinpath(
            "10-annotate", "region_db", "minimap",
            "{sample}.minimap.wd",
            "{sample}.{path_id}.{region_db}.aln.paf.gz"
        )
    benchmark: DIR_RSRC.joinpath(
            "10-annotate", "region_db", "minimap",
            "{sample}.{path_id}.{region_db}.aln.mm2.rsrc"
        )
    wildcard_constraints:
        sample=CONSTRAINT_ALL_SAMPLES
    conda:
        DIR_ENVS.joinpath("biotools", "align_tools.yaml")
    threads: CPU_LOW
    resources:
        mem_mb=lambda wildcards, attempt: 16384 * attempt,
        time_hrs=lambda wildcards, attempt: attempt * attempt
    shell:
        "minimap2 -x asm20 -t {threads} -N 5 -p 0.95 -L -c --eqx --MD "
        "{input.fasta} {input.region_db} | gzip > {output.paf}"


rule normalize_paf_align_region_db:
    input:
        paf = rule.minimap_align_region_db.output.paf,
    output:
        tsv = DIR_PROC.joinpath(
            "10-annotate", "region_db", "minimap",
            "{sample}.minimap.wd",
            "{sample}.{path_id}.{region_db}.aln.norm-paf.tsv.gz"
        )
    conda:
        DIR_ENVS.joinpath("biotools", "align_tools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    params:
        script=find_script("normalize_paf")
    shell:
        "{params.script} --input {input.paf} --output {output.tsv}"


rule run_all_minimap_region_db:
    input:
        tsv = expand(
            rules.normalize_paf_align_region_db.output.tsv,
            match_sample_path_id,
            sample=SAMPLES,
            path_id=PATH_IDS,
            region_db=MINIMAP_REGION_DB_NAMES
        )
