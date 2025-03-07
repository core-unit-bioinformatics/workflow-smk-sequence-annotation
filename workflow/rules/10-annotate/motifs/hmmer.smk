
rule hmmer_motif_search:
    """NB: the reported hits can EITHER
    be thresholded on the E-value [-E] OR
    on the score [-T], but not both in the same run.
    The implementation here only thresholds on
    the E-value (if specified for the motif) and
    then later labels hits above the score threshold
    (if specified for the motif) as high-quality
    """
    input:
        fasta = rules.check_input_sequences.output.norm_file,
        motif = DIR_GLOBAL_REF.joinpath("{motif}.fasta")
    output:
        txt = DIR_PROC.joinpath(
            "10-annotate", "motifs", "hmmer",
            "{sample}.hmmer.wd",
            "{sample}.{path_id}.{motif}.hmmer-out.txt"
        ),
        table = DIR_PROC.joinpath(
            "10-annotate", "motifs", "hmmer",
            "{sample}.hmmer.wd",
            "{sample}.{path_id}.{motif}.hmmer-table.txt"
        ),
    benchmark:
        DIR_RSRC.joinpath(
            "10-annotate", "motifs", "hmmer",
            "{sample}.{path_id}.{motif}.hmmer.rsrc"
        )
    wildcard_constraints:
        sample=CONSTRAINT_ALL_SAMPLES
    conda:
        DIR_ENVS.joinpath("biotools", "motif_tools.yaml")
    threads: lambda wildcards: get_num_threads_hmmer(wildcards.motif)
    resources:
        mem_mb = lambda wildcards, attempt: attempt * get_mem_mb_hmmer(wildcards.motif),
        time_hrs = lambda wildcards, attempt: attempt * hmmer_scaling("time", wildcards.motif)
    params:
        evalue_t = lambda wildcards: hmmer_threshold_value("evalue_t", wildcards.motif),
        alphabet = "--dna"
    shell:
        "nhmmer --cpu {threads} {params.alphabet} "
        "-o {output.txt} --tblout {output.table} "
        "-E {params.evalue_t} "
        "{input.motif} {input.fasta}"


rule normalize_hmmer_output_table:
    input:
        txt_table = rules.hmmer_motif_search.output.table
    output:
        tsv = DIR_RES.joinpath(
            "annotations", "motifs", "hmmer",
            "{sample}",
            "{sample}.{path_id}.{motif}.hmmer-tblout-norm.tsv.gz"
        ),
        tmp_bed = temp(
            DIR_RES.joinpath(
                "annotations", "motifs", "hmmer", "tmp",
                "{sample}.{path_id}.{motif}.hmmer-tblout-norm.bed"
        )),
        bed = DIR_RES.joinpath(
            "annotations", "motifs", "hmmer",
            "{sample}",
            "{sample}.{path_id}.{motif}.hmmer-tblout-norm.bed.gz"
        )
    conda:
        DIR_ENVS.joinpath("scripts", "pyseq.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    params:
        script=find_script("norm_hmmer_table"),
        score_t = lambda wildcards: (
            f"-score-t {hmmer_threshold_value('score', wildcards.motif)}"
            if hmmer_threshold_value('score', wildcards.motif) > 0 else ""
        )
    shell:
        "{params.script} --hmmer-table {input.txt_table} --add-metadata "
        "{params.score_t} --output-table {output.tsv} --output-bedlike {output.tmp_bed}"
            " && "
        "bgzip -c {output.tmp_bed} > {output.bed}"
            " && "
        "tabix -p bed {output.bed}"


rule compress_raw_hmmer_output:
    input:
        table = rules.hmmer_motif_search.output.table,
        text = rules.hmmer_motif_search.output.txt
    output:
        table = DIR_RES.joinpath(
            "annotations", "motifs", "hmmer", "{sample}", "raw",
            "{sample}.{path_id}.{motif}.hmmer-tblout.txt.gz"
        ),
        text = DIR_RES.joinpath(
            "annotations", "motifs", "hmmer", "{sample}", "raw",
            "{sample}.{path_id}.{motif}.hmmer-out.txt.gz"
        )
    shell:
        "gzip -c {input.table} > {output.table}"
            " && "
        "gzip -c {input.text} > {output.text}"


rule subset_hmmer_high_quality_hits:
    input:
        bedlike = rules.normalize_hmmer_output_table.output.bed
    output:
        bedlike = temp(
            DIR_RES.joinpath(
                "annotations", "motifs", "hmmer", "tmp",
                "{sample}.{path_id}.{motif}.hmmer-tblout-norm.hiq.bed"
            )
        )
    run:
        import pandas as pd

        df = pd.read_csv(input.bedlike, sep="\t", header=0)
        if "high_quality_hit" in df.columns:
            df = df.loc[df["high_quality_hit"] > 0, :].copy()
        df.to_csv(output.bedlike, sep="\t", header=0, index=False)
    # END OF RUN BLOCK


rule compress_subset_hmmer_output:
    input:
        bedlike = rules.subset_hmmer_high_quality_hits.output.bedlike
    output:
        bed = DIR_RES.joinpath(
            "annotations", "motifs", "hmmer",
            "{sample}",
            "{sample}.{path_id}.{motif}.hmmer-tblout-norm.hiq.bed.gz"
        )
    conda:
        DIR_ENVS.joinpath("biotools", "align_tools.yaml")
    shell:
        "bgzip -c {input.bedlike} > {output.bedlike}"
            " && "
        "tabix -p bed {output.bed}"


rule run_all_hmmer_motif_searches:
    input:
        tables = expand(
            rules.normalize_hmmer_output_table.output.tsv,
            match_sample_path_id,
            sample=SAMPLES,
            path_id=PATH_IDS,
            motif=HMMER_MOTIF_NAMES
        ),
        bedlike = expand(
            rules.normalize_hmmer_output_table.output.bed,
            match_sample_path_id,
            sample=SAMPLES,
            path_id=PATH_IDS,
            motif=HMMER_MOTIF_NAMES
        ),
        raw_out = expand(
            rules.compress_raw_hmmer_output.output.text,
            match_sample_path_id,
            sample=SAMPLES,
            path_id=PATH_IDS,
            motif=HMMER_MOTIF_NAMES
        ),
        hiq = expand(
            rules.compress_subset_hmmer_output.output.bed,
            match_sample_path_id,
            sample=SAMPLES,
            path_id=PATH_IDS,
            motif=HMMER_MOTIF_NAMES
        )
