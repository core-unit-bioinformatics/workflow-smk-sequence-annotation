
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
            "{sample}.{path_id}.{motif}.hmmer-tblout-norm.tsv.gz"
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
        "{params.script} --hmmer-table {input.table} --add-metadata "
        "{params.score_t} --output-table {output.tsv}"


rule compress_raw_hmmer_output:
    input:
        table = rules.hmmer_motif_search.output.table,
        text = rules.hmmer_motif_search.output.txt
    output:
        table = DIR_RES.joinpath(
            "annotations", "motifs", "hmmer", "raw",
            "{sample}.{path_id}.{motif}.hmmer-tblout.txt.gz"
        ),
        text = DIR_RES.joinpath(
            "annotations", "motifs", "hmmer", "raw",
            "{sample}.{path_id}.{motif}.hmmer-out.txt.gz"
        )
    shell:
        "gzip {input.table} > {output.table}"
            " && "
        "gzip {input.text} > {output.text}"


rule run_all_hmmer_motif_searches:
    input:
        tables = expand(
            rules.normalize_hmmer_output_table.output.tsv,
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
        )
