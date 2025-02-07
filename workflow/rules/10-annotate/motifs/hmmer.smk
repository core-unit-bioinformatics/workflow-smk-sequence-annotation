
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


rule run_all_hmmer_motif_searches:
    input:
        tables = expand(
            rules.hmmer_motif_search.output.table,
            match_sample_path_id,
            sample=SAMPLES,
            path_id=PATH_IDS,
            motif=HMMER_MOTIF_NAMES
        )
