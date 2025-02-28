
rule repeatmasker_default_run:
    """The output file listing of this rule
    omit the ".cat/.cat.gz" file because
    this seems only to be gzip-compressed
    if it reaches a certain size. Pulling
    that out of the 'proc/' hierarchy is done
    in a subsequent (dynamic) python rule.
    """
    input:
        setup_ok = rules.setup_repeatmasker.input.ok,
        fasta = rules.check_input_sequences.output.norm_file
    output:
        masked_fasta = DIR_PROC.joinpath(
            "10-annotate", "repeats", "repeatmasker",
            "{sample}.repeatmasker.wd",
            "{sample}.{path_id}.fasta.masked"
        ),
        summary = DIR_PROC.joinpath(
            "10-annotate", "repeats", "repeatmasker",
            "{sample}.repeatmasker.wd",
            "{sample}.{path_id}.fasta.tbl"
        ),
        table = DIR_PROC.joinpath(
            "10-annotate", "repeats", "repeatmasker",
            "{sample}.repeatmasker.wd",
            "{sample}.{path_id}.fasta.out"
        ),
    log:
        DIR_LOG.joinpath("10-annotate", "repeats", "repeatmasker",
            "{sample}.{path_id}.repeatmasker.log"
        )
    benchmark:
        DIR_RSRC.joinpath("10-annotate", "repeats", "repeatmasker",
            "{sample}.{path_id}.repeatmasker.rsrc"
        )
    wildcard_constraints:
        sample=CONSTRAINT_ALL_SAMPLES
    conda:
        DIR_ENVS.joinpath("biotools", "motif_tools.yaml")
    threads: lambda wildcards, attempt: CPU_LOW if attempt < 2 else CPU_MEDIUM
    resources:
        mem_mb = lambda wildcards, attempt, input: attempt * get_repeatmasker_run_memory_mb(input.size_mb),
        time_hrs = lambda wildcards, attempt, input: attempt * get_repeatmasker_run_time_hrs(input.size_mb)
    params:
        species=REPEATMASKER_SPECIES,
        outdir=lambda wildcards, output: pathlib.Path(output.summary).parent
    shell:
        "RepeatMasker -pa {threads} -dir {params.outdir} "
        "-xsmall -gff -s -species {params.species} "
        "{input.fasta} &> {log}"


rule run_all_repeatmasker_default:
    """Why the shell call?
    Failed RepeatMasker runs do not clean up after themselves
    (presumably to keep debugging information intact) and
    there is no switch to change that behavior. Hence,
    this trigger rule also performs the cleanup operation.
    Obviously, this can only be kicked off after all jobs
    have completed (restarted until completion),
    which makes RepeatMasker a typical candidate
    for leaving behind garbage in case the pipeline is
    interrupted in some way. Extremely annoying!!!
    """
    input:
        checks = expand(
            rules.repeatmasker_default_run.output.summary,
            match_sample_path_id,
            sample=SAMPLES,
            path_id=PATH_IDS
        )
    shell:
        "rm -rf RM_*"


