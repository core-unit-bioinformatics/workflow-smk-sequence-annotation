
rule repeatmasker_default_run:
    input:
        setup_ok = rules.setup_repeatmasker.input.ok,
        fasta = rules.check_input_sequences.output.norm_file
    output:
        check = DIR_PROC.joinpath(
            "10-annotate", "repeats", "repeatmasker",
            "{sample}.repeatmasker.wd",
            "{sample}.{path_id}.ok"
        )
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
    threads: lambda wildcards, attempt: CPU_LOW if attempt < 2 else CPU_MEDIUM
    resources:
        mem_mb = lambda wildcards, attempt, input: attempt * get_repeatmasker_run_memory_mb(input.size_mb),
        time_hrs = lambda wildcards, attempt, input: attempt * get_repeatmasker_run_time_hrs(input.size_mb)
    params:
        species=REPEATMASKER_SPECIES,
        outdir=lambda wildcards, output: pathlib.Path(output.check).parent
    shell:
        "RepeatMasker -pa {threads} -s -dir {params.outdir} "
        "-species {params.species} {input.fasta} &> {log}"
            " && "
        "touch {output.check}"


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
            rules.repeatmasker_default_run.output.check,
            match_sample_path_id,
            sample=SAMPLES,
            path_id=PATH_IDS
        )
    shell:
        "rm -rf RM_*"


