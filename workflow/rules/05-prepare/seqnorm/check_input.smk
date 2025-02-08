
rule check_input_sequences:
    """This preprocessing rule also exists
    because RepeatMasker cannot work on
    compressed input files.
    """
    input:
        input_file = lambda wildcards: SAMPLE_INPUT[wildcards.sample]["by_id"][wildcards.path_id]
    output:
        norm_file = DIR_PROC.joinpath(
            "05-prepare", "seqnorm",
            "{sample}", "{sample}.{path_id}.fasta"
        )
    conda:
        DIR_ENVS.joinpath("scripts", "pyseq.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt,
        time_hrs=lambda wildcards, attempt: attempt * attempt
    params:
        check_names="--no-name-check" if SKIP_SEQUENCE_HEADER_NAME_CHECK else "",
        script=find_script("norm_seq_input")
    shell:
        "{params.script} {params.check_names} "
        "--input-file {input.input_file} "
        "--output-file {output.norm_file}"

