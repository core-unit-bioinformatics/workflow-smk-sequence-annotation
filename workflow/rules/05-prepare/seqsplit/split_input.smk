"""
Convenience sub-module to split
larger (whole-genome) input files
into individual sequences and to
create a suitable sample sheet file
to restart the workflow with.
"""


rule split_input_by_sequence:
    input:
        fasta = lambda wildcards: SAMPLE_INPUT[wildcards.sample]["all_files"]
    output:
        tsv = DIR_PROC.joinpath(
            "05-prepare", "seqsplit",
            "{sample}", "{sample}.splits.tsv"
        )
    conda:
        DIR_ENVS.joinpath("biotools", "motif_tools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt,
        time_hrs=lambda wildcards, attempt: attempt
    params:
        out_dir = lambda output: pl.Path(output.tsv).parent,
        script = find_script("splitfa.py")
    shell:
        "{params.script} --sample {wildcards.sample} "
        "--input {input.fasta} --output {params.out_dir}"


rule merge_all_split_sample_sheets:
    input:
        tables = expand(
            rules.split_input_by_sequence.output.tsv,
            sample=SAMPLES
        )
    output:
        tsv = DIR_PROC.joinpath(
            "05-prepare", "seqsplit", "samples.splitseq.tsv"
        )
    run:
        import pandas as pd

        merge = []
        for tsv_file in input.tables:
            df = pd.read_csv(tsv_file, sep="\t", header=0)
        merge = pd.concat(merge, axis=0, ignore_index=False)

        merge.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK
