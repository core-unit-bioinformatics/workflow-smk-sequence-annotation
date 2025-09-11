"""This module implements two different alignment-based
labeling strategies. The first one takes a 'region database'
(= a FASTA file w/ more than one entry) and aligns those
regions to the sample sequence.

The second strategy does the inverse approach and aligns
the sample sequences to a labeled reference sequence, i.e.
this one requires a reference FASTA file and a corresponding
BED file.
"""

### first set of rules: region database approach

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
            "{sample}.{path_id}.{region_db}.region-db-aln.paf.gz"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "10-annotate", "region_db", "minimap",
            "{sample}.{path_id}.{region_db}.region-db-aln.mm2.rsrc"
        )
    wildcard_constraints:
        sample=CONSTRAINT_ALL_SAMPLES
    conda:
        DIR_ENVS.joinpath("biotools", "align_tools.yaml")
    threads: CPU_LOW
    resources:
        mem_mb=lambda wildcards, attempt: 24576 * attempt,
        time_hrs=lambda wildcards, attempt: attempt * attempt
    shell:
        "minimap2 -x asm20 -t {threads} -N 100 -p 0.9 -L -c --eqx --MD "
        "{input.fasta} {input.region_db} | gzip > {output.paf}"


rule normalize_paf_align_region_db:
    input:
        paf = rules.minimap_align_region_db.output.paf,
    output:
        tsv = DIR_RES.joinpath(
            "annotations", "regions", "minimap",
            "{sample}",
            "{sample}.{path_id}.{region_db}.mm2-region-db-aln.norm-paf.tsv.gz"
        )
    conda:
        DIR_ENVS.joinpath("biotools", "align_tools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    params:
        script=find_script("normalize_paf")
    shell:
        "{params.script} --input {input.paf} --output {output.tsv}"


rule dump_region_db_bed:
    """The script called in this rule makes use of PyRanges
    to simply cluster/merge overlapping alignment intervals
    with identical query name (= region from the region database)
    and strand. Computed statistics give a rough idea
    of the original alignment quality (pct. id. matches and so on)
    but may not be accurate if the CIGAR string is not part of
    the PAF file.

    This script is fairly "dumb" in the sense that it largely
    implements a format conversion w/ some descriptive statistics.
    """
    input:
        norm_paf = rules.normalize_paf_align_region_db.output.tsv,
    output:
        tmp_bed = temp(
            DIR_PROC.joinpath(
                "10-annotate", "region_db", "annotation",
                "{sample}.{path_id}.{region_db}.regiondb-labeled.bed"
        )),
        bed = DIR_RES.joinpath(
            "annotations", "regions", "minimap",
            "{sample}",
            "{sample}.{path_id}.{region_db}.mm2-region-db.bed.gz"
        )
    conda:
        DIR_ENVS.joinpath("scripts", "pyregions.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    params:
        script=find_script("dump_region_db_bed")
    shell:
        "{params.script} --input-paf {input.norm_paf} --output-bed {output.tmp_bed}"
            " && "
        "bgzip -c {output.tmp_bed} > {output.bed}"
            " && "
        "tabix -p bed {output.bed}"


rule run_all_minimap_region_db:
    input:
        tsv = expand(
            rules.dump_region_db_bed.output.bed,
            match_sample_path_id,
            sample=SAMPLES,
            path_id=PATH_IDS,
            region_db=MINIMAP_REGION_DB_NAMES
        )


### second set of rules: labeled reference approach


rule minimap_align_labeled_reference:
    input:
        fasta = rules.check_input_sequences.output.norm_file,
        label_ref = lambda wildcards: DIR_GLOBAL_REF.joinpath(
            MINIMAP_LABELED_REFERENCES[wildcards.labelref]["sequence"]
        )
    output:
        paf = DIR_PROC.joinpath(
            "10-annotate", "labeled_ref", "minimap",
            "{sample}.minimap.wd",
            "{sample}.{path_id}.{labelref}.label-ref-aln.paf.gz"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "10-annotate", "labeled_ref", "minimap",
            "{sample}.{path_id}.{labelref}.label-ref-aln.mm2.rsrc"
        )
    wildcard_constraints:
        sample=CONSTRAINT_ALL_SAMPLES
    conda:
        DIR_ENVS.joinpath("biotools", "align_tools.yaml")
    threads: CPU_LOW
    resources:
        mem_mb=lambda wildcards, attempt: 24576 * attempt,
        time_hrs=lambda wildcards, attempt: attempt * attempt
    shell:
        "minimap2 -x asm20 -t {threads} --secondary=no -L -c --eqx --MD "
        "{input.label_ref} {input.fasta} | gzip > {output.paf}"


rule normalize_paf_align_labeled_ref:
    input:
        paf = rules.minimap_align_labeled_reference.output.paf,
    output:
        tsv = DIR_PROC.joinpath(
            "10-annotate", "labeled_ref", "minimap",
            "{sample}.minimap.wd",
            "{sample}.{path_id}.{labelref}.label-ref-aln.norm-paf.tsv.gz"
        )
    conda:
        DIR_ENVS.joinpath("biotools", "align_tools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    params:
        script=find_script("normalize_paf")
    shell:
        "{params.script} --input {input.paf} --output {output.tsv}"


rule trim_paf_align_labeled_ref:
    """From rustybam's cli help:

    >>>
    This is a function for lifting over coordinates from a reference (<BED>) to a query using a PAF file
    [...]
    The returned file is a PAF file that is trimmed to the regions in the bed file.
    Even the cigar in the returned PAF file is trimmed so it can be used downstream
    <<<

    So we use rustybam here to lift the regions annotated
    in the reference (= the target in the PAF) over to the
    query. The output is, however, not a lifted BED but the
    trimmed PAF which still needs to be reduced to a
    BED-like format.
    """
    input:
        paf = rules.minimap_align_labeled_reference.output.paf,
        labels = lambda wildcards: DIR_GLOBAL_REF.joinpath(
            MINIMAP_LABELED_REFERENCES[wildcards.labelref]["labels"]
        )
    output:
        trimmed_paf = DIR_PROC.joinpath(
            "10-annotate", "labeled_ref", "minimap",
            "{sample}.minimap.wd",
            "{sample}.{path_id}.{labelref}.label-ref-aln.trimmed.paf.gz"
        )
    conda:
        DIR_ENVS.joinpath("biotools", "align_tools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt
    shell:
        "rustybam liftover --bed {input.labels} {input.paf} | gzip > {output.trimmed_paf}"


rule normalize_trimmed_paf_align_labeled_ref:
    input:
        trimmed_paf = rules.trim_paf_align_labeled_ref.output.trimmed_paf,
    output:
        tsv = DIR_PROC.joinpath(
            "10-annotate", "labeled_ref", "minimap",
            "{sample}.minimap.wd",
            "{sample}.{path_id}.{labelref}.label-ref-aln.trimmed.norm-paf.tsv.gz"
        )
    conda:
        DIR_ENVS.joinpath("biotools", "align_tools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    params:
        script=find_script("normalize_paf")
    shell:
        "{params.script} --input {input.trimmed_paf} --output {output.tsv}"


rule dump_labeled_ref_bed:
    input:
        norm_paf = rules.normalize_trimmed_paf_align_labeled_ref.output.tsv,
        labels = lambda wildcards: DIR_GLOBAL_REF.joinpath(
            MINIMAP_LABELED_REFERENCES[wildcards.labelref]["labels"]
        )
    output:
        tmp_bed = temp(
            DIR_PROC.joinpath(
                "10-annotate", "labeled_ref", "annotation",
                "{sample}.{path_id}.{labelref}.labeled-ref.bed"
        )),
        bed = DIR_RES.joinpath(
            "annotations", "regions", "minimap",
            "{sample}",
            "{sample}.{path_id}.{labelref}.mm2-label-ref.bed.gz"
        )
    conda:
        DIR_ENVS.joinpath("scripts", "pyregions.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt
    params:
        script=find_script("dump_labeled_ref_bed")
    shell:
        "{params.script} --input-paf {input.norm_paf} --input-bed {input.labels} "
        "--output-bed {output.tmp_bed}"
            " && "
        "bgzip -c {output.tmp_bed} > {output.bed}"
            " && "
        "tabix -p bed {output.bed}"


rule run_all_minimap_labeled_ref:
    input:
        tsv = expand(
            rules.dump_labeled_ref_bed.output.bed,
            match_sample_path_id,
            sample=SAMPLES,
            path_id=PATH_IDS,
            labelref=MINIMAP_LABELED_REFERENCE_NAMES
        )
