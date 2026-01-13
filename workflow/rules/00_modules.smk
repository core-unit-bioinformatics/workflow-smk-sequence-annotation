"""
Use this module to list all includes
required for your pipeline - do not
add your pipeline-specific modules
to "commons/00_commons.smk"
"""

include: "00-settings/workflow.smk"

include: "00-settings/samples/pyutils.smk"
include: "00-settings/samples/sample_table.smk"

include: "00-settings/tools/repeatmasker.smk"
include: "00-settings/tools/hmmer.smk"
include: "00-settings/tools/minimap.smk"

include: "05-prepare/seqsplit/split_input.smk"
include: "05-prepare/deployment/repeatmasker.smk"
include: "05-prepare/seqnorm/check_input.smk"

include: "10-annotate/motifs/pyutils.smk"
include: "10-annotate/motifs/hmmer.smk"
include: "10-annotate/repeats/pyutils.smk"
include: "10-annotate/repeats/repeatmasker.smk"
include: "10-annotate/regions/minimap.smk"

include: "20-combine/pyutils.smk"
include: "20-combine/combine_all.smk"
include: "20-combine/combine_group.smk"
include: "20-combine/annotate.smk"
