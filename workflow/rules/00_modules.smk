"""
Use this module to list all includes
required for your pipeline - do not
add your pipeline-specific modules
to "commons/00_commons.smk"
"""

include: "00-settings/tools/repeatmasker.smk"
include: "05-prepare/deployment/repeatmasker.smk"
