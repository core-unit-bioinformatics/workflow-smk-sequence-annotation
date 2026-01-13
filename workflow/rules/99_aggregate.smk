"""
Use this module to extend the default
workflow output (a list of target files)
per sub-module.
The WORKFLOW_OUTPUT list is referenced
in the main Snakefile
"""

WORKFLOW_OUTPUT = []
# Example for extending the output
# with output from another module
# (remember to include that module
# in 00_modules.smk):
# WORKFLOW_OUTPUT.extend(MODULE_OUTPUT)

if SPLIT_INPUT_BY_SEQUENCE:
    # is this option is set, only that one rule
    # must be executed and the workflow restarted
    # afterwards with the new sample sheet.
    WORKFLOW_OUTPUT.append(rules.merge_all_split_sample_sheets.output.tsv)
