
RUN_HMMER = config.get("run_hmmer", False)

USE_HMMER_IN_LABEL_COMBINATION = config.get("use_hmmer_in_label_combination", False)

HMMER_MOTIF_PARAMS = config.get("hmmer_motif_params", None)

if HMMER_MOTIF_PARAMS is not None:
    assert isinstance(HMMER_MOTIF_PARAMS, dict)
    HMMER_MOTIF_NAMES = sorted(HMMER_MOTIF_PARAMS.keys())
else:
    HMMER_MOTIF_NAMES = []
