
RUN_HMMER = config.get("run_hmmer", False)

HMMER_MOTIF_PARAMS = config.get("hmmer_motif_params", None)
assert isinstance(HMMER_MOTIF_PARAMS, dict)

if HMMER_MOTIF_PARAMS is not None:

    HMMER_MOTIF_NAMES = sorted(HMMER_MOTIF_PARAMS.keys())
