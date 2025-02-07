import pathlib


def get_num_threads_hmmer(motif_name):
    """This separate helper function exists
    to encapsulate the THREADS value set for
    a default run, i.e. CPU_LOW
    """
    num_threads = int(min(CPU_MAX, CPU_LOW * hmmer_scaling("cpu", motif_name)))
    return num_threads


def get_mem_mb_hmmer(motif_name):
    """Empirically, HMMER uses about 3g
    per CPU thread on average
    """

    num_threads = get_num_threads_hmmer(motif_name)
    mem_scale_factor = hmmer_scaling("mem", motif_name)

    mem_per_thread = 3072 * mem_scale_factor
    total_mem = int(mem_per_thread * num_threads)
    return total_mem


def hmmer_scaling(resource, motif_name):

    assert resource in ["cpu", "mem", "time"]
    try:
        scaling_factor = int(HMMER_MOTIF_PARAMS[motif_name][f"scale_{resource}"])
    except KeyError:
        scaling_factor = 1
    return scaling_factor


def hmmer_threshold_value(threshold, motif_name):

    # defaults as per HMMER cli interface
    _DEFAULT_HMMER_EVALUE = "10"
    _DEFAULT_HMMER_SCORE = 0

    assert threshold in ["evalue_t", "evalue", "score", "score_t"]

    t_value = None
    if threshold in ["evalue", "evalue_t"]:
        # NB: only the E-value threshold is used in the HMMER call
        t_value = HMMER_MOTIF_PARAMS[motif_name].get("evalue_t", _DEFAULT_HMMER_EVALUE)

    if threshold in ["score", "score_t"]:
        t_value = HMMER_MOTIF_PARAMS[motif_name].get("score_t", _DEFAULT_HMMER_EVALUE)

    assert t_value is not None

    return t_value
