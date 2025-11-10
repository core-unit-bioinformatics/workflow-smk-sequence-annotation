
def get_repeatmasker_run_memory_mb(input_size_mb, compressed=False):

    threshold_tiny = 1
    threshold_small = 100
    threshold_normal = 1000

    if compressed:
        # ~gzip compressed FASTA vs uncompressed
        # compressed fasta has smaller file size,
        # but total sequence is the relevant factor
        compression_scaling = 4
    else:
        compression_scaling = 1

    if input_size_mb < threshold_tiny * compression_scaling:
        mem_mb = 16384
    elif input_size_mb < threshold_small * compression_scaling:
        mem_mb = 24576
    elif input_size_mb < threshold_normal * compression_scaling:
        mem_mb = 32768
    else:
        mem_mb = 229376
    return mem_mb


def get_repeatmasker_run_time_hrs(input_size_mb, compressed=False):

    threshold_tiny = 1
    threshold_small = 100
    threshold_normal = 1000

    if compressed:
        # ~gzip compressed FASTA vs uncompressed
        # compressed fasta has smaller file size,
        # but total sequence is the relevant factor
        compression_scaling = 4
    else:
        compression_scaling = 1

    if input_size_mb < threshold_tiny * compression_scaling:
        time_hrs = 0
    elif input_size_mb < threshold_small * compression_scaling:
        time_hrs = 1
    elif input_size_mb < threshold_normal * compression_scaling:
        time_hrs = 71
    else:
        time_hrs = 84
    return time_hrs
