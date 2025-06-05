
def select_combination_input(wildcards):

    if COMBINE_LABELS_BY_GROUP:
        _MULTI_ANNOTATION_INPUT = _MULTI_ANNOTATION_INPUT_GRP[wildcards.group_prefix]
    else:
        _MULTI_ANNOTATION_INPUT = _MULTI_ANNOTATION_INPUT_ALL
    return _MULTI_ANNOTATION_INPUT


def select_combination_labels(wildcards):

    if COMBINE_LABELS_BY_GROUP:
        _MULTI_ANNOTATION_LABELS = _MULTI_ANNOTATION_LABELS_GRP[wildcards.group_prefix]
    else:
        _MULTI_ANNOTATION_LABELS = _MULTI_ANNOTATION_LABELS_ALL
    return _MULTI_ANNOTATION_LABELS
