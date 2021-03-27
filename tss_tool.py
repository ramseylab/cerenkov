def min_tss_dist(dist_seq):
    """
    Find the TSS distance in `dist_seq` with minimum absolute value.

    We prefer negative TSS distances because we prefer upstream SNPs.
    That is to say, if the minimum absolute value in `dist_seq` is `n`,
        and both `n` and `-n` are present in `dist_seq`, choose `-n`
    The preference is performed by adding 0.1 to the distances and then rank their absolute value ascendingly.
    E.g. |-5 + 0.1| < |5 + 0.1|, so -5 is picked.

    :param dist_seq:
    :return:
    """
    amended_dist = [abs(i + 0.1) for i in dist_seq]
    # `index` only returns index of the first instance even if there are multiple min values
    min_index = amended_dist.index(min(amended_dist))

    return dist_seq[min_index]
