from IO import BedBase, BedBaseCI, IncompatabilityError


def compare_pvalue_ci(reference_pvalue_ci: BedBaseCI,
                      comparison_pvalue_ci: BedBaseCI) -> BedBase:
    """
    Compares the smallest pvalue of the reference dataset with the largest
    from the comparison dataset. This is essential knowledge required for
    determining pseudopeaks.

    Args:
        reference_pvalue_ci: BedBaseCI object containing the reference
            dataset.
        comparison_pvalue_ci: BedBaseCI object containing the comparison
            dataset.

    Returns:
        A BedBase object containing the comparison dataset.
    """
    if not reference_pvalue_ci.has_same_positions(comparison_pvalue_ci):
        raise IncompatabilityError(
            "Reference and comparison files must be over the same region.")
    # For a base in the comparison dataset to be in contention to be a
    # psuedopeak, one part of the criteria is to have an overlapping (or
    # better) confidence interval with the reference dataset
    pvalues_to_beat = reference_pvalue_ci.get("LOWER_SCORE").to_numpy()
    contender_pvalues = comparison_pvalue_ci.get("UPPER_SCORE").to_numpy()
    is_significant = (contender_pvalues > pvalues_to_beat)

    compared_pvalues = BedBase(
        CHR=reference_pvalue_ci.get("CHR"),
        BASE=reference_pvalue_ci.get("BASE"),
        SCORE=is_significant.astype(int)
    )
    return compared_pvalues


def determine_psuedopeaks(comparison_pvalues: BedBase,
                          compared_pvalues: BedBase,
                          reference_labelled_peaks: BedBase,
                          cutoff: float) -> BedBase:
    """
    Determines the pseudopeaks from the comparison dataset. For more
    information on pseudopeaks, go to the documentation.

    Args:
        comparison_pvalues: BedBase object containing the comparison dataset.
        compared_pvalues: BedBase object containing the comparison dataset.
        reference_labelled_peaks: BedBase object containing the reference
            dataset.
        cutoff: The cutoff for determining whether a peak is a pseudopeak.

    Returns:
        A BedBase object containing the pseudopeaks.
    """
    if not comparison_pvalues.has_same_positions(compared_pvalues):
        raise IncompatabilityError(
            "All BedBase files must be over the same region.")
    if not comparison_pvalues.has_same_positions(reference_labelled_peaks):
        raise IncompatabilityError(
            "All BedBase files must be over the same region.")
    if not compared_pvalues.has_same_positions(reference_labelled_peaks):
        raise IncompatabilityError(
            "All BedBase files must be over the same region.")

    pvalue = comparison_pvalues.get("SCORE").to_numpy()
    passed_ci_comparison = compared_pvalues.get("SCORE").to_numpy()
    peak_type = reference_labelled_peaks.get("SCORE").to_numpy()

    # criteria for psuedopeak
    is_pseudopeak = (
        (passed_ci_comparison == 1) & (peak_type == 2)
    ) | (
        (pvalue > cutoff) & (peak_type == 1)
    )
    pseudopeaks = BedBase(
        CHR=comparison_pvalues.get("CHR"),
        BASE=comparison_pvalues.get("BASE"),
        SCORE=is_pseudopeak.astype(int)
    )
    return pseudopeaks
