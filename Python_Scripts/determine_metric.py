from IO import BedBase, IncompatabilityError


def calculate_metric(reference_labelled_peaks: BedBase,
                     psuedopeaks: BedBase,
                     include_merged_peaks: bool = True) -> float:
    """
    Calculates the metric that compares the number of peaks in the reference
    dataset to the number of peaks in the pseudopeaks dataset.

    Args:
        reference_labelled_peaks: BedBase object containing the reference
            dataset.
        psuedopeaks: BedBase object containing the pseudopeaks dataset.
        include_merged_peaks: Whether to include merged peaks in the
            calculation. If False, only single peaks will be counted.

    Returns:
        The metric value.
    """
    if not reference_labelled_peaks.has_same_positions(psuedopeaks):
        raise IncompatabilityError(
            "Reference peaks must be over the same region as pseudopeaks.")
    peaks_in_reference = reference_labelled_peaks.get("SCORE")
    if include_merged_peaks:
        number_of_reference_peaks = (peaks_in_reference > 1).sum()
    else:
        number_of_reference_peaks = (peaks_in_reference == 1).sum()

    number_of_pseudopeaks = (psuedopeaks.get("SCORE") == 1).sum()
    metric = number_of_pseudopeaks / number_of_reference_peaks
    if metric > 1:
        print("Peaks in comparison dataset is greater than in reference "
              "dataset. For a better result, consider switching the order of ",
              "each dataset.")

    return metric
