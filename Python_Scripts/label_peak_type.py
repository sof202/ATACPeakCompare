from extract_region import extract_bedbase_region
from IO import Bed, BedGraph, BedBase, IncompatabilityError
import numpy as np
import pandas as pd


def convert_narrow_peak_to_bedbase(peak_data: Bed,
                                   chromosome: str,
                                   start: int,
                                   end: int) -> BedBase:
    """Converts peaks in a selected region of narrow peak data into bedbase
    format

    Args
        peak_data (pandas.DataFrame): Narrow peak data in BED3+7 format.
        chromosome (str): Chromosome to extract.
        start (int): Start of region.
        end (int): End of region.

    Returns:
        A pandas.DataFrame in bedbase format of the selected region where the
        score column is 0 if no peak is at that base, and 1 if there is a peak.
    """
    peak_data = BedGraph(
        CHR=peak_data.get("CHR"),
        START=peak_data.get("START"),
        END=peak_data.get("END"),
        SCORE=pd.Series(np.zeros(len(peak_data.get("START"))))
    )
    peak_data = extract_bedbase_region(peak_data, chromosome, start, end)

    # extract_bedbase_region() will return NaN values in the SCORE column for
    # all bases that do not exist in the narrow peak file. This is exactly the
    # bases that are not within a peak. A trick here is to convert to Boolean
    # values to convert all non-NaN scores into 1 and all NaN scores to 0.
    score = peak_data.get("SCORE").fillna(0).astype(bool).astype(int)
    peak_data.get()["SCORE"] = score
    return peak_data


def label_peak_type(unmerged_peaks: BedBase,
                    merged_peaks: BedBase) -> BedBase:
    """Using two bedbase files from narrow peak files the types of peak are
    determined.

    Args:
        unmerged_peaks (pandas DataFrame): A set of peaks before merging peaks
        in bedbase format.
        merged_peaks (pandas DataFrame): A set of peaks after merging peaks in
        bedbase format.

    Returns:
        A pandas DataFrame in bedbase format with a score column indicating the
        peak type: 0 indicates not a peak, 1 indicates a peak before merging,
        2 indicates a peak only after mering. Returns None if the input
        data frames do not align (different bases).
    """
    if not unmerged_peaks.has_same_positions(merged_peaks):
        raise IncompatabilityError(
            "Merged and unmerged peaks have incompatible regions")

    unmerged_peaks = unmerged_peaks.get()
    merged_peaks = merged_peaks.get()
    labelled_peaks = pd.merge(
        unmerged_peaks,
        merged_peaks,
        on=["CHR", "BASE"],
        how="inner"
    )

    # Note: Because MACS deletes small peaks AFTER merging, the peaks in the
    # merged peaks file is necessarily a superset of the peaks in the unmerged
    # peaks file.
    labelled_peaks["SCORE"] = labelled_peaks["SCORE_x"] + \
        labelled_peaks["SCORE_y"]
    labelled_peaks = BedBase(
        labelled_peaks["CHR"],
        labelled_peaks["BASE"],
        labelled_peaks["SCORE"]
    )
    return labelled_peaks
