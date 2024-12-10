import argparse
from create_confidence_intervals import generate_pvalue_ci
from determine_metric import calculate_metric
from determine_psuedo_peaks import compare_pvalue_ci, determine_psuedopeaks
from extract_region import extract_bedbase_region
from label_peak_type import label_peak_type, convert_narrow_peak_to_bedbase
from IO import BedGraph, Bed
import sys


def main(args: argparse.Namespace) -> None:
    chromosome = args.chromosome
    try:
        start = int(args.start)
        end = int(args.end)
    except ValueError:
        print("start and end must be integers")
        sys.exit(1)

    try:
        cutoff = float(args.cutoff)
    except ValueError:
        print("cutoff must be a floating point number")
        sys.exit(1)

    reference_merged_peaks = convert_narrow_peak_to_bedbase(
        Bed.read_from_file(args.reference_merged_peaks_file),
        chromosome,
        start,
        end
    )
    reference_unmerged_peaks = convert_narrow_peak_to_bedbase(
        Bed.read_from_file(args.reference_unmerged_peaks_file),
        chromosome,
        start,
        end
    )
    reference_bias_track = extract_bedbase_region(
        BedGraph.read_from_file(args.reference_bias_track_file),
        chromosome,
        start,
        end
    )
    reference_coverage_track = extract_bedbase_region(
        BedGraph.read_from_file(args.reference_coverage_track_file),
        chromosome,
        start,
        end
    )
    comparison_bias_track = extract_bedbase_region(
        BedGraph.read_from_file(args.comparison_bias_track_file),
        chromosome,
        start,
        end
    )
    comparison_coverage_track = extract_bedbase_region(
        BedGraph.read_from_file(args.comparison_coverage_track_file),
        chromosome,
        start,
        end
    )
    comparison_pvalue_track = extract_bedbase_region(
        BedGraph.read_from_file(args.comparison_pvalue_file),
        chromosome,
        start,
        end
    )
    reference_labelled_peaks = label_peak_type(
        reference_merged_peaks,
        reference_unmerged_peaks
    )
    reference_pvalue_ci = generate_pvalue_ci(
        reference_bias_track,
        reference_coverage_track
    )
    comparison_pvalue_ci = generate_pvalue_ci(
        comparison_bias_track,
        comparison_coverage_track
    )
    compared_pvalues = compare_pvalue_ci(
        reference_pvalue_ci,
        comparison_pvalue_ci
    )
    pseudopeaks = determine_psuedopeaks(
        comparison_pvalue_track,
        compared_pvalues,
        reference_labelled_peaks,
        cutoff
    )
    if args.unmerged:
        metric = calculate_metric(
            reference_labelled_peaks,
            pseudopeaks,
            include_merged_peaks=False
        )
    else:
        metric = calculate_metric(
            reference_labelled_peaks,
            pseudopeaks,
            include_merged_peaks=True
        )
    print(f"The metric for these datasets over the range {start} to {end} ",
          f"for chromosome {chromosome} is: {metric}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="PeakCompare",
        description="Determine how likely a peak is to be in two datasets."
    )
    parser.add_argument(
        "--unmerged",
        action="store_true",
        help=("Set this if you want to discount peaks that are a result of "
              "merging when calculating the metric.")
    )
    parser.add_argument(
        "chromosome",
        help="The chromosome of the region you wish to inspect."
    )
    parser.add_argument(
        "start",
        help=("The base pair position at the start of the region you wish to "
              "inspect")
    )
    parser.add_argument(
        "end",
        help=("The base pair position at the end of the region you wish to "
              "inspect")
    )
    parser.add_argument(
        "reference_merged_peaks_file",
        help=("The narrow peak file from reference dataset where peaks are"
              "merged.")
    )
    parser.add_argument(
        "reference_unmerged_peaks_file",
        help=("The narrow peak file from reference dataset where peaks are"
              "not merged.")
    )
    parser.add_argument(
        "reference_bias_track_file",
        help="The bias track file for the reference dataset."
    )
    parser.add_argument(
        "reference_coverage_track_file",
        help="The coverage track (pileup) file for the reference dataset."
    )
    parser.add_argument(
        "comparison_bias_track_file",
        help="The bias track file for the comparison dataset."
    )
    parser.add_argument(
        "comparison_coverage_track_file",
        help="The coverage track (pileup) file for the comparison dataset."
    )
    parser.add_argument(
        "comparison_pvalue_file",
        help="The pvalues for the comparison dataset."
    )
    parser.add_argument(
        "cutoff",
        help="The cutoff used to call peaks in the reference dataset."
    )
    args = parser.parse_args()
    main(args)
