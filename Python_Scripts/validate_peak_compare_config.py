import argparse
import sys
from config_file_functions import (
    validate_eol_format,
    get_config_variables,
    is_positive_integer,
    is_larger,
    is_positive_float,
    is_variable_missing,
    path_does_not_exist
)


def validate_variable_existence(config_variables: dict) -> None:
    variable_missing = False
    # Always required
    variables_to_check = [
        "CHROMOSOME",
        "START",
        "END",
        "CUTOFF",
        "REFERENCE_DATASET_OUTPUT_DIRECTORY",
        "REFERENCE_SAMPLE_NAME",
        "COMPARISON_DATASET_OUTPUT_DIRECTORY",
        "COMPARISON_SAMPLE_NAME",
        "REFERENCE_MERGED_PEAK_FILE",
        "REFERENCE_UNMERGED_PEAK_FILE",
        "REFERENCE_BIAS_TRACK_FILE",
        "REFERENCE_COVERAGE_TRACK_FILE",
        "COMPARISON_BIAS_TRACK_FILE",
        "COMPARISON_COVERAGE_TRACK_FILE",
        "COMPARISON_PVALUE_FILE",
        "LOG_DIRECTORY"
    ]
    for variable in variables_to_check:
        variable_missing = variable_missing or is_variable_missing(
            variable, config_variables)
    # If a required variable is missing, we cannot continue, so exit now.
    if variable_missing:
        sys.exit(1)


def all_variables_correct(config_variables: dict) -> bool:
    variables_correct = True
    positive_integer_variables = [
        "START",
        "END"
    ]
    for variable in positive_integer_variables:
        variables_correct = variables_correct and is_positive_integer(
            config_variables[variable],
            variable
        )
    variables_correct = variables_correct and is_positive_float(
        config_variables["CUTOFF"],
        "CUTOFF"
    )
    variables_correct = variables_correct and is_larger(
        config_variables["START"],
        config_variables["END"]
    )
    return variables_correct


def any_file_paths_missing(config_variables: dict) -> bool:
    file_missing = False
    paths_to_check = [
        "REFERENCE_DATASET_OUTPUT_DIRECTORY",
        "COMPARISON_DATASET_OUTPUT_DIRECTORY"
    ]
    for path in paths_to_check:
        file_missing = file_missing or path_does_not_exist(
            config_variables[path], path, dir=True
        )
    return file_missing


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="ChromCompare config file checker",
        description="Checks whether config file for ChromCompare is valid."
    )
    parser.add_argument('file_path')
    args = parser.parse_args()
    validate_eol_format(args.file_path)
    config_variables = get_config_variables(args.file_path)
    validate_variable_existence(config_variables)
    config_malformed = False
    config_malformed = not all_variables_correct(config_variables)
    config_malformed = config_malformed or any_file_paths_missing(
        config_variables)
    if config_malformed:
        sys.exit(1)
    sys.exit(0)
