import argparse
import sys
from config_file_functions import (
    validate_eol_format,
    get_config_variables,
    is_positive_integer,
    is_larger,
    is_variable_missing,
    path_exists
)


def validate_variable_existence(config_variables: dict) -> None:
    variable_missing = False
    # Always required
    variables_to_check = [
        "DEBUG_MODE",
        "OUTPUT_DIRECTORY",
        "INPUT_FILE",
        "SAMPLE_NAME",
        "FILE_TYPE",
        "BUILD_MODEL",
        "GENOME_SIZE",
        "SMALL_LOCAL_SIZE",
        "LARGE_LOCAL_SIZE"
    ]
    if config_variables["BUILD_MODEL"] == "1":
        variables_to_check.append("MFOLD_LOWER")
        variables_to_check.append("MFOLD_UPPER")
    else:
        if not is_variable_missing("CONTROL_FILE", config_variables, True):
            variables_to_check.append("NUMBER_OF_CONTROL_READS")
            variables_to_check.append("CONTROL_FRAGMENT_LENGTH")
        variables_to_check.append("READ_LENGTH")
        variables_to_check.append("FRAGMENT_LENGTH")
        variables_to_check.append("NUMBER_OF_READS")
    for variable in variables_to_check:
        variable_missing = variable_missing or is_variable_missing(
            variable, config_variables)

    cutoff_missing = is_variable_missing(
        "CUTOFF",
        config_variables,
        quiet=True
    )
    average_peak_missing = is_variable_missing(
        "AVERAGE_PEAK_LENGTH",
        config_variables,
        quiet=True
    )
    if cutoff_missing and average_peak_missing:
        variable_missing = True
        print(
            "You must have one of CUTOFF or AVERAGE_PEAK_LENGTH defined in",
            "the config file."
        )

    # If a required variable is missing, we cannot continue, so exit now.
    if variable_missing:
        sys.exit(1)


def any_variables_incorrect(config_variables: dict) -> bool:
    variable_incorrect = False
    positive_integer_variables = [
        "GENOME_SIZE",
        "SMALL_LOCAL_SIZE",
        "LARGE_LOCAL_SIZE"
    ]
    ascending_variables = [
        ("SMALL_LOCAL_SIZE", "LARGE_LOCAL_SIZE")
    ]
    if config_variables["BUILD_MODEL"] == 1:
        positive_integer_variables.append("MFOLD_LOWER")
        positive_integer_variables.append("MFOLD_UPPER")
        ascending_variables.append(("MFOLD_LOWER", "MFOLD_UPPER"))
    else:
        if not is_variable_missing("CONTROL_FILE", config_variables, True):
            positive_integer_variables.append("NUMBER_OF_CONTROL_READS")
            positive_integer_variables.append("CONTROL_FRAGMENT_LENGTH")
        positive_integer_variables.append("READ_LENGTH")
        positive_integer_variables.append("FRAGMENT_LENGTH")
        positive_integer_variables.append("NUMBER_OF_READS")
    if is_variable_missing("CUTOFF", config_variables, quiet=True):
        positive_integer_variables.append("AVERAGE_PEAK_LENGTH")

    for variable in positive_integer_variables:
        variable_incorrect = variable_incorrect or is_positive_integer(
            config_variables[variable],
            variable
        )
    for pair in ascending_variables:
        variable_pair = tuple([config_variables[key] for key in pair])
        variable_incorrect = variable_incorrect or is_larger(*variable_pair)
    return variable_incorrect


def any_file_paths_missing(config_variables: dict) -> bool:
    file_missing = False
    paths_to_check = [
        "INPUT_FILE"
    ]
    if not is_variable_missing("CONTROL_FILE", config_variables):
        paths_to_check.append("CONTROL_FILE")
    for path in paths_to_check:
        file_missing = file_missing or path_exists(
            config_variables[path], path
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
    config_malformed = any_variables_incorrect(config_variables)
    config_malformed = config_malformed or any_file_paths_missing(
        config_variables)
    if config_malformed:
        sys.exit(1)
    sys.exit(0)
