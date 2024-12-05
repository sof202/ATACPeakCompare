import argparse
import sys
from pathlib import Path


def validate_eol_format(file_path: str) -> None:
    try:
        with open(file_path, "rb") as config_file:
            content: bytes = config_file.read()
    except IOError:
        print(
            "Could not read file",
            file_path,
            "Please ensure that this file exists and has read permissions."
        )
        sys.exit(1)
    if b'\r' in content:
        print(
            "Carriage return character found in file.",
            "Please ensure that config file uses Linux EOL characters only."
        )
        sys.exit(1)


def get_config_variables(file_path: str) -> dict:
    config_variables: dict = {}
    with open(file_path, "r") as config_file:
        for line in config_file:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            variable, value = line.split("=", 1)
            value = value.strip('"')
            if value == "":
                print(
                    "Could not parse config file.",
                    f"No value was given for {variable}."
                )
                sys.exit(1)
            config_variables[variable] = value
    return config_variables


def is_variable_missing(variable: str, config_variables: dict) -> bool:
    if variable not in config_variables:
        print(
            f"{variable} is missing from config file.",
            "Please check the example config for what is required."
        )
        return True
    return False


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
    if config_variables["BUILD_MODEL"] == 1:
        variables_to_check.append("MFOLD_LOWER")
        variables_to_check.append("MFOLD_UPPER")
    else:
        if "CONTROL_FILE" in config_variables:
            variables_to_check.append("NUMBER_OF_CONTROL_READS")
            variables_to_check.append("CONTROL_FRAGMENT_LENGTH")
        variables_to_check.append("READ_LENGTH")
        variables_to_check.append("FRAGMENT_LENGTH")
        variables_to_check.append("NUMBER_OF_READS")
    for variable in variables_to_check:
        variable_missing = variable_missing or is_variable_missing(
            variable, config_variables)

    cutoff_missing = is_variable_missing("CUTOFF", config_variables)
    average_peak_missing = is_variable_missing(
        "AVERAGE_PEAK_LENGTH", config_variables)
    if cutoff_missing and average_peak_missing:
        print(
            "You must have one of CUTOFF or AVERAGE_PEAK_LENGTH defined in",
            "the config file."
        )

    # If a required variable is missing, we cannot continue, so exit now.
    if variable_missing:
        sys.exit(1)


def is_positive_integer(x: str) -> bool:
    message = f"{x} must be a positive integer."
    if not x.isdigit():
        print(message)
        return False
    if int(x) < 0:
        print(message)
        return False
    return True


def is_larger(x: str, y: str) -> bool:
    if not (x.isdigit() and y.isdigit()):
        print(f"{x} and {y} must be integers.")
        return False
    if int(x) >= int(y):
        print(f"{y} must be strictly greater than {x}")
        return False
    return True


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
        if not is_variable_missing("CONTROL_FILE", config_variables):
            positive_integer_variables.append("NUMBER_OF_CONTROL_READS")
            positive_integer_variables.append("CONTROL_FRAGMENT_LENGTH")
        positive_integer_variables.append("READ_LENGTH")
        positive_integer_variables.append("FRAGMENT_LENGTH")
        positive_integer_variables.append("NUMBER_OF_READS")
    if is_variable_missing("CUTOFF", config_variables):
        positive_integer_variables.append("AVERAGE_PEAK_LENGTH")

    for variable in positive_integer_variables:
        variable_incorrect = variable_incorrect or is_positive_integer(
            config_variables[variable]
        )
    for pair in ascending_variables:
        variable_pair = tuple([config_variables[key] for key in pair])
        variable_incorrect = variable_incorrect or is_larger(*variable_pair)
    return variable_incorrect


def path_exists(file_path: str, variable_name: str) -> bool:
    if Path(file_path).is_file():
        print(
            f"The path given for {variable_name} doesn't exist.",
            f"Please check this path: {file_path}"
        )
        return False
    return True


def any_file_paths_missing(config_variables: dict) -> bool:
    file_missing = False
    paths_to_check = [
        "OUTPUT_DIRECTORY",
        "INPUT_FILE",
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
    config_malformed = config_malformed or any_variables_incorrect(
        config_variables)
    config_malformed = config_malformed or any_file_paths_missing(
        config_variables)
    if config_malformed:
        sys.exit(1)
    sys.exit(0)
