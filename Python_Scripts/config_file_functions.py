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


def is_variable_missing(variable: str,
                        config_variables: dict,
                        quiet: bool = False) -> bool:
    if variable not in config_variables:
        if not quiet:
            print(
                f"{variable} is missing from config file.",
                "Please check the example config for what is required."
            )
        return True
    return False


def is_positive_float(x: str, variable_name: str) -> bool:
    message = f"{variable_name} must be a positive floating point number."
    try:
        if float(x) > 0:
            return True
        else:
            print(message)
            return False
    except ValueError:
        print(message)
        return False


def is_positive_integer(x: str, variable_name: str) -> bool:
    message = f"{variable_name} must be a positive integer."
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


def path_does_not_exist(file_path: str, variable_name: str) -> bool:
    if not Path(file_path).is_file():
        print(
            f"The path given for {variable_name} doesn't exist.",
            f"Please check this path: {file_path}"
        )
        return True
    return False
