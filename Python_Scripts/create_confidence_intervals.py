import numpy as np
import pandas as pd
from IO import BedBase
from scipy.stats import poisson
from typing import NamedTuple


class ConfidenceInterval(NamedTuple):
    lower: np.ndarray
    upper: np.ndarray


def caculate_exact_lambda_ci(lambdas: np.ndarray,
                             reads: np.ndarray) -> ConfidenceInterval:
    lower = np.exp(np.log(lambdas) - np.sqrt(1 / reads * lambdas))
    upper = np.exp(np.log(lambdas) + np.sqrt(1 / reads * lambdas))
    return ConfidenceInterval(lower=lower, upper=upper)


def generate_bias_track_ci(bias_bedbase: BedBase,
                           coverage_bedbase: BedBase
                           ) -> pd.DataFrame:
    bias_bedbase = bias_bedbase.get()
    coverage_bedbase = coverage_bedbase.get()
    lambda_ci = caculate_exact_lambda_ci(
        lambdas=bias_bedbase["SCORE"].to_numpy(),
        reads=coverage_bedbase["SCORE"].to_numpy()
    )
    bias_bedbase_ci = pd.DataFrame({
        "CHR": coverage_bedbase["CHR"],
        "BASE": coverage_bedbase["BASE"],
        "LOWER": pd.Series(lambda_ci[0]),
        "UPPER": pd.Series(lambda_ci[1])
    })
    return bias_bedbase_ci


def calculate_pavlue(reads: np.ndarray, lambdas: np.ndarray) -> np.ndarray:
    return poisson.cdf(reads, lambdas)


def generate_pvalue_ci(bias_bedbase: BedBase,
                       coverage_bedbase: BedBase) -> pd.DataFrame:
    bias_bedbase_ci = generate_bias_track_ci(bias_bedbase, coverage_bedbase)
    coverage_bedbase = coverage_bedbase.get()
    bias_bedbase_ci = bias_bedbase_ci.get()

    lower_pvalue = calculate_pavlue(
        coverage_bedbase["SCORE"],
        bias_bedbase_ci["LOWER"]
    )
    upper_pvalue = calculate_pavlue(
        coverage_bedbase["SCORE"],
        bias_bedbase_ci["UPPER"]
    )
    pvalues_bedbase_ci = pd.DataFrame({
        "CHR": coverage_bedbase["CHR"],
        "BASE": coverage_bedbase["BASE"],
        "LOWER": pd.Series(lower_pvalue),
        "UPPER": pd.Series(upper_pvalue)
    })
    return pvalues_bedbase_ci
