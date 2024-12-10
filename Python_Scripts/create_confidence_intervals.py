import numpy as np
import pandas as pd
from IO import BedBase, BedBaseCI, IncompatabilityError
from scipy.stats import poisson
from typing import NamedTuple


class ConfidenceInterval(NamedTuple):
    lower: np.ndarray
    upper: np.ndarray


def caculate_exact_lambda_ci(lambdas: np.ndarray,
                             reads: np.ndarray) -> ConfidenceInterval:
    # zero reads results in divide by zero. Adding a small number here results
    # in huge confidence intervals, adding a large number here is dishonest.
    # A middle ground here is adding 1 or 2 to the read count.
    if any(reads == 0):
        reads += 1
    lower = np.exp(np.log(lambdas) - np.sqrt(1 / reads * lambdas))
    upper = np.exp(np.log(lambdas) + np.sqrt(1 / reads * lambdas))
    return ConfidenceInterval(lower=lower, upper=upper)


def generate_bias_track_ci(bias_bedbase: BedBase,
                           coverage_bedbase: BedBase
                           ) -> BedBaseCI:
    lambda_ci = caculate_exact_lambda_ci(
        lambdas=bias_bedbase.get("SCORE").to_numpy(),
        reads=coverage_bedbase.get("SCORE").to_numpy()
    )
    bias_bedbase_ci = BedBaseCI(
        CHR=coverage_bedbase.get("CHR"),
        BASE=coverage_bedbase.get("BASE"),
        LOWER_SCORE=pd.Series(lambda_ci[0]),
        UPPER_SCORE=pd.Series(lambda_ci[1])
    )
    return bias_bedbase_ci


def calculate_pavlue(reads: np.ndarray, lambdas: np.ndarray) -> np.ndarray:
    return poisson.cdf(reads, lambdas)


def generate_pvalue_ci(bias_bedbase: BedBase,
                       coverage_bedbase: BedBase) -> BedBaseCI:
    if not bias_bedbase.has_same_positions(coverage_bedbase):
        raise IncompatabilityError(
            "Bias track and coverage track are over different regions.")

    bias_bedbase_ci = generate_bias_track_ci(bias_bedbase, coverage_bedbase)

    # A higher lambda in the poisson distribution will cause the same number
    # of reads to generate a lower pvalue. To stay consistent with naming, we
    # switch the order of upper and lower below.
    lower_pvalue = calculate_pavlue(
        coverage_bedbase.get("SCORE"),
        bias_bedbase_ci.get("UPPER_SCORE")
    )
    upper_pvalue = calculate_pavlue(
        coverage_bedbase.get("SCORE"),
        bias_bedbase_ci.get("LOWER_SCORE")
    )
    pvalues_bedbase_ci = BedBaseCI(
        CHR=coverage_bedbase.get("CHR"),
        BASE=coverage_bedbase.get("BASE"),
        LOWER_SCORE=pd.Series(np.nan_to_num(lower_pvalue)),
        UPPER_SCORE=pd.Series(np.nan_to_num(upper_pvalue))
    )
    return pvalues_bedbase_ci
