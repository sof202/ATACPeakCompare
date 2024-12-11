import numpy as np
import pandas as pd
from IO import BedBase, BedBaseCI, IncompatabilityError
from scipy.stats import poisson, norm
from typing import NamedTuple


class ConfidenceInterval(NamedTuple):
    def __init__(self, lower: np.ndarray, upper: np.ndarray):
        self.lower = lower
        self.upper = upper


def calculate_lambda_ci(lambdas: np.ndarray,
                        significance: float,
                        window_size: int = 50) -> ConfidenceInterval:
    variances = np.zeros_like(lambdas)
    for i in range(len(lambdas)):
        start = max(0, i - window_size // 2)
        end = min(len(lambdas), i + window_size // 2 + 1)
        variances[i] = np.var(lambdas[start:end])
    standard_error = np.sqrt(variances / window_size)
    z_a = norm.ppf(significance)
    lower = lambdas - z_a * standard_error
    upper = lambdas + z_a * standard_error

    # Poisson distribution doesn't take kindly to non-positive lambdas
    lower = np.clip(lower, a_min=np.min(lambdas), a_max=None)
    return ConfidenceInterval(lower=lower, upper=upper)


def generate_bias_track_ci(bias_bedbase: BedBase,
                           significance: float) -> BedBaseCI:
    lambda_ci = calculate_lambda_ci(
        lambdas=bias_bedbase.get("SCORE").to_numpy(),
        significance=significance
    )
    bias_bedbase_ci = BedBaseCI(
        CHR=bias_bedbase.get("CHR"),
        BASE=bias_bedbase.get("BASE"),
        LOWER_SCORE=pd.Series(lambda_ci[0]),
        UPPER_SCORE=pd.Series(lambda_ci[1])
    )
    return bias_bedbase_ci


def calculate_pavlue(reads: np.ndarray, lambdas: np.ndarray) -> np.ndarray:
    return poisson.cdf(reads, lambdas)


def generate_pvalue_ci(bias_bedbase: BedBase,
                       coverage_bedbase: BedBase,
                       significance: float,
                       window_size: int = 50) -> BedBaseCI:
    if not bias_bedbase.has_same_positions(coverage_bedbase):
        raise IncompatabilityError(
            "Bias track and coverage track are over different regions.")

    bias_bedbase_ci = generate_bias_track_ci(
        bias_bedbase,
        significance,
        window_size
    )

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
