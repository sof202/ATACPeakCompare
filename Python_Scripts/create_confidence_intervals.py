import numpy as np
import pandas as pd
from IO import BedBase, BedBaseCI, IncompatabilityError
from scipy.stats import poisson, norm
from typing import NamedTuple


class ConfidenceInterval(NamedTuple):
    """
    Represents a set of confidence intervals.
    """
    lower: np.ndarray
    upper: np.ndarray


def calculate_lambda_ci(lambdas: np.ndarray,
                        significance: float = 0.95,
                        window_size: int = 50) -> ConfidenceInterval:
    """
    Calculates confidence intervals for lambda values, considering variance in
    surrounding values, handling edge cases correctly.

    Args:
        lambdas: A NumPy array of lambda values.
        significance: The significance level for the confidence interval
        (e.g., 0.95 for a 95% CI).
        window_size: The size of the sliding window to calculate variance.

    Returns:
        A ConfidenceInterval object containing the lower and upper bounds of
        the confidence interval.
    """
    variances = np.zeros_like(lambdas)
    sample_sizes = np.zeros_like(lambdas)
    for i in range(len(lambdas)):
        start = max(0, i - window_size // 2)
        end = min(len(lambdas), i + window_size // 2 + 1)
        window = lambdas[start:end]
        variances[i] = np.var(window)
        sample_sizes[i] = len(window)

    standard_errors = np.sqrt(variances / sample_sizes)
    z_a = norm.ppf(significance)
    lower = lambdas - z_a * standard_errors
    upper = lambdas + z_a * standard_errors

    # Poisson distribution doesn't take kindly to non-positive lambdas
    lower = np.clip(lower, a_min=np.min(lambdas), a_max=None)
    return ConfidenceInterval(lower=lower, upper=upper)


def generate_bias_track_ci(bias_bedbase: BedBase,
                           significance: float = 0.95,
                           window_size: int = 50) -> BedBaseCI:
    lambda_ci = calculate_lambda_ci(
        lambdas=bias_bedbase.get("SCORE").to_numpy(),
        significance=significance,
        window_size=window_size
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
                       significance: float = 0.95,
                       window_size: int = 50) -> BedBaseCI:
    """
    Generates a confidence interval for the pvalue of the coverage track
    given the bias track.

    Args:
        bias_bedbase: A BedBase object containing the bias track.
        coverage_bedbase: A BedBase object containing the coverage track.
        significance: The significance level for the confidence interval
        (e.g., 0.95 for a 95% CI).
        window_size: The size of the sliding window to calculate variance.

    Returns:
        A BedBaseCI object containing the lower and upper bounds of the
        confidence interval.
    """
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
