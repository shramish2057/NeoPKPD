"""
OpenPKPD NCA Bioequivalence - FDA/EMA compliant bioequivalence analysis.

This module provides bioequivalence assessment tools including:
- 90% confidence interval calculation
- TOST (Two One-Sided Tests) analysis
- Geometric mean ratio calculation
- Within-subject CV estimation
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

from .._core import _require_julia, _to_julia_float_vector


# ============================================================================
# Geometric Mean and Ratio
# ============================================================================

def geometric_mean(values: List[float]) -> float:
    """
    Calculate geometric mean.

    GM = exp(mean(log(values)))

    Args:
        values: List of positive values

    Returns:
        Geometric mean

    Example:
        >>> gm = geometric_mean([10, 20, 30])
        >>> print(f"Geometric mean: {gm:.2f}")
    """
    jl = _require_julia()
    v = _to_julia_float_vector(jl, values)
    return float(jl.OpenPKPDCore.geometric_mean(v))


def geometric_mean_ratio(
    test_values: List[float],
    reference_values: List[float]
) -> float:
    """
    Calculate geometric mean ratio (GMR) of test to reference.

    GMR = exp(mean(log(test)) - mean(log(reference)))

    Args:
        test_values: Test formulation values (Cmax or AUC)
        reference_values: Reference formulation values

    Returns:
        Geometric mean ratio (test/reference)

    Example:
        >>> gmr = geometric_mean_ratio(test_cmax, reference_cmax)
        >>> print(f"GMR: {gmr:.4f}")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, test_values)
    r = _to_julia_float_vector(jl, reference_values)
    return float(jl.OpenPKPDCore.geometric_mean_ratio(t, r))


def within_subject_cv(
    test_values: List[float],
    reference_values: List[float]
) -> float:
    """
    Estimate within-subject coefficient of variation from crossover data.

    CV_intra = sqrt(exp(MSE) - 1) * 100%

    Where MSE is from ANOVA on log-transformed data.

    Args:
        test_values: Test values (paired with reference)
        reference_values: Reference values (paired)

    Returns:
        Within-subject CV (%)

    Example:
        >>> cv = within_subject_cv(test_cmax, reference_cmax)
        >>> print(f"Intra-subject CV: {cv:.1f}%")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, test_values)
    r = _to_julia_float_vector(jl, reference_values)
    return float(jl.OpenPKPDCore.within_subject_cv(t, r))


# ============================================================================
# 90% Confidence Interval
# ============================================================================

def bioequivalence_90ci(
    test_values: List[float],
    reference_values: List[float],
    log_transform: bool = True
) -> Dict[str, float]:
    """
    Calculate 90% confidence interval for geometric mean ratio.

    Uses the standard two-sequence, two-period crossover analysis.

    Args:
        test_values: Test formulation values (paired with reference)
        reference_values: Reference formulation values
        log_transform: Apply log transformation (default: True)

    Returns:
        Dict with keys:
        - gmr: Geometric mean ratio
        - ci_lower: Lower bound of 90% CI
        - ci_upper: Upper bound of 90% CI
        - cv_intra: Intra-subject CV (%)
        - n: Number of subjects

    Example:
        >>> result = bioequivalence_90ci(test_cmax, reference_cmax)
        >>> print(f"GMR: {result['gmr']:.4f}")
        >>> print(f"90% CI: ({result['ci_lower']:.4f}, {result['ci_upper']:.4f})")
        >>> if result['ci_lower'] >= 0.80 and result['ci_upper'] <= 1.25:
        ...     print("Bioequivalent!")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, test_values)
    r = _to_julia_float_vector(jl, reference_values)
    result = jl.OpenPKPDCore.bioequivalence_90ci(t, r, log_transform=log_transform)

    return {
        "gmr": float(result.gmr),
        "ci_lower": float(result.ci_lower),
        "ci_upper": float(result.ci_upper),
        "cv_intra": float(result.cv_intra),
        "n": int(result.n),
    }


# ============================================================================
# TOST Analysis
# ============================================================================

def tost_analysis(
    test_values: List[float],
    reference_values: List[float],
    theta_lower: float = 0.80,
    theta_upper: float = 1.25,
    alpha: float = 0.05
) -> Dict[str, Any]:
    """
    Perform Two One-Sided Tests (TOST) procedure for bioequivalence.

    Tests:
    - H01: muT/muR <= theta_lower (lower bound)
    - H02: muT/muR >= theta_upper (upper bound)

    BE is concluded if both null hypotheses are rejected.

    Args:
        test_values: Test formulation values
        reference_values: Reference formulation values
        theta_lower: Lower equivalence bound (default: 0.80)
        theta_upper: Upper equivalence bound (default: 1.25)
        alpha: Significance level (default: 0.05)

    Returns:
        Dict with keys:
        - parameter: Parameter analyzed
        - t_lower: T-statistic for lower bound test
        - t_upper: T-statistic for upper bound test
        - p_lower: P-value for lower bound test
        - p_upper: P-value for upper bound test
        - reject_lower: Whether lower bound H0 is rejected
        - reject_upper: Whether upper bound H0 is rejected
        - conclusion: "bioequivalent" or "not_bioequivalent"

    Example:
        >>> result = tost_analysis(test_auc, reference_auc)
        >>> print(f"Conclusion: {result['conclusion']}")
        >>> print(f"p-values: lower={result['p_lower']:.4f}, upper={result['p_upper']:.4f}")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, test_values)
    r = _to_julia_float_vector(jl, reference_values)
    result = jl.OpenPKPDCore.tost_analysis(
        t, r,
        theta_lower=theta_lower,
        theta_upper=theta_upper,
        alpha=alpha
    )

    return {
        "parameter": str(result.parameter),
        "t_lower": float(result.t_lower),
        "t_upper": float(result.t_upper),
        "p_lower": float(result.p_lower),
        "p_upper": float(result.p_upper),
        "reject_lower": bool(result.reject_lower),
        "reject_upper": bool(result.reject_upper),
        "conclusion": str(result.be_conclusion),
    }


# ============================================================================
# BE Conclusion
# ============================================================================

def be_conclusion(
    ci_lower: float,
    ci_upper: float,
    theta_lower: float = 0.80,
    theta_upper: float = 1.25
) -> str:
    """
    Determine bioequivalence conclusion from confidence interval.

    Args:
        ci_lower: Lower bound of 90% CI
        ci_upper: Upper bound of 90% CI
        theta_lower: Lower equivalence bound (default: 0.80)
        theta_upper: Upper equivalence bound (default: 1.25)

    Returns:
        "bioequivalent", "not_bioequivalent", or "inconclusive"

    Example:
        >>> result = bioequivalence_90ci(test, reference)
        >>> conclusion = be_conclusion(result['ci_lower'], result['ci_upper'])
        >>> print(f"BE conclusion: {conclusion}")
    """
    jl = _require_julia()
    result = jl.OpenPKPDCore.be_conclusion(
        ci_lower, ci_upper,
        theta_lower=theta_lower,
        theta_upper=theta_upper
    )
    return str(result)


# ============================================================================
# Complete BE Analysis
# ============================================================================

@dataclass
class BioequivalenceResult:
    """
    Complete bioequivalence analysis result.

    Attributes:
        parameter: Parameter analyzed (e.g., "cmax", "auc")
        n_test: Number of test subjects
        n_reference: Number of reference subjects
        gmr: Geometric mean ratio
        ci_lower: Lower bound of 90% CI
        ci_upper: Upper bound of 90% CI
        cv_intra: Intra-subject CV (%)
        conclusion: BE conclusion
        be_limits: BE acceptance limits used
    """
    parameter: str
    n_test: int
    n_reference: int
    gmr: float
    ci_lower: float
    ci_upper: float
    cv_intra: float
    conclusion: str
    be_limits: Tuple[float, float]


def run_bioequivalence(
    parameter: str,
    test_values: List[float],
    reference_values: List[float],
    be_limits: Tuple[float, float] = (0.80, 1.25)
) -> BioequivalenceResult:
    """
    Run complete bioequivalence analysis for a parameter.

    Args:
        parameter: Parameter name (e.g., "cmax", "auc")
        test_values: Test formulation values
        reference_values: Reference formulation values
        be_limits: BE acceptance limits (default: (0.80, 1.25))

    Returns:
        BioequivalenceResult with complete analysis

    Example:
        >>> result = run_bioequivalence("cmax", test_cmax, ref_cmax)
        >>> print(f"Parameter: {result.parameter}")
        >>> print(f"GMR: {result.gmr:.4f} ({result.ci_lower:.4f}, {result.ci_upper:.4f})")
        >>> print(f"Conclusion: {result.conclusion}")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, test_values)
    r = _to_julia_float_vector(jl, reference_values)
    param_sym = jl.Symbol(parameter)

    result = jl.OpenPKPDCore.create_be_result(
        param_sym, t, r,
        be_limits=be_limits
    )

    return BioequivalenceResult(
        parameter=str(result.parameter),
        n_test=int(result.n_test),
        n_reference=int(result.n_reference),
        gmr=float(result.gmr),
        ci_lower=float(result.ci_lower),
        ci_upper=float(result.ci_upper),
        cv_intra=float(result.cv_intra),
        conclusion=str(result.be_conclusion),
        be_limits=be_limits,
    )
