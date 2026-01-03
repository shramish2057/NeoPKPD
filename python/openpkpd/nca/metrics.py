"""
OpenPKPD NCA Metrics - Individual NCA metric functions.

This module provides Python bindings to the Julia NCA metric functions.
All functions call the corresponding Julia implementations for FDA/EMA
compliant calculations.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

from .._core import _require_julia, _to_julia_float_vector


# ============================================================================
# Primary Exposure Metrics
# ============================================================================

def nca_cmax(concentrations: List[float]) -> float:
    """
    Find maximum observed concentration (Cmax).

    Args:
        concentrations: List of concentration values

    Returns:
        Maximum concentration

    Example:
        >>> cmax = nca_cmax([0.0, 1.2, 2.0, 1.5, 0.8])
        >>> print(f"Cmax: {cmax}")
    """
    jl = _require_julia()
    c = _to_julia_float_vector(jl, concentrations)
    return float(jl.OpenPKPDCore.nca_cmax(c))


def nca_tmax(times: List[float], concentrations: List[float]) -> float:
    """
    Find time of maximum observed concentration (Tmax).

    Args:
        times: List of time points
        concentrations: List of concentration values

    Returns:
        Time of maximum concentration

    Example:
        >>> tmax = nca_tmax([0.0, 1.0, 2.0, 4.0], [0.0, 1.5, 2.0, 1.0])
        >>> print(f"Tmax: {tmax}")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    return float(jl.OpenPKPDCore.nca_tmax(t, c))


def nca_cmin(concentrations: List[float]) -> float:
    """
    Find minimum observed concentration (Cmin).

    Args:
        concentrations: List of concentration values

    Returns:
        Minimum concentration

    Example:
        >>> cmin = nca_cmin([2.5, 4.2, 3.8, 2.5])
        >>> print(f"Cmin: {cmin}")
    """
    jl = _require_julia()
    c = _to_julia_float_vector(jl, concentrations)
    return float(jl.OpenPKPDCore.nca_cmin(c))


def nca_clast(
    times: List[float],
    concentrations: List[float],
    lloq: float = 0.0
) -> float:
    """
    Find last measurable concentration (Clast).

    Args:
        times: List of time points
        concentrations: List of concentration values
        lloq: Lower limit of quantification

    Returns:
        Last measurable concentration
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    return float(jl.OpenPKPDCore.nca_clast(t, c, lloq=lloq))


def nca_tlast(
    times: List[float],
    concentrations: List[float],
    lloq: float = 0.0
) -> float:
    """
    Find time of last measurable concentration (Tlast).

    Args:
        times: List of time points
        concentrations: List of concentration values
        lloq: Lower limit of quantification

    Returns:
        Time of last measurable concentration
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    return float(jl.OpenPKPDCore.nca_tlast(t, c, lloq=lloq))


def nca_cavg(
    times: List[float],
    concentrations: List[float],
    tau: float,
    method: str = "lin_log_mixed"
) -> float:
    """
    Calculate average concentration over a dosing interval (Cavg).

    Cavg = AUC0-tau / tau

    Args:
        times: List of time points
        concentrations: List of concentration values
        tau: Dosing interval
        method: AUC calculation method ("linear", "log_linear", "lin_log_mixed")

    Returns:
        Average concentration over dosing interval
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    config = _make_nca_config(jl, method=method)
    return float(jl.OpenPKPDCore.nca_cavg(t, c, tau, config))


# ============================================================================
# Lambda-z Estimation
# ============================================================================

def estimate_lambda_z(
    times: List[float],
    concentrations: List[float],
    min_points: int = 3,
    r2_threshold: float = 0.9,
    method: str = "lin_log_mixed"
) -> Dict[str, Any]:
    """
    Estimate terminal elimination rate constant (lambda_z).

    Uses log-linear regression on the terminal phase to estimate lambda_z.

    Args:
        times: List of time points
        concentrations: List of concentration values
        min_points: Minimum points for regression (default: 3)
        r2_threshold: Minimum R-squared threshold (default: 0.9)
        method: AUC calculation method

    Returns:
        Dict with keys:
        - lambda_z: Terminal elimination rate constant (or None)
        - t_half: Terminal half-life (or None)
        - r_squared: R-squared of regression (or None)
        - n_points: Number of points used
        - quality_flag: Quality assessment ("good", "warning", "insufficient")
        - warnings: List of quality warnings

    Example:
        >>> result = estimate_lambda_z([0, 1, 2, 4, 8, 12, 24], [0, 2, 1.8, 1.2, 0.6, 0.3, 0.075])
        >>> print(f"Lambda_z: {result['lambda_z']}")
        >>> print(f"t1/2: {result['t_half']}")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    config = _make_nca_config(jl, method=method, lambda_z_min_points=min_points,
                              lambda_z_r2_threshold=r2_threshold)

    result = jl.OpenPKPDCore.estimate_lambda_z(t, c, config)

    return {
        "lambda_z": _maybe_float(result.lambda_z),
        "t_half": _maybe_float(result.t_half),
        "r_squared": _maybe_float(result.r_squared),
        "adjusted_r_squared": _maybe_float(result.adjusted_r_squared),
        "intercept": _maybe_float(result.intercept),
        "n_points": int(result.n_points),
        "start_time": float(result.start_time),
        "end_time": float(result.end_time),
        "points_used": list(result.points_used),
        "quality_flag": str(result.quality_flag),
        "warnings": list(result.warnings),
    }


def nca_half_life(lambda_z: float) -> float:
    """
    Calculate terminal half-life from lambda_z.

    t1/2 = ln(2) / lambda_z

    Args:
        lambda_z: Terminal elimination rate constant

    Returns:
        Terminal half-life

    Example:
        >>> t_half = nca_half_life(0.1)
        >>> print(f"Half-life: {t_half:.2f} hours")
    """
    jl = _require_julia()
    return float(jl.OpenPKPDCore.nca_half_life(lambda_z))


# ============================================================================
# AUC Calculations
# ============================================================================

def auc_0_t(
    times: List[float],
    concentrations: List[float],
    method: str = "lin_log_mixed"
) -> float:
    """
    Calculate AUC from time 0 to last measurable concentration (AUC0-t).

    Args:
        times: List of time points
        concentrations: List of concentration values
        method: Calculation method ("linear", "log_linear", "lin_log_mixed")

    Returns:
        AUC from 0 to last measurable concentration

    Example:
        >>> auc = auc_0_t([0, 1, 2, 4, 8], [0, 2, 1.5, 1.0, 0.5])
        >>> print(f"AUC0-t: {auc:.2f}")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    config = _make_nca_config(jl, method=method)
    return float(jl.OpenPKPDCore.auc_0_t(t, c, config))


def auc_0_inf(
    times: List[float],
    concentrations: List[float],
    lambda_z: float,
    clast: float,
    method: str = "lin_log_mixed"
) -> Tuple[float, float]:
    """
    Calculate AUC extrapolated to infinity (AUC0-inf).

    AUC0-inf = AUC0-t + Clast/lambda_z

    Args:
        times: List of time points
        concentrations: List of concentration values
        lambda_z: Terminal elimination rate constant
        clast: Last measurable concentration
        method: Calculation method

    Returns:
        Tuple of (AUC0-inf, extrapolation percentage)

    Example:
        >>> auc_inf, extra_pct = auc_0_inf(t, c, lambda_z=0.1, clast=0.5)
        >>> print(f"AUC0-inf: {auc_inf:.2f}, Extrapolation: {extra_pct:.1f}%")
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    config = _make_nca_config(jl, method=method)
    result = jl.OpenPKPDCore.auc_0_inf(t, c, lambda_z, clast, config)
    return (float(result[0]), float(result[1]))


def auc_0_tau(
    times: List[float],
    concentrations: List[float],
    tau: float,
    method: str = "lin_log_mixed"
) -> float:
    """
    Calculate AUC over a dosing interval (AUC0-tau).

    Args:
        times: List of time points
        concentrations: List of concentration values
        tau: Dosing interval
        method: Calculation method

    Returns:
        AUC over the dosing interval
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    config = _make_nca_config(jl, method=method)
    return float(jl.OpenPKPDCore.auc_0_tau(t, c, tau, config))


def auc_partial(
    times: List[float],
    concentrations: List[float],
    t_start: float,
    t_end: float,
    method: str = "lin_log_mixed"
) -> float:
    """
    Calculate partial AUC between two time points.

    Args:
        times: List of time points
        concentrations: List of concentration values
        t_start: Start time
        t_end: End time
        method: Calculation method

    Returns:
        Partial AUC
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    config = _make_nca_config(jl, method=method)
    return float(jl.OpenPKPDCore.auc_partial(t, c, t_start, t_end, config))


def aumc_0_t(
    times: List[float],
    concentrations: List[float],
    method: str = "lin_log_mixed"
) -> float:
    """
    Calculate AUMC from time 0 to last measurable concentration.

    AUMC = Area Under the (first) Moment Curve = integral of t*C(t)

    Args:
        times: List of time points
        concentrations: List of concentration values
        method: Calculation method

    Returns:
        AUMC from 0 to last measurable concentration
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    config = _make_nca_config(jl, method=method)
    return float(jl.OpenPKPDCore.aumc_0_t(t, c, config))


def aumc_0_inf(
    times: List[float],
    concentrations: List[float],
    lambda_z: float,
    clast: float,
    tlast: float,
    method: str = "lin_log_mixed"
) -> float:
    """
    Calculate AUMC extrapolated to infinity.

    Args:
        times: List of time points
        concentrations: List of concentration values
        lambda_z: Terminal elimination rate constant
        clast: Last measurable concentration
        tlast: Time of last measurable concentration
        method: Calculation method

    Returns:
        AUMC extrapolated to infinity
    """
    jl = _require_julia()
    t = _to_julia_float_vector(jl, times)
    c = _to_julia_float_vector(jl, concentrations)
    config = _make_nca_config(jl, method=method)
    return float(jl.OpenPKPDCore.aumc_0_inf(t, c, lambda_z, clast, tlast, config))


# ============================================================================
# PK Parameters
# ============================================================================

def nca_mrt(
    aumc_0_inf: float,
    auc_0_inf: float,
    route: str = "extravascular",
    t_inf: float = 0.0
) -> float:
    """
    Calculate Mean Residence Time (MRT).

    For extravascular: MRT = AUMC0-inf / AUC0-inf
    For IV infusion: MRT = AUMC0-inf / AUC0-inf - Tinf/2

    Args:
        aumc_0_inf: AUMC extrapolated to infinity
        auc_0_inf: AUC extrapolated to infinity
        route: Administration route ("extravascular", "iv_bolus", "iv_infusion")
        t_inf: Infusion duration (for iv_infusion)

    Returns:
        Mean residence time
    """
    jl = _require_julia()
    route_sym = jl.Symbol(route)
    return float(jl.OpenPKPDCore.nca_mrt(aumc_0_inf, auc_0_inf, route=route_sym, t_inf=t_inf))


def nca_cl_f(dose: float, auc_0_inf: float) -> float:
    """
    Calculate apparent clearance (CL/F).

    CL/F = Dose / AUC0-inf

    Args:
        dose: Administered dose
        auc_0_inf: AUC extrapolated to infinity

    Returns:
        Apparent clearance
    """
    jl = _require_julia()
    return float(jl.OpenPKPDCore.nca_cl_f(dose, auc_0_inf))


def nca_vz_f(dose: float, lambda_z: float, auc_0_inf: float) -> float:
    """
    Calculate apparent volume of distribution (Vz/F).

    Vz/F = Dose / (lambda_z * AUC0-inf)

    Args:
        dose: Administered dose
        lambda_z: Terminal elimination rate constant
        auc_0_inf: AUC extrapolated to infinity

    Returns:
        Apparent volume of distribution
    """
    jl = _require_julia()
    return float(jl.OpenPKPDCore.nca_vz_f(dose, lambda_z, auc_0_inf))


def nca_vss(cl: float, mrt: float) -> float:
    """
    Calculate volume of distribution at steady state (Vss).

    Vss = CL * MRT

    Args:
        cl: Clearance
        mrt: Mean residence time

    Returns:
        Volume at steady state
    """
    jl = _require_julia()
    return float(jl.OpenPKPDCore.nca_vss(cl, mrt))


def nca_bioavailability(
    auc_test: float,
    dose_test: float,
    auc_reference: float,
    dose_reference: float
) -> float:
    """
    Calculate relative bioavailability (F).

    F = (AUC_test / Dose_test) / (AUC_ref / Dose_ref)

    Args:
        auc_test: AUC of test formulation
        dose_test: Dose of test formulation
        auc_reference: AUC of reference formulation
        dose_reference: Dose of reference formulation

    Returns:
        Relative bioavailability
    """
    jl = _require_julia()
    return float(jl.OpenPKPDCore.nca_bioavailability(auc_test, dose_test, auc_reference, dose_reference))


# ============================================================================
# Multiple Dose Metrics
# ============================================================================

def nca_accumulation_index(auc_ss: float, auc_sd: float) -> float:
    """
    Calculate accumulation index (Rac).

    Rac = AUC_ss / AUC_sd

    Args:
        auc_ss: AUC0-tau at steady state
        auc_sd: AUC0-inf from single dose

    Returns:
        Accumulation index
    """
    jl = _require_julia()
    return float(jl.OpenPKPDCore.nca_accumulation_index(auc_ss, auc_sd))


def nca_ptf(cmax: float, cmin: float, cavg: float) -> float:
    """
    Calculate Peak-Trough Fluctuation (PTF) percentage.

    PTF = 100 * (Cmax - Cmin) / Cavg

    Args:
        cmax: Maximum concentration in dosing interval
        cmin: Minimum (trough) concentration
        cavg: Average concentration over dosing interval

    Returns:
        Peak-trough fluctuation (%)
    """
    jl = _require_julia()
    return float(jl.OpenPKPDCore.nca_ptf(cmax, cmin, cavg))


def nca_swing(cmax: float, cmin: float) -> float:
    """
    Calculate Swing percentage.

    Swing = 100 * (Cmax - Cmin) / Cmin

    Args:
        cmax: Maximum concentration in dosing interval
        cmin: Minimum (trough) concentration

    Returns:
        Swing (%)
    """
    jl = _require_julia()
    return float(jl.OpenPKPDCore.nca_swing(cmax, cmin))


def nca_linearity_index(doses: List[float], aucs: List[float]) -> Dict[str, Any]:
    """
    Assess dose proportionality using power model.

    AUC = alpha * Dose^beta
    Linear if beta is approximately 1.0.

    Args:
        doses: List of dose levels
        aucs: List of corresponding AUC values

    Returns:
        Dict with keys: beta, r_squared, is_linear
    """
    jl = _require_julia()
    d = _to_julia_float_vector(jl, doses)
    a = _to_julia_float_vector(jl, aucs)
    result = jl.OpenPKPDCore.nca_linearity_index(d, a)
    return {
        "beta": float(result.beta),
        "r_squared": float(result.r_squared),
        "is_linear": bool(result.is_linear),
    }


def nca_time_to_steady_state(lambda_z: float, fraction: float = 0.90) -> float:
    """
    Estimate time to reach a fraction of steady state.

    t_ss = -ln(1 - fraction) / lambda_z

    Args:
        lambda_z: Terminal elimination rate constant
        fraction: Fraction of steady state (default: 0.90)

    Returns:
        Time to reach specified fraction of steady state
    """
    jl = _require_julia()
    return float(jl.OpenPKPDCore.nca_time_to_steady_state(lambda_z, fraction=fraction))


# ============================================================================
# Helper Functions
# ============================================================================

def _make_nca_config(
    jl: Any,
    method: str = "lin_log_mixed",
    lambda_z_min_points: int = 3,
    lambda_z_r2_threshold: float = 0.9,
    extrapolation_max_pct: float = 20.0,
    lloq: Optional[float] = None,
) -> Any:
    """Create Julia NCAConfig object."""
    # Select method type
    if method == "linear":
        method_obj = jl.OpenPKPDCore.LinearMethod()
    elif method == "log_linear":
        method_obj = jl.OpenPKPDCore.LogLinearMethod()
    else:
        method_obj = jl.OpenPKPDCore.LinLogMixedMethod()

    if lloq is not None:
        return jl.OpenPKPDCore.NCAConfig(
            method=method_obj,
            lambda_z_min_points=lambda_z_min_points,
            lambda_z_r2_threshold=lambda_z_r2_threshold,
            extrapolation_max_pct=extrapolation_max_pct,
            lloq=lloq,
        )
    else:
        return jl.OpenPKPDCore.NCAConfig(
            method=method_obj,
            lambda_z_min_points=lambda_z_min_points,
            lambda_z_r2_threshold=lambda_z_r2_threshold,
            extrapolation_max_pct=extrapolation_max_pct,
        )


def _maybe_float(val: Any) -> Optional[float]:
    """Convert Julia value to float, returning None for Julia nothing."""
    if val is None:
        return None
    try:
        return float(val)
    except (TypeError, ValueError):
        return None
