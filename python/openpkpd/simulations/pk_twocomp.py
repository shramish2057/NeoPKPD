"""
Two-Compartment PK Simulations

This module provides simulation functions for two-compartment PK models:
- IV bolus (TwoCompIVBolus)
- Oral first-order absorption (TwoCompOral)
"""

from typing import Any, Dict, List

from .._core import _require_julia, _simresult_to_py, _to_julia_vector


def simulate_pk_twocomp_iv_bolus(
    cl: float,
    v1: float,
    q: float,
    v2: float,
    doses: List[Dict[str, float]],
    t0: float,
    t1: float,
    saveat: List[float],
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a two-compartment IV bolus PK simulation.

    Model structure:
    - Central compartment: receives IV bolus, connected to peripheral
    - Peripheral compartment: equilibrates with central
    - Elimination from central compartment

    Model equations:
        dA_central/dt = -k10*A_central - k12*A_central + k21*A_peripheral
        dA_peripheral/dt = k12*A_central - k21*A_peripheral
        C = A_central / V1

    Micro-rate constants:
        k10 = CL/V1 (elimination)
        k12 = Q/V1 (central to peripheral)
        k21 = Q/V2 (peripheral to central)

    Args:
        cl: Clearance from central compartment (volume/time)
        v1: Volume of central compartment
        q: Inter-compartmental clearance
        v2: Volume of peripheral compartment
        doses: List of dose events, each a dict with 'time' and 'amount'
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys:
        - t: List of time points
        - states: Dict of state variables (A_central, A_peripheral)
        - observations: Dict of observables (conc)
        - metadata: Dict of run metadata

    Example:
        >>> result = openpkpd.simulate_pk_twocomp_iv_bolus(
        ...     cl=10.0, v1=50.0, q=5.0, v2=100.0,
        ...     doses=[{"time": 0.0, "amount": 500.0}],
        ...     t0=0.0, t1=48.0,
        ...     saveat=list(range(49))
        ... )
        >>> print(f"Cmax: {max(result['observations']['conc'])}")
    """
    jl = _require_julia()

    DoseEvent = jl.OpenPKPDCore.DoseEvent
    ModelSpec = jl.OpenPKPDCore.ModelSpec
    TwoCompIVBolus = jl.OpenPKPDCore.TwoCompIVBolus
    TwoCompIVBolusParams = jl.OpenPKPDCore.TwoCompIVBolusParams
    SimGrid = jl.OpenPKPDCore.SimGrid
    SolverSpec = jl.OpenPKPDCore.SolverSpec

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    spec = ModelSpec(
        TwoCompIVBolus(), "py_twocomp_iv",
        TwoCompIVBolusParams(float(cl), float(v1), float(q), float(v2)),
        doses_vec
    )
    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.OpenPKPDCore.simulate(spec, grid, solver)
    return _simresult_to_py(res)


def simulate_pk_twocomp_oral(
    ka: float,
    cl: float,
    v1: float,
    q: float,
    v2: float,
    doses: List[Dict[str, float]],
    t0: float,
    t1: float,
    saveat: List[float],
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a two-compartment oral first-order absorption PK simulation.

    Model structure:
    - Gut compartment: first-order absorption into central
    - Central compartment: connected to peripheral, elimination
    - Peripheral compartment: equilibrates with central

    Args:
        ka: Absorption rate constant (1/time)
        cl: Clearance from central compartment (volume/time)
        v1: Volume of central compartment
        q: Inter-compartmental clearance
        v2: Volume of peripheral compartment
        doses: List of dose events
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys:
        - t: List of time points
        - states: Dict of state variables (A_gut, A_central, A_peripheral)
        - observations: Dict of observables (conc)
        - metadata: Dict of run metadata

    Example:
        >>> result = openpkpd.simulate_pk_twocomp_oral(
        ...     ka=1.0, cl=10.0, v1=50.0, q=5.0, v2=100.0,
        ...     doses=[{"time": 0.0, "amount": 500.0}],
        ...     t0=0.0, t1=48.0,
        ...     saveat=list(range(49))
        ... )
    """
    jl = _require_julia()

    DoseEvent = jl.OpenPKPDCore.DoseEvent
    ModelSpec = jl.OpenPKPDCore.ModelSpec
    TwoCompOral = jl.OpenPKPDCore.TwoCompOral
    TwoCompOralParams = jl.OpenPKPDCore.TwoCompOralParams
    SimGrid = jl.OpenPKPDCore.SimGrid
    SolverSpec = jl.OpenPKPDCore.SolverSpec

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    spec = ModelSpec(
        TwoCompOral(), "py_twocomp_oral",
        TwoCompOralParams(float(ka), float(cl), float(v1), float(q), float(v2)),
        doses_vec
    )
    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.OpenPKPDCore.simulate(spec, grid, solver)
    return _simresult_to_py(res)
