"""
Three-Compartment PK Simulations

This module provides simulation functions for three-compartment PK models:
- IV bolus (ThreeCompIVBolus)
"""

from typing import Any, Dict, List

from .._core import _require_julia, _simresult_to_py, _to_julia_vector


def simulate_pk_threecomp_iv_bolus(
    cl: float,
    v1: float,
    q2: float,
    v2: float,
    q3: float,
    v3: float,
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
    Run a three-compartment IV bolus PK simulation (mammillary model).

    Model structure:
    - Central compartment: receives IV bolus, connected to two peripheral compartments
    - Shallow peripheral (periph1): rapid equilibration with central
    - Deep peripheral (periph2): slow equilibration with central
    - Elimination from central compartment

    The concentration-time profile shows tri-exponential decay:
    - Alpha phase: rapid initial decline (distribution to shallow peripheral)
    - Beta phase: intermediate decline (distribution to deep peripheral)
    - Gamma phase: terminal elimination phase

    Args:
        cl: Clearance from central compartment (volume/time)
        v1: Volume of central compartment
        q2: Inter-compartmental clearance to shallow peripheral
        v2: Volume of shallow peripheral compartment
        q3: Inter-compartmental clearance to deep peripheral
        v3: Volume of deep peripheral compartment
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
        - states: Dict of state variables (A_central, A_periph1, A_periph2)
        - observations: Dict of observables (conc)
        - metadata: Dict of run metadata

    Example:
        >>> result = openpkpd.simulate_pk_threecomp_iv_bolus(
        ...     cl=10.0, v1=50.0, q2=10.0, v2=80.0, q3=2.0, v3=200.0,
        ...     doses=[{"time": 0.0, "amount": 1000.0}],
        ...     t0=0.0, t1=72.0,
        ...     saveat=list(range(73))
        ... )
    """
    jl = _require_julia()

    DoseEvent = jl.OpenPKPDCore.DoseEvent
    ModelSpec = jl.OpenPKPDCore.ModelSpec
    ThreeCompIVBolus = jl.OpenPKPDCore.ThreeCompIVBolus
    ThreeCompIVBolusParams = jl.OpenPKPDCore.ThreeCompIVBolusParams
    SimGrid = jl.OpenPKPDCore.SimGrid
    SolverSpec = jl.OpenPKPDCore.SolverSpec

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    spec = ModelSpec(
        ThreeCompIVBolus(), "py_threecomp_iv",
        ThreeCompIVBolusParams(float(cl), float(v1), float(q2), float(v2), float(q3), float(v3)),
        doses_vec
    )
    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.OpenPKPDCore.simulate(spec, grid, solver)
    return _simresult_to_py(res)
