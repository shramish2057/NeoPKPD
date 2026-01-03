"""
PK-PD Simulations

This module provides simulation functions for coupled PK-PD models:
- Direct Emax (DirectEmax)
- Sigmoid Emax / Hill equation (SigmoidEmax)
- Indirect response turnover (IndirectResponseTurnover)
- Biophase equilibration / effect compartment (BiophaseEquilibration)
"""

from typing import Any, Dict, List, Optional

from .._core import _require_julia, _simresult_to_py, _to_julia_vector


def simulate_pkpd_direct_emax(
    cl: float,
    v: float,
    doses: List[Dict[str, float]],
    e0: float,
    emax: float,
    ec50: float,
    t0: float,
    t1: float,
    saveat: List[float],
    pk_kind: str = "OneCompIVBolus",
    ka: Optional[float] = None,
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a coupled PK-PD simulation with direct Emax effect model.

    PD equation:
        Effect = E0 + (Emax * C) / (EC50 + C)

    Args:
        cl: Clearance (volume/time)
        v: Volume of distribution
        doses: List of dose events
        e0: Baseline effect
        emax: Maximum effect
        ec50: Concentration at 50% of Emax
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        pk_kind: PK model type ("OneCompIVBolus" or "OneCompOralFirstOrder")
        ka: Absorption rate constant (required if pk_kind is oral)
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys:
        - t: List of time points
        - states: Dict of state variables
        - observations: Dict of observables (conc, effect)
        - metadata: Dict of run metadata

    Example:
        >>> result = openpkpd.simulate_pkpd_direct_emax(
        ...     cl=1.0, v=10.0,
        ...     doses=[{"time": 0.0, "amount": 100.0}],
        ...     e0=0.0, emax=100.0, ec50=5.0,
        ...     t0=0.0, t1=24.0,
        ...     saveat=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
        ... )
        >>> print(f"Max effect: {max(result['observations']['effect'])}")
    """
    jl = _require_julia()

    DoseEvent = jl.OpenPKPDCore.DoseEvent
    ModelSpec = jl.OpenPKPDCore.ModelSpec
    SimGrid = jl.OpenPKPDCore.SimGrid
    SolverSpec = jl.OpenPKPDCore.SolverSpec
    PDSpec = jl.OpenPKPDCore.PDSpec
    DirectEmax = jl.OpenPKPDCore.DirectEmax
    DirectEmaxParams = jl.OpenPKPDCore.DirectEmaxParams

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    if pk_kind == "OneCompIVBolus":
        OneCompIVBolus = jl.OpenPKPDCore.OneCompIVBolus
        OneCompIVBolusParams = jl.OpenPKPDCore.OneCompIVBolusParams
        pk_spec = ModelSpec(OneCompIVBolus(), "py_pkpd_pk", OneCompIVBolusParams(float(cl), float(v)), doses_vec)
    elif pk_kind == "OneCompOralFirstOrder":
        if ka is None:
            raise ValueError("ka required for OneCompOralFirstOrder PK model")
        OneCompOralFirstOrder = jl.OpenPKPDCore.OneCompOralFirstOrder
        OneCompOralFirstOrderParams = jl.OpenPKPDCore.OneCompOralFirstOrderParams
        pk_spec = ModelSpec(
            OneCompOralFirstOrder(), "py_pkpd_pk",
            OneCompOralFirstOrderParams(float(ka), float(cl), float(v)), doses_vec
        )
    else:
        raise ValueError(f"Unsupported pk_kind: {pk_kind}")

    pd_spec = PDSpec(
        DirectEmax(), "py_pkpd_pd",
        DirectEmaxParams(float(e0), float(emax), float(ec50)),
        jl.Symbol("conc"), jl.Symbol("effect")
    )

    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.OpenPKPDCore.simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)
    return _simresult_to_py(res)


def simulate_pkpd_indirect_response(
    cl: float,
    v: float,
    doses: List[Dict[str, float]],
    kin: float,
    kout: float,
    r0: float,
    imax: float,
    ic50: float,
    t0: float,
    t1: float,
    saveat: List[float],
    pk_kind: str = "OneCompIVBolus",
    ka: Optional[float] = None,
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a coupled PK-PD simulation with indirect response turnover model.

    PD equations:
        I(C) = (Imax * C) / (IC50 + C)
        dR/dt = Kin - Kout * (1 - I(C)) * R

    Args:
        cl: Clearance (volume/time)
        v: Volume of distribution
        doses: List of dose events
        kin: Zero-order production rate
        kout: First-order elimination rate
        r0: Baseline response (should equal Kin/Kout for steady state)
        imax: Maximum inhibition (0 to 1)
        ic50: Concentration at 50% of Imax
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        pk_kind: PK model type ("OneCompIVBolus" or "OneCompOralFirstOrder")
        ka: Absorption rate constant (required if pk_kind is oral)
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys:
        - t: List of time points
        - states: Dict of state variables (including R for response)
        - observations: Dict of observables (conc, response)
        - metadata: Dict of run metadata

    Example:
        >>> result = openpkpd.simulate_pkpd_indirect_response(
        ...     cl=1.0, v=10.0,
        ...     doses=[{"time": 0.0, "amount": 100.0}],
        ...     kin=10.0, kout=0.1, r0=100.0, imax=0.9, ic50=5.0,
        ...     t0=0.0, t1=120.0,
        ...     saveat=list(range(121))
        ... )
        >>> print(f"Min response: {min(result['observations']['response'])}")
    """
    jl = _require_julia()

    DoseEvent = jl.OpenPKPDCore.DoseEvent
    ModelSpec = jl.OpenPKPDCore.ModelSpec
    SimGrid = jl.OpenPKPDCore.SimGrid
    SolverSpec = jl.OpenPKPDCore.SolverSpec
    PDSpec = jl.OpenPKPDCore.PDSpec
    IndirectResponseTurnover = jl.OpenPKPDCore.IndirectResponseTurnover
    IndirectResponseTurnoverParams = jl.OpenPKPDCore.IndirectResponseTurnoverParams

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    if pk_kind == "OneCompIVBolus":
        OneCompIVBolus = jl.OpenPKPDCore.OneCompIVBolus
        OneCompIVBolusParams = jl.OpenPKPDCore.OneCompIVBolusParams
        pk_spec = ModelSpec(OneCompIVBolus(), "py_pkpd_pk", OneCompIVBolusParams(float(cl), float(v)), doses_vec)
    elif pk_kind == "OneCompOralFirstOrder":
        if ka is None:
            raise ValueError("ka required for OneCompOralFirstOrder PK model")
        OneCompOralFirstOrder = jl.OpenPKPDCore.OneCompOralFirstOrder
        OneCompOralFirstOrderParams = jl.OpenPKPDCore.OneCompOralFirstOrderParams
        pk_spec = ModelSpec(
            OneCompOralFirstOrder(), "py_pkpd_pk",
            OneCompOralFirstOrderParams(float(ka), float(cl), float(v)), doses_vec
        )
    else:
        raise ValueError(f"Unsupported pk_kind: {pk_kind}")

    pd_spec = PDSpec(
        IndirectResponseTurnover(), "py_pkpd_pd",
        IndirectResponseTurnoverParams(float(kin), float(kout), float(r0), float(imax), float(ic50)),
        jl.Symbol("conc"), jl.Symbol("response")
    )

    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.OpenPKPDCore.simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)
    return _simresult_to_py(res)


def simulate_pkpd_sigmoid_emax(
    cl: float,
    v: float,
    doses: List[Dict[str, float]],
    e0: float,
    emax: float,
    ec50: float,
    gamma: float,
    t0: float,
    t1: float,
    saveat: List[float],
    pk_kind: str = "OneCompIVBolus",
    ka: Optional[float] = None,
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a coupled PK-PD simulation with Sigmoid Emax (Hill equation) effect model.

    PD equation (Hill equation):
        Effect = E0 + (Emax * C^gamma) / (EC50^gamma + C^gamma)

    The gamma (Hill coefficient) controls the steepness of the response:
    - gamma = 1: Standard hyperbolic Emax model
    - gamma > 1: Steeper, more "switch-like" response (threshold effect)
    - gamma < 1: More gradual response

    Common gamma values: 0.5-5 for most drugs, 3-6 for neuromuscular blockers

    Args:
        cl: Clearance (volume/time)
        v: Volume of distribution
        doses: List of dose events
        e0: Baseline effect
        emax: Maximum effect (can be negative for inhibitory effects)
        ec50: Concentration at 50% of Emax
        gamma: Hill coefficient (steepness parameter)
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        pk_kind: PK model type ("OneCompIVBolus" or "OneCompOralFirstOrder")
        ka: Absorption rate constant (required if pk_kind is oral)
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys:
        - t: List of time points
        - states: Dict of state variables
        - observations: Dict of observables (conc, effect)
        - metadata: Dict of run metadata

    Example:
        >>> result = openpkpd.simulate_pkpd_sigmoid_emax(
        ...     cl=1.0, v=10.0,
        ...     doses=[{"time": 0.0, "amount": 100.0}],
        ...     e0=0.0, emax=100.0, ec50=5.0, gamma=2.0,
        ...     t0=0.0, t1=24.0,
        ...     saveat=list(range(25))
        ... )
        >>> # Effect curve will be steeper than standard Emax
    """
    jl = _require_julia()

    DoseEvent = jl.OpenPKPDCore.DoseEvent
    ModelSpec = jl.OpenPKPDCore.ModelSpec
    SimGrid = jl.OpenPKPDCore.SimGrid
    SolverSpec = jl.OpenPKPDCore.SolverSpec
    PDSpec = jl.OpenPKPDCore.PDSpec
    SigmoidEmax = jl.OpenPKPDCore.SigmoidEmax
    SigmoidEmaxParams = jl.OpenPKPDCore.SigmoidEmaxParams

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    if pk_kind == "OneCompIVBolus":
        OneCompIVBolus = jl.OpenPKPDCore.OneCompIVBolus
        OneCompIVBolusParams = jl.OpenPKPDCore.OneCompIVBolusParams
        pk_spec = ModelSpec(OneCompIVBolus(), "py_pkpd_pk", OneCompIVBolusParams(float(cl), float(v)), doses_vec)
    elif pk_kind == "OneCompOralFirstOrder":
        if ka is None:
            raise ValueError("ka required for OneCompOralFirstOrder PK model")
        OneCompOralFirstOrder = jl.OpenPKPDCore.OneCompOralFirstOrder
        OneCompOralFirstOrderParams = jl.OpenPKPDCore.OneCompOralFirstOrderParams
        pk_spec = ModelSpec(
            OneCompOralFirstOrder(), "py_pkpd_pk",
            OneCompOralFirstOrderParams(float(ka), float(cl), float(v)), doses_vec
        )
    else:
        raise ValueError(f"Unsupported pk_kind: {pk_kind}")

    pd_spec = PDSpec(
        SigmoidEmax(), "py_pkpd_pd",
        SigmoidEmaxParams(float(e0), float(emax), float(ec50), float(gamma)),
        jl.Symbol("conc"), jl.Symbol("effect")
    )

    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.OpenPKPDCore.simulate_pkpd(pk_spec, pd_spec, grid, solver)
    return _simresult_to_py(res)


def simulate_pkpd_biophase_equilibration(
    cl: float,
    v: float,
    doses: List[Dict[str, float]],
    ke0: float,
    e0: float,
    emax: float,
    ec50: float,
    t0: float,
    t1: float,
    saveat: List[float],
    pk_kind: str = "OneCompIVBolus",
    ka: Optional[float] = None,
    alg: str = "Tsit5",
    reltol: float = 1e-10,
    abstol: float = 1e-12,
    maxiters: int = 10**7,
) -> Dict[str, Any]:
    """
    Run a coupled PK-PD simulation with biophase equilibration (effect compartment) model.

    This model introduces a hypothetical effect site compartment to account for
    temporal delay between plasma concentration changes and pharmacodynamic effects.

    PD equations:
        dCe/dt = ke0 * (Cp - Ce)  # Effect compartment equilibration
        Effect = E0 + (Emax * Ce) / (EC50 + Ce)  # Emax effect from effect site

    The equilibration half-life t1/2,ke0 = ln(2)/ke0 indicates:
    - Small ke0 (long t1/2): Slow equilibration, significant hysteresis
    - Large ke0 (short t1/2): Fast equilibration, near-direct effect

    Common applications: Anesthetics, CNS-active drugs, neuromuscular blockers

    Note: This function uses quasi-steady state approximation (Ce = Cp).
    For true effect compartment dynamics with hysteresis, use the full ODE approach.

    Args:
        cl: Clearance (volume/time)
        v: Volume of distribution
        doses: List of dose events
        ke0: Effect site equilibration rate constant (1/time)
        e0: Baseline effect
        emax: Maximum effect
        ec50: Effect site concentration at 50% of Emax
        t0: Simulation start time
        t1: Simulation end time
        saveat: List of time points for output
        pk_kind: PK model type ("OneCompIVBolus" or "OneCompOralFirstOrder")
        ka: Absorption rate constant (required if pk_kind is oral)
        alg: ODE solver algorithm (default: "Tsit5")
        reltol: Relative tolerance (default: 1e-10)
        abstol: Absolute tolerance (default: 1e-12)
        maxiters: Maximum solver iterations (default: 10^7)

    Returns:
        Dict with keys:
        - t: List of time points
        - states: Dict of state variables
        - observations: Dict of observables (conc, effect)
        - metadata: Dict of run metadata

    Example:
        >>> result = openpkpd.simulate_pkpd_biophase_equilibration(
        ...     cl=1.0, v=10.0,
        ...     doses=[{"time": 0.0, "amount": 100.0}],
        ...     ke0=0.5, e0=0.0, emax=100.0, ec50=5.0,
        ...     t0=0.0, t1=24.0,
        ...     saveat=list(range(25))
        ... )
        >>> # Effect will lag behind concentration due to equilibration delay
    """
    jl = _require_julia()

    DoseEvent = jl.OpenPKPDCore.DoseEvent
    ModelSpec = jl.OpenPKPDCore.ModelSpec
    SimGrid = jl.OpenPKPDCore.SimGrid
    SolverSpec = jl.OpenPKPDCore.SolverSpec
    PDSpec = jl.OpenPKPDCore.PDSpec
    BiophaseEquilibration = jl.OpenPKPDCore.BiophaseEquilibration
    BiophaseEquilibrationParams = jl.OpenPKPDCore.BiophaseEquilibrationParams

    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    if pk_kind == "OneCompIVBolus":
        OneCompIVBolus = jl.OpenPKPDCore.OneCompIVBolus
        OneCompIVBolusParams = jl.OpenPKPDCore.OneCompIVBolusParams
        pk_spec = ModelSpec(OneCompIVBolus(), "py_pkpd_pk", OneCompIVBolusParams(float(cl), float(v)), doses_vec)
    elif pk_kind == "OneCompOralFirstOrder":
        if ka is None:
            raise ValueError("ka required for OneCompOralFirstOrder PK model")
        OneCompOralFirstOrder = jl.OpenPKPDCore.OneCompOralFirstOrder
        OneCompOralFirstOrderParams = jl.OpenPKPDCore.OneCompOralFirstOrderParams
        pk_spec = ModelSpec(
            OneCompOralFirstOrder(), "py_pkpd_pk",
            OneCompOralFirstOrderParams(float(ka), float(cl), float(v)), doses_vec
        )
    else:
        raise ValueError(f"Unsupported pk_kind: {pk_kind}")

    pd_spec = PDSpec(
        BiophaseEquilibration(), "py_pkpd_pd",
        BiophaseEquilibrationParams(float(ke0), float(e0), float(emax), float(ec50)),
        jl.Symbol("conc"), jl.Symbol("effect")
    )

    grid = SimGrid(float(t0), float(t1), [float(x) for x in saveat])
    solver = SolverSpec(jl.Symbol(alg), float(reltol), float(abstol), int(maxiters))

    res = jl.OpenPKPDCore.simulate_pkpd(pk_spec, pd_spec, grid, solver)
    return _simresult_to_py(res)
