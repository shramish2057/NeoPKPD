using Pkg
Pkg.activate("core/OpenPKPDCore")
Pkg.instantiate()

using OpenPKPDCore
using JSON

function write(path, artifact)
    open(path, "w") do io
        JSON.print(io, artifact)
    end
end

function gen_pk_iv_bolus()
    pk = ModelSpec(
        OneCompIVBolus(),
        "golden_pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [
            DoseEvent(0.0, 60.0),
            DoseEvent(0.0, 40.0),   # duplicate time to lock semantics summing
            DoseEvent(12.0, 25.0),
        ],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(pk, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res)
end

function gen_pk_oral()
    pk = ModelSpec(
        OneCompOralFirstOrder(),
        "golden_pk_oral",
        OneCompOralFirstOrderParams(1.2, 5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 50.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(pk, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res)
end

function gen_pk_then_pd_direct_emax()
    pk = ModelSpec(
        OneCompIVBolus(),
        "golden_pk_then_pd",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    pd = PDSpec(
        DirectEmax(),
        "golden_emax",
        DirectEmaxParams(10.0, 40.0, 0.8),
        :conc,
        :effect,
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_pkpd(pk, pd, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res, pd_spec=pd)
end

function gen_coupled_pkpd_turnover_oral()
    pk = ModelSpec(
        OneCompOralFirstOrder(),
        "golden_coupled_oral",
        OneCompOralFirstOrderParams(1.2, 5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 50.0)],
    )

    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseTurnover(),
        "golden_turnover",
        IndirectResponseTurnoverParams(Kin, Kout, Rss, 0.8, 0.5),
        :conc,
        :response,
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res, pd_spec=pd)
end

function gen_population_pk_iv()
    base = ModelSpec(
        OneCompIVBolus(),
        "golden_pop_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    iiv = IIVSpec(LogNormalIIV(), Dict(:CL => 0.2, :V => 0.1), UInt64(7777), 5)
    pop = PopulationSpec(base, iiv, nothing, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_population(pop, grid, solver)

    return serialize_population_execution(population_spec=pop, grid=grid, solver=solver, result=res)
end

function gen_sensitivity_single_iv()
    spec = ModelSpec(
        OneCompIVBolus(),
        "golden_sens_single",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    plan = PerturbationPlan("CL_up_10pct", [Perturbation(RelativePerturbation(), :CL, 0.10)])

    res = run_sensitivity(spec, grid, solver; plan = plan, observation = :conc)

    return serialize_sensitivity_execution(model_spec = spec, grid = grid, solver = solver, result = res)
end

function gen_sensitivity_population_iv()
    base = ModelSpec(
        OneCompIVBolus(),
        "golden_sens_pop_base",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    iiv = IIVSpec(LogNormalIIV(), Dict(:CL => 0.2, :V => 0.1), UInt64(7777), 5)
    pop = PopulationSpec(base, iiv, nothing, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    plan = PerturbationPlan("CL_up_10pct", [Perturbation(RelativePerturbation(), :CL, 0.10)])

    res = run_population_sensitivity(pop, grid, solver; plan = plan, observation = :conc, probs = [0.05, 0.95])

    return serialize_population_sensitivity_execution(population_spec = pop, grid = grid, solver = solver, result = res)
end

function gen_population_iov_iv()
    base = ModelSpec(
        OneCompIVBolus(),
        "golden_iov_pop_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 100.0)],
    )

    iov = IOVSpec(LogNormalIIV(), Dict(:CL => 0.3), UInt64(1234), OccasionDefinition(:dose_times))
    pop = PopulationSpec(base, nothing, iov, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_population(pop, grid, solver)

    return serialize_population_execution(population_spec=pop, grid=grid, solver=solver, result=res)
end

function gen_population_iov_pkpd()
    pk = ModelSpec(
        OneCompIVBolus(),
        "golden_iov_pkpd",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 100.0)],
    )

    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseTurnover(),
        "golden_turnover",
        IndirectResponseTurnoverParams(Kin, Kout, Rss, 0.8, 0.5),
        :conc,
        :response,
    )

    iov = IOVSpec(LogNormalIIV(), Dict(:CL => 0.3), UInt64(2222), OccasionDefinition(:dose_times))
    pop = PopulationSpec(pk, nothing, iov, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_population(pop, grid, solver; pd_spec=pd)

    return serialize_population_execution(
        population_spec=pop,
        grid=grid,
        solver=solver,
        result=res,
        pd_spec=pd,
    )
end


function gen_population_time_varying_covariate_iv()
    base = ModelSpec(
        OneCompIVBolus(),
        "golden_tv_cov_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    cm = CovariateModel(
        "cl_tv",
        [CovariateEffect(LinearCovariate(), :CL, :CLMULT, 1.0, 1.0)],
    )

    tv = TimeVaryingCovariates(Dict(
        :CLMULT => TimeCovariateSeries(StepTimeCovariate(), [0.0, 10.0], [1.0, 2.0]),
    ))

    covs = [IndividualCovariates(Dict{Symbol,Float64}(), tv)]
    pop = PopulationSpec(base, nothing, nothing, cm, covs)

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_population(pop, grid, solver)

    return serialize_population_execution(population_spec=pop, grid=grid, solver=solver, result=res)
end

# =====================================================
# New PK Models
# =====================================================

function gen_pk_twocomp_iv()
    pk = ModelSpec(
        TwoCompIVBolus(),
        "golden_pk_twocomp_iv",
        TwoCompIVBolusParams(10.0, 50.0, 5.0, 100.0),  # CL=10, V1=50, Q=5, V2=100
        [DoseEvent(0.0, 500.0), DoseEvent(12.0, 300.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:0.5:48.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(pk, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res)
end

function gen_pk_twocomp_oral()
    pk = ModelSpec(
        TwoCompOral(),
        "golden_pk_twocomp_oral",
        TwoCompOralParams(1.0, 10.0, 50.0, 5.0, 100.0),  # Ka=1, CL=10, V1=50, Q=5, V2=100
        [DoseEvent(0.0, 500.0), DoseEvent(12.0, 300.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:0.5:48.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(pk, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res)
end

function gen_pk_threecomp_iv()
    pk = ModelSpec(
        ThreeCompIVBolus(),
        "golden_pk_threecomp_iv",
        ThreeCompIVBolusParams(10.0, 50.0, 10.0, 80.0, 2.0, 200.0),
        [DoseEvent(0.0, 1000.0)],
    )

    grid = SimGrid(0.0, 72.0, collect(0.0:1.0:72.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(pk, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res)
end

function gen_pk_transit_absorption()
    pk = ModelSpec(
        TransitAbsorption(),
        "golden_pk_transit",
        TransitAbsorptionParams(5, 2.0, 1.0, 10.0, 50.0),  # N=5, Ktr=2, Ka=1, CL=10, V=50
        [DoseEvent(0.0, 500.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(pk, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res)
end

function gen_pk_michaelis_menten()
    pk = ModelSpec(
        MichaelisMentenElimination(),
        "golden_pk_mm",
        MichaelisMentenEliminationParams(100.0, 5.0, 50.0),  # Vmax=100, Km=5, V=50
        [DoseEvent(0.0, 500.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:0.5:48.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(pk, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res)
end

# =====================================================
# New PD Models
# =====================================================

function gen_pkpd_sigmoid_emax()
    pk = ModelSpec(
        OneCompIVBolus(),
        "golden_pk_for_semax",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    pd = PDSpec(
        SigmoidEmax(),
        "golden_semax",
        SigmoidEmaxParams(10.0, 40.0, 0.8, 2.0),  # E0=10, Emax=40, EC50=0.8, gamma=2
        :conc,
        :effect,
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_pkpd(pk, pd, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res, pd_spec=pd)
end

function gen_pkpd_biophase_equilibration()
    pk = ModelSpec(
        OneCompIVBolus(),
        "golden_pk_for_biophase",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    pd = PDSpec(
        BiophaseEquilibration(),
        "golden_biophase",
        BiophaseEquilibrationParams(0.5, 10.0, 40.0, 0.8),  # ke0=0.5, E0=10, Emax=40, EC50=0.8
        :conc,
        :effect,
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_pkpd(pk, pd, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res, pd_spec=pd)
end

# =====================================================
# Integration: New PK + PD combinations
# =====================================================

function gen_twocomp_with_sigmoid_emax()
    pk = ModelSpec(
        TwoCompIVBolus(),
        "golden_twocomp_semax",
        TwoCompIVBolusParams(10.0, 50.0, 5.0, 100.0),
        [DoseEvent(0.0, 500.0)],
    )

    pd = PDSpec(
        SigmoidEmax(),
        "golden_semax_2c",
        SigmoidEmaxParams(0.0, 100.0, 2.0, 2.0),
        :conc,
        :effect,
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:1.0:48.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_pkpd(pk, pd, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res, pd_spec=pd)
end

function gen_transit_with_turnover()
    pk = ModelSpec(
        TransitAbsorption(),
        "golden_transit_turnover",
        TransitAbsorptionParams(5, 2.0, 1.0, 10.0, 50.0),
        [DoseEvent(0.0, 500.0)],
    )

    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseTurnover(),
        "golden_turnover_transit",
        IndirectResponseTurnoverParams(Kin, Kout, Rss, 0.8, 1.0),
        :conc,
        :response,
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res, pd_spec=pd)
end

function main()
    mkpath("validation/golden")

    artifacts = Dict(
        # Original goldens
        "pk_iv_bolus.json" => gen_pk_iv_bolus(),
        "pk_oral.json" => gen_pk_oral(),
        "pk_then_pd_direct_emax.json" => gen_pk_then_pd_direct_emax(),
        "pkpd_coupled_turnover_oral.json" => gen_coupled_pkpd_turnover_oral(),
        "population_pk_iv.json" => gen_population_pk_iv(),
        "sensitivity_single_iv.json" => gen_sensitivity_single_iv(),
        "sensitivity_population_iv.json" => gen_sensitivity_population_iv(),
        "population_iov_iv.json" => gen_population_iov_iv(),
        "population_iov_pkpd.json" => gen_population_iov_pkpd(),
        "population_time_varying_cov_iv.json" => gen_population_time_varying_covariate_iv(),

        # New PK model goldens
        "pk_twocomp_iv.json" => gen_pk_twocomp_iv(),
        "pk_twocomp_oral.json" => gen_pk_twocomp_oral(),
        "pk_threecomp_iv.json" => gen_pk_threecomp_iv(),
        "pk_transit_absorption.json" => gen_pk_transit_absorption(),
        "pk_michaelis_menten.json" => gen_pk_michaelis_menten(),

        # New PD model goldens
        "pkpd_sigmoid_emax.json" => gen_pkpd_sigmoid_emax(),
        "pkpd_biophase_equilibration.json" => gen_pkpd_biophase_equilibration(),

        # Integration goldens
        "pkpd_twocomp_sigmoid_emax.json" => gen_twocomp_with_sigmoid_emax(),
        "pkpd_transit_turnover.json" => gen_transit_with_turnover(),
    )

    for (fname, art) in artifacts
        path = joinpath("validation/golden", fname)
        write(path, art)
        println("Wrote: $(path)")
    end
end

main()
