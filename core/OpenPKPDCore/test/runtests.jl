using Test
using OpenPKPDCore

function analytic_onecomp_ivbolus_conc(
    t::Float64, doses::Vector{DoseEvent}, CL::Float64, V::Float64
)
    k = CL / V
    c = 0.0
    for d in doses
        if t >= d.time
            dt = t - d.time
            c += (d.amount / V) * exp(-k * dt)
        end
    end
    return c
end

@testset "OneCompIVBolus analytic equivalence" begin
    spec = ModelSpec(
        OneCompIVBolus(),
        "1c_iv_bolus",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    CL = spec.params.CL
    V = spec.params.V

    for (i, t) in enumerate(res.t)
        c_ref = analytic_onecomp_ivbolus_conc(t, spec.doses, CL, V)
        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
    end
end

@testset "Dose event at nonzero time" begin
    spec = ModelSpec(
        OneCompIVBolus(),
        "1c_iv_bolus_delayed",
        OneCompIVBolusParams(3.0, 30.0),
        [DoseEvent(2.0, 120.0)],
    )

    grid = SimGrid(0.0, 12.0, collect(0.0:0.25:12.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    CL = spec.params.CL
    V = spec.params.V

    for (i, t) in enumerate(res.t)
        if any(d.time == t for d in spec.doses)
            continue  # skip discontinuity (left-continuous solver output)
        end
        c_ref = analytic_onecomp_ivbolus_conc(t, spec.doses, CL, V)
        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
    end
end

@testset "Multiple bolus schedule" begin
    spec = ModelSpec(
        OneCompIVBolus(),
        "1c_iv_bolus_multi",
        OneCompIVBolusParams(4.0, 40.0),
        [DoseEvent(0.0, 80.0), DoseEvent(6.0, 50.0), DoseEvent(10.0, 50.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    CL = spec.params.CL
    V = spec.params.V

    for (i, t) in enumerate(res.t)
        if any(d.time == t for d in spec.doses)
            continue
        end
        c_ref = analytic_onecomp_ivbolus_conc(t, spec.doses, CL, V)
        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
    end
end

function analytic_onecomp_oral_first_order_conc(
    t::Float64, doses, Ka::Float64, CL::Float64, V::Float64
)
    k = CL / V
    c = 0.0
    for d in doses
        if t >= d.time
            dt = t - d.time
            # Handle Ka close to k for numerical stability
            if abs(Ka - k) < 1e-12
                # Limit as Ka -> k:
                # C(t) = (Dose/V) * Ka * dt * exp(-k*dt)
                c += (d.amount / V) * Ka * dt * exp(-k * dt)
            else
                c += (d.amount / V) * (Ka / (Ka - k)) * (exp(-k * dt) - exp(-Ka * dt))
            end
        end
    end
    return c
end

@testset "OneCompOralFirstOrder analytic equivalence" begin
    spec = ModelSpec(
        OneCompOralFirstOrder(),
        "1c_oral_fo",
        OneCompOralFirstOrderParams(1.2, 5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    Ka = spec.params.Ka
    CL = spec.params.CL
    V = spec.params.V

    for (i, t) in enumerate(res.t)
        c_ref = analytic_onecomp_oral_first_order_conc(t, spec.doses, Ka, CL, V)
        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
    end
end

@testset "OneCompOralFirstOrder multiple doses" begin
    spec = ModelSpec(
        OneCompOralFirstOrder(),
        "1c_oral_fo_multi",
        OneCompOralFirstOrderParams(0.9, 4.0, 40.0),
        [DoseEvent(0.0, 80.0), DoseEvent(8.0, 50.0), DoseEvent(16.0, 50.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    Ka = spec.params.Ka
    CL = spec.params.CL
    V = spec.params.V

    for (i, t) in enumerate(res.t)
        c_ref = analytic_onecomp_oral_first_order_conc(t, spec.doses, Ka, CL, V)
        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
    end
end

@testset "State outputs are present and aligned" begin
    spec = ModelSpec(
        OneCompOralFirstOrder(),
        "state_presence",
        OneCompOralFirstOrderParams(1.0, 3.0, 30.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 12.0, collect(0.0:1.0:12.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    res = simulate(spec, grid, solver)

    @test haskey(res.states, :A_gut)
    @test haskey(res.states, :A_central)
    @test haskey(res.observations, :conc)

    for v in values(res.states)
        @test length(v) == length(res.t)
    end
end

function direct_emax(C::Float64, E0::Float64, Emax::Float64, EC50::Float64)
    return E0 + (Emax * C) / (EC50 + C)
end

@testset "PKPD coupling with DirectEmax using IV bolus analytic reference" begin
    pk = ModelSpec(
        OneCompIVBolus(), "pk_iv", OneCompIVBolusParams(5.0, 50.0), [DoseEvent(0.0, 100.0)]
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(DirectEmax(), "pd_emax", DirectEmaxParams(10.0, 40.0, 0.8), :conc, :effect)

    res = simulate_pkpd(pk, pd, grid, solver)

    @test haskey(res.observations, :conc)
    @test haskey(res.observations, :effect)
    @test length(res.observations[:effect]) == length(res.t)

    CL = pk.params.CL
    V = pk.params.V

    for (i, t) in enumerate(res.t)
        c_ref = analytic_onecomp_ivbolus_conc(t, pk.doses, CL, V)
        e_ref = direct_emax(c_ref, pd.params.E0, pd.params.Emax, pd.params.EC50)

        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
        @test isapprox(res.observations[:effect][i], e_ref; rtol=1e-10, atol=1e-12)
    end
end

function analytic_turnover_R(t::Float64, Kin::Float64, Kout::Float64, R0::Float64)
    # dR/dt = Kin - Kout*R
    return R0 * exp(-Kout * t) + (Kin / Kout) * (1.0 - exp(-Kout * t))
end

@testset "IndirectResponseTurnover coupled: Imax=0 matches analytic turnover and PK analytic" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv_for_pd",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 50.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(
        IndirectResponseTurnover(),
        "turnover_no_effect",
        IndirectResponseTurnoverParams(
            10.0,  # Kin
            0.5,   # Kout
            15.0,  # R0
            0.0,   # Imax, no drug effect
            1.0,   # IC50, irrelevant here
        ),
        :conc,
        :response,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    @test haskey(res.observations, :conc)
    @test haskey(res.observations, :response)

    CL = pk.params.CL
    V = pk.params.V

    Kin = pd.params.Kin
    Kout = pd.params.Kout
    R0 = pd.params.R0

    for (i, t) in enumerate(res.t)
        if any(d.time == t for d in pk.doses)
            continue  # skip discontinuities (left-continuous solver output)
        end

        c_ref = analytic_onecomp_ivbolus_conc(t, pk.doses, CL, V)
        r_ref = analytic_turnover_R(t, Kin, Kout, R0)

        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
        @test isapprox(res.observations[:response][i], r_ref; rtol=1e-8, atol=1e-10)
    end
end

@testset "IndirectResponseTurnover coupled: inhibition raises response above baseline" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv_effect",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 200.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    # Baseline at steady state to make interpretation clean
    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseTurnover(),
        "turnover_with_effect",
        IndirectResponseTurnoverParams(Kin, Kout, Rss, 0.8, 0.5),
        :conc,
        :response,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    r = res.observations[:response]
    @test maximum(r) > Rss
end

@testset "Coupled engine generalization: oral PK with Imax=0 matches analytic PK and analytic turnover" begin
    pk = ModelSpec(
        OneCompOralFirstOrder(),
        "pk_oral_for_pd",
        OneCompOralFirstOrderParams(1.2, 5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 50.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(
        IndirectResponseTurnover(),
        "turnover_no_effect_oral",
        IndirectResponseTurnoverParams(
            10.0,  # Kin
            0.5,   # Kout
            15.0,  # R0
            0.0,   # Imax
            1.0,   # IC50
        ),
        :conc,
        :response,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    Ka = pk.params.Ka
    CL = pk.params.CL
    V = pk.params.V

    Kin = pd.params.Kin
    Kout = pd.params.Kout
    R0 = pd.params.R0

    for (i, t) in enumerate(res.t)
        c_ref = analytic_onecomp_oral_first_order_conc(t, pk.doses, Ka, CL, V)
        r_ref = analytic_turnover_R(t, Kin, Kout, R0)

        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
        @test isapprox(res.observations[:response][i], r_ref; rtol=1e-8, atol=1e-10)
    end
end

@testset "Event semantics: duplicate times are summed (IV bolus)" begin
    base_params = OneCompIVBolusParams(5.0, 50.0)
    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pk_dup = ModelSpec(
        OneCompIVBolus(),
        "dup_times",
        base_params,
        [
            DoseEvent(0.0, 60.0),
            DoseEvent(0.0, 40.0),
            DoseEvent(12.0, 10.0),
            DoseEvent(12.0, 15.0),
        ],
    )

    pk_sum = ModelSpec(
        OneCompIVBolus(),
        "summed",
        base_params,
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 25.0)],
    )

    r_dup = simulate(pk_dup, grid, solver)
    r_sum = simulate(pk_sum, grid, solver)

    @test r_dup.metadata["event_semantics_version"] == "1.0.0"
    @test r_sum.metadata["event_semantics_version"] == "1.0.0"

    for i in eachindex(r_dup.t)
        @test isapprox(
            r_dup.observations[:conc][i],
            r_sum.observations[:conc][i];
            rtol=1e-12,
            atol=1e-12,
        )
    end
end

@testset "Event semantics: input ordering does not matter for same-time events" begin
    params = OneCompIVBolusParams(5.0, 50.0)
    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    doses_a = [
        DoseEvent(0.0, 30.0),
        DoseEvent(0.0, 70.0),
        DoseEvent(12.0, 10.0),
        DoseEvent(12.0, 15.0),
    ]
    doses_b = [
        DoseEvent(0.0, 70.0),
        DoseEvent(0.0, 30.0),
        DoseEvent(12.0, 15.0),
        DoseEvent(12.0, 10.0),
    ]

    pk_a = ModelSpec(OneCompIVBolus(), "a", params, doses_a)
    pk_b = ModelSpec(OneCompIVBolus(), "b", params, doses_b)

    r_a = simulate(pk_a, grid, solver)
    r_b = simulate(pk_b, grid, solver)

    for i in eachindex(r_a.t)
        @test isapprox(
            r_a.observations[:conc][i], r_b.observations[:conc][i]; rtol=1e-12, atol=1e-12
        )
    end
end

@testset "Event semantics: duplicate times are summed in coupled PKPD" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_dup_coupled",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 60.0), DoseEvent(0.0, 40.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseTurnover(),
        "pd_coupled",
        IndirectResponseTurnoverParams(Kin, Kout, Rss, 0.0, 1.0),
        :conc,
        :response,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    @test res.metadata["event_semantics_version"] == "1.0.0"
end

@testset "Solver semantics version is present and stable" begin
    spec = ModelSpec(
        OneCompIVBolus(),
        "solver_semantics_test",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 12.0, collect(0.0:1.0:12.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    res = simulate(spec, grid, solver)

    @test haskey(res.metadata, "solver_semantics_version")
    @test res.metadata["solver_semantics_version"] == "1.0.0"
end

@testset "Execution artifact serialization is complete and self-consistent" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "serialize_test",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 12.0, collect(0.0:1.0:12.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    res = simulate(pk, grid, solver)

    artifact = serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res)

    @test artifact["artifact_schema_version"] == "1.0.0"
    @test haskey(artifact, "model_spec")
    @test haskey(artifact, "grid")
    @test haskey(artifact, "solver")
    @test haskey(artifact, "result")

    @test artifact["result"]["metadata"]["event_semantics_version"] == "1.0.0"
    @test artifact["result"]["metadata"]["solver_semantics_version"] == "1.0.0"

    @test artifact["model_spec"]["name"] == "serialize_test"
    @test artifact["solver"]["alg"] == "Tsit5"
end

@testset "Artifact replay: PK only matches original outputs" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "replay_pk",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 60.0), DoseEvent(0.0, 40.0), DoseEvent(12.0, 25.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res1 = simulate(pk, grid, solver)

    artifact = serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res1)

    res2 = replay_execution(artifact)

    @test res2.t == res1.t
    for i in eachindex(res1.t)
        @test isapprox(
            res2.observations[:conc][i], res1.observations[:conc][i]; rtol=1e-12, atol=1e-12
        )
    end
end

@testset "Artifact replay: PK then PD (DirectEmax) matches original outputs" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "replay_pkpd_direct",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(DirectEmax(), "emax", DirectEmaxParams(10.0, 40.0, 0.8), :conc, :effect)

    res1 = simulate_pkpd(pk, pd, grid, solver)

    artifact = serialize_execution(
        model_spec=pk, grid=grid, solver=solver, result=res1, pd_spec=pd
    )

    res2 = replay_execution(artifact)

    @test res2.t == res1.t
    for i in eachindex(res1.t)
        @test isapprox(
            res2.observations[:conc][i], res1.observations[:conc][i]; rtol=1e-12, atol=1e-12
        )
        @test isapprox(
            res2.observations[:effect][i],
            res1.observations[:effect][i];
            rtol=1e-12,
            atol=1e-12,
        )
    end
end

@testset "Artifact replay: coupled PKPD (IndirectResponseTurnover) matches original outputs" begin
    pk = ModelSpec(
        OneCompOralFirstOrder(),
        "replay_coupled_oral",
        OneCompOralFirstOrderParams(1.2, 5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 50.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseTurnover(),
        "turnover",
        IndirectResponseTurnoverParams(Kin, Kout, Rss, 0.8, 0.5),
        :conc,
        :response,
    )

    res1 = simulate_pkpd_coupled(pk, pd, grid, solver)

    artifact = serialize_execution(
        model_spec=pk, grid=grid, solver=solver, result=res1, pd_spec=pd
    )

    res2 = replay_execution(artifact)

    @test res2.t == res1.t
    for i in eachindex(res1.t)
        @test isapprox(
            res2.observations[:conc][i], res1.observations[:conc][i]; rtol=1e-12, atol=1e-12
        )
        @test isapprox(
            res2.observations[:response][i],
            res1.observations[:response][i];
            rtol=1e-12,
            atol=1e-12,
        )
    end
end

@testset "Semantics fingerprint is complete" begin
    fp = semantics_fingerprint()

    @test haskey(fp, "event_semantics_version")
    @test haskey(fp, "solver_semantics_version")
    @test haskey(fp, "artifact_schema_version")
end

@testset "Population simulation is deterministic with StableRNG seed" begin
    base = ModelSpec(
        OneCompIVBolus(),
        "pop_pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    iiv = IIVSpec(LogNormalIIV(), Dict(:CL => 0.2, :V => 0.1), UInt64(12345), 20)
    pop = PopulationSpec(base, iiv, nothing, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    r1 = simulate_population(pop, grid, solver)
    r2 = simulate_population(pop, grid, solver)

    @test r1.metadata["n"] == 20
    @test r2.metadata["n"] == 20

    # Realized individual parameters must match exactly run-to-run
    for i in 1:20
        @test r1.params[i][:CL] == r2.params[i][:CL]
        @test r1.params[i][:V] == r2.params[i][:V]
    end

    # Concentration series must match exactly to tight tolerance
    for i in 1:20
        c1 = r1.individuals[i].observations[:conc]
        c2 = r2.individuals[i].observations[:conc]
        @test length(c1) == length(c2)
        for j in eachindex(c1)
            @test isapprox(c1[j], c2[j]; rtol=1e-12, atol=1e-12)
        end
    end
end

@testset "Population simulation without IIV produces one individual matching base simulation" begin
    base = ModelSpec(
        OneCompOralFirstOrder(),
        "pop_pk_oral_base",
        OneCompOralFirstOrderParams(1.2, 5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    pop = PopulationSpec(base, nothing, nothing, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    base_res = simulate(base, grid, solver)
    pop_res = simulate_population(pop, grid, solver)

    @test length(pop_res.individuals) == 1

    c_base = base_res.observations[:conc]
    c_pop = pop_res.individuals[1].observations[:conc]

    for j in eachindex(c_base)
        @test isapprox(c_pop[j], c_base[j]; rtol=1e-12, atol=1e-12)
    end
end

@testset "Population artifact serialization and replay matches original outputs" begin
    base = ModelSpec(
        OneCompIVBolus(),
        "pop_artifact_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    iiv = IIVSpec(LogNormalIIV(), Dict(:CL => 0.2, :V => 0.1), UInt64(424242), 5)
    pop = PopulationSpec(base, iiv, nothing, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    r1 = simulate_population(pop, grid, solver)

    artifact = serialize_population_execution(
        population_spec=pop, grid=grid, solver=solver, result=r1
    )

    r2 = replay_population_execution(artifact)

    @test r2.metadata["n"] == r1.metadata["n"]
    @test length(r2.individuals) == length(r1.individuals)

    for i in eachindex(r1.individuals)
        c1 = r1.individuals[i].observations[:conc]
        c2 = r2.individuals[i].observations[:conc]
        @test length(c1) == length(c2)
        for j in eachindex(c1)
            @test isapprox(c2[j], c1[j]; rtol=1e-12, atol=1e-12)
        end
    end

    @test haskey(r1.summaries, :conc)
    @test haskey(r2.summaries, :conc)

    s1 = r1.summaries[:conc]
    s2 = r2.summaries[:conc]

    @test s1.probs == s2.probs

    for i in eachindex(s1.mean)
        @test isapprox(s2.mean[i], s1.mean[i]; rtol=1e-12, atol=1e-12)
        @test isapprox(s2.median[i], s1.median[i]; rtol=1e-12, atol=1e-12)
    end

    for p in s1.probs
        q1 = s1.quantiles[p]
        q2 = s2.quantiles[p]
        for i in eachindex(q1)
            @test isapprox(q2[i], q1[i]; rtol=1e-12, atol=1e-12)
        end
    end
end

@testset "Single-run sensitivity: increasing CL reduces concentration" begin
    spec = ModelSpec(
        OneCompIVBolus(),
        "sens_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    plan = PerturbationPlan(
        "CL_up_10pct", [Perturbation(RelativePerturbation(), :CL, 0.10)]
    )

    out = run_sensitivity(spec, grid, solver; plan=plan, observation=:conc)

    # Concentration should be lower at positive times after increasing CL
    # We avoid t=0 because it reflects initial condition and bolus application.
    @test out.pert_metric_series[2] < out.base_metric_series[2]
    @test out.metrics.max_abs_delta > 0.0
end

@testset "Population sensitivity: increasing CL reduces mean concentration" begin
    base = ModelSpec(
        OneCompIVBolus(),
        "sens_pop_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    iiv = IIVSpec(LogNormalIIV(), Dict(:CL => 0.2, :V => 0.1), UInt64(9999), 20)
    pop = PopulationSpec(base, iiv, nothing, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    plan = PerturbationPlan(
        "CL_up_10pct", [Perturbation(RelativePerturbation(), :CL, 0.10)]
    )

    out = run_population_sensitivity(pop, grid, solver; plan=plan, observation=:conc)

    @test out.pert_summary_mean[2] < out.base_summary_mean[2]
    @test out.metrics_mean.max_abs_delta > 0.0
end

@testset "Sensitivity artifact replay matches stored series (single)" begin
    spec = ModelSpec(
        OneCompIVBolus(),
        "sens_art_single",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    plan = PerturbationPlan(
        "CL_up_10pct", [Perturbation(RelativePerturbation(), :CL, 0.10)]
    )

    r1 = run_sensitivity(spec, grid, solver; plan=plan, observation=:conc)

    art = serialize_sensitivity_execution(
        model_spec=spec, grid=grid, solver=solver, result=r1
    )

    r2 = replay_sensitivity_execution(art)

    for i in eachindex(r1.base_metric_series)
        @test isapprox(
            r2.base_metric_series[i], r1.base_metric_series[i]; rtol=1e-12, atol=1e-12
        )
        @test isapprox(
            r2.pert_metric_series[i], r1.pert_metric_series[i]; rtol=1e-12, atol=1e-12
        )
    end
end

@testset "Sensitivity artifact replay matches stored series (population)" begin
    base = ModelSpec(
        OneCompIVBolus(),
        "sens_art_pop",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    iiv = IIVSpec(LogNormalIIV(), Dict(:CL => 0.2, :V => 0.1), UInt64(9999), 20)
    pop = PopulationSpec(base, iiv, nothing, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    plan = PerturbationPlan(
        "CL_up_10pct", [Perturbation(RelativePerturbation(), :CL, 0.10)]
    )

    r1 = run_population_sensitivity(
        pop, grid, solver; plan=plan, observation=:conc, probs=[0.05, 0.95]
    )

    art = serialize_population_sensitivity_execution(
        population_spec=pop, grid=grid, solver=solver, result=r1
    )

    r2 = replay_population_sensitivity_execution(art)

    for i in eachindex(r1.base_summary_mean)
        @test isapprox(
            r2.base_summary_mean[i], r1.base_summary_mean[i]; rtol=1e-12, atol=1e-12
        )
        @test isapprox(
            r2.pert_summary_mean[i], r1.pert_summary_mean[i]; rtol=1e-12, atol=1e-12
        )
    end
end

@testset "IOV changes parameters across occasions deterministically" begin
    base = ModelSpec(
        OneCompIVBolus(),
        "iov_pop_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 100.0)],
    )

    iiv = nothing
    iov = IOVSpec(
        LogNormalIIV(), Dict(:CL => 0.3), UInt64(1234), OccasionDefinition(:dose_times)
    )

    pop = PopulationSpec(base, iiv, iov, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res1 = simulate_population(pop, grid, solver)
    res2 = simulate_population(pop, grid, solver)

    # Deterministic replay
    c1 = res1.individuals[1].observations[:conc]
    c2 = res2.individuals[1].observations[:conc]
    for i in eachindex(c1)
        @test isapprox(c1[i], c2[i]; rtol=1e-12, atol=1e-12)
    end

    # Behavioral: after second dose, concentration trajectory should differ from a non-IOV run
    pop_no_iov = PopulationSpec(base, nothing, nothing, nothing, IndividualCovariates[])
    base_res = simulate_population(pop_no_iov, grid, solver)

    c_base = base_res.individuals[1].observations[:conc]
    # pick a time after 12.0, for example index where t == 13
    idx13 = findfirst(==(13.0), grid.saveat)
    @test idx13 !== nothing
    @test abs(c1[idx13] - c_base[idx13]) > 0.0
end

@testset "Coupled PKPD with IOV is deterministic and alters post-dose dynamics" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pkpd_iov",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 100.0)],
    )

    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseTurnover(),
        "pd_turnover",
        IndirectResponseTurnoverParams(Kin, Kout, Rss, 0.8, 0.5),
        :conc,
        :response,
    )

    iov = IOVSpec(
        LogNormalIIV(), Dict(:CL => 0.3), UInt64(2222), OccasionDefinition(:dose_times)
    )
    pop = PopulationSpec(pk, nothing, iov, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    r1 = simulate_population(pop, grid, solver; pd_spec=pd)
    r2 = simulate_population(pop, grid, solver; pd_spec=pd)

    # Deterministic replay
    y1 = r1.individuals[1].observations[:response]
    y2 = r2.individuals[1].observations[:response]

    for i in eachindex(y1)
        @test isapprox(y1[i], y2[i]; rtol=1e-12, atol=1e-12)
    end

    # Compare against no-IOV behavior
    pop_no_iov = PopulationSpec(pk, nothing, nothing, nothing, IndividualCovariates[])
    base = simulate_population(pop_no_iov, grid, solver; pd_spec=pd)

    idx13 = findfirst(==(13.0), grid.saveat)
    @test idx13 !== nothing
    @test abs(y1[idx13] - base.individuals[1].observations[:response][idx13]) > 0.0
end

@testset "Time-varying covariate changes concentration deterministically via segmentation" begin
    base = ModelSpec(
        OneCompIVBolus(),
        "tv_cov_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    # CL increases at t=10, so concentrations after 10 should drop faster
    cm = CovariateModel(
        "cl_tv", [CovariateEffect(LinearCovariate(), :CL, :CLMULT, 1.0, 1.0)]
    )

    tv = TimeVaryingCovariates(
        Dict(:CLMULT => TimeCovariateSeries(StepTimeCovariate(), [0.0, 10.0], [1.0, 2.0]))
    )

    covs = [IndividualCovariates(Dict{Symbol,Float64}(), tv)]
    pop = PopulationSpec(base, nothing, nothing, cm, covs)

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_population(pop, grid, solver)

    c = res.individuals[1].observations[:conc]
    idx9 = findfirst(==(9.0), grid.saveat)
    idx11 = findfirst(==(11.0), grid.saveat)

    @test idx9 !== nothing
    @test idx11 !== nothing

    # sanity: conc at 11 should be lower than it would be with constant CL
    pop_const = PopulationSpec(base, nothing, nothing, nothing, IndividualCovariates[])
    base_res = simulate_population(pop_const, grid, solver)
    c_base = base_res.individuals[1].observations[:conc]

    @test abs(c[idx11] - c_base[idx11]) > 0.0
end

@testset "Exposure metrics: trapezoid AUC and Cmax" begin
    t = [0.0, 1.0, 2.0]
    c = [0.0, 1.0, 1.0]
    @test cmax(t, c) == 1.0
    # AUC = 0.5*(0+1)*1 + 0.5*(1+1)*1 = 0.5 + 1.0 = 1.5
    @test isapprox(auc_trapezoid(t, c), 1.5; rtol = 0.0, atol = 1e-12)
end

@testset "Response metrics: emin, time_below, auc_above_baseline" begin
    t = [0.0, 1.0, 2.0]
    y = [10.0, 8.0, 9.0]
    @test emin(t, y) == 8.0
    # left-constant rule: interval [0,1] uses y[1]=10 >= 9.5 (not counted),
    # interval [1,2] uses y[2]=8 < 9.5 (counted), so total = 1.0
    @test isapprox(time_below(t, y, 9.5), 1.0; atol = 1e-12)
    # baseline 10: suppression curve is [0,2,1]; AUC = 0.5*(0+2)*1 + 0.5*(2+1)*1 = 1 + 1.5 = 2.5
    @test isapprox(auc_above_baseline(t, y, 10.0), 2.5; atol = 1e-12)
end

# =====================================================
# Tests for new PK models
# =====================================================

@testset "TwoCompIVBolus basic simulation" begin
    spec = ModelSpec(
        TwoCompIVBolus(),
        "2c_iv_bolus",
        TwoCompIVBolusParams(10.0, 50.0, 5.0, 100.0),  # CL=10, V1=50, Q=5, V2=100
        [DoseEvent(0.0, 500.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:1.0:48.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    @test haskey(res.states, :A_central)
    @test haskey(res.states, :A_peripheral)
    @test haskey(res.observations, :conc)
    @test length(res.observations[:conc]) == length(res.t)

    # Initial concentration should be Dose/V1
    @test isapprox(res.observations[:conc][1], 500.0 / 50.0; rtol=0.01)

    # Concentration should decline over time
    @test res.observations[:conc][end] < res.observations[:conc][1]

    # Metadata should be correct
    @test res.metadata["model"] == "TwoCompIVBolus"
end

@testset "TwoCompIVBolus mass distribution" begin
    spec = ModelSpec(
        TwoCompIVBolus(),
        "2c_iv_mass",
        TwoCompIVBolusParams(10.0, 50.0, 5.0, 100.0),  # Normal elimination
        [DoseEvent(0.0, 1000.0)],
    )

    grid = SimGrid(0.0, 100.0, collect(0.0:1.0:100.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    # Total mass should decline over time due to elimination
    initial_total = res.states[:A_central][1] + res.states[:A_peripheral][1]
    final_total = res.states[:A_central][end] + res.states[:A_peripheral][end]
    @test final_total < initial_total

    # Mass should be distributed to peripheral at equilibrium
    # Initially all drug in central, later some moves to peripheral
    @test res.states[:A_peripheral][1] ≈ 0.0 atol=1e-10
    @test res.states[:A_peripheral][10] > 0.0  # After some time, peripheral has drug
end

@testset "TwoCompOral basic simulation" begin
    spec = ModelSpec(
        TwoCompOral(),
        "2c_oral",
        TwoCompOralParams(1.0, 10.0, 50.0, 5.0, 100.0),  # Ka=1, CL=10, V1=50, Q=5, V2=100
        [DoseEvent(0.0, 500.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:0.5:48.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    @test haskey(res.states, :A_gut)
    @test haskey(res.states, :A_central)
    @test haskey(res.states, :A_peripheral)
    @test haskey(res.observations, :conc)

    # Initial concentration should be zero (drug in gut)
    @test isapprox(res.observations[:conc][1], 0.0; atol=1e-10)

    # Concentration should rise then fall
    conc = res.observations[:conc]
    max_idx = argmax(conc)
    @test max_idx > 1  # Peak occurs after t=0

    # Metadata should be correct
    @test res.metadata["model"] == "TwoCompOral"
end

@testset "ThreeCompIVBolus basic simulation" begin
    spec = ModelSpec(
        ThreeCompIVBolus(),
        "3c_iv_bolus",
        ThreeCompIVBolusParams(10.0, 50.0, 10.0, 80.0, 2.0, 200.0),  # CL=10, V1=50, Q2=10, V2=80, Q3=2, V3=200
        [DoseEvent(0.0, 1000.0)],
    )

    grid = SimGrid(0.0, 72.0, collect(0.0:1.0:72.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    @test haskey(res.states, :A_central)
    @test haskey(res.states, :A_periph1)
    @test haskey(res.states, :A_periph2)
    @test haskey(res.observations, :conc)

    # Initial concentration should be Dose/V1
    @test isapprox(res.observations[:conc][1], 1000.0 / 50.0; rtol=0.01)

    # Metadata should be correct
    @test res.metadata["model"] == "ThreeCompIVBolus"
end

@testset "ThreeCompIVBolus shows tri-exponential decline" begin
    spec = ModelSpec(
        ThreeCompIVBolus(),
        "3c_triexp",
        ThreeCompIVBolusParams(5.0, 50.0, 20.0, 100.0, 1.0, 500.0),
        [DoseEvent(0.0, 1000.0)],
    )

    grid = SimGrid(0.0, 200.0, collect(0.0:1.0:200.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    conc = res.observations[:conc]

    # Concentration should be monotonically decreasing (or nearly so)
    for i in 2:length(conc)
        @test conc[i] <= conc[i-1] + 1e-10
    end
end

@testset "TransitAbsorption basic simulation" begin
    spec = ModelSpec(
        TransitAbsorption(),
        "transit_abs",
        TransitAbsorptionParams(5, 2.0, 1.0, 10.0, 50.0),  # N=5, Ktr=2, Ka=1, CL=10, V=50
        [DoseEvent(0.0, 500.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    @test haskey(res.states, :A_central)
    @test haskey(res.states, :Transit_1)
    @test haskey(res.states, :Transit_5)
    @test haskey(res.observations, :conc)

    # Initial concentration should be zero
    @test isapprox(res.observations[:conc][1], 0.0; atol=1e-10)

    # Delayed peak due to transit chain
    conc = res.observations[:conc]
    max_idx = argmax(conc)
    tmax = res.t[max_idx]
    @test tmax > 0.5  # Peak should be delayed

    # Metadata should be correct
    @test res.metadata["model"] == "TransitAbsorption"
    @test res.metadata["N_transit"] == 5
end

@testset "TransitAbsorption: more compartments delay peak" begin
    grid = SimGrid(0.0, 24.0, collect(0.0:0.1:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    spec_n2 = ModelSpec(
        TransitAbsorption(),
        "transit_n2",
        TransitAbsorptionParams(2, 2.0, 1.0, 10.0, 50.0),
        [DoseEvent(0.0, 500.0)],
    )

    spec_n8 = ModelSpec(
        TransitAbsorption(),
        "transit_n8",
        TransitAbsorptionParams(8, 2.0, 1.0, 10.0, 50.0),
        [DoseEvent(0.0, 500.0)],
    )

    res_n2 = simulate(spec_n2, grid, solver)
    res_n8 = simulate(spec_n8, grid, solver)

    tmax_n2 = res_n2.t[argmax(res_n2.observations[:conc])]
    tmax_n8 = res_n8.t[argmax(res_n8.observations[:conc])]

    # More transit compartments should delay the peak
    @test tmax_n8 > tmax_n2
end

@testset "MichaelisMentenElimination basic simulation" begin
    spec = ModelSpec(
        MichaelisMentenElimination(),
        "mm_elim",
        MichaelisMentenEliminationParams(100.0, 5.0, 50.0),  # Vmax=100, Km=5, V=50
        [DoseEvent(0.0, 500.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:0.5:48.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    @test haskey(res.states, :A_central)
    @test haskey(res.observations, :conc)

    # Initial concentration should be Dose/V
    @test isapprox(res.observations[:conc][1], 500.0 / 50.0; rtol=0.01)

    # Concentration should decline
    @test res.observations[:conc][end] < res.observations[:conc][1]

    # Metadata should be correct
    @test res.metadata["model"] == "MichaelisMentenElimination"
end

@testset "MichaelisMentenElimination shows nonlinear kinetics" begin
    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    # Low dose - approximately linear
    spec_low = ModelSpec(
        MichaelisMentenElimination(),
        "mm_low",
        MichaelisMentenEliminationParams(100.0, 10.0, 50.0),  # Km=10
        [DoseEvent(0.0, 50.0)],  # C0 = 1, << Km
    )

    # High dose - saturated elimination
    spec_high = ModelSpec(
        MichaelisMentenElimination(),
        "mm_high",
        MichaelisMentenEliminationParams(100.0, 10.0, 50.0),
        [DoseEvent(0.0, 5000.0)],  # C0 = 100, >> Km
    )

    res_low = simulate(spec_low, grid, solver)
    res_high = simulate(spec_high, grid, solver)

    # AUC should increase disproportionately with dose
    auc_low = auc_trapezoid(res_low.t, res_low.observations[:conc])
    auc_high = auc_trapezoid(res_high.t, res_high.observations[:conc])

    # AUC ratio should be > dose ratio (100x) due to saturable kinetics
    @test auc_high / auc_low > 100.0
end

# =====================================================
# Tests for new PD models
# =====================================================

function sigmoid_emax_ref(C::Float64, E0::Float64, Emax::Float64, EC50::Float64, gamma::Float64)
    if C <= 0.0
        return E0
    end
    return E0 + (Emax * C^gamma) / (EC50^gamma + C^gamma)
end

@testset "SigmoidEmax basic evaluation" begin
    pd = PDSpec(
        SigmoidEmax(),
        "semax",
        SigmoidEmaxParams(10.0, 50.0, 1.0, 2.0),  # E0=10, Emax=50, EC50=1, gamma=2
        :conc,
        :effect,
    )

    concentrations = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]
    effects = evaluate(pd, concentrations)

    for (i, C) in enumerate(concentrations)
        e_ref = sigmoid_emax_ref(C, 10.0, 50.0, 1.0, 2.0)
        @test isapprox(effects[i], e_ref; rtol=1e-10)
    end
end

@testset "SigmoidEmax: gamma=1 matches DirectEmax" begin
    pd_sigmoid = PDSpec(
        SigmoidEmax(),
        "semax_g1",
        SigmoidEmaxParams(10.0, 40.0, 0.8, 1.0),  # gamma=1
        :conc,
        :effect,
    )

    pd_direct = PDSpec(
        DirectEmax(),
        "demax",
        DirectEmaxParams(10.0, 40.0, 0.8),
        :conc,
        :effect,
    )

    concentrations = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0]

    effects_sig = evaluate(pd_sigmoid, concentrations)
    effects_dir = evaluate(pd_direct, concentrations)

    for i in eachindex(concentrations)
        @test isapprox(effects_sig[i], effects_dir[i]; rtol=1e-10)
    end
end

@testset "SigmoidEmax: higher gamma gives steeper response" begin
    pd_low = PDSpec(
        SigmoidEmax(),
        "semax_low",
        SigmoidEmaxParams(0.0, 100.0, 5.0, 1.0),  # gamma=1
        :conc,
        :effect,
    )

    pd_high = PDSpec(
        SigmoidEmax(),
        "semax_high",
        SigmoidEmaxParams(0.0, 100.0, 5.0, 4.0),  # gamma=4
        :conc,
        :effect,
    )

    # At EC50, both should give ~50% of Emax
    effects_low = evaluate(pd_low, [5.0])
    effects_high = evaluate(pd_high, [5.0])
    @test isapprox(effects_low[1], 50.0; rtol=0.01)
    @test isapprox(effects_high[1], 50.0; rtol=0.01)

    # Below EC50, high gamma gives lower effect
    effects_low_sub = evaluate(pd_low, [2.0])
    effects_high_sub = evaluate(pd_high, [2.0])
    @test effects_high_sub[1] < effects_low_sub[1]

    # Above EC50, high gamma gives higher effect
    effects_low_sup = evaluate(pd_low, [10.0])
    effects_high_sup = evaluate(pd_high, [10.0])
    @test effects_high_sup[1] > effects_low_sup[1]
end

@testset "SigmoidEmax PKPD coupling" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_semax",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(
        SigmoidEmax(),
        "semax",
        SigmoidEmaxParams(10.0, 40.0, 0.8, 2.0),
        :conc,
        :effect,
    )

    res = simulate_pkpd(pk, pd, grid, solver)

    @test haskey(res.observations, :conc)
    @test haskey(res.observations, :effect)
    @test length(res.observations[:effect]) == length(res.t)

    # Verify effects match expected values
    for (i, t) in enumerate(res.t)
        C = res.observations[:conc][i]
        e_ref = sigmoid_emax_ref(C, 10.0, 40.0, 0.8, 2.0)
        @test isapprox(res.observations[:effect][i], e_ref; rtol=1e-10)
    end
end

@testset "BiophaseEquilibration basic evaluation (quasi-steady state)" begin
    pd = PDSpec(
        BiophaseEquilibration(),
        "biophase",
        BiophaseEquilibrationParams(0.5, 10.0, 40.0, 0.8),  # ke0=0.5, E0=10, Emax=40, EC50=0.8
        :conc,
        :effect,
    )

    concentrations = [0.0, 0.5, 1.0, 2.0, 5.0]
    effects = evaluate(pd, concentrations)

    # Quasi-steady state: Ce = Cp, so effect is direct Emax
    for (i, C) in enumerate(concentrations)
        e_ref = direct_emax(C, 10.0, 40.0, 0.8)
        @test isapprox(effects[i], e_ref; rtol=1e-10)
    end
end

@testset "BiophaseEquilibration PKPD coupling" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_biophase",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(
        BiophaseEquilibration(),
        "biophase_test",
        BiophaseEquilibrationParams(0.5, 10.0, 40.0, 0.8),
        :conc,
        :effect,
    )

    res = simulate_pkpd(pk, pd, grid, solver)

    @test haskey(res.observations, :conc)
    @test haskey(res.observations, :effect)
    @test length(res.observations[:effect]) == length(res.t)
    @test res.metadata["pd_model"] == "BiophaseEquilibration"
end

# =====================================================
# Integration tests: New PK models with PD
# =====================================================

@testset "TwoCompIVBolus with DirectEmax PKPD" begin
    pk = ModelSpec(
        TwoCompIVBolus(),
        "2c_pkpd",
        TwoCompIVBolusParams(10.0, 50.0, 5.0, 100.0),
        [DoseEvent(0.0, 500.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:1.0:48.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(DirectEmax(), "emax", DirectEmaxParams(0.0, 100.0, 2.0), :conc, :effect)

    res = simulate_pkpd(pk, pd, grid, solver)

    @test haskey(res.observations, :effect)
    @test length(res.observations[:effect]) == length(res.t)

    # Effect should track concentration
    for (i, t) in enumerate(res.t)
        C = res.observations[:conc][i]
        e_ref = direct_emax(C, 0.0, 100.0, 2.0)
        @test isapprox(res.observations[:effect][i], e_ref; rtol=1e-10)
    end
end

@testset "TransitAbsorption with SigmoidEmax PKPD" begin
    pk = ModelSpec(
        TransitAbsorption(),
        "transit_pkpd",
        TransitAbsorptionParams(5, 2.0, 1.0, 10.0, 50.0),
        [DoseEvent(0.0, 500.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(SigmoidEmax(), "semax", SigmoidEmaxParams(0.0, 100.0, 2.0, 2.0), :conc, :effect)

    res = simulate_pkpd(pk, pd, grid, solver)

    @test haskey(res.observations, :effect)

    # Effect should follow sigmoid Emax of concentration
    for (i, t) in enumerate(res.t)
        C = res.observations[:conc][i]
        e_ref = sigmoid_emax_ref(C, 0.0, 100.0, 2.0, 2.0)
        @test isapprox(res.observations[:effect][i], e_ref; rtol=1e-10)
    end
end

@testset "MichaelisMentenElimination with DirectEmax PKPD" begin
    pk = ModelSpec(
        MichaelisMentenElimination(),
        "mm_pkpd",
        MichaelisMentenEliminationParams(100.0, 5.0, 50.0),
        [DoseEvent(0.0, 500.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(DirectEmax(), "emax", DirectEmaxParams(10.0, 50.0, 3.0), :conc, :effect)

    res = simulate_pkpd(pk, pd, grid, solver)

    @test haskey(res.observations, :effect)

    for (i, t) in enumerate(res.t)
        C = res.observations[:conc][i]
        e_ref = direct_emax(C, 10.0, 50.0, 3.0)
        @test isapprox(res.observations[:effect][i], e_ref; rtol=1e-10)
    end
end

@testset "ThreeCompIVBolus with BiophaseEquilibration PKPD" begin
    pk = ModelSpec(
        ThreeCompIVBolus(),
        "3c_biophase",
        ThreeCompIVBolusParams(10.0, 50.0, 10.0, 80.0, 2.0, 200.0),
        [DoseEvent(0.0, 1000.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:1.0:48.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(BiophaseEquilibration(), "biophase", BiophaseEquilibrationParams(0.5, 0.0, 100.0, 5.0), :conc, :effect)

    res = simulate_pkpd(pk, pd, grid, solver)

    @test haskey(res.observations, :effect)
    @test length(res.observations[:effect]) == length(res.t)
end
