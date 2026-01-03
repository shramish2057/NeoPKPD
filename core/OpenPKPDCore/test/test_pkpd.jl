# Coupled PK-PD model tests

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
