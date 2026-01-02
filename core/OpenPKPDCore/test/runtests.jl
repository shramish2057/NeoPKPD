using Test
using OpenPKPDCore

function analytic_onecomp_ivbolus_conc(t::Float64, doses, CL::Float64, V::Float64)
    k = CL / V
    c = 0.0
    for d in doses
        if t >= d.time
            c += (d.amount / V) * exp(-k * (t - d.time))
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
        @test isapprox(res.conc[i], c_ref; rtol=1e-8, atol=1e-10)
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
        @test isapprox(res.conc[i], c_ref; rtol=1e-8, atol=1e-10)
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
        @test isapprox(res.conc[i], c_ref; rtol=1e-8, atol=1e-10)
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
        @test isapprox(res.conc[i], c_ref; rtol=1e-8, atol=1e-10)
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
        @test isapprox(res.conc[i], c_ref; rtol=1e-8, atol=1e-10)
    end
end
