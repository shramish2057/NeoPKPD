using DifferentialEquations
using SciMLBase

export simulate

function validate(grid::SimGrid)
    if grid.t0 < 0.0
        error("t0 must be >= 0, got $(grid.t0)")
    end
    if !(grid.t1 > grid.t0)
        error("t1 must be > t0")
    end
    if isempty(grid.saveat)
        error("saveat must not be empty")
    end
    if any(t -> t < grid.t0 || t > grid.t1, grid.saveat)
        error("All saveat values must be within [t0, t1]")
    end
    if !issorted(grid.saveat)
        error("saveat must be sorted ascending")
    end
    return nothing
end

function validate(solver::SolverSpec)
    if !(solver.reltol > 0.0)
        error("Expected positive value for reltol, got $(solver.reltol)")
    end
    if !(solver.abstol > 0.0)
        error("Expected positive value for abstol, got $(solver.abstol)")
    end
    if solver.maxiters < 1
        error("maxiters must be >= 1")
    end
    return nothing
end

function _solver_alg(alg::Symbol)
    if alg == :Tsit5
        return Tsit5()
    elseif alg == :Rosenbrock23
        return Rosenbrock23()
    else
        error("Unsupported solver alg: $(alg). Supported: :Tsit5, :Rosenbrock23")
    end
end

function _dose_callback(doses::Vector{DoseEvent}, t0::Float64, t1::Float64)
    dose_times = Float64[]
    dose_amounts = Float64[]

    for d in doses
        if d.time > t0 && d.time <= t1
            push!(dose_times, d.time)
            push!(dose_amounts, d.amount)
        end
    end

    if isempty(dose_times)
        return nothing
    end

    function affect!(integrator)
        idx = findfirst(==(integrator.t), dose_times)
        if idx === nothing
            error("Internal error: dose time not found for t=$(integrator.t)")
        end
        integrator.u[1] += dose_amounts[idx]
    end

    return PresetTimeCallback(dose_times, affect!)
end

function simulate(
    spec::ModelSpec{OneCompIVBolus,OneCompIVBolusParams}, grid::SimGrid, solver::SolverSpec
)
    validate(spec)
    validate(grid)
    validate(solver)

    CL = spec.params.CL
    V = spec.params.V

    A0 = 0.0
    for d in spec.doses
        if d.time == grid.t0
            A0 += d.amount
        end
    end

    p = (CL=CL, V=V)
    u0 = [A0]
    tspan = (grid.t0, grid.t1)

    prob = ODEProblem(_ode_onecomp_ivbolus!, u0, tspan, p)
    cb = _dose_callback(spec.doses, grid.t0, grid.t1)

    sol = solve(
        prob,
        _solver_alg(solver.alg);
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
    )

    A = [u[1] for u in sol.u]
    C = [a / V for a in A]

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.0",
        "model" => "OneCompIVBolus",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "dose_schedule" => [(d.time, d.amount) for d in spec.doses],
        "deterministic_output_grid" => true,
    )

    return SimResult(Vector{Float64}(sol.t), A, C, metadata)
end

function simulate(
    spec::ModelSpec{OneCompOralFirstOrder,OneCompOralFirstOrderParams},
    grid::SimGrid,
    solver::SolverSpec,
)
    validate(spec)
    validate(grid)
    validate(solver)

    Ka = spec.params.Ka
    CL = spec.params.CL
    V = spec.params.V

    # Oral bolus doses add to Agut
    Agut0 = 0.0
    for d in spec.doses
        if d.time == grid.t0
            Agut0 += d.amount
        end
    end

    p = (Ka=Ka, CL=CL, V=V)
    u0 = [Agut0, 0.0]
    tspan = (grid.t0, grid.t1)

    prob = ODEProblem(_ode_onecomp_oral_first_order!, u0, tspan, p)

    # Dose callback adds to gut compartment
    dose_times = Float64[]
    dose_amounts = Float64[]
    for d in spec.doses
        if d.time > grid.t0 && d.time <= grid.t1
            push!(dose_times, d.time)
            push!(dose_amounts, d.amount)
        end
    end

    cb = nothing
    if !isempty(dose_times)
        function affect!(integrator)
            idx = findfirst(==(integrator.t), dose_times)
            if idx === nothing
                error("Internal error: dose time not found for t=$(integrator.t)")
            end
            integrator.u[1] += dose_amounts[idx]
        end
        cb = PresetTimeCallback(dose_times, affect!)
    end

    sol = solve(
        prob,
        _solver_alg(solver.alg);
        reltol=solver.reltol,
        abstol=solver.abstol,
        maxiters=solver.maxiters,
        saveat=grid.saveat,
        callback=cb,
    )

    Agut = [u[1] for u in sol.u]
    Acent = [u[2] for u in sol.u]
    C = [a / V for a in Acent]

    metadata = Dict{String,Any}(
        "engine_version" => "0.1.0",
        "model" => "OneCompOralFirstOrder",
        "solver_alg" => String(solver.alg),
        "reltol" => solver.reltol,
        "abstol" => solver.abstol,
        "dose_schedule" => [(d.time, d.amount) for d in spec.doses],
        "deterministic_output_grid" => true,
    )

    # For now SimResult.amount stores the central amount for multi-state models
    # We keep this convention explicit and will formalize output structs later.
    return SimResult(Vector{Float64}(sol.t), Acent, C, metadata)
end
