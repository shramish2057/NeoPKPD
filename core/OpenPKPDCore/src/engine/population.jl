using StableRNGs
using Distributions

export PopulationResult, simulate_population

"""
Population simulation output.

individuals:
- Vector of SimResult, one per individual

params:
- Vector of Dicts describing realized individual parameters

metadata:
- includes seed, n, iiv kind, and omega map
"""
struct PopulationResult
    individuals::Vector{SimResult}
    params::Vector{Dict{Symbol, Float64}}
    metadata::Dict{String, Any}
end

function validate(iiv::IIVSpec{LogNormalIIV})
    if iiv.n < 1
        error("IIVSpec n must be >= 1")
    end
    for (k, ω) in iiv.omegas
        _require_positive("omega for $(k)", ω)
    end
    return nothing
end

function _realize_params_log_normal(
    base_params,
    omegas::Dict{Symbol, Float64},
    rng,
)
    # base_params is a typed struct, we produce a new typed struct with per-individual values
    T = typeof(base_params)
    fn = fieldnames(T)

    vals = Dict{Symbol, Float64}()

    for f in fn
        θ = Float64(getfield(base_params, f))
        if haskey(omegas, f)
            ω = omegas[f]
            η = rand(rng, Normal(0.0, ω))
            vals[f] = θ * exp(η)
        else
            vals[f] = θ
        end
    end

    # Reconstruct typed params in declared field order
    new_params = T((vals[f] for f in fn)...)

    return new_params, vals
end

"""
simulate_population runs deterministic population simulations for supported PK models.

Current support:
- OneCompIVBolusParams: CL, V
- OneCompOralFirstOrderParams: Ka, CL, V

IIV is log-normal, per-parameter independent (diagonal omega) in v1.

Covariates are accepted but not used yet, kept for forward compatibility.
"""
function simulate_population(
    pop::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec,
)
    base = pop.base_model_spec

    if pop.iiv === nothing
        n = 1
        seed = UInt64(0)
        iiv_kind = "none"
        omegas = Dict{Symbol, Float64}()
        rng = StableRNG(0)
    else
        validate(pop.iiv)
        n = pop.iiv.n
        seed = pop.iiv.seed
        iiv_kind = string(typeof(pop.iiv.kind))
        omegas = pop.iiv.omegas
        rng = StableRNG(seed)
    end

    # Covariates vector must be empty or length n
    if !isempty(pop.covariates) && length(pop.covariates) != n
        error("covariates must be empty or length n")
    end

    individuals = Vector{SimResult}(undef, n)
    realized_params = Vector{Dict{Symbol, Float64}}(undef, n)

    for i in 1:n
        # Each individual gets deterministic stream from the main StableRNG by consuming in order.
        # This is reproducible and stable across platforms.
        if pop.iiv === nothing
            spec_i = base
            individuals[i] = simulate(spec_i, grid, solver)

            # capture base params as realized
            T = typeof(base.params)
            fn = fieldnames(T)
            d = Dict{Symbol, Float64}()
            for f in fn
                d[f] = Float64(getfield(base.params, f))
            end
            realized_params[i] = d
        else
            new_params, d = _realize_params_log_normal(base.params, omegas, rng)
            realized_params[i] = d

            spec_i = ModelSpec(base.kind, base.name * "_i$(i)", new_params, base.doses)
            individuals[i] = simulate(spec_i, grid, solver)
        end
    end

    metadata = Dict{String, Any}(
        "n" => n,
        "seed" => seed,
        "iiv_kind" => iiv_kind,
        "omegas" => Dict(String(k) => v for (k, v) in omegas),
        "engine_version" => "0.1.0",
        "event_semantics_version" => EVENT_SEMANTICS_VERSION,
        "solver_semantics_version" => SOLVER_SEMANTICS_VERSION,
    )

    return PopulationResult(individuals, realized_params, metadata)
end
