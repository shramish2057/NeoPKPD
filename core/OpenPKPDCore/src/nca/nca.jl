# NCA Module
# Non-Compartmental Analysis following FDA/EMA guidance

export run_nca, run_population_nca, summarize_population_nca
export nca_from_simresult, round_nca_result

# Include NCA components
include("specs.jl")
include("lambda_z.jl")
include("auc.jl")
include("exposure_metrics.jl")
include("pk_parameters.jl")
include("multiple_dose.jl")
include("bioequivalence.jl")

# =============================================================================
# Main NCA Entry Point
# =============================================================================

"""
    run_nca(t, c, dose; config=NCAConfig(), dosing_type=:single, tau=nothing, route=:extravascular)

Perform complete Non-Compartmental Analysis on concentration-time data.

# Arguments
- `t::Vector{Float64}`: Time points (sorted, ascending)
- `c::Vector{Float64}`: Concentration values
- `dose::Float64`: Administered dose

# Keyword Arguments
- `config::NCAConfig`: NCA configuration (default: standard FDA/EMA settings)
- `dosing_type::Symbol`: `:single`, `:multiple`, or `:steady_state` (default: :single)
- `tau::Union{Float64,Nothing}`: Dosing interval for multiple dose (required if dosing_type != :single)
- `route::Symbol`: Administration route - `:extravascular`, `:iv_bolus`, or `:iv_infusion`
- `t_inf::Float64`: Infusion duration (required for :iv_infusion)

# Returns
- `NCAResult`: Complete NCA results

# Example
```julia
t = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
c = [0.0, 1.2, 2.0, 1.8, 1.2, 0.6, 0.3, 0.075]
dose = 100.0

result = run_nca(t, c, dose)
println("Cmax: \$(result.cmax)")
println("AUC0-inf: \$(result.auc_0_inf)")
println("t1/2: \$(result.t_half)")
```
"""
function run_nca(
    t::Vector{Float64},
    c::Vector{Float64},
    dose::Float64;
    config::NCAConfig = NCAConfig(),
    dosing_type::Symbol = :single,
    tau::Union{Float64,Nothing} = nothing,
    route::Symbol = :extravascular,
    t_inf::Float64 = 0.0
)
    # Validate inputs
    n = length(t)
    @assert n == length(c) "Time and concentration vectors must have same length"
    @assert n >= 3 "Need at least 3 data points for NCA"
    @assert dose > 0.0 "Dose must be positive"
    @assert issorted(t) "Time points must be sorted ascending"

    if dosing_type in [:multiple, :steady_state]
        @assert tau !== nothing "Dosing interval (tau) required for multiple dose analysis"
        @assert tau > 0.0 "Dosing interval must be positive"
    end

    warnings = String[]
    quality_flags = Symbol[]

    # Apply BLQ handling
    c_clean = _apply_blq_handling(c, config)

    # ==========================================================================
    # Primary Exposure Metrics
    # ==========================================================================

    cmax, tmax_idx = find_cmax(t, c_clean)
    tmax = t[tmax_idx]

    clast, tlast, tlast_idx = find_clast(t, c_clean; lloq=config.lloq !== nothing ? config.lloq : 0.0)

    cmin = nothing
    cavg = nothing

    if dosing_type in [:multiple, :steady_state] && tau !== nothing
        cmin = nca_cmin(c_clean[findall(ti -> 0.0 <= ti <= tau, t)])
        cavg = nca_cavg(t, c_clean, tau, config)
    end

    # ==========================================================================
    # AUC Calculations
    # ==========================================================================

    auc_0t = auc_0_t(t, c_clean, config)
    aumc_0t = aumc_0_t(t, c_clean, config)

    # ==========================================================================
    # Lambda-z Estimation
    # ==========================================================================

    lambda_z_result = estimate_lambda_z(t, c_clean, config; tmax_idx=tmax_idx)

    # Extract lambda_z values
    lambda_z = lambda_z_result.lambda_z
    t_half = lambda_z_result.t_half

    if lambda_z_result.quality_flag == :insufficient
        push!(quality_flags, :lambda_z_insufficient)
        append!(warnings, lambda_z_result.warnings)
    elseif lambda_z_result.quality_flag == :warning
        push!(quality_flags, :lambda_z_warning)
        append!(warnings, lambda_z_result.warnings)
    end

    # ==========================================================================
    # AUC Extrapolation (if lambda_z available)
    # ==========================================================================

    auc_inf = nothing
    auc_extra_pct = nothing
    aumc_inf = nothing

    if lambda_z !== nothing && clast > 0.0
        auc_inf, auc_extra_pct = auc_0_inf(t, c_clean, lambda_z, clast, config)
        aumc_inf = aumc_0_inf(t, c_clean, lambda_z, clast, tlast, config)

        if auc_extra_pct > config.extrapolation_max_pct
            push!(quality_flags, :high_extrapolation)
            push!(warnings, "AUC extrapolation $(round(auc_extra_pct, digits=1))% exceeds $(config.extrapolation_max_pct)% threshold")
        end
    end

    # ==========================================================================
    # Multiple Dose / Steady State AUC
    # ==========================================================================

    auc_tau = nothing
    if tau !== nothing
        auc_tau = auc_0_tau(t, c_clean, tau, config)
    end

    # ==========================================================================
    # PK Parameters
    # ==========================================================================

    mrt = nothing
    cl_f = nothing
    vz_f = nothing
    vss = nothing

    if auc_inf !== nothing && aumc_inf !== nothing
        mrt = nca_mrt(aumc_inf, auc_inf; route=route, t_inf=t_inf)

        cl_f = nca_cl_f(dose, auc_inf)

        if lambda_z !== nothing
            vz_f = nca_vz_f(dose, lambda_z, auc_inf)
        end

        if route == :iv_bolus && mrt !== nothing && cl_f !== nothing
            vss = nca_vss(cl_f, mrt)
        end
    end

    # For IV bolus, CL and Vz (not CL/F and Vz/F)
    # The naming convention changes but calculation is same
    # Users should interpret based on route

    # ==========================================================================
    # Multiple Dose Metrics
    # ==========================================================================

    accumulation_index = nothing
    ptf = nothing
    swing = nothing

    if dosing_type in [:multiple, :steady_state] && tau !== nothing && cmin !== nothing
        if auc_inf !== nothing && auc_tau !== nothing
            accumulation_index = nca_accumulation_index(auc_tau, auc_inf)
        end

        if cavg !== nothing && cmin > 0.0
            ptf = nca_ptf(cmax, cmin, cavg)
            swing = nca_swing(cmax, cmin)
        end
    end

    # ==========================================================================
    # Dose-Normalized Metrics
    # ==========================================================================

    cmax_dn = nca_dose_normalized_cmax(cmax, dose)
    auc_dn = auc_inf !== nothing ? nca_dose_normalized_auc(auc_inf, dose) : nothing

    # ==========================================================================
    # Build Result
    # ==========================================================================

    metadata = Dict{String,Any}(
        "dose" => dose,
        "route" => String(route),
        "dosing_type" => String(dosing_type),
        "n_points" => n,
        "config_method" => string(typeof(config.method))
    )

    if tau !== nothing
        metadata["tau"] = tau
    end

    return NCAResult(
        # Primary exposure
        cmax, tmax, cmin, clast, tlast, cavg,
        # AUC
        auc_0t, auc_inf, auc_extra_pct, auc_tau, aumc_0t, aumc_inf,
        # Terminal phase
        lambda_z_result,
        # PK parameters
        t_half, mrt, cl_f, vz_f, vss,
        # Multiple dose
        accumulation_index, ptf, swing,
        # Dose-normalized
        cmax_dn, auc_dn,
        # Quality
        quality_flags, warnings, metadata
    )
end

# =============================================================================
# Population NCA
# =============================================================================

"""
    run_population_nca(pop_result, dose; config=NCAConfig(), observation=:conc)

Perform NCA analysis on each individual in a population simulation result.

# Arguments
- `pop_result::PopulationResult`: Population simulation result
- `dose::Float64`: Administered dose

# Keyword Arguments
- `config::NCAConfig`: NCA configuration
- `observation::Symbol`: Observation to analyze (default: :conc)
- `dosing_type::Symbol`: Dosing type (default: :single)
- `tau::Union{Float64,Nothing}`: Dosing interval for multiple dose
- `route::Symbol`: Administration route (default: :extravascular)

# Returns
- `Vector{NCAResult}`: NCA results for each individual
"""
function run_population_nca(
    pop_result::PopulationResult,
    dose::Float64;
    config::NCAConfig = NCAConfig(),
    observation::Symbol = :conc,
    dosing_type::Symbol = :single,
    tau::Union{Float64,Nothing} = nothing,
    route::Symbol = :extravascular
)
    n_individuals = length(pop_result.individuals)
    results = Vector{NCAResult}(undef, n_individuals)

    for i in 1:n_individuals
        ind = pop_result.individuals[i]

        t = ind.t
        c = ind.observations[observation]

        results[i] = run_nca(
            t, c, dose;
            config=config,
            dosing_type=dosing_type,
            tau=tau,
            route=route
        )
    end

    return results
end

"""
    summarize_population_nca(nca_results; parameters=[:cmax, :auc_0_inf, :t_half])

Summarize NCA results across a population.

# Arguments
- `nca_results::Vector{NCAResult}`: NCA results from population
- `parameters::Vector{Symbol}`: Parameters to summarize

# Returns
- `Dict{Symbol, NamedTuple}`: Summary statistics for each parameter
"""
function summarize_population_nca(
    nca_results::Vector{NCAResult};
    parameters::Vector{Symbol} = [:cmax, :auc_0_inf, :t_half, :cl_f]
)
    n = length(nca_results)
    summaries = Dict{Symbol, NamedTuple}()

    for param in parameters
        values = Float64[]

        for res in nca_results
            val = getfield(res, param)
            if val !== nothing
                push!(values, val)
            end
        end

        if !isempty(values)
            n_valid = length(values)
            mean_val = sum(values) / n_valid
            sorted_vals = sort(values)
            median_val = n_valid % 2 == 0 ?
                (sorted_vals[n_valid÷2] + sorted_vals[n_valid÷2+1]) / 2 :
                sorted_vals[(n_valid+1)÷2]
            sd_val = sqrt(sum((values .- mean_val).^2) / (n_valid - 1))
            cv_val = sd_val / mean_val * 100.0
            gm_val = exp(sum(log.(values)) / n_valid)

            summaries[param] = (
                n = n_valid,
                mean = mean_val,
                median = median_val,
                sd = sd_val,
                cv_pct = cv_val,
                geometric_mean = gm_val,
                min = minimum(values),
                max = maximum(values)
            )
        end
    end

    return summaries
end

# =============================================================================
# Helper Functions
# =============================================================================

"""
    _apply_blq_handling(c, config)

Apply BLQ handling to concentration data.
"""
function _apply_blq_handling(c::Vector{Float64}, config::NCAConfig)
    if config.lloq === nothing
        return c
    end

    c_clean = copy(c)
    lloq = config.lloq

    for i in eachindex(c_clean)
        if c_clean[i] < lloq
            if config.blq_handling isa BLQZero
                c_clean[i] = 0.0
            elseif config.blq_handling isa BLQLLOQHalf
                c_clean[i] = lloq / 2.0
            elseif config.blq_handling isa BLQMissing
                c_clean[i] = NaN  # Will be handled during calculations
            end
        end
    end

    return c_clean
end

"""
    nca_from_simresult(result, dose; config=NCAConfig(), route=:extravascular)

Convenience function to run NCA directly on a SimResult.

# Arguments
- `result::SimResult`: Simulation result
- `dose::Float64`: Administered dose

# Returns
- `NCAResult`: NCA analysis results
"""
function nca_from_simresult(
    result::SimResult,
    dose::Float64;
    config::NCAConfig = NCAConfig(),
    route::Symbol = :extravascular,
    observation::Symbol = :conc
)
    t = result.t
    c = result.observations[observation]

    return run_nca(t, c, dose; config=config, route=route)
end

# =============================================================================
# Regulatory Rounding
# =============================================================================

"""
    round_nca_result(value, significant_digits)

Round NCA result to specified significant digits for regulatory reporting.
"""
function round_nca_result(value::Float64, significant_digits::Int)
    if value == 0.0
        return 0.0
    end

    d = ceil(Int, log10(abs(value)))
    factor = 10.0^(significant_digits - d)

    return round(value * factor) / factor
end
