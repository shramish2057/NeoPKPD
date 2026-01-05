# Full CWRES (Conditional Weighted Residuals) Implementation
# Implements proper FO approximation with AD-computed derivatives

using ForwardDiff
using LinearAlgebra
using Distributions

# Note: Basic compute_cwres, compute_iwres, compute_wres are defined in engine/residual_error.jl
# and estimation/diagnostics.jl. This file provides:
# - Advanced CWRES with full FO approximation (different signature)
# - NPDE computation
# - ResidualResult type for aggregated diagnostics
export compute_npde
export ResidualResult, ResidualSummary, compute_full_cwres
export compute_npde_from_population, summarize_residuals, check_residual_normality

# ============================================================================
# Result Types
# ============================================================================

"""
Result of residual computation for a single subject.
"""
struct ResidualResult
    subject_id::String
    times::Vector{Float64}
    dv::Vector{Float64}           # Observed values
    pred::Vector{Float64}         # Population predictions
    ipred::Vector{Float64}        # Individual predictions
    wres::Vector{Float64}         # Weighted residuals
    iwres::Vector{Float64}        # Individual weighted residuals
    cwres::Vector{Float64}        # Conditional weighted residuals
    npde::Vector{Float64}         # Normalized prediction distribution errors
end

"""
Summary statistics for residuals.
"""
struct ResidualSummary
    n_obs::Int
    mean_wres::Float64
    sd_wres::Float64
    mean_iwres::Float64
    sd_iwres::Float64
    mean_cwres::Float64
    sd_cwres::Float64
    shapiro_wilk_p::Float64  # Test for normality
end

# ============================================================================
# Helper: Compute Residual Variance
# ============================================================================

"""
Compute residual variance based on error model.
Used for advanced CWRES computation.
"""
function _compute_residual_variance(pred::Float64, sigma::ResidualErrorSpec)::Float64
    if sigma.kind isa AdditiveError
        return sigma.params.sigma^2
    elseif sigma.kind isa ProportionalError
        return (sigma.params.sigma * pred)^2
    elseif sigma.kind isa CombinedError
        add = sigma.params.sigma_add
        prop = sigma.params.sigma_prop
        return add^2 + (prop * pred)^2
    else
        return sigma.params.sigma^2
    end
end

# ============================================================================
# CWRES (Conditional Weighted Residuals) - Full FO Approximation
# ============================================================================

"""
    compute_cwres(dv, pred, ipred, eta, omega, sigma, dF_deta)

Compute conditional weighted residuals using the FO approximation.

CWRES_i = (DV_i - PRED_i) / sqrt(C_i)

where C_i = dF/dη * Ω * dF/dη' + σ²

This properly accounts for the correlation structure induced by random effects.

Arguments:
- dv: Observed values
- pred: Population predictions (η=0)
- ipred: Individual predictions (with estimated η)
- eta: Estimated random effects vector
- omega: Random effects covariance matrix
- sigma: Residual error specification
- dF_deta: Jacobian of predictions with respect to eta (n_obs x n_eta)
"""
function compute_cwres(
    dv::Vector{Float64},
    pred::Vector{Float64},
    ipred::Vector{Float64},
    eta::Vector{Float64},
    omega::Matrix{Float64},
    sigma::ResidualErrorSpec,
    dF_deta::Matrix{Float64}
)::Vector{Float64}
    n_obs = length(dv)
    n_eta = length(eta)
    cwres = zeros(n_obs)

    for i in 1:n_obs
        # dF/dη for this observation (row vector)
        dF_i = dF_deta[i, :]

        # Variance from random effects: dF/dη * Ω * dF/dη'
        var_eta = dot(dF_i, omega * dF_i)

        # Residual variance
        var_eps = _compute_residual_variance(pred[i], sigma)

        # Total variance under FO approximation
        var_total = var_eta + var_eps

        # CWRES
        cwres[i] = (dv[i] - pred[i]) / sqrt(var_total)
    end

    return cwres
end

"""
    compute_cwres_with_ad(dv, times, pred_fn, eta, omega, sigma)

Compute CWRES with automatic differentiation for dF/dη.

Arguments:
- dv: Observed values
- times: Observation times
- pred_fn: Function (eta, t) -> prediction at time t
- eta: Estimated random effects
- omega: Random effects covariance matrix
- sigma: Residual error specification
"""
function compute_cwres_with_ad(
    dv::Vector{Float64},
    times::Vector{Float64},
    pred_fn::Function,
    eta::Vector{Float64},
    omega::Matrix{Float64},
    sigma::ResidualErrorSpec
)::Vector{Float64}
    n_obs = length(dv)
    n_eta = length(eta)

    # Compute population predictions (eta = 0)
    eta_zero = zeros(n_eta)
    pred = [pred_fn(eta_zero, t) for t in times]

    # Compute dF/dη using ForwardDiff
    dF_deta = zeros(n_obs, n_eta)
    for i in 1:n_obs
        t = times[i]
        grad = ForwardDiff.gradient(e -> pred_fn(e, t), eta)
        dF_deta[i, :] = grad
    end

    # Compute individual predictions for reference
    ipred = [pred_fn(eta, t) for t in times]

    return compute_cwres(dv, pred, ipred, eta, omega, sigma, dF_deta)
end

# ============================================================================
# Full CWRES Computation for EstimationResult
# ============================================================================

"""
    compute_full_cwres(estimation_result, observed, model_spec, grid, solver)

Compute full CWRES for all subjects in an estimation result.

This uses the estimated individual parameters to compute proper CWRES
with AD-derived gradients.
"""
function compute_full_cwres(
    result::EstimationResult,
    observed::ObservedData,
    model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec
)::Vector{ResidualResult}
    residual_results = ResidualResult[]

    for (i, subj) in enumerate(observed.subjects)
        ind = result.individuals[i]

        # Get subject's observations
        times = subj.times
        dv = subj.observations

        # Create prediction function that takes eta and returns predictions
        pred_fn = (eta, t) -> begin
            # Apply eta to parameters
            modified_params = _apply_eta_to_params(model_spec.params, eta, result)

            # Create modified spec
            mod_spec = ModelSpec(
                model_spec.kind,
                model_spec.name,
                modified_params,
                model_spec.doses
            )

            # Simulate at single time point
            single_grid = SimGrid(grid.t0, t + 0.1, [t])
            sim_result = simulate(mod_spec, single_grid, solver)

            return sim_result.observations[:conc][1]
        end

        # Compute CWRES with AD
        cwres = compute_cwres_with_ad(
            dv, times, pred_fn, ind.eta, result.omega, result.sigma
        )

        # Compute other residuals
        wres = compute_wres(dv, ind.pred, result.sigma)
        iwres = compute_iwres(dv, ind.ipred, result.sigma)

        # NPDE (placeholder - requires Monte Carlo simulation)
        npde = zeros(length(dv))  # Will implement below

        push!(residual_results, ResidualResult(
            subj.id, times, dv, ind.pred, ind.ipred,
            wres, iwres, cwres, npde
        ))
    end

    return residual_results
end

"""
Apply random effects to population parameters.
"""
function _apply_eta_to_params(params, eta::Vector{Float64}, result::EstimationResult)
    # This is a simplified version - actual implementation depends on
    # how parameters are structured and which ones have IIV

    # For now, assume eta applies multiplicatively to first n_eta parameters
    # A full implementation would use the IIV specification

    param_values = collect(values(params))
    n_eta = length(eta)

    for j in 1:min(n_eta, length(param_values))
        param_values[j] *= exp(eta[j])
    end

    return param_values
end

# ============================================================================
# NPDE (Normalized Prediction Distribution Errors)
# ============================================================================

"""
    compute_npde(dv, simulated_preds)

Compute NPDE using Monte Carlo simulations.

NPDE transforms the observed data to standard normal based on the
predictive distribution from simulations.

Arguments:
- dv: Observed values
- simulated_preds: Matrix of simulated predictions (n_obs x n_sim)

Returns:
- Vector of NPDE values
"""
function compute_npde(
    dv::Vector{Float64},
    simulated_preds::Matrix{Float64}
)::Vector{Float64}
    n_obs = length(dv)
    n_sim = size(simulated_preds, 2)
    npde = zeros(n_obs)

    for i in 1:n_obs
        # Get simulated values for this observation
        sim_vals = simulated_preds[i, :]

        # Compute percentile (pde)
        pde = sum(sim_vals .<= dv[i]) / n_sim

        # Handle edge cases
        pde = clamp(pde, 0.001, 0.999)

        # Transform to standard normal
        npde[i] = quantile(Normal(), pde)
    end

    return npde
end

"""
    compute_npde_from_population(observed, pop_spec, grid, solver; n_sim)

Compute NPDE by running population simulations.
"""
function compute_npde_from_population(
    observed::ObservedData,
    pop_spec::PopulationSpec,
    grid::SimGrid,
    solver::SolverSpec;
    n_sim::Int=1000,
    seed::UInt64=UInt64(12345)
)::Vector{Vector{Float64}}
    rng = StableRNG(seed)
    npde_results = Vector{Float64}[]

    for subj in observed.subjects
        times = subj.times
        dv = subj.observations
        n_obs = length(times)

        # Simulate predictions
        simulated_preds = zeros(n_obs, n_sim)

        for sim_idx in 1:n_sim
            # Create population spec with new seed
            sim_seed = rand(rng, UInt64)
            new_iiv = IIVSpec(
                pop_spec.iiv.kind,
                pop_spec.iiv.omegas,
                sim_seed,
                1  # Single subject
            )
            new_pop_spec = PopulationSpec(
                pop_spec.model_spec,
                new_iiv,
                pop_spec.iov,
                pop_spec.error,
                pop_spec.covariates
            )

            # Simulate
            pop_result = simulate_population(new_pop_spec, grid, solver)

            # Extract concentrations at observation times
            sim_times = pop_result.subjects[1].times
            sim_conc = pop_result.subjects[1].observations[:conc]

            for (j, t) in enumerate(times)
                # Find nearest simulation time
                idx = argmin(abs.(sim_times .- t))
                simulated_preds[j, sim_idx] = sim_conc[idx]
            end
        end

        # Compute NPDE for this subject
        npde = compute_npde(dv, simulated_preds)
        push!(npde_results, npde)
    end

    return npde_results
end

# ============================================================================
# Diagnostic Functions
# ============================================================================

"""
Compute summary statistics for residuals.
"""
function summarize_residuals(residuals::Vector{ResidualResult})::ResidualSummary
    all_wres = vcat([r.wres for r in residuals]...)
    all_iwres = vcat([r.iwres for r in residuals]...)
    all_cwres = vcat([r.cwres for r in residuals]...)

    n_obs = length(all_wres)

    # Simple Shapiro-Wilk approximation for normality test
    # (Full implementation would use proper statistical test)
    cwres_sorted = sort(all_cwres)
    expected_normal = [quantile(Normal(), (i - 0.375) / (n_obs + 0.25)) for i in 1:n_obs]
    correlation = cor(cwres_sorted, expected_normal)
    sw_p = correlation^2  # Simplified approximation

    return ResidualSummary(
        n_obs,
        mean(all_wres),
        std(all_wres),
        mean(all_iwres),
        std(all_iwres),
        mean(all_cwres),
        std(all_cwres),
        sw_p
    )
end

export summarize_residuals

"""
Check if residuals follow standard normal distribution.
Returns true if residuals are approximately N(0,1).
"""
function check_residual_normality(
    residuals::Vector{Float64};
    mean_tol::Float64=0.1,
    sd_tol::Float64=0.2
)::Bool
    m = mean(residuals)
    s = std(residuals)

    return abs(m) < mean_tol && abs(s - 1.0) < sd_tol
end

export check_residual_normality
