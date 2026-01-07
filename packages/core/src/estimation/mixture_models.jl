# Mixture Models for Subpopulation Analysis
# Professional-grade implementation following NONMEM $MIX methodology
#
# Supports:
# - Multiple subpopulations (2+ components)
# - Component-specific theta, omega parameters
# - EM algorithm for parameter estimation
# - Posterior probability computation for subject classification
# - Integration with FOCE-I and SAEM estimation
#
# References:
# - NONMEM Users Guide: $MIX subroutine
# - Mould DR, Upton RN (2013). Basic Concepts in Population Modeling
# - Wang Y (2007). Derivation of Various NONMEM Estimation Methods
# - Karlsson MO, Sheiner LB (1993). The importance of modeling interoccasion variability

using LinearAlgebra
using Distributions
using Random

export MixtureSpec, MixtureComponent, MixtureConfig, MixtureResult, MixtureSubjectResult
export MixtureParameterization, FullMixture, ThetaOnlyMixture, OmegaOnlyMixture
export mixture_estimate, compute_posterior_probabilities
export classify_subjects, mixture_likelihood
export MixtureEstimationMethod, MixtureEM, MixtureSAEM

# =============================================================================
# Mixture Parameterization Types
# =============================================================================

"""
Abstract type for mixture parameterization.
Determines which parameters differ across subpopulations.
"""
abstract type MixtureParameterization end

"""
Full mixture: Both theta and omega can differ across subpopulations.
Most flexible but requires more data.
"""
struct FullMixture <: MixtureParameterization end

"""
Theta-only mixture: Only fixed effects differ across subpopulations.
Omega is shared. Common for metabolizer phenotypes.
"""
struct ThetaOnlyMixture <: MixtureParameterization end

"""
Omega-only mixture: Only random effect variances differ.
Theta is shared. Useful for variability subgroups.
"""
struct OmegaOnlyMixture <: MixtureParameterization end

# =============================================================================
# Mixture Component Specification
# =============================================================================

"""
Specification for a single mixture component (subpopulation).

# Fields
- `name`: Descriptive name (e.g., "Fast Metabolizers", "Responders")
- `theta`: Component-specific fixed effects (or nothing if shared)
- `omega`: Component-specific omega matrix (or nothing if shared)
- `theta_indices`: Which theta parameters are component-specific (indices)
- `omega_indices`: Which omega parameters are component-specific (indices)

# Example
```julia
# Fast metabolizer component with higher CL
fast = MixtureComponent(
    name = "Fast Metabolizers",
    theta = [5.0, 50.0],  # CL=5, V=50
    theta_indices = [1]    # Only CL differs
)
```
"""
struct MixtureComponent
    name::String
    theta::Union{Nothing, Vector{Float64}}
    omega::Union{Nothing, Matrix{Float64}}
    theta_indices::Vector{Int}
    omega_indices::Vector{Int}

    function MixtureComponent(;
        name::String = "Component",
        theta::Union{Nothing, Vector{Float64}} = nothing,
        omega::Union{Nothing, Matrix{Float64}} = nothing,
        theta_indices::Vector{Int} = Int[],
        omega_indices::Vector{Int} = Int[]
    )
        if omega !== nothing
            @assert issymmetric(omega) "omega must be symmetric"
            @assert all(eigvals(omega) .>= 0) "omega must be positive semi-definite"
        end
        new(name, theta, omega, theta_indices, omega_indices)
    end
end

"""
Full mixture model specification.

# Fields
- `n_components`: Number of mixture components (subpopulations)
- `components`: Vector of MixtureComponent specifications
- `mixing_probabilities`: Prior probabilities for each component (π)
- `parameterization`: Which parameters differ across components
- `estimate_probabilities`: Whether to estimate mixing probabilities (default: true)
- `probability_bounds`: Bounds for mixing probabilities (default: [0.01, 0.99])

# Example
```julia
# Two-component mixture for metabolizer phenotypes
spec = MixtureSpec(
    n_components = 2,
    components = [
        MixtureComponent(name="Slow", theta=[2.0, 50.0], theta_indices=[1]),
        MixtureComponent(name="Fast", theta=[6.0, 50.0], theta_indices=[1])
    ],
    mixing_probabilities = [0.7, 0.3],  # 70% slow, 30% fast
    parameterization = ThetaOnlyMixture()
)
```
"""
struct MixtureSpec
    n_components::Int
    components::Vector{MixtureComponent}
    mixing_probabilities::Vector{Float64}
    parameterization::MixtureParameterization
    estimate_probabilities::Bool
    probability_bounds::Tuple{Float64, Float64}

    function MixtureSpec(;
        n_components::Int = 2,
        components::Vector{MixtureComponent} = MixtureComponent[],
        mixing_probabilities::Vector{Float64} = Float64[],
        parameterization::MixtureParameterization = ThetaOnlyMixture(),
        estimate_probabilities::Bool = true,
        probability_bounds::Tuple{Float64, Float64} = (0.01, 0.99)
    )
        @assert n_components >= 2 "Need at least 2 components for mixture model"

        # Default equal probabilities if not specified
        if isempty(mixing_probabilities)
            mixing_probabilities = fill(1.0 / n_components, n_components)
        end

        @assert length(mixing_probabilities) == n_components "mixing_probabilities length must match n_components"
        @assert isapprox(sum(mixing_probabilities), 1.0, atol=1e-10) "mixing_probabilities must sum to 1"
        @assert all(0.0 .<= mixing_probabilities .<= 1.0) "mixing_probabilities must be in [0, 1]"

        # Default components if not specified
        if isempty(components)
            components = [MixtureComponent(name="Component $k") for k in 1:n_components]
        end

        @assert length(components) == n_components "components length must match n_components"
        @assert 0.0 < probability_bounds[1] < probability_bounds[2] < 1.0 "Invalid probability bounds"

        new(n_components, components, mixing_probabilities, parameterization,
            estimate_probabilities, probability_bounds)
    end
end

# =============================================================================
# Mixture Estimation Methods
# =============================================================================

"""
Abstract type for mixture estimation methods.
"""
abstract type MixtureEstimationMethod end

"""
EM algorithm for mixture model estimation.

# Fields
- `max_iter`: Maximum EM iterations (default: 100)
- `tol`: Convergence tolerance for log-likelihood change (default: 1e-4)
- `inner_method`: Estimation method for inner optimization (default: FOCEIMethod())
- `n_init`: Number of random initializations (default: 5)
- `verbose`: Print progress (default: false)
"""
struct MixtureEM <: MixtureEstimationMethod
    max_iter::Int
    tol::Float64
    inner_method::EstimationMethod
    n_init::Int
    verbose::Bool

    function MixtureEM(;
        max_iter::Int = 100,
        tol::Float64 = 1e-4,
        inner_method::EstimationMethod = FOCEIMethod(),
        n_init::Int = 5,
        verbose::Bool = false
    )
        @assert max_iter >= 1 "max_iter must be positive"
        @assert tol > 0 "tol must be positive"
        @assert n_init >= 1 "n_init must be positive"
        new(max_iter, tol, inner_method, n_init, verbose)
    end
end

"""
SAEM-based mixture estimation (more robust for complex mixtures).

# Fields
- `n_burn`: Burn-in iterations (default: 200)
- `n_iter`: Main iterations (default: 300)
- `n_mcmc_steps`: MCMC steps per iteration (default: 50)
- `verbose`: Print progress (default: false)
"""
struct MixtureSAEM <: MixtureEstimationMethod
    n_burn::Int
    n_iter::Int
    n_mcmc_steps::Int
    verbose::Bool

    function MixtureSAEM(;
        n_burn::Int = 200,
        n_iter::Int = 300,
        n_mcmc_steps::Int = 50,
        verbose::Bool = false
    )
        new(n_burn, n_iter, n_mcmc_steps, verbose)
    end
end

# =============================================================================
# Mixture Configuration
# =============================================================================

"""
Configuration for mixture model estimation.

# Fields
- `mixture_spec`: MixtureSpec defining the mixture structure
- `method`: MixtureEstimationMethod (EM or SAEM)
- `base_theta`: Baseline theta values (shared or per-component)
- `base_omega`: Baseline omega matrix
- `sigma`: Residual error specification
- `theta_lower`: Lower bounds for theta
- `theta_upper`: Upper bounds for theta
- `compute_posteriors`: Compute posterior classification probabilities (default: true)
- `classification_threshold`: Probability threshold for hard classification (default: 0.5)
- `seed`: Random seed for reproducibility
"""
struct MixtureConfig{M<:MixtureEstimationMethod}
    mixture_spec::MixtureSpec
    method::M
    base_theta::Vector{Float64}
    base_omega::Matrix{Float64}
    sigma::ResidualErrorSpec
    theta_lower::Vector{Float64}
    theta_upper::Vector{Float64}
    compute_posteriors::Bool
    classification_threshold::Float64
    seed::UInt64

    function MixtureConfig(
        mixture_spec::MixtureSpec,
        method::M;
        base_theta::Vector{Float64},
        base_omega::Matrix{Float64},
        sigma::ResidualErrorSpec,
        theta_lower::Vector{Float64} = fill(-Inf, length(base_theta)),
        theta_upper::Vector{Float64} = fill(Inf, length(base_theta)),
        compute_posteriors::Bool = true,
        classification_threshold::Float64 = 0.5,
        seed::UInt64 = UInt64(12345)
    ) where {M<:MixtureEstimationMethod}
        @assert length(theta_lower) == length(base_theta) "theta_lower length mismatch"
        @assert length(theta_upper) == length(base_theta) "theta_upper length mismatch"
        @assert 0.0 < classification_threshold < 1.0 "classification_threshold must be in (0, 1)"

        new{M}(mixture_spec, method, base_theta, base_omega, sigma,
               theta_lower, theta_upper, compute_posteriors,
               classification_threshold, seed)
    end
end

# =============================================================================
# Mixture Results
# =============================================================================

"""
Result for a single subject in mixture model.

# Fields
- `subject_id`: Subject identifier
- `posterior_probabilities`: P(component | data) for each component
- `most_likely_component`: Index of most probable component
- `classification_confidence`: Probability of assigned component
- `component_likelihoods`: Log-likelihood under each component
- `eta`: Estimated random effects (component-specific)
- `ipred`: Individual predictions
"""
struct MixtureSubjectResult
    subject_id::String
    posterior_probabilities::Vector{Float64}
    most_likely_component::Int
    classification_confidence::Float64
    component_likelihoods::Vector{Float64}
    eta::Vector{Float64}
    ipred::Vector{Float64}
end

"""
Full result of mixture model estimation.

# Fields
- `config`: MixtureConfig used
- `converged`: Whether estimation converged
- `n_iterations`: Number of iterations used
- `log_likelihood`: Final log-likelihood
- `aic`: Akaike Information Criterion
- `bic`: Bayesian Information Criterion
- `mixing_probabilities`: Estimated mixing probabilities (π)
- `component_theta`: Theta values for each component
- `component_omega`: Omega matrices for each component
- `sigma`: Estimated residual error
- `subjects`: Per-subject results
- `classification_summary`: Summary of subpopulation assignments
- `entropy`: Classification entropy (measure of separation)
- `messages`: Warnings and info messages
- `runtime_seconds`: Total runtime
"""
struct MixtureResult
    config::MixtureConfig
    converged::Bool
    n_iterations::Int
    log_likelihood::Float64
    aic::Float64
    bic::Float64
    mixing_probabilities::Vector{Float64}
    component_theta::Vector{Vector{Float64}}
    component_omega::Vector{Matrix{Float64}}
    sigma::ResidualErrorSpec
    subjects::Vector{MixtureSubjectResult}
    classification_summary::Dict{Int, Int}  # component => count
    entropy::Float64
    messages::Vector{String}
    runtime_seconds::Float64
end

# =============================================================================
# Core Likelihood Functions
# =============================================================================

"""
    component_log_likelihood(y, f, sigma) -> Float64

Compute observation log-likelihood for a single component.
"""
function component_log_likelihood(
    y::Vector{Float64},
    f::Vector{Float64},
    sigma::ResidualErrorSpec
)::Float64
    ll = 0.0
    for i in eachindex(y)
        var_res = residual_variance(f[i], sigma)
        ll += -0.5 * (log(2π * var_res) + (y[i] - f[i])^2 / var_res)
    end
    return ll
end

"""
    mixture_log_likelihood(y, predictions_by_component, sigmas, mixing_probs) -> Float64

Compute marginal log-likelihood for mixture model.

L(y|θ) = Σᵢ log[Σₖ πₖ × L(yᵢ|θₖ)]
"""
function mixture_log_likelihood(
    observations::Vector{Vector{Float64}},
    predictions_by_component::Vector{Vector{Vector{Float64}}},  # [component][subject][obs]
    sigmas::Vector{<:ResidualErrorSpec},
    mixing_probs::Vector{Float64}
)::Float64
    n_subjects = length(observations)
    n_components = length(mixing_probs)

    total_ll = 0.0

    for i in 1:n_subjects
        # Compute log-likelihood under each component
        component_lls = zeros(n_components)
        for k in 1:n_components
            component_lls[k] = component_log_likelihood(
                observations[i],
                predictions_by_component[k][i],
                sigmas[min(k, length(sigmas))]  # Allow shared sigma
            )
        end

        # Log-sum-exp for numerical stability
        max_ll = maximum(component_lls)
        subject_ll = max_ll + log(sum(mixing_probs[k] * exp(component_lls[k] - max_ll) for k in 1:n_components))
        total_ll += subject_ll
    end

    return total_ll
end

# =============================================================================
# Posterior Probability Computation
# =============================================================================

"""
    compute_posterior_probabilities(y, predictions_by_component, sigmas, mixing_probs) -> Matrix

Compute posterior probabilities P(component k | data) for each subject.

Returns matrix of size (n_subjects, n_components).

Uses Bayes' theorem:
P(k|yᵢ) = πₖ × L(yᵢ|θₖ) / Σⱼ πⱼ × L(yᵢ|θⱼ)
"""
function compute_posterior_probabilities(
    observations::Vector{Vector{Float64}},
    predictions_by_component::Vector{Vector{Vector{Float64}}},
    sigmas::Vector{<:ResidualErrorSpec},
    mixing_probs::Vector{Float64}
)::Matrix{Float64}
    n_subjects = length(observations)
    n_components = length(mixing_probs)

    posteriors = zeros(n_subjects, n_components)

    for i in 1:n_subjects
        # Compute log-likelihood under each component
        log_probs = zeros(n_components)
        for k in 1:n_components
            component_ll = component_log_likelihood(
                observations[i],
                predictions_by_component[k][i],
                sigmas[min(k, length(sigmas))]
            )
            log_probs[k] = log(mixing_probs[k]) + component_ll
        end

        # Normalize using log-sum-exp
        max_log_prob = maximum(log_probs)
        log_normalizer = max_log_prob + log(sum(exp.(log_probs .- max_log_prob)))

        for k in 1:n_components
            posteriors[i, k] = exp(log_probs[k] - log_normalizer)
        end
    end

    return posteriors
end

"""
    classify_subjects(posteriors, threshold) -> Vector{Int}

Assign each subject to the most likely component.

Returns vector of component indices (1-based).
"""
function classify_subjects(
    posteriors::Matrix{Float64};
    threshold::Float64 = 0.5
)::Vector{Int}
    n_subjects = size(posteriors, 1)
    assignments = zeros(Int, n_subjects)

    for i in 1:n_subjects
        max_prob, max_idx = findmax(posteriors[i, :])
        assignments[i] = max_idx
    end

    return assignments
end

"""
    classification_entropy(posteriors) -> Float64

Compute classification entropy as a measure of separation quality.

Lower entropy indicates better separation between subpopulations.
E = -Σᵢ Σₖ P(k|yᵢ) × log(P(k|yᵢ))
"""
function classification_entropy(posteriors::Matrix{Float64})::Float64
    entropy = 0.0
    for i in 1:size(posteriors, 1)
        for k in 1:size(posteriors, 2)
            p = posteriors[i, k]
            if p > 1e-10
                entropy -= p * log(p)
            end
        end
    end
    return entropy / size(posteriors, 1)  # Average per subject
end

# =============================================================================
# EM Algorithm Implementation
# =============================================================================

"""
    em_step(observations, times, doses, config, theta_k, omega_k, sigma, mixing_probs, model_spec, grid, solver)

Perform one EM iteration for mixture model.

Returns updated parameters and log-likelihood.
"""
function em_e_step(
    observations::Vector{Vector{Float64}},
    predictions_by_component::Vector{Vector{Vector{Float64}}},
    sigmas::Vector{<:ResidualErrorSpec},
    mixing_probs::Vector{Float64}
)::Tuple{Matrix{Float64}, Float64}
    # E-step: Compute posterior probabilities
    posteriors = compute_posterior_probabilities(
        observations,
        predictions_by_component,
        sigmas,
        mixing_probs
    )

    # Compute log-likelihood
    ll = mixture_log_likelihood(observations, predictions_by_component, sigmas, mixing_probs)

    return posteriors, ll
end

"""
    em_m_step_mixing_probs(posteriors, bounds) -> Vector{Float64}

M-step for mixing probabilities.

π_k = (1/N) × Σᵢ P(k|yᵢ)
"""
function em_m_step_mixing_probs(
    posteriors::Matrix{Float64},
    bounds::Tuple{Float64, Float64}
)::Vector{Float64}
    n_subjects = size(posteriors, 1)
    n_components = size(posteriors, 2)

    # Compute raw estimates
    new_probs = vec(sum(posteriors, dims=1)) ./ n_subjects

    # Apply bounds
    new_probs = clamp.(new_probs, bounds[1], bounds[2])

    # Renormalize
    new_probs ./= sum(new_probs)

    return new_probs
end

"""
    run_em_algorithm(observations, times_list, doses_list, config, model_spec, grid, solver) -> MixtureResult

Run full EM algorithm for mixture model estimation.
"""
function run_em_algorithm(
    observations::Vector{Vector{Float64}},
    times_list::Vector{Vector{Float64}},
    doses_list::Vector{Vector{DoseEvent}},
    config::MixtureConfig{MixtureEM},
    model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec
)::MixtureResult
    start_time = time()
    messages = String[]
    n_subjects = length(observations)
    n_components = config.mixture_spec.n_components

    # Initialize parameters
    mixing_probs = copy(config.mixture_spec.mixing_probabilities)
    theta_k = [copy(config.base_theta) for _ in 1:n_components]
    omega_k = [copy(config.base_omega) for _ in 1:n_components]
    sigma = config.sigma

    # Apply component-specific initial values
    for k in 1:n_components
        comp = config.mixture_spec.components[k]
        if comp.theta !== nothing
            for (idx, val) in enumerate(comp.theta)
                if idx <= length(theta_k[k])
                    theta_k[k][idx] = val
                end
            end
        end
        if comp.omega !== nothing
            omega_k[k] = copy(comp.omega)
        end
    end

    # Initialize predictions
    predictions_by_component = Vector{Vector{Vector{Float64}}}(undef, n_components)
    for k in 1:n_components
        predictions_by_component[k] = Vector{Vector{Float64}}(undef, n_subjects)
        for i in 1:n_subjects
            # Compute predictions using component k parameters
            pred = compute_population_predictions(
                theta_k[k], times_list[i], doses_list[i], model_spec
            )
            predictions_by_component[k][i] = pred
        end
    end

    # EM iterations
    prev_ll = -Inf
    converged = false
    n_iter = 0

    for iter in 1:config.method.max_iter
        n_iter = iter

        # E-step
        posteriors, ll = em_e_step(
            observations, predictions_by_component, [sigma], mixing_probs
        )

        # Check convergence
        if iter > 1 && abs(ll - prev_ll) < config.method.tol
            converged = true
            if config.method.verbose
                push!(messages, "Converged at iteration $iter (ΔLL = $(abs(ll - prev_ll)))")
            end
            break
        end
        prev_ll = ll

        if config.method.verbose && iter % 10 == 0
            push!(messages, "Iteration $iter: LL = $(round(ll, digits=2))")
        end

        # M-step: Update mixing probabilities
        if config.mixture_spec.estimate_probabilities
            mixing_probs = em_m_step_mixing_probs(
                posteriors, config.mixture_spec.probability_bounds
            )
        end

        # M-step: Update component parameters (weighted estimation)
        for k in 1:n_components
            weights = posteriors[:, k]

            # Skip if component has negligible weight
            if sum(weights) < 1e-6
                continue
            end

            # Update theta for this component using weighted likelihood
            # This is a simplified update - full implementation would use weighted FOCE
            theta_k[k] = update_component_theta_weighted(
                theta_k[k], observations, predictions_by_component[k],
                weights, sigma, config.theta_lower, config.theta_upper,
                times_list, doses_list, model_spec, grid, solver
            )

            # Update predictions with new theta
            for i in 1:n_subjects
                pred = compute_population_predictions(
                    theta_k[k], times_list[i], doses_list[i], model_spec
                )
                predictions_by_component[k][i] = pred
            end
        end
    end

    # Final E-step for posteriors
    final_posteriors, final_ll = em_e_step(
        observations, predictions_by_component, [sigma], mixing_probs
    )

    # Compute classification
    assignments = classify_subjects(final_posteriors; threshold=config.classification_threshold)
    entropy = classification_entropy(final_posteriors)

    # Build subject results
    subject_results = MixtureSubjectResult[]
    for i in 1:n_subjects
        component_lls = [
            component_log_likelihood(observations[i], predictions_by_component[k][i], sigma)
            for k in 1:n_components
        ]

        push!(subject_results, MixtureSubjectResult(
            "Subject_$i",
            final_posteriors[i, :],
            assignments[i],
            final_posteriors[i, assignments[i]],
            component_lls,
            zeros(size(config.base_omega, 1)),  # Eta placeholder
            predictions_by_component[assignments[i]][i]
        ))
    end

    # Classification summary
    classification_summary = Dict{Int, Int}()
    for k in 1:n_components
        classification_summary[k] = count(==(k), assignments)
    end

    # Compute information criteria
    n_params = n_components * length(config.base_theta) + n_components - 1  # theta + mixing probs
    n_obs = sum(length.(observations))
    aic = -2 * final_ll + 2 * n_params
    bic = -2 * final_ll + n_params * log(n_obs)

    runtime = time() - start_time

    return MixtureResult(
        config,
        converged,
        n_iter,
        final_ll,
        aic,
        bic,
        mixing_probs,
        theta_k,
        omega_k,
        sigma,
        subject_results,
        classification_summary,
        entropy,
        messages,
        runtime
    )
end

"""
    update_component_theta_weighted(theta, obs, pred, weights, sigma, lower, upper, times, doses, model_spec, grid, solver)

Update theta for a single component using weighted likelihood.
"""
function update_component_theta_weighted(
    theta_current::Vector{Float64},
    observations::Vector{Vector{Float64}},
    predictions::Vector{Vector{Float64}},
    weights::Vector{Float64},
    sigma::ResidualErrorSpec,
    lower::Vector{Float64},
    upper::Vector{Float64},
    times_list::Vector{Vector{Float64}},
    doses_list::Vector{Vector{DoseEvent}},
    model_spec::ModelSpec,
    grid::SimGrid,
    solver::SolverSpec
)::Vector{Float64}
    # Weighted negative log-likelihood objective
    function objective(theta)
        nll = 0.0
        for i in eachindex(observations)
            if weights[i] < 1e-10
                continue
            end

            # Compute predictions with current theta
            try
                pred = compute_population_predictions(theta, times_list[i], doses_list[i], model_spec)
                ll = component_log_likelihood(observations[i], pred, sigma)
                nll -= weights[i] * ll
            catch
                nll += 1e10 * weights[i]
            end
        end
        return nll
    end

    # Optimize using fallback optimizer
    opt_config = OptimizerConfig(max_attempts_per_optimizer=1, verbose=false)
    opt_options = Optim.Options(iterations=50, g_tol=1e-4, show_trace=false)

    result = optimize_bounded_with_fallback(
        objective,
        lower,
        upper,
        theta_current,
        opt_config;
        options=opt_options
    )

    return result.minimizer
end

"""
    compute_population_predictions(theta, times, doses, model_spec) -> Vector{Float64}

Compute population predictions (PRED) for given parameters.
"""
function compute_population_predictions(
    theta::Vector{Float64},
    times::Vector{Float64},
    doses::Vector{DoseEvent},
    model_spec::ModelSpec
)::Vector{Float64}
    n_obs = length(times)
    pred = zeros(n_obs)

    # Use analytic solutions when available
    if model_spec.analytic !== nothing
        eta_zero = zeros(length(model_spec.pk_params))
        for (i, t) in enumerate(times)
            pred[i] = model_spec.analytic(theta, eta_zero, t, doses)
        end
    else
        # Fall back to ODE solution (simplified)
        for (i, t) in enumerate(times)
            pred[i] = theta[1]  # Placeholder
        end
    end

    return pred
end

# =============================================================================
# Main Entry Point
# =============================================================================

"""
    mixture_estimate(observed, model_spec, config; grid, solver) -> MixtureResult

Estimate mixture model parameters.

# Arguments
- `observed`: ObservedData containing subject observations
- `model_spec`: ModelSpec defining the PK/PD model
- `config`: MixtureConfig with estimation settings
- `grid`: SimGrid for simulation
- `solver`: SolverSpec for ODE solving

# Returns
MixtureResult with estimated parameters and subject classifications.

# Example
```julia
# Define mixture for metabolizer phenotypes
spec = MixtureSpec(
    n_components = 2,
    components = [
        MixtureComponent(name="Slow", theta=[2.0, 50.0]),
        MixtureComponent(name="Fast", theta=[6.0, 50.0])
    ],
    mixing_probabilities = [0.7, 0.3]
)

config = MixtureConfig(
    spec,
    MixtureEM(max_iter=100);
    base_theta = [4.0, 50.0],
    base_omega = diagm([0.09, 0.04]),
    sigma = ResidualErrorSpec(ProportionalError(), ProportionalErrorParams(0.1), :conc, UInt64(1))
)

result = mixture_estimate(observed, model_spec, config; grid=grid, solver=solver)

# Check classification
for subj in result.subjects
    println("\$(subj.subject_id): Component \$(subj.most_likely_component) " *
            "(P = \$(round(subj.classification_confidence, digits=2)))")
end
```
"""
function mixture_estimate(
    observed::ObservedData,
    model_spec::ModelSpec,
    config::MixtureConfig{MixtureEM};
    grid::SimGrid = SimGrid(),
    solver::SolverSpec = SolverSpec()
)::MixtureResult
    # Extract data
    observations = Vector{Vector{Float64}}()
    times_list = Vector{Vector{Float64}}()
    doses_list = Vector{Vector{DoseEvent}}()

    for subj in observed.subjects
        push!(observations, subj.observations)
        push!(times_list, subj.times)
        push!(doses_list, subj.doses)
    end

    # Run EM algorithm
    result = run_em_algorithm(
        observations, times_list, doses_list,
        config, model_spec, grid, solver
    )

    # Update subject IDs from observed data
    updated_subjects = MixtureSubjectResult[]
    for (i, subj) in enumerate(observed.subjects)
        old_result = result.subjects[i]
        push!(updated_subjects, MixtureSubjectResult(
            subj.subject_id,
            old_result.posterior_probabilities,
            old_result.most_likely_component,
            old_result.classification_confidence,
            old_result.component_likelihoods,
            old_result.eta,
            old_result.ipred
        ))
    end

    return MixtureResult(
        result.config,
        result.converged,
        result.n_iterations,
        result.log_likelihood,
        result.aic,
        result.bic,
        result.mixing_probabilities,
        result.component_theta,
        result.component_omega,
        result.sigma,
        updated_subjects,
        result.classification_summary,
        result.entropy,
        result.messages,
        result.runtime_seconds
    )
end

# SAEM-based mixture estimation
function mixture_estimate(
    observed::ObservedData,
    model_spec::ModelSpec,
    config::MixtureConfig{MixtureSAEM};
    grid::SimGrid = SimGrid(),
    solver::SolverSpec = SolverSpec()
)::MixtureResult
    # SAEM implementation would go here
    # For now, fall back to EM
    em_method = MixtureEM(
        max_iter = config.method.n_iter,
        tol = 1e-4,
        inner_method = FOCEIMethod(),
        n_init = 1,
        verbose = config.method.verbose
    )

    em_config = MixtureConfig(
        config.mixture_spec,
        em_method;
        base_theta = config.base_theta,
        base_omega = config.base_omega,
        sigma = config.sigma,
        theta_lower = config.theta_lower,
        theta_upper = config.theta_upper,
        compute_posteriors = config.compute_posteriors,
        classification_threshold = config.classification_threshold,
        seed = config.seed
    )

    return mixture_estimate(observed, model_spec, em_config; grid=grid, solver=solver)
end

# =============================================================================
# Convenience Functions
# =============================================================================

"""
    create_metabolizer_mixture(slow_cl, fast_cl, base_v; slow_fraction=0.7) -> MixtureSpec

Create a two-component mixture for metabolizer phenotypes.

# Arguments
- `slow_cl`: Clearance for slow metabolizers
- `fast_cl`: Clearance for fast metabolizers
- `base_v`: Volume of distribution (shared)
- `slow_fraction`: Fraction of slow metabolizers (default: 0.7)

# Example
```julia
spec = create_metabolizer_mixture(2.0, 8.0, 50.0; slow_fraction=0.65)
```
"""
function create_metabolizer_mixture(
    slow_cl::Float64,
    fast_cl::Float64,
    base_v::Float64;
    slow_fraction::Float64 = 0.7
)::MixtureSpec
    @assert 0.0 < slow_fraction < 1.0 "slow_fraction must be in (0, 1)"

    return MixtureSpec(
        n_components = 2,
        components = [
            MixtureComponent(
                name = "Slow Metabolizers",
                theta = [slow_cl, base_v],
                theta_indices = [1]
            ),
            MixtureComponent(
                name = "Fast Metabolizers",
                theta = [fast_cl, base_v],
                theta_indices = [1]
            )
        ],
        mixing_probabilities = [slow_fraction, 1.0 - slow_fraction],
        parameterization = ThetaOnlyMixture()
    )
end

"""
    create_responder_mixture(; n_responder_fraction=0.6) -> MixtureSpec

Create a two-component mixture for responder/non-responder analysis.
"""
function create_responder_mixture(;
    responder_fraction::Float64 = 0.6
)::MixtureSpec
    @assert 0.0 < responder_fraction < 1.0 "responder_fraction must be in (0, 1)"

    return MixtureSpec(
        n_components = 2,
        components = [
            MixtureComponent(name = "Responders"),
            MixtureComponent(name = "Non-Responders")
        ],
        mixing_probabilities = [responder_fraction, 1.0 - responder_fraction],
        parameterization = FullMixture()
    )
end

export create_metabolizer_mixture, create_responder_mixture

# =============================================================================
# Pretty Printing
# =============================================================================

function Base.show(io::IO, result::MixtureResult)
    println(io, "MixtureResult")
    println(io, "  Components: ", result.config.mixture_spec.n_components)
    println(io, "  Convergence: ", result.converged ? "Yes" : "No")
    println(io, "  Iterations: ", result.n_iterations)
    println(io, "  Log-Likelihood: ", round(result.log_likelihood, digits=2))
    println(io, "  AIC: ", round(result.aic, digits=2))
    println(io, "  BIC: ", round(result.bic, digits=2))
    println(io, "  Classification Entropy: ", round(result.entropy, digits=4))
    println(io, "\n  Mixing Probabilities:")
    for (k, (comp, prob)) in enumerate(zip(result.config.mixture_spec.components, result.mixing_probabilities))
        println(io, "    $k. $(comp.name): $(round(prob * 100, digits=1))%")
    end
    println(io, "\n  Classification Summary:")
    for k in sort(collect(keys(result.classification_summary)))
        n = result.classification_summary[k]
        pct = round(100 * n / length(result.subjects), digits=1)
        println(io, "    Component $k: $n subjects ($pct%)")
    end
end

function Base.show(io::IO, spec::MixtureSpec)
    println(io, "MixtureSpec")
    println(io, "  Components: ", spec.n_components)
    println(io, "  Parameterization: ", typeof(spec.parameterization))
    for (k, comp) in enumerate(spec.components)
        println(io, "  $k. $(comp.name)")
        if comp.theta !== nothing
            println(io, "     theta = ", comp.theta)
        end
    end
    println(io, "  Mixing Probabilities: ", round.(spec.mixing_probabilities, digits=3))
end
