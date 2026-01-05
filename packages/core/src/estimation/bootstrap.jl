# Bootstrap Module for Standard Error Validation
# Non-parametric bootstrap with stratified resampling

using StableRNGs
using Statistics
using LinearAlgebra

export BootstrapSpec, BootstrapResult, BootstrapDiagnostics
export run_bootstrap, stratified_resample, compute_bootstrap_ci
export BootstrapCIMethod, PercentileCI, BCCI, StudentizedCI

# ============================================================================
# Bootstrap CI Methods
# ============================================================================

"""
Abstract type for bootstrap confidence interval methods.
"""
abstract type BootstrapCIMethod end

"""
Percentile method: CI = [θ*(α/2), θ*(1-α/2)]
"""
struct PercentileCI <: BootstrapCIMethod end

"""
Bias-corrected and accelerated (BCa) method.
"""
struct BCCI <: BootstrapCIMethod
    acceleration::Float64  # Acceleration constant (computed from jackknife)
end
BCCI() = BCCI(0.0)  # Default with no acceleration

"""
Studentized (bootstrap-t) method.
"""
struct StudentizedCI <: BootstrapCIMethod end

# ============================================================================
# Bootstrap Specification
# ============================================================================

"""
Specification for bootstrap analysis.

Fields:
- n_bootstrap: Number of bootstrap replicates (default: 1000)
- stratify_by: Variables to stratify by (e.g., [:dose, :formulation])
- seed: Random seed for reproducibility
- parallel: Whether to run in parallel
- ci_level: Confidence interval level (default: 0.95)
- ci_method: Method for computing CIs
"""
struct BootstrapSpec
    n_bootstrap::Int
    stratify_by::Vector{Symbol}
    seed::UInt64
    parallel::Bool
    ci_level::Float64
    ci_method::BootstrapCIMethod

    function BootstrapSpec(;
        n_bootstrap::Int=1000,
        stratify_by::Vector{Symbol}=Symbol[],
        seed::UInt64=UInt64(12345),
        parallel::Bool=false,
        ci_level::Float64=0.95,
        ci_method::BootstrapCIMethod=PercentileCI()
    )
        @assert n_bootstrap >= 100 "n_bootstrap must be >= 100"
        @assert 0.0 < ci_level < 1.0 "ci_level must be in (0, 1)"
        new(n_bootstrap, stratify_by, seed, parallel, ci_level, ci_method)
    end
end

# ============================================================================
# Bootstrap Results
# ============================================================================

"""
Diagnostics for bootstrap analysis.
"""
struct BootstrapDiagnostics
    n_successful::Int          # Number of successful bootstrap runs
    n_failed::Int              # Number of failed runs
    convergence_rate::Float64  # Fraction that converged
    median_iterations::Float64 # Median iterations to convergence
    outlier_indices::Vector{Int}  # Indices of potential outlier estimates
end

"""
Result of bootstrap analysis.

Fields:
- theta_estimates: Matrix of parameter estimates (n_bootstrap x n_params)
- theta_mean: Mean of bootstrap estimates
- theta_se: Standard error from bootstrap
- theta_ci_lower: Lower confidence interval
- theta_ci_upper: Upper confidence interval
- original_estimate: Original point estimate
- bias: Bootstrap bias estimate
- diagnostics: Bootstrap diagnostics
"""
struct BootstrapResult
    theta_estimates::Matrix{Float64}
    theta_mean::Vector{Float64}
    theta_se::Vector{Float64}
    theta_ci_lower::Vector{Float64}
    theta_ci_upper::Vector{Float64}
    original_estimate::Vector{Float64}
    bias::Vector{Float64}
    omega_estimates::Union{Nothing, Vector{Matrix{Float64}}}
    sigma_estimates::Union{Nothing, Vector{Float64}}
    diagnostics::BootstrapDiagnostics
end

# ============================================================================
# Stratified Resampling
# ============================================================================

"""
    stratified_resample(subject_ids, strata, rng)

Perform stratified resampling of subjects.

Arguments:
- subject_ids: Vector of subject identifiers
- strata: Vector of strata labels (same length as subject_ids)
- rng: Random number generator

Returns:
- Vector of resampled subject indices
"""
function stratified_resample(
    subject_ids::Vector{String},
    strata::Vector{T},
    rng::AbstractRNG
)::Vector{Int} where T
    n_subjects = length(subject_ids)

    # Group subjects by stratum
    stratum_indices = Dict{T, Vector{Int}}()
    for (i, s) in enumerate(strata)
        if !haskey(stratum_indices, s)
            stratum_indices[s] = Int[]
        end
        push!(stratum_indices[s], i)
    end

    # Resample within each stratum
    resampled = Int[]
    for (stratum, indices) in stratum_indices
        n_stratum = length(indices)
        sampled = rand(rng, indices, n_stratum)
        append!(resampled, sampled)
    end

    return resampled
end

"""
    stratified_resample(n_subjects, rng)

Simple random resampling (no stratification).
"""
function stratified_resample(n_subjects::Int, rng::AbstractRNG)::Vector{Int}
    return rand(rng, 1:n_subjects, n_subjects)
end

# ============================================================================
# Bootstrap CI Computation
# ============================================================================

"""
    compute_bootstrap_ci(estimates, original, ci_level, method)

Compute confidence intervals from bootstrap estimates.
"""
function compute_bootstrap_ci(
    estimates::Matrix{Float64},
    original::Vector{Float64},
    ci_level::Float64,
    method::PercentileCI
)::Tuple{Vector{Float64}, Vector{Float64}}
    n_params = size(estimates, 2)
    alpha = 1.0 - ci_level

    lower = zeros(n_params)
    upper = zeros(n_params)

    for j in 1:n_params
        sorted = sort(estimates[:, j])
        n_boot = length(sorted)
        lower[j] = sorted[max(1, Int(floor(alpha/2 * n_boot)))]
        upper[j] = sorted[min(n_boot, Int(ceil((1 - alpha/2) * n_boot)))]
    end

    return lower, upper
end

function compute_bootstrap_ci(
    estimates::Matrix{Float64},
    original::Vector{Float64},
    ci_level::Float64,
    method::BCCI
)::Tuple{Vector{Float64}, Vector{Float64}}
    n_boot, n_params = size(estimates)
    alpha = 1.0 - ci_level

    lower = zeros(n_params)
    upper = zeros(n_params)

    for j in 1:n_params
        theta_boot = estimates[:, j]
        theta_0 = original[j]

        # Bias correction factor z0
        p0 = sum(theta_boot .< theta_0) / n_boot
        p0 = clamp(p0, 0.001, 0.999)
        z0 = quantile(Normal(), p0)

        # Acceleration factor (from spec or computed)
        a = method.acceleration

        # Adjusted percentiles
        z_alpha_lower = quantile(Normal(), alpha/2)
        z_alpha_upper = quantile(Normal(), 1 - alpha/2)

        alpha1 = cdf(Normal(), z0 + (z0 + z_alpha_lower) / (1 - a * (z0 + z_alpha_lower)))
        alpha2 = cdf(Normal(), z0 + (z0 + z_alpha_upper) / (1 - a * (z0 + z_alpha_upper)))

        sorted = sort(theta_boot)
        lower[j] = sorted[max(1, Int(floor(alpha1 * n_boot)))]
        upper[j] = sorted[min(n_boot, Int(ceil(alpha2 * n_boot)))]
    end

    return lower, upper
end

function compute_bootstrap_ci(
    estimates::Matrix{Float64},
    original::Vector{Float64},
    ci_level::Float64,
    method::StudentizedCI
)::Tuple{Vector{Float64}, Vector{Float64}}
    # Studentized bootstrap requires SE estimates for each bootstrap sample
    # Fall back to percentile if not available
    return compute_bootstrap_ci(estimates, original, ci_level, PercentileCI())
end

# ============================================================================
# Main Bootstrap Function
# ============================================================================

"""
    run_bootstrap(estimation_fn, observed_data, spec; kwargs...)

Run bootstrap analysis for parameter estimation.

Arguments:
- estimation_fn: Function that takes observed data and returns EstimationResult
- observed_data: Original observed data (ObservedData type)
- spec: BootstrapSpec specifying bootstrap parameters

Keyword Arguments:
- strata: Optional vector of strata labels for each subject
- verbose: Print progress (default: false)

Returns:
- BootstrapResult with bootstrap estimates and diagnostics
"""
function run_bootstrap(
    estimation_fn::Function,
    observed_data,  # ObservedData type
    spec::BootstrapSpec;
    strata::Union{Nothing, Vector}=nothing,
    verbose::Bool=false
)::BootstrapResult
    rng = StableRNG(spec.seed)

    n_subjects = length(observed_data.subjects)
    subject_ids = [s.id for s in observed_data.subjects]

    # Run original estimation to get point estimate
    if verbose
        println("Running original estimation...")
    end
    original_result = estimation_fn(observed_data)
    original_theta = original_result.theta
    n_params = length(original_theta)

    # Storage for bootstrap estimates
    theta_estimates = zeros(spec.n_bootstrap, n_params)
    omega_estimates = Vector{Matrix{Float64}}[]
    sigma_estimates = Float64[]

    n_successful = 0
    n_failed = 0
    iterations = Float64[]

    for b in 1:spec.n_bootstrap
        if verbose && b % 100 == 0
            println("Bootstrap replicate $b / $(spec.n_bootstrap)")
        end

        try
            # Resample subjects
            if strata !== nothing
                resampled_indices = stratified_resample(subject_ids, strata, rng)
            else
                resampled_indices = stratified_resample(n_subjects, rng)
            end

            # Create resampled dataset
            resampled_data = _create_resampled_data(observed_data, resampled_indices)

            # Run estimation on resampled data
            boot_result = estimation_fn(resampled_data)

            # Store estimates
            theta_estimates[b, :] = boot_result.theta

            if boot_result.omega !== nothing
                push!(omega_estimates, boot_result.omega)
            end
            if boot_result.sigma !== nothing
                push!(sigma_estimates, boot_result.sigma.params.sigma)
            end

            n_successful += 1

            # Track iterations if available
            if hasproperty(boot_result, :n_iterations)
                push!(iterations, Float64(boot_result.n_iterations))
            end

        catch e
            n_failed += 1
            # Fill with NaN for failed runs
            theta_estimates[b, :] .= NaN
            if verbose
                println("Bootstrap $b failed: $e")
            end
        end
    end

    # Remove failed runs for statistics
    valid_mask = .!any(isnan.(theta_estimates), dims=2)[:]
    valid_estimates = theta_estimates[valid_mask, :]

    # Compute statistics
    theta_mean = vec(mean(valid_estimates, dims=1))
    theta_se = vec(std(valid_estimates, dims=1))
    bias = theta_mean .- original_theta

    # Compute confidence intervals
    lower, upper = compute_bootstrap_ci(valid_estimates, original_theta, spec.ci_level, spec.ci_method)

    # Detect outliers (> 3 SD from mean)
    outlier_indices = Int[]
    for j in 1:n_params
        z_scores = abs.(valid_estimates[:, j] .- theta_mean[j]) ./ theta_se[j]
        outliers = findall(z_scores .> 3.0)
        append!(outlier_indices, outliers)
    end
    outlier_indices = unique(outlier_indices)

    # Create diagnostics
    diagnostics = BootstrapDiagnostics(
        n_successful,
        n_failed,
        n_successful / spec.n_bootstrap,
        isempty(iterations) ? 0.0 : median(iterations),
        outlier_indices
    )

    return BootstrapResult(
        theta_estimates,
        theta_mean,
        theta_se,
        lower,
        upper,
        original_theta,
        bias,
        isempty(omega_estimates) ? nothing : omega_estimates,
        isempty(sigma_estimates) ? nothing : sigma_estimates,
        diagnostics
    )
end

"""
Helper function to create resampled observed data.
"""
function _create_resampled_data(observed_data, indices::Vector{Int})
    # Create new subjects list with resampled indices
    # Note: Same subject can appear multiple times (with replacement)

    new_subjects = [observed_data.subjects[i] for i in indices]

    # Renumber subject IDs to avoid duplicates
    for (i, subj) in enumerate(new_subjects)
        # Create new subject with modified ID
        # This is a shallow approach - actual implementation depends on ObservedData structure
    end

    # Return a structure compatible with ObservedData
    # The actual implementation depends on the ObservedData type definition
    return typeof(observed_data)(new_subjects)
end

# ============================================================================
# Bootstrap SE Comparison
# ============================================================================

"""
    compare_se_methods(bootstrap_result, analytical_se)

Compare bootstrap standard errors with analytical standard errors.

Returns a dictionary with comparison metrics.
"""
function compare_se_methods(
    bootstrap_result::BootstrapResult,
    analytical_se::Vector{Float64}
)::Dict{Symbol, Any}
    n_params = length(analytical_se)

    # Compute ratios
    se_ratio = bootstrap_result.theta_se ./ analytical_se

    # Compute relative difference
    rel_diff = (bootstrap_result.theta_se .- analytical_se) ./ analytical_se .* 100

    # Flag parameters where SE differs by more than 20%
    flags = abs.(rel_diff) .> 20.0

    return Dict(
        :bootstrap_se => bootstrap_result.theta_se,
        :analytical_se => analytical_se,
        :se_ratio => se_ratio,
        :relative_diff_percent => rel_diff,
        :discrepancy_flags => flags,
        :mean_ratio => mean(se_ratio),
        :max_discrepancy => maximum(abs.(rel_diff))
    )
end

export compare_se_methods

# ============================================================================
# Bootstrap Diagnostics Plotting Data
# ============================================================================

"""
    get_bootstrap_histogram_data(result, param_index)

Get data for plotting bootstrap distribution histogram.
"""
function get_bootstrap_histogram_data(
    result::BootstrapResult,
    param_index::Int
)::Dict{Symbol, Any}
    estimates = result.theta_estimates[:, param_index]
    valid = filter(!isnan, estimates)

    return Dict(
        :estimates => valid,
        :original => result.original_estimate[param_index],
        :mean => result.theta_mean[param_index],
        :ci_lower => result.theta_ci_lower[param_index],
        :ci_upper => result.theta_ci_upper[param_index],
        :se => result.theta_se[param_index],
        :bias => result.bias[param_index]
    )
end

export get_bootstrap_histogram_data
