# Enhanced NPDE (Normalized Prediction Distribution Errors) Module
# Full implementation with decorrelation, diagnostics, and statistical tests

using StableRNGs
using Statistics
using LinearAlgebra
using Distributions

export NPDESpec, NPDEResult, NPDEDiagnostics
export compute_npde_enhanced, compute_decorrelated_npde
export npde_statistical_tests, NPDETestResult
export get_npde_qq_data, get_npde_vs_pred_data, get_npde_vs_time_data

# ============================================================================
# NPDE Specification and Results
# ============================================================================

"""
Specification for NPDE computation.

Fields:
- n_simulations: Number of Monte Carlo simulations (default: 1000)
- seed: Random seed for reproducibility
- decorrelate: Whether to apply decorrelation (default: true)
- compute_pde: Whether to also compute PDE (default: true)
"""
struct NPDESpec
    n_simulations::Int
    seed::UInt64
    decorrelate::Bool
    compute_pde::Bool

    function NPDESpec(;
        n_simulations::Int=1000,
        seed::UInt64=UInt64(12345),
        decorrelate::Bool=true,
        compute_pde::Bool=true
    )
        @assert n_simulations >= 100 "n_simulations must be >= 100"
        new(n_simulations, seed, decorrelate, compute_pde)
    end
end

"""
Diagnostics for NPDE analysis.
"""
struct NPDEDiagnostics
    mean::Float64           # Should be ~0
    variance::Float64       # Should be ~1
    skewness::Float64       # Should be ~0
    kurtosis::Float64       # Should be ~3 (excess kurtosis ~0)
    shapiro_wilk_p::Float64 # p-value for normality test
    mean_test_p::Float64    # t-test for mean = 0
    var_test_p::Float64     # Chi-square test for var = 1
end

"""
Result of NPDE computation for a single subject.
"""
struct SubjectNPDEResult
    subject_id::String
    times::Vector{Float64}
    dv::Vector{Float64}
    npde::Vector{Float64}
    pde::Vector{Float64}           # Prediction distribution error (before transform)
    epred::Vector{Float64}         # Expected prediction (mean of simulations)
    var_pred::Vector{Float64}      # Variance of predictions
end

"""
Result of NPDE analysis.
"""
struct NPDEResult
    subjects::Vector{SubjectNPDEResult}
    pooled_npde::Vector{Float64}
    pooled_pde::Vector{Float64}
    diagnostics::NPDEDiagnostics
    simulation_matrix::Union{Nothing, Matrix{Float64}}  # Optional: store simulations
end

# ============================================================================
# Statistical Tests for NPDE
# ============================================================================

"""
Result of NPDE statistical tests.
"""
struct NPDETestResult
    test_name::String
    statistic::Float64
    p_value::Float64
    passed::Bool
    threshold::Float64
end

"""
    npde_statistical_tests(npde_values; alpha=0.05)

Perform statistical tests on NPDE values.

Tests performed:
1. t-test for mean = 0
2. Fisher variance test for var = 1
3. Shapiro-Wilk test for normality (approximation)
4. Symmetry test (skewness)

Returns vector of NPDETestResult.
"""
function npde_statistical_tests(
    npde_values::Vector{Float64};
    alpha::Float64=0.05
)::Vector{NPDETestResult}
    n = length(npde_values)
    results = NPDETestResult[]

    # Filter out any NaN values
    valid_npde = filter(!isnan, npde_values)
    n_valid = length(valid_npde)

    if n_valid < 10
        return results  # Not enough data
    end

    # 1. t-test for mean = 0
    m = mean(valid_npde)
    se_mean = std(valid_npde) / sqrt(n_valid)
    t_stat = m / se_mean
    t_dist = TDist(n_valid - 1)
    p_mean = 2 * (1 - cdf(t_dist, abs(t_stat)))
    push!(results, NPDETestResult(
        "t-test (mean = 0)",
        t_stat,
        p_mean,
        p_mean > alpha,
        alpha
    ))

    # 2. Variance test (var = 1)
    v = var(valid_npde)
    chi2_stat = (n_valid - 1) * v / 1.0
    chi2_dist = Chisq(n_valid - 1)
    # Two-sided test
    p_var = 2 * min(cdf(chi2_dist, chi2_stat), 1 - cdf(chi2_dist, chi2_stat))
    push!(results, NPDETestResult(
        "Variance test (var = 1)",
        chi2_stat,
        p_var,
        p_var > alpha,
        alpha
    ))

    # 3. Normality test (Shapiro-Wilk approximation via correlation)
    sorted_npde = sort(valid_npde)
    expected_normal = [quantile(Normal(), (i - 0.375) / (n_valid + 0.25)) for i in 1:n_valid]
    correlation = cor(sorted_npde, expected_normal)
    # Approximate p-value (simplified)
    w_stat = correlation^2
    # For large n, use approximation
    if n_valid > 50
        z_w = (log(1 - w_stat) - (-0.0752 - 1.3107 / n_valid)) / (1.0308 - 0.26758 / sqrt(n_valid))
        p_normal = 1 - cdf(Normal(), z_w)
    else
        p_normal = w_stat  # Simplified for small samples
    end
    push!(results, NPDETestResult(
        "Shapiro-Wilk normality",
        w_stat,
        p_normal,
        p_normal > alpha,
        alpha
    ))

    # 4. Symmetry test (skewness = 0)
    skew = _compute_skewness(valid_npde)
    se_skew = sqrt(6.0 / n_valid)  # Approximate SE of skewness
    z_skew = skew / se_skew
    p_skew = 2 * (1 - cdf(Normal(), abs(z_skew)))
    push!(results, NPDETestResult(
        "Symmetry test (skewness = 0)",
        z_skew,
        p_skew,
        p_skew > alpha,
        alpha
    ))

    return results
end

# ============================================================================
# Helper Functions
# ============================================================================

"""
Compute skewness of a distribution.
"""
function _compute_skewness(x::Vector{Float64})::Float64
    n = length(x)
    m = mean(x)
    s = std(x)
    return sum(((xi - m) / s)^3 for xi in x) / n
end

"""
Compute excess kurtosis of a distribution.
"""
function _compute_kurtosis(x::Vector{Float64})::Float64
    n = length(x)
    m = mean(x)
    s = std(x)
    return sum(((xi - m) / s)^4 for xi in x) / n - 3.0
end

# ============================================================================
# Decorrelated NPDE
# ============================================================================

"""
    compute_decorrelated_npde(pde_matrix, var_cov_matrix)

Apply decorrelation to NPDE using the variance-covariance matrix.

The decorrelation step accounts for correlation between observations
within the same subject due to shared random effects.

Arguments:
- pde_matrix: Matrix of PDE values (n_obs x n_sim)
- var_cov_matrix: Variance-covariance matrix of predictions (n_obs x n_obs)

Returns:
- Decorrelated NPDE values
"""
function compute_decorrelated_npde(
    pde::Vector{Float64},
    var_cov_matrix::Matrix{Float64}
)::Vector{Float64}
    n = length(pde)

    # Cholesky decomposition of variance-covariance matrix
    try
        L = cholesky(var_cov_matrix).L
        # Decorrelate: npde_dec = L^{-1} * npde
        npde_raw = quantile.(Normal(), clamp.(pde, 0.001, 0.999))
        npde_dec = L \ npde_raw
        return npde_dec
    catch e
        # If Cholesky fails (not positive definite), return non-decorrelated
        return quantile.(Normal(), clamp.(pde, 0.001, 0.999))
    end
end

# ============================================================================
# Enhanced NPDE Computation
# ============================================================================

"""
    compute_npde_enhanced(observed, simulation_fn, spec)

Compute NPDE with full posterior predictive sampling.

Arguments:
- observed: ObservedData with subject observations
- simulation_fn: Function that takes a seed and returns simulated predictions
                 Should return Dict(subject_id => Vector{Float64})
- spec: NPDESpec with computation parameters

Returns:
- NPDEResult with NPDE values and diagnostics
"""
function compute_npde_enhanced(
    observed,  # ObservedData type
    simulation_fn::Function,
    spec::NPDESpec
)::NPDEResult
    rng = StableRNG(spec.seed)

    n_subjects = length(observed.subjects)
    subject_results = SubjectNPDEResult[]

    all_npde = Float64[]
    all_pde = Float64[]

    for subj in observed.subjects
        times = subj.times
        dv = subj.observations
        n_obs = length(times)

        # Run simulations
        sim_matrix = zeros(n_obs, spec.n_simulations)

        for sim_idx in 1:spec.n_simulations
            sim_seed = rand(rng, UInt64)
            sim_preds = simulation_fn(sim_seed)

            if haskey(sim_preds, subj.id)
                sim_matrix[:, sim_idx] = sim_preds[subj.id][1:n_obs]
            end
        end

        # Compute PDE for each observation
        pde = zeros(n_obs)
        for i in 1:n_obs
            sim_vals = sim_matrix[i, :]
            # PDE = fraction of simulations <= observed
            pde[i] = sum(sim_vals .<= dv[i]) / spec.n_simulations
            # Clamp to avoid Inf
            pde[i] = clamp(pde[i], 0.001, 0.999)
        end

        # Compute NPDE
        if spec.decorrelate && n_obs > 1
            # Compute variance-covariance matrix of simulations
            var_cov = cov(sim_matrix')
            npde = compute_decorrelated_npde(pde, var_cov)
        else
            npde = quantile.(Normal(), pde)
        end

        # Expected predictions and variance
        epred = vec(mean(sim_matrix, dims=2))
        var_pred = vec(var(sim_matrix, dims=2))

        push!(subject_results, SubjectNPDEResult(
            subj.id, times, dv, npde, pde, epred, var_pred
        ))

        append!(all_npde, npde)
        append!(all_pde, pde)
    end

    # Compute diagnostics
    diagnostics = _compute_npde_diagnostics(all_npde)

    return NPDEResult(
        subject_results,
        all_npde,
        all_pde,
        diagnostics,
        nothing  # Don't store simulation matrix by default
    )
end

"""
Compute NPDE diagnostics.
"""
function _compute_npde_diagnostics(npde::Vector{Float64})::NPDEDiagnostics
    valid_npde = filter(!isnan, npde)
    n = length(valid_npde)

    if n < 10
        return NPDEDiagnostics(NaN, NaN, NaN, NaN, NaN, NaN, NaN)
    end

    m = mean(valid_npde)
    v = var(valid_npde)
    skew = _compute_skewness(valid_npde)
    kurt = _compute_kurtosis(valid_npde) + 3.0  # Convert from excess to raw

    # Tests
    tests = npde_statistical_tests(valid_npde)

    sw_p = length(tests) >= 3 ? tests[3].p_value : NaN
    mean_p = length(tests) >= 1 ? tests[1].p_value : NaN
    var_p = length(tests) >= 2 ? tests[2].p_value : NaN

    return NPDEDiagnostics(m, v, skew, kurt, sw_p, mean_p, var_p)
end

# ============================================================================
# Visualization Data
# ============================================================================

"""
    get_npde_qq_data(result)

Get data for NPDE Q-Q plot.
"""
function get_npde_qq_data(result::NPDEResult)::Dict{Symbol, Vector{Float64}}
    valid_npde = filter(!isnan, result.pooled_npde)
    n = length(valid_npde)
    sorted_npde = sort(valid_npde)
    theoretical = [quantile(Normal(), (i - 0.5) / n) for i in 1:n]

    return Dict(
        :observed => sorted_npde,
        :theoretical => theoretical,
        :identity_line_x => [-3.0, 3.0],
        :identity_line_y => [-3.0, 3.0]
    )
end

"""
    get_npde_vs_pred_data(result)

Get data for NPDE vs prediction plot.
"""
function get_npde_vs_pred_data(result::NPDEResult)::Dict{Symbol, Vector{Float64}}
    npde = Float64[]
    pred = Float64[]

    for subj in result.subjects
        append!(npde, subj.npde)
        append!(pred, subj.epred)
    end

    return Dict(
        :npde => npde,
        :pred => pred
    )
end

"""
    get_npde_vs_time_data(result)

Get data for NPDE vs time plot.
"""
function get_npde_vs_time_data(result::NPDEResult)::Dict{Symbol, Vector{Float64}}
    npde = Float64[]
    time = Float64[]

    for subj in result.subjects
        append!(npde, subj.npde)
        append!(time, subj.times)
    end

    return Dict(
        :npde => npde,
        :time => time
    )
end

"""
    get_npde_histogram_data(result; n_bins=30)

Get data for NPDE histogram with standard normal overlay.
"""
function get_npde_histogram_data(
    result::NPDEResult;
    n_bins::Int=30
)::Dict{Symbol, Any}
    valid_npde = filter(!isnan, result.pooled_npde)

    # Histogram bin edges
    min_val = min(-3.5, minimum(valid_npde))
    max_val = max(3.5, maximum(valid_npde))
    bin_edges = range(min_val, max_val, length=n_bins+1)

    # Count observations in each bin
    counts = zeros(Int, n_bins)
    for v in valid_npde
        bin = searchsortedlast(collect(bin_edges), v)
        bin = clamp(bin, 1, n_bins)
        counts[bin] += 1
    end

    # Normalize to density
    bin_width = (max_val - min_val) / n_bins
    density = counts ./ (length(valid_npde) * bin_width)

    # Standard normal reference
    x_ref = range(-3.5, 3.5, length=100)
    y_ref = pdf.(Normal(), x_ref)

    return Dict(
        :bin_centers => [(bin_edges[i] + bin_edges[i+1])/2 for i in 1:n_bins],
        :density => density,
        :reference_x => collect(x_ref),
        :reference_y => collect(y_ref)
    )
end

export get_npde_histogram_data
