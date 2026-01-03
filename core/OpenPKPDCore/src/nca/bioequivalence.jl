# Bioequivalence Analysis
# FDA/EMA compliant bioequivalence assessment

export bioequivalence_90ci, tost_analysis, be_conclusion, create_be_result
export geometric_mean_ratio, geometric_mean, within_subject_cv

# =============================================================================
# Geometric Mean Ratio
# =============================================================================

"""
    geometric_mean_ratio(test_values, reference_values)

Calculate geometric mean ratio (GMR) of test to reference.

GMR = exp(mean(log(test)) - mean(log(reference)))

# Arguments
- `test_values::Vector{Float64}`: Test formulation values (Cmax or AUC)
- `reference_values::Vector{Float64}`: Reference formulation values

# Returns
- `Float64`: Geometric mean ratio (test/reference)
"""
function geometric_mean_ratio(test_values::Vector{Float64}, reference_values::Vector{Float64})
    @assert !isempty(test_values) "Test values cannot be empty"
    @assert !isempty(reference_values) "Reference values cannot be empty"
    @assert all(v -> v > 0, test_values) "All test values must be positive"
    @assert all(v -> v > 0, reference_values) "All reference values must be positive"

    log_test_mean = sum(log.(test_values)) / length(test_values)
    log_ref_mean = sum(log.(reference_values)) / length(reference_values)

    return exp(log_test_mean - log_ref_mean)
end

"""
    geometric_mean(values)

Calculate geometric mean.

GM = exp(mean(log(values)))

# Arguments
- `values::Vector{Float64}`: Values to average

# Returns
- `Float64`: Geometric mean
"""
function geometric_mean(values::Vector{Float64})
    @assert !isempty(values) "Values cannot be empty"
    @assert all(v -> v > 0, values) "All values must be positive"

    return exp(sum(log.(values)) / length(values))
end

# =============================================================================
# Within-Subject CV
# =============================================================================

"""
    within_subject_cv(test_values, reference_values)

Estimate within-subject coefficient of variation from crossover data.

CV_intra = sqrt(exp(MSE) - 1) × 100%

Where MSE is from ANOVA on log-transformed data.

# Arguments
- `test_values::Vector{Float64}`: Test values (paired)
- `reference_values::Vector{Float64}`: Reference values (paired)

# Returns
- `Float64`: Within-subject CV (%)
"""
function within_subject_cv(test_values::Vector{Float64}, reference_values::Vector{Float64})
    @assert length(test_values) == length(reference_values) "Vectors must have same length (paired data)"
    n = length(test_values)
    @assert n >= 2 "Need at least 2 subjects"

    # Log-transform
    log_test = log.(test_values)
    log_ref = log.(reference_values)

    # Calculate within-subject variance
    # For 2-period crossover: σ²_w = Σ(d_i - mean(d))² / (2(n-1))
    # where d_i = log(test_i) - log(ref_i)

    diffs = log_test .- log_ref
    mean_diff = sum(diffs) / n
    ss_within = sum((diffs .- mean_diff).^2)

    # MSE for crossover
    mse = ss_within / (n - 1)

    # CV = sqrt(exp(MSE) - 1) × 100
    cv = sqrt(exp(mse) - 1.0) * 100.0

    return cv
end

# =============================================================================
# 90% Confidence Interval
# =============================================================================

"""
    bioequivalence_90ci(test_values, reference_values; log_transform=true)

Calculate 90% confidence interval for geometric mean ratio.

Uses the standard two-sequence, two-period crossover analysis.

# Arguments
- `test_values::Vector{Float64}`: Test formulation values (paired with reference)
- `reference_values::Vector{Float64}`: Reference formulation values
- `log_transform::Bool`: Apply log transformation (default: true)

# Returns
- `NamedTuple`: (gmr, ci_lower, ci_upper, cv_intra, n)
"""
function bioequivalence_90ci(
    test_values::Vector{Float64},
    reference_values::Vector{Float64};
    log_transform::Bool = true
)
    @assert length(test_values) == length(reference_values) "Vectors must have same length (paired data)"
    n = length(test_values)
    @assert n >= 3 "Need at least 3 subjects for CI calculation"

    if log_transform
        log_test = log.(test_values)
        log_ref = log.(reference_values)
    else
        log_test = test_values
        log_ref = reference_values
    end

    # Calculate difference for each subject
    diffs = log_test .- log_ref
    mean_diff = sum(diffs) / n

    # Standard error of the mean difference
    ss_within = sum((diffs .- mean_diff).^2)
    mse = ss_within / (n - 1)
    se = sqrt(2.0 * mse / n)  # SE for crossover design

    # t-critical value for 90% CI (two-sided)
    # df = n - 1 for simple analysis
    t_crit = _t_critical(n - 1, 0.10)

    # CI on log scale
    log_ci_lower = mean_diff - t_crit * se
    log_ci_upper = mean_diff + t_crit * se

    if log_transform
        # Back-transform to ratio scale
        gmr = exp(mean_diff)
        ci_lower = exp(log_ci_lower)
        ci_upper = exp(log_ci_upper)
        cv_intra = sqrt(exp(mse) - 1.0) * 100.0
    else
        gmr = mean_diff
        ci_lower = log_ci_lower
        ci_upper = log_ci_upper
        cv_intra = sqrt(mse) * 100.0 / abs(mean_diff)  # Approximate CV
    end

    return (gmr=gmr, ci_lower=ci_lower, ci_upper=ci_upper, cv_intra=cv_intra, n=n)
end

"""
    _t_critical(df, alpha)

Calculate critical t-value for two-sided confidence interval.

Uses approximation for t-distribution quantile.
"""
function _t_critical(df::Int, alpha::Float64)
    # Use Cornish-Fisher expansion for approximation
    # For high df, approaches normal distribution

    p = 1.0 - alpha / 2.0  # Upper tail probability

    # Normal quantile approximation
    z = _normal_quantile(p)

    if df >= 30
        return z
    end

    # t-distribution correction for small df
    # Using Wilson-Hilferty approximation
    g1 = 1.0 / df
    g2 = 1.0 / (df^2)

    t = z + (z^3 + z) * g1 / 4.0 +
        (5.0 * z^5 + 16.0 * z^3 + 3.0 * z) * g2 / 96.0

    return t
end

"""
    _normal_quantile(p)

Calculate normal distribution quantile using Abramowitz-Stegun approximation.
"""
function _normal_quantile(p::Float64)
    @assert 0.0 < p < 1.0 "Probability must be between 0 and 1"

    if p == 0.5
        return 0.0
    end

    if p > 0.5
        sign = 1.0
        p_adj = 1.0 - p
    else
        sign = -1.0
        p_adj = p
    end

    t = sqrt(-2.0 * log(p_adj))

    # Coefficients for rational approximation
    c0 = 2.515517
    c1 = 0.802853
    c2 = 0.010328
    d1 = 1.432788
    d2 = 0.189269
    d3 = 0.001308

    z = t - (c0 + c1*t + c2*t^2) / (1.0 + d1*t + d2*t^2 + d3*t^3)

    return sign * z
end

# =============================================================================
# TOST Analysis
# =============================================================================

"""
    tost_analysis(test_values, reference_values; theta_lower=0.80, theta_upper=1.25, alpha=0.05)

Perform Two One-Sided Tests (TOST) procedure for bioequivalence.

Tests:
- H01: μT/μR ≤ θL (lower bound)
- H02: μT/μR ≥ θU (upper bound)

BE is concluded if both null hypotheses are rejected.

# Arguments
- `test_values::Vector{Float64}`: Test formulation values
- `reference_values::Vector{Float64}`: Reference formulation values
- `theta_lower::Float64`: Lower equivalence bound (default: 0.80)
- `theta_upper::Float64`: Upper equivalence bound (default: 1.25)
- `alpha::Float64`: Significance level (default: 0.05)

# Returns
- `TOSTResult`: Complete TOST analysis result
"""
function tost_analysis(
    test_values::Vector{Float64},
    reference_values::Vector{Float64};
    theta_lower::Float64 = 0.80,
    theta_upper::Float64 = 1.25,
    alpha::Float64 = 0.05
)
    @assert length(test_values) == length(reference_values) "Vectors must have same length"
    n = length(test_values)
    @assert n >= 3 "Need at least 3 subjects"

    # Log-transform
    log_test = log.(test_values)
    log_ref = log.(reference_values)

    # Calculate difference
    diffs = log_test .- log_ref
    mean_diff = sum(diffs) / n

    # Standard error
    ss_within = sum((diffs .- mean_diff).^2)
    mse = ss_within / (n - 1)
    se = sqrt(2.0 * mse / n)

    # Log of bounds
    log_theta_lower = log(theta_lower)
    log_theta_upper = log(theta_upper)

    # T-statistics for TOST
    t_lower = (mean_diff - log_theta_lower) / se  # Test against lower bound
    t_upper = (log_theta_upper - mean_diff) / se  # Test against upper bound

    # Degrees of freedom
    df = n - 1

    # P-values (one-sided)
    p_lower = _t_cdf_upper(t_lower, df)
    p_upper = _t_cdf_upper(t_upper, df)

    # Rejection decisions
    reject_lower = p_lower < alpha
    reject_upper = p_upper < alpha

    # BE conclusion
    if reject_lower && reject_upper
        conclusion = :bioequivalent
    else
        conclusion = :not_bioequivalent
    end

    return TOSTResult(
        :generic,
        t_lower,
        t_upper,
        p_lower,
        p_upper,
        reject_lower,
        reject_upper,
        conclusion
    )
end

"""
    _t_cdf_upper(t, df)

Calculate upper tail probability of t-distribution.
Uses numerical approximation.
"""
function _t_cdf_upper(t::Float64, df::Int)
    # For large df, use normal approximation
    if df >= 30
        return _normal_cdf_upper(t)
    end

    # Beta function approximation for t-distribution CDF
    x = df / (df + t^2)

    if t >= 0
        return 0.5 * _incomplete_beta(df/2.0, 0.5, x)
    else
        return 1.0 - 0.5 * _incomplete_beta(df/2.0, 0.5, x)
    end
end

"""
    _normal_cdf_upper(z)

Calculate upper tail probability of standard normal distribution.
"""
function _normal_cdf_upper(z::Float64)
    # Approximation using error function relation
    return 0.5 * (1.0 - erf(z / sqrt(2.0)))
end

"""
    _incomplete_beta(a, b, x)

Regularized incomplete beta function approximation.
"""
function _incomplete_beta(a::Float64, b::Float64, x::Float64)
    if x <= 0.0
        return 0.0
    elseif x >= 1.0
        return 1.0
    end

    # Use continued fraction approximation
    # Simplified version for t-distribution case

    # For t-dist: a = df/2, b = 0.5
    # Use series expansion for small x
    if x < (a + 1.0) / (a + b + 2.0)
        return _beta_cf(a, b, x) * x^a * (1-x)^b / (a * _beta(a, b))
    else
        return 1.0 - _beta_cf(b, a, 1-x) * (1-x)^b * x^a / (b * _beta(a, b))
    end
end

"""
    _lgamma(x)

Log gamma function approximation using Stirling's formula.
"""
function _lgamma(x::Float64)
    # Lanczos approximation for log gamma
    # This is accurate to ~15 digits for x > 0

    if x <= 0.0
        error("lgamma requires positive argument")
    end

    # Coefficients for Lanczos approximation
    g = 7.0
    c = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    ]

    if x < 0.5
        # Use reflection formula
        return log(pi / sin(pi * x)) - _lgamma(1.0 - x)
    end

    x -= 1.0

    a = c[1]
    for i in 2:9
        a += c[i] / (x + Float64(i - 1))
    end

    t = x + g + 0.5
    return 0.5 * log(2.0 * pi) + (x + 0.5) * log(t) - t + log(a)
end

"""
    _beta(a, b)

Beta function approximation.
"""
function _beta(a::Float64, b::Float64)
    return exp(_lgamma(a) + _lgamma(b) - _lgamma(a + b))
end

"""
    _beta_cf(a, b, x)

Continued fraction for incomplete beta.
"""
function _beta_cf(a::Float64, b::Float64, x::Float64)
    max_iter = 100
    eps = 1e-10

    qab = a + b
    qap = a + 1.0
    qam = a - 1.0

    c = 1.0
    d = 1.0 - qab * x / qap

    if abs(d) < eps
        d = eps
    end
    d = 1.0 / d
    h = d

    for m in 1:max_iter
        m2 = 2 * m

        # Even step
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < eps
            d = eps
        end
        c = 1.0 + aa / c
        if abs(c) < eps
            c = eps
        end
        d = 1.0 / d
        h *= d * c

        # Odd step
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < eps
            d = eps
        end
        c = 1.0 + aa / c
        if abs(c) < eps
            c = eps
        end
        d = 1.0 / d
        del = d * c
        h *= del

        if abs(del - 1.0) < eps
            break
        end
    end

    return h
end

# =============================================================================
# BE Conclusion
# =============================================================================

"""
    be_conclusion(ci_lower, ci_upper; theta_lower=0.80, theta_upper=1.25)

Determine bioequivalence conclusion from confidence interval.

# Arguments
- `ci_lower::Float64`: Lower bound of 90% CI
- `ci_upper::Float64`: Upper bound of 90% CI
- `theta_lower::Float64`: Lower equivalence bound (default: 0.80)
- `theta_upper::Float64`: Upper equivalence bound (default: 1.25)

# Returns
- `Symbol`: :bioequivalent, :not_bioequivalent, or :inconclusive
"""
function be_conclusion(
    ci_lower::Float64,
    ci_upper::Float64;
    theta_lower::Float64 = 0.80,
    theta_upper::Float64 = 1.25
)
    if ci_lower >= theta_lower && ci_upper <= theta_upper
        return :bioequivalent
    elseif ci_upper < theta_lower || ci_lower > theta_upper
        return :not_bioequivalent
    else
        return :inconclusive
    end
end

"""
    create_be_result(parameter, test_values, reference_values; be_limits=(0.80, 1.25))

Create complete BioequivalenceResult.

# Arguments
- `parameter::Symbol`: Parameter analyzed (:cmax, :auc, etc.)
- `test_values::Vector{Float64}`: Test formulation values
- `reference_values::Vector{Float64}`: Reference formulation values
- `be_limits::Tuple{Float64,Float64}`: BE acceptance limits

# Returns
- `BioequivalenceResult`: Complete BE analysis result
"""
function create_be_result(
    parameter::Symbol,
    test_values::Vector{Float64},
    reference_values::Vector{Float64};
    be_limits::Tuple{Float64,Float64} = (0.80, 1.25)
)
    ci_result = bioequivalence_90ci(test_values, reference_values)

    conclusion = be_conclusion(
        ci_result.ci_lower,
        ci_result.ci_upper;
        theta_lower=be_limits[1],
        theta_upper=be_limits[2]
    )

    return BioequivalenceResult(
        parameter,
        length(test_values),
        length(reference_values),
        ci_result.gmr,
        ci_result.ci_lower,
        ci_result.ci_upper,
        ci_result.cv_intra,
        conclusion,
        be_limits
    )
end
