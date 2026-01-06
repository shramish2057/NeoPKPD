# Adaptive Trial Simulation
# Implements interim analyses, stopping rules, and conditional power

using StableRNGs
using SpecialFunctions: erf, erfinv

export simulate_adaptive_trial, AdaptiveTrialResult, InterimResult
export check_efficacy_stop, check_futility_stop, compute_conditional_power
export get_efficacy_boundary, get_alpha_spending

# ============================================================================
# Result Types
# ============================================================================

"""
    InterimResult

Results from a single interim analysis.

# Fields
- `analysis_number`: Which interim analysis (1, 2, ...)
- `information_fraction`: Fraction of total planned information
- `n_enrolled`: Number of subjects enrolled at this interim
- `primary_effect`: Estimated treatment effect
- `primary_se`: Standard error of effect estimate
- `primary_z_stat`: Z-statistic for primary endpoint
- `p_value`: P-value for primary endpoint
- `efficacy_boundary`: Efficacy stopping boundary at this interim
- `futility_boundary`: Futility stopping boundary
- `conditional_power`: Conditional power estimate
- `stop_for_efficacy`: Whether stopped for efficacy
- `stop_for_futility`: Whether stopped for futility
"""
struct InterimResult
    analysis_number::Int
    information_fraction::Float64
    n_enrolled::Int
    primary_effect::Float64
    primary_se::Float64
    primary_z_stat::Float64
    p_value::Float64
    efficacy_boundary::Float64
    futility_boundary::Union{Nothing, Float64}
    conditional_power::Float64
    stop_for_efficacy::Bool
    stop_for_futility::Bool
end

"""
    AdaptiveTrialResult

Complete results from adaptive trial simulation.

# Fields
- `final_result`: Final endpoint analysis
- `interim_results`: Vector of InterimResult from each interim
- `stopped_early`: Whether trial stopped before planned completion
- `stop_reason`: Reason for early stopping (:efficacy, :futility, :none)
- `final_n`: Final number of subjects enrolled
- `initial_n`: Planned sample size
- `seed`: Random seed used
"""
struct AdaptiveTrialResult
    final_result::Dict{Symbol, Any}
    interim_results::Vector{InterimResult}
    stopped_early::Bool
    stop_reason::Union{Nothing, Symbol}
    final_n::Int
    initial_n::Int
    seed::UInt64
end

# ============================================================================
# Alpha Spending Functions
# ============================================================================

"""
    get_alpha_spending(spending_function::Symbol, info_fraction::Float64, alpha::Float64)

Compute cumulative alpha spent at a given information fraction.

# Arguments
- `spending_function`: Type of spending function (:obrien_fleming, :pocock, :haybittle_peto)
- `info_fraction`: Information fraction (0 to 1)
- `alpha`: Overall alpha level (typically 0.025 one-sided)

# Returns
- Cumulative alpha spent at this information fraction
"""
function get_alpha_spending(spending_function::Symbol, info_fraction::Float64, alpha::Float64)
    t = info_fraction

    if spending_function == :obrien_fleming
        # O'Brien-Fleming: α(t) = 1 - Φ(z_α/√t) for one-sided alpha
        # Very conservative early, spends most alpha at end
        z_alpha = sqrt(2) * erfinv(1 - 2*alpha)
        return 1 - 0.5 * (1 + erf(z_alpha / sqrt(t) / sqrt(2)))

    elseif spending_function == :pocock
        # Pocock: α(t) = α * log(1 + (e-1)*t)
        # Constant boundaries, more aggressive early
        return alpha * log(1 + (exp(1) - 1) * t)

    elseif spending_function == :haybittle_peto
        # Haybittle-Peto: Use z=3 for interim, full alpha at final
        # Very conservative for early analyses
        if t < 1.0
            return 0.001  # Approximately 2*(1-Φ(3))
        else
            return alpha
        end

    elseif spending_function == :linear
        # Linear spending (Kim-DeMets family with rho=1)
        return alpha * t

    elseif spending_function == :none
        # No alpha spending - all at final
        if t >= 1.0
            return alpha
        else
            return 0.0
        end

    else
        # Default to O'Brien-Fleming
        z_alpha = sqrt(2) * erfinv(1 - 2*alpha)
        return 1 - 0.5 * (1 + erf(z_alpha / sqrt(t) / sqrt(2)))
    end
end

"""
    get_efficacy_boundary(spending_function::Symbol, info_fraction::Float64, alpha::Float64)

Get the efficacy stopping boundary (critical z-value) at a given information fraction.

# Returns
- Z-value threshold for efficacy stopping
"""
function get_efficacy_boundary(spending_function::Symbol, info_fraction::Float64, alpha::Float64)
    # Compute incremental alpha for this look
    alpha_spent = get_alpha_spending(spending_function, info_fraction, alpha)

    # Convert to z-value
    if alpha_spent <= 0
        return Inf  # No stopping at this interim
    elseif alpha_spent >= 1
        return 0.0
    else
        # Z-value corresponding to this alpha (one-sided)
        return sqrt(2) * erfinv(1 - 2*alpha_spent)
    end
end

# ============================================================================
# Stopping Rules
# ============================================================================

"""
    check_efficacy_stop(z_stat::Float64, boundary::Float64)

Check if trial should stop for efficacy.

# Arguments
- `z_stat`: Observed z-statistic
- `boundary`: Efficacy boundary (z-value)

# Returns
- Boolean: true if should stop for efficacy
"""
function check_efficacy_stop(z_stat::Float64, boundary::Float64)::Bool
    return z_stat > boundary
end

"""
    check_futility_stop(conditional_power::Float64, threshold::Float64)

Check if trial should stop for futility based on conditional power.

# Arguments
- `conditional_power`: Estimated conditional power
- `threshold`: Futility threshold (typically 0.10 to 0.20)

# Returns
- Boolean: true if should stop for futility
"""
function check_futility_stop(conditional_power::Float64, threshold::Float64)::Bool
    return conditional_power < threshold
end

"""
    compute_conditional_power(effect::Float64, se::Float64, info_frac::Float64,
                               target_effect::Float64; alpha=0.025)

Compute conditional power given interim results.

Conditional power is the probability of reaching significance at the end
of the trial, given the current data and assuming a specific true effect.

# Arguments
- `effect`: Observed treatment effect at interim
- `se`: Standard error of effect at interim
- `info_frac`: Information fraction at interim
- `target_effect`: Assumed true effect for power calculation
- `alpha`: Significance level (one-sided)

# Returns
- Conditional power (0 to 1)
"""
function compute_conditional_power(
    effect::Float64,
    se::Float64,
    info_frac::Float64,
    target_effect::Float64;
    alpha::Float64 = 0.025
)::Float64

    if info_frac >= 1.0 || info_frac <= 0.0 || se <= 0
        return NaN
    end

    # Information at interim
    info_interim = 1 / se^2

    # Projected information at final
    info_final = info_interim / info_frac

    # Information remaining
    info_remaining = info_final - info_interim

    # Z-statistic at interim
    z_interim = effect / se

    # Z-value required for final significance
    z_crit = sqrt(2) * erfinv(1 - 2*alpha)

    # Under the alternative (true effect = target_effect):
    # Final Z ~ N(target_effect * sqrt(info_final), 1)
    # Conditional on interim data, remaining contribution is:
    # target_effect * sqrt(info_remaining) + noise

    # Predicted final z-statistic (weighted combination)
    # z_final = (z_interim * sqrt(info_interim) + z_remaining * sqrt(info_remaining)) / sqrt(info_final)

    # Expected z_remaining under alternative
    expected_z_remaining = target_effect * sqrt(info_remaining)

    # Predicted final z
    predicted_z_final = (z_interim * sqrt(info_interim) + expected_z_remaining) / sqrt(info_final)

    # Variance of remaining contribution (under alternative)
    var_remaining = 1.0

    # Conditional power
    # P(Z_final > z_crit | interim data, theta = target_effect)
    # = P(Z_remaining > required_remaining)
    # where required_remaining = (z_crit * sqrt(info_final) - z_interim * sqrt(info_interim)) / sqrt(info_remaining)

    required_remaining = (z_crit * sqrt(info_final) - z_interim * sqrt(info_interim)) / sqrt(info_remaining)

    # Under H1: Z_remaining ~ N(target_effect * sqrt(info_remaining), 1)
    mean_remaining = target_effect * sqrt(info_remaining)

    # P(Z_remaining > required_remaining) = 1 - Φ((required_remaining - mean_remaining) / 1)
    cp = 1 - 0.5 * (1 + erf((required_remaining - mean_remaining) / sqrt(2)))

    return max(0.0, min(1.0, cp))
end

# ============================================================================
# Adaptive Trial Simulation
# ============================================================================

"""
    simulate_adaptive_trial(trial_spec::TrialSpec, design::AdaptiveDesign;
                            grid=nothing, solver=nothing)

Simulate an adaptive clinical trial with interim analyses and stopping rules.

# Arguments
- `trial_spec`: Trial specification
- `design`: AdaptiveDesign with interim analysis schedule and rules

# Returns
- `AdaptiveTrialResult` with complete results including interim analyses

# Example
```julia
base_design = ParallelDesign(2)
adaptive = AdaptiveDesign(base_design;
    interim_analyses = [0.5],  # One interim at 50%
    alpha_spending = :obrien_fleming
)
result = simulate_adaptive_trial(trial_spec, adaptive)
```
"""
function simulate_adaptive_trial(
    trial_spec::TrialSpec,
    design::AdaptiveDesign;
    grid = nothing,
    solver = nothing,
    target_effect::Float64 = 0.0,  # For conditional power (0 = use observed)
    futility_threshold::Float64 = 0.10
)::AdaptiveTrialResult

    rng = StableRNG(trial_spec.seed)

    # Validate inputs
    if isempty(trial_spec.arms)
        error("Trial must have at least one arm")
    end

    # Total planned sample size
    planned_n = sum(arm.n_subjects for arm in trial_spec.arms)

    # Interim analysis schedule
    interim_fractions = sort(design.interim_analyses)
    push!(interim_fractions, 1.0)  # Add final analysis

    # Get alpha level (default 0.025 one-sided)
    alpha = get(design.adaptation_rules, :alpha, 0.025)

    # Futility settings
    futility_enabled = get(design.adaptation_rules, :futility_enabled, true)

    # Generate full virtual population upfront
    population = generate_virtual_population(trial_spec.virtual_population, planned_n)

    # Cumulative data storage
    cumulative_results = Dict{String, Vector{Any}}()
    cumulative_values = Dict{String, Dict{Symbol, Vector{Float64}}}()

    for arm in trial_spec.arms
        cumulative_results[arm.name] = []
        cumulative_values[arm.name] = Dict{Symbol, Vector{Float64}}()
        for ep in trial_spec.endpoints
            cumulative_values[arm.name][ep.name] = Float64[]
        end
    end

    # Track interim results
    interim_results = InterimResult[]
    stopped_early = false
    stop_reason = nothing
    current_n = 0

    # Process each interim analysis point
    for (analysis_num, target_frac) in enumerate(interim_fractions)
        target_n = round(Int, target_frac * planned_n)
        n_to_enroll = target_n - current_n

        if n_to_enroll <= 0
            continue
        end

        # Enroll subjects up to this target
        for arm in trial_spec.arms
            arm_target = round(Int, arm.n_subjects * target_frac)
            arm_current = length(cumulative_results[arm.name])
            arm_n_to_enroll = arm_target - arm_current

            for i in 1:arm_n_to_enroll
                # Get subject from population
                subj_idx = arm_current + i
                if subj_idx > length(population)
                    break
                end
                subject = population[subj_idx]

                # Simulate individual
                result = simulate_individual(
                    arm, subject, trial_spec.duration_days,
                    trial_spec.pk_sampling_times;
                    rng = rng, grid = grid, solver = solver
                )

                push!(cumulative_results[arm.name], result)

                # Calculate endpoints
                for ep in trial_spec.endpoints
                    value = calculate_individual_endpoint(result, ep, subject)
                    push!(cumulative_values[arm.name][ep.name], value)
                end
            end
        end

        current_n = sum(length(cumulative_results[arm.name]) for arm in trial_spec.arms)

        # Skip interim analysis at final (just complete enrollment)
        if target_frac >= 1.0
            continue
        end

        # Perform interim analysis on primary endpoint
        if !isempty(trial_spec.endpoints)
            primary_ep = trial_spec.endpoints[1]

            # Get arm values
            arm_names = [arm.name for arm in trial_spec.arms]
            if length(arm_names) >= 2
                v1 = cumulative_values[arm_names[1]][primary_ep.name]
                v2 = cumulative_values[arm_names[2]][primary_ep.name]

                # Filter NaN
                v1_valid = filter(!isnan, v1)
                v2_valid = filter(!isnan, v2)

                if length(v1_valid) >= 2 && length(v2_valid) >= 2
                    # Compute effect and SE
                    n1, n2 = length(v1_valid), length(v2_valid)
                    mean1, mean2 = sum(v1_valid)/n1, sum(v2_valid)/n2
                    var1 = sum((v1_valid .- mean1).^2) / (n1 - 1)
                    var2 = sum((v2_valid .- mean2).^2) / (n2 - 1)

                    effect = mean2 - mean1
                    se = sqrt(var1/n1 + var2/n2)
                    z_stat = effect / se

                    # Two-tailed p-value
                    p_value = 2 * (1 - 0.5 * (1 + erf(abs(z_stat) / sqrt(2))))

                    # Get boundaries
                    info_frac = target_frac
                    efficacy_boundary = get_efficacy_boundary(design.alpha_spending, info_frac, alpha)

                    # Conditional power
                    cp_target = target_effect > 0 ? target_effect : effect
                    cp = compute_conditional_power(effect, se, info_frac, cp_target; alpha = alpha)

                    # Check stopping rules
                    stop_eff = check_efficacy_stop(abs(z_stat), efficacy_boundary)
                    stop_fut = futility_enabled && check_futility_stop(cp, futility_threshold)

                    # Record interim result
                    push!(interim_results, InterimResult(
                        analysis_num,
                        info_frac,
                        current_n,
                        effect,
                        se,
                        z_stat,
                        p_value,
                        efficacy_boundary,
                        futility_enabled ? futility_threshold : nothing,
                        cp,
                        stop_eff,
                        stop_fut
                    ))

                    # Execute stopping
                    if stop_eff
                        stopped_early = true
                        stop_reason = :efficacy
                        break
                    elseif stop_fut
                        stopped_early = true
                        stop_reason = :futility
                        break
                    end
                end
            end
        end
    end

    # Build final results
    arm_results = Dict{String, ArmResult}()
    for arm in trial_spec.arms
        n_enrolled = length(cumulative_results[arm.name])

        # Build summary stats
        summary_stats = Dict{Symbol, Dict{Symbol, Float64}}()
        for ep in trial_spec.endpoints
            values = cumulative_values[arm.name][ep.name]
            valid = filter(!isnan, values)
            if !isempty(valid)
                mean_val = sum(valid) / length(valid)
                var_val = length(valid) > 1 ? sum((valid .- mean_val).^2) / (length(valid) - 1) : 0.0
                summary_stats[ep.name] = Dict{Symbol, Float64}(
                    :mean => mean_val,
                    :sd => sqrt(var_val),
                    :n => Float64(length(valid))
                )
            end
        end

        arm_results[arm.name] = ArmResult(
            arm.name,
            n_enrolled,
            n_enrolled,  # Assume no dropout for now
            0,
            cumulative_results[arm.name],
            cumulative_values[arm.name],
            summary_stats
        )
    end

    # Final endpoint analysis
    final_analyses = Dict{Symbol, Dict{Symbol, Any}}()
    for ep in trial_spec.endpoints
        final_analyses[ep.name] = analyze_trial_endpoint(arm_results, ep, design.base_design)
    end

    return AdaptiveTrialResult(
        final_analyses,
        interim_results,
        stopped_early,
        stop_reason,
        current_n,
        planned_n,
        trial_spec.seed
    )
end
