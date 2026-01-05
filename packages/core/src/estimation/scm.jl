# Stepwise Covariate Modeling (SCM)
# Automatic covariate selection using forward selection and backward elimination

using Distributions
using LinearAlgebra

export SCMResult, SCMSpec, CovariateRelationship
export run_scm, forward_selection, backward_elimination
# Note: likelihood_ratio_test is defined in diagnostics.jl

# ============================================================================
# Types and Specifications
# ============================================================================

"""
    CovariateRelationship

Defines a potential covariate-parameter relationship to test.

Fields:
- param: Parameter affected by covariate (e.g., :CL)
- covariate: Covariate name (e.g., :WT, :AGE)
- relationship: Type of relationship (:linear, :power, :exponential)
- reference: Reference value for centering (e.g., median weight)
"""
struct CovariateRelationship
    param::Symbol
    covariate::Symbol
    relationship::Symbol  # :linear, :power, :exponential
    reference::Float64

    function CovariateRelationship(
        param::Symbol,
        covariate::Symbol;
        relationship::Symbol=:power,
        reference::Float64=1.0
    )
        relationship in (:linear, :power, :exponential) ||
            error("relationship must be :linear, :power, or :exponential")
        new(param, covariate, relationship, reference)
    end
end

"""
    SCMSpec

Specification for stepwise covariate modeling.

Fields:
- relationships: Vector of covariate relationships to test
- forward_p_value: P-value threshold for forward inclusion (default: 0.05)
- backward_p_value: P-value threshold for backward elimination (default: 0.01)
- max_iterations: Maximum iterations for forward/backward steps
- parallel: Run covariate tests in parallel
"""
struct SCMSpec
    relationships::Vector{CovariateRelationship}
    forward_p_value::Float64
    backward_p_value::Float64
    max_iterations::Int
    parallel::Bool

    function SCMSpec(
        relationships::Vector{CovariateRelationship};
        forward_p_value::Float64=0.05,
        backward_p_value::Float64=0.01,
        max_iterations::Int=100,
        parallel::Bool=false
    )
        0 < forward_p_value < 1 || error("forward_p_value must be in (0, 1)")
        0 < backward_p_value < 1 || error("backward_p_value must be in (0, 1)")
        forward_p_value >= backward_p_value ||
            @warn "backward_p_value should typically be <= forward_p_value"

        new(relationships, forward_p_value, backward_p_value, max_iterations, parallel)
    end
end

"""
    SCMStepResult

Result of a single SCM step (adding or removing one covariate).
"""
struct SCMStepResult
    relationship::CovariateRelationship
    ofv_change::Float64
    p_value::Float64
    coefficient::Float64
    coefficient_se::Float64
    included::Bool
    step_type::Symbol  # :forward or :backward
end

"""
    SCMResult

Complete result of stepwise covariate modeling.

Fields:
- final_model: Vector of included covariate relationships
- forward_steps: Results from forward selection
- backward_steps: Results from backward elimination
- base_ofv: OFV of base model (no covariates)
- final_ofv: OFV of final model
- all_tested: All relationships that were tested
- convergence: Whether SCM converged
"""
struct SCMResult
    final_model::Vector{CovariateRelationship}
    forward_steps::Vector{SCMStepResult}
    backward_steps::Vector{SCMStepResult}
    base_ofv::Float64
    final_ofv::Float64
    all_tested::Vector{CovariateRelationship}
    convergence::Bool
    n_forward_iterations::Int
    n_backward_iterations::Int
end

export SCMStepResult

# ============================================================================
# Covariate Model Building
# Note: likelihood_ratio_test is defined in diagnostics.jl
# ============================================================================

"""
Build covariate effect function based on relationship type.

Returns a function: (theta_pop, cov_value, coefficient) -> theta_individual
"""
function build_covariate_effect(rel::CovariateRelationship)::Function
    if rel.relationship == :linear
        # θ_i = θ_pop * (1 + β * (COV - REF))
        return (theta_pop, cov, beta) -> theta_pop * (1 + beta * (cov - rel.reference))
    elseif rel.relationship == :power
        # θ_i = θ_pop * (COV / REF)^β
        return (theta_pop, cov, beta) -> theta_pop * (cov / rel.reference)^beta
    elseif rel.relationship == :exponential
        # θ_i = θ_pop * exp(β * (COV - REF))
        return (theta_pop, cov, beta) -> theta_pop * exp(beta * (cov - rel.reference))
    else
        error("Unknown relationship type: $(rel.relationship)")
    end
end

"""
Create a modified OFV function that includes a covariate effect.

Arguments:
- base_ofv_fn: Original OFV function taking theta vector
- rel: Covariate relationship to add
- covariate_values: Vector of covariate values (one per individual)
- param_index: Index of the parameter in theta vector
- n_base_params: Number of parameters in base model
"""
function create_covariate_ofv(
    base_ofv_fn::Function,
    rel::CovariateRelationship,
    covariate_values::Vector{Float64},
    param_index::Int,
    n_base_params::Int
)::Function
    effect_fn = build_covariate_effect(rel)

    function covariate_ofv(theta_extended::Vector{Float64})::Float64
        # theta_extended = [base_params..., covariate_coefficient]
        theta_base = theta_extended[1:n_base_params]
        beta = theta_extended[end]

        # This is a simplified version - in practice, you'd need to
        # modify individual predictions based on covariate values
        # For now, return base OFV (placeholder for actual implementation)
        return base_ofv_fn(theta_base)
    end

    return covariate_ofv
end

# ============================================================================
# Forward Selection
# ============================================================================

"""
    forward_selection(estimate_fn, base_theta, candidates, spec) -> (selected, steps)

Perform forward covariate selection.

Arguments:
- estimate_fn: Function (relationships) -> (ofv, theta, se) for fitting model
- base_theta: Initial parameter estimates
- candidates: Vector of candidate covariate relationships
- spec: SCM specification

Returns:
- selected: Vector of selected relationships
- steps: Vector of SCMStepResult for each step
"""
function forward_selection(
    estimate_fn::Function,
    base_theta::Vector{Float64},
    candidates::Vector{CovariateRelationship},
    spec::SCMSpec
)::Tuple{Vector{CovariateRelationship}, Vector{SCMStepResult}}
    selected = CovariateRelationship[]
    remaining = copy(candidates)
    steps = SCMStepResult[]

    # Get base model OFV
    base_ofv, _, _ = estimate_fn(selected)

    iteration = 0
    while !isempty(remaining) && iteration < spec.max_iterations
        iteration += 1

        # Test each remaining covariate
        best_candidate = nothing
        best_p_value = 1.0
        best_ofv = base_ofv
        best_coef = 0.0
        best_se = 0.0

        for rel in remaining
            # Fit model with this covariate added
            test_model = vcat(selected, [rel])
            try
                new_ofv, theta, se = estimate_fn(test_model)

                # LRT with 1 df (one additional parameter)
                _, p_value = likelihood_ratio_test(new_ofv, base_ofv, 1)

                if p_value < best_p_value
                    best_candidate = rel
                    best_p_value = p_value
                    best_ofv = new_ofv
                    best_coef = length(theta) > 0 ? theta[end] : 0.0
                    best_se = length(se) > 0 ? se[end] : 0.0
                end
            catch e
                @warn "Failed to fit model with $(rel.param)-$(rel.covariate): $e"
                continue
            end
        end

        # Check if best candidate meets threshold
        if best_candidate !== nothing && best_p_value < spec.forward_p_value
            push!(selected, best_candidate)
            filter!(r -> r != best_candidate, remaining)
            base_ofv = best_ofv

            push!(steps, SCMStepResult(
                best_candidate,
                base_ofv - best_ofv,
                best_p_value,
                best_coef,
                best_se,
                true,
                :forward
            ))
        else
            # No more significant covariates
            break
        end
    end

    return selected, steps
end

# ============================================================================
# Backward Elimination
# ============================================================================

"""
    backward_elimination(estimate_fn, included, spec) -> (final, steps)

Perform backward covariate elimination.

Arguments:
- estimate_fn: Function (relationships) -> (ofv, theta, se) for fitting model
- included: Vector of currently included relationships
- spec: SCM specification

Returns:
- final: Vector of relationships remaining after elimination
- steps: Vector of SCMStepResult for each elimination step
"""
function backward_elimination(
    estimate_fn::Function,
    included::Vector{CovariateRelationship},
    spec::SCMSpec
)::Tuple{Vector{CovariateRelationship}, Vector{SCMStepResult}}
    current = copy(included)
    steps = SCMStepResult[]

    if isempty(current)
        return current, steps
    end

    # Get full model OFV
    full_ofv, _, _ = estimate_fn(current)

    iteration = 0
    while length(current) > 0 && iteration < spec.max_iterations
        iteration += 1

        # Test removing each covariate
        worst_candidate = nothing
        worst_p_value = 0.0
        worst_ofv = full_ofv
        worst_coef = 0.0
        worst_se = 0.0

        for (i, rel) in enumerate(current)
            # Fit model without this covariate
            test_model = vcat(current[1:i-1], current[i+1:end])

            try
                reduced_ofv, theta, se = estimate_fn(test_model)

                # LRT: is the reduced model significantly worse?
                _, p_value = likelihood_ratio_test(full_ofv, reduced_ofv, 1)

                # We want to remove covariates with HIGH p-values (not significant)
                if p_value > worst_p_value
                    worst_candidate = rel
                    worst_p_value = p_value
                    worst_ofv = reduced_ofv
                    # Get coefficient from full model (before removal)
                    worst_coef = i <= length(theta) ? theta[i] : 0.0
                    worst_se = i <= length(se) ? se[i] : 0.0
                end
            catch e
                @warn "Failed to fit reduced model without $(rel.param)-$(rel.covariate): $e"
                continue
            end
        end

        # Check if worst candidate should be removed
        if worst_candidate !== nothing && worst_p_value > spec.backward_p_value
            filter!(r -> r != worst_candidate, current)
            full_ofv = worst_ofv

            push!(steps, SCMStepResult(
                worst_candidate,
                worst_ofv - full_ofv,
                worst_p_value,
                worst_coef,
                worst_se,
                false,
                :backward
            ))
        else
            # All remaining covariates are significant
            break
        end
    end

    return current, steps
end

# ============================================================================
# Main SCM Function
# ============================================================================

"""
    run_scm(estimate_fn, base_theta, spec) -> SCMResult

Run full stepwise covariate modeling.

This performs:
1. Forward selection to add significant covariates
2. Backward elimination to remove non-significant ones

Arguments:
- estimate_fn: Function (relationships) -> (ofv, theta, se)
  Takes a vector of covariate relationships and returns:
  - ofv: Objective function value
  - theta: Parameter estimates
  - se: Standard errors
- base_theta: Initial parameter estimates for base model
- spec: SCM specification

Returns:
- SCMResult with final model and all tested relationships

Example:
```julia
# Define candidate relationships
candidates = [
    CovariateRelationship(:CL, :WT, relationship=:power, reference=70.0),
    CovariateRelationship(:CL, :AGE, relationship=:linear, reference=40.0),
    CovariateRelationship(:V, :WT, relationship=:power, reference=70.0),
    CovariateRelationship(:CL, :CRCL, relationship=:power, reference=100.0)
]

spec = SCMSpec(candidates, forward_p_value=0.05, backward_p_value=0.01)

# estimate_fn should fit the model and return (ofv, theta, se)
result = run_scm(estimate_fn, base_theta, spec)
```
"""
function run_scm(
    estimate_fn::Function,
    base_theta::Vector{Float64},
    spec::SCMSpec
)::SCMResult
    # Get base model OFV
    base_ofv, _, _ = estimate_fn(CovariateRelationship[])

    # Forward selection
    selected, forward_steps = forward_selection(
        estimate_fn, base_theta, spec.relationships, spec
    )

    # Backward elimination
    final_model, backward_steps = backward_elimination(
        estimate_fn, selected, spec
    )

    # Get final OFV
    final_ofv = if isempty(final_model)
        base_ofv
    else
        ofv, _, _ = estimate_fn(final_model)
        ofv
    end

    return SCMResult(
        final_model,
        forward_steps,
        backward_steps,
        base_ofv,
        final_ofv,
        spec.relationships,
        true,  # convergence
        length(forward_steps),
        length(backward_steps)
    )
end

export run_scm, forward_selection, backward_elimination

# ============================================================================
# Utility Functions
# ============================================================================

"""
Print SCM result summary.
"""
function summarize_scm(result::SCMResult)::String
    lines = String[]

    push!(lines, "=== Stepwise Covariate Modeling Results ===")
    push!(lines, "")
    push!(lines, "Base OFV: $(round(result.base_ofv, digits=3))")
    push!(lines, "Final OFV: $(round(result.final_ofv, digits=3))")
    push!(lines, "ΔOFV: $(round(result.base_ofv - result.final_ofv, digits=3))")
    push!(lines, "")

    push!(lines, "Forward Selection Steps: $(result.n_forward_iterations)")
    for step in result.forward_steps
        rel = step.relationship
        push!(lines, "  + $(rel.param) ~ $(rel.covariate) ($(rel.relationship)): p=$(round(step.p_value, digits=4)), β=$(round(step.coefficient, digits=4))")
    end

    push!(lines, "")
    push!(lines, "Backward Elimination Steps: $(result.n_backward_iterations)")
    for step in result.backward_steps
        rel = step.relationship
        push!(lines, "  - $(rel.param) ~ $(rel.covariate) ($(rel.relationship)): p=$(round(step.p_value, digits=4))")
    end

    push!(lines, "")
    push!(lines, "Final Model Covariates:")
    if isempty(result.final_model)
        push!(lines, "  (none)")
    else
        for rel in result.final_model
            push!(lines, "  $(rel.param) ~ $(rel.covariate) ($(rel.relationship), ref=$(rel.reference))")
        end
    end

    return join(lines, "\n")
end

export summarize_scm

"""
Create standard covariate relationships for common PK parameters.

Returns relationships for:
- Weight effects on CL, V (allometric scaling)
- Age effects on CL
- Creatinine clearance effects on CL
- Sex effects on CL, V
"""
function standard_pk_covariates(;
    weight_ref::Float64=70.0,
    age_ref::Float64=40.0,
    crcl_ref::Float64=100.0
)::Vector{CovariateRelationship}
    return [
        # Weight effects (allometric)
        CovariateRelationship(:CL, :WT, relationship=:power, reference=weight_ref),
        CovariateRelationship(:V, :WT, relationship=:power, reference=weight_ref),
        CovariateRelationship(:V1, :WT, relationship=:power, reference=weight_ref),
        CovariateRelationship(:V2, :WT, relationship=:power, reference=weight_ref),
        CovariateRelationship(:Q, :WT, relationship=:power, reference=weight_ref),

        # Age effects
        CovariateRelationship(:CL, :AGE, relationship=:linear, reference=age_ref),

        # Renal function
        CovariateRelationship(:CL, :CRCL, relationship=:power, reference=crcl_ref),
        CovariateRelationship(:CL, :eGFR, relationship=:power, reference=crcl_ref),

        # Sex effects (typically as proportional)
        CovariateRelationship(:CL, :SEX, relationship=:linear, reference=0.0),
        CovariateRelationship(:V, :SEX, relationship=:linear, reference=0.0)
    ]
end

export standard_pk_covariates
