# Lambda-z Terminal Slope Estimation
# FDA/EMA compliant terminal elimination rate constant estimation

export estimate_lambda_z

"""
    estimate_lambda_z(t, c, config; exclude_indices=Int[], tmax_idx=nothing)

Estimate the terminal elimination rate constant (λz) using log-linear regression.

# Algorithm
1. Identify Cmax and its time (Tmax)
2. Consider only points after Tmax with positive concentrations
3. Iterate over different starting points (minimum points to last)
4. For each candidate range, perform log-linear regression
5. Select range with highest adjusted R² meeting threshold
6. Validate span covers minimum half-lives requirement

# Arguments
- `t::Vector{Float64}`: Time points (sorted, ascending)
- `c::Vector{Float64}`: Concentration values
- `config::NCAConfig`: NCA configuration

# Keyword Arguments
- `exclude_indices::Vector{Int}`: Indices to exclude from regression
- `tmax_idx::Union{Int,Nothing}`: Pre-computed Tmax index (optional)

# Returns
- `LambdaZResult`: Complete lambda_z estimation result
"""
function estimate_lambda_z(
    t::Vector{Float64},
    c::Vector{Float64},
    config::NCAConfig;
    exclude_indices::Vector{Int} = Int[],
    tmax_idx::Union{Int,Nothing} = nothing
)
    n = length(t)
    @assert n == length(c) "Time and concentration vectors must have same length"

    warnings = String[]

    # Find Tmax if not provided
    if tmax_idx === nothing
        tmax_idx = argmax(c)
    end

    # Get indices after Tmax with positive concentrations (terminal phase candidates)
    terminal_indices = Int[]
    for i in (tmax_idx + 1):n
        if c[i] > 0.0 && !(i in exclude_indices)
            push!(terminal_indices, i)
        end
    end

    # Also include Tmax point if it has positive concentration
    if c[tmax_idx] > 0.0 && !(tmax_idx in exclude_indices)
        pushfirst!(terminal_indices, tmax_idx)
    end

    n_terminal = length(terminal_indices)

    # Check minimum points requirement
    if n_terminal < config.lambda_z_min_points
        push!(warnings, "Insufficient points for lambda_z estimation (found $n_terminal, need $(config.lambda_z_min_points))")
        return LambdaZResult(
            warnings=warnings,
            quality_flag=:insufficient
        )
    end

    # Try different starting points and find best fit
    best_result = nothing
    best_adj_r2 = -Inf

    # Iterate from using all terminal points down to minimum points
    for start_offset in 0:(n_terminal - config.lambda_z_min_points)
        idx_range = terminal_indices[(start_offset + 1):end]
        n_points = length(idx_range)

        if n_points < config.lambda_z_min_points
            continue
        end

        # Extract time and log(concentration) for regression
        t_reg = [t[i] for i in idx_range]
        log_c_reg = [log(c[i]) for i in idx_range]

        # Perform linear regression: log(C) = intercept - lambda_z * t
        slope, intercept, r2, adj_r2 = _linear_regression(t_reg, log_c_reg)

        # Lambda_z is negative of slope (since C = C0 * exp(-lambda_z * t))
        lambda_z = -slope

        # Skip if lambda_z is not positive (not monotonically decreasing)
        if lambda_z <= 0.0
            continue
        end

        # Calculate half-life
        t_half = log(2.0) / lambda_z

        # Check if time span covers minimum half-lives
        time_span = t_reg[end] - t_reg[1]
        half_lives_covered = time_span / t_half

        # Check R² threshold
        if r2 < config.lambda_z_r2_threshold
            continue
        end

        # Select best by adjusted R²
        if adj_r2 > best_adj_r2
            best_adj_r2 = adj_r2

            # Determine quality flag
            quality_flag = :good
            result_warnings = String[]

            if half_lives_covered < config.lambda_z_span_half_lives
                quality_flag = :warning
                push!(result_warnings, "Time span covers only $(round(half_lives_covered, digits=2)) half-lives (recommended: $(config.lambda_z_span_half_lives))")
            end

            if r2 < 0.95
                if quality_flag == :good
                    quality_flag = :warning
                end
                push!(result_warnings, "R² = $(round(r2, digits=4)) is below 0.95")
            end

            best_result = LambdaZResult(
                lambda_z=lambda_z,
                t_half=t_half,
                r_squared=r2,
                adjusted_r_squared=adj_r2,
                intercept=intercept,
                n_points=n_points,
                start_time=t_reg[1],
                end_time=t_reg[end],
                points_used=idx_range,
                quality_flag=quality_flag,
                warnings=result_warnings
            )
        end
    end

    # Return best result or insufficient
    if best_result === nothing
        push!(warnings, "No valid lambda_z regression found meeting R² threshold of $(config.lambda_z_r2_threshold)")
        return LambdaZResult(
            warnings=warnings,
            quality_flag=:insufficient
        )
    end

    return best_result
end

"""
    _linear_regression(x, y)

Perform simple linear regression: y = a + b*x

Returns: (slope, intercept, r_squared, adjusted_r_squared)
"""
function _linear_regression(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    @assert n == length(y)
    @assert n >= 2

    # Calculate means
    x_mean = sum(x) / n
    y_mean = sum(y) / n

    # Calculate sums of squares
    ss_xx = 0.0
    ss_yy = 0.0
    ss_xy = 0.0

    for i in 1:n
        dx = x[i] - x_mean
        dy = y[i] - y_mean
        ss_xx += dx * dx
        ss_yy += dy * dy
        ss_xy += dx * dy
    end

    # Calculate slope and intercept
    slope = ss_xy / ss_xx
    intercept = y_mean - slope * x_mean

    # Calculate R²
    ss_res = 0.0
    for i in 1:n
        y_pred = intercept + slope * x[i]
        ss_res += (y[i] - y_pred)^2
    end

    r_squared = 1.0 - ss_res / ss_yy

    # Adjusted R² (accounts for degrees of freedom)
    # adj_r2 = 1 - (1 - r2) * (n - 1) / (n - p - 1), where p = 1 (one predictor)
    adjusted_r_squared = 1.0 - (1.0 - r_squared) * (n - 1) / (n - 2)

    return (slope, intercept, r_squared, adjusted_r_squared)
end

"""
    _find_best_lambda_z_range(t, c, tmax_idx, config)

Helper to find the optimal range for lambda_z estimation.
Returns indices of points to use.
"""
function _find_best_lambda_z_range(
    t::Vector{Float64},
    c::Vector{Float64},
    tmax_idx::Int,
    config::NCAConfig
)
    # Get all valid terminal phase points
    valid_indices = Int[]
    for i in tmax_idx:length(t)
        if c[i] > 0.0
            push!(valid_indices, i)
        end
    end

    return valid_indices
end
