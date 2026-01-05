# Parameter Estimation Module for NLME Models
# Entry point for all estimation functionality

# Types and configuration
include("estimation_types.jl")

# Core estimation algorithm
include("estimate.jl")

# Individual method implementations
include("laplacian.jl")
include("foce.jl")  # FOCE-I with full Laplacian correction
include("saem.jl")

# Gradients and optimization utilities
include("gradients.jl")

# Standard errors and covariance matrix
include("standard_errors.jl")

# Covariance step (R and S matrices, sandwich estimators)
include("covariance.jl")

# Diagnostics (CWRES, IWRES, etc.)
include("diagnostics.jl")

# Bayesian estimation (HMC/NUTS)
include("bayesian.jl")

# Stepwise Covariate Modeling (SCM)
include("scm.jl")

# Bootstrap for SE validation
include("bootstrap.jl")

# Model Averaging (AIC/BIC)
include("model_averaging.jl")
