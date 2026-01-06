# Trial Module
# Clinical trial simulation for OpenPKPD

# Include all trial submodules
include("specs.jl")
include("designs.jl")
include("regimens.jl")
include("population.jl")
include("events.jl")
include("endpoints.jl")
include("crossover.jl")  # Crossover-specific analysis
include("power.jl")
include("simulation.jl")  # Real subject exposure simulation
include("execution.jl")
include("adaptive.jl")    # Adaptive trials with interim analyses
include("escalation.jl")
