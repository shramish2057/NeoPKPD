# Specifications are pure data. No solver logic and no hidden defaults.

export ModelKind,
    OneCompIVBolus, OneCompOralFirstOrder, OneCompIVBolusParams, OneCompOralFirstOrderParams
export DoseEvent, ModelSpec
export SolverSpec, SimGrid, SimResult

abstract type ModelKind end

struct OneCompIVBolus <: ModelKind end
struct OneCompOralFirstOrder <: ModelKind end

struct DoseEvent
    time::Float64
    amount::Float64
end

struct OneCompIVBolusParams
    CL::Float64
    V::Float64
end

struct OneCompOralFirstOrderParams
    Ka::Float64
    CL::Float64
    V::Float64
end

struct ModelSpec{K<:ModelKind,P}
    kind::K
    name::String
    params::P
    doses::Vector{DoseEvent}
end

struct SolverSpec
    alg::Symbol
    reltol::Float64
    abstol::Float64
    maxiters::Int
end

struct SimGrid
    t0::Float64
    t1::Float64
    saveat::Vector{Float64}
end

struct SimResult
    t::Vector{Float64}
    amount::Vector{Float64}
    conc::Vector{Float64}
    metadata::Dict{String,Any}
end
