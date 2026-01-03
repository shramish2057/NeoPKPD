module OpenPKPDCLI

using ArgParse
using JSON
using OpenPKPDCore

# ============================================================================
# Utility Functions
# ============================================================================

function _die(msg::String)
    println(stderr, "Error: ", msg)
    exit(1)
end

function _info(msg::String)
    println(stderr, msg)
end

function _write_json(path::String, data::Any)
    open(path, "w") do io
        JSON.print(io, data, 2)
    end
end

# ============================================================================
# Version Command
# ============================================================================

const VERSION_HELP = """
Display version information for OpenPKPD.

Shows the OpenPKPD version along with semantic versions for event handling,
solver behavior, and artifact schema.

Example:
  openpkpd version
"""

function cmd_version()
    println("OpenPKPD " * OpenPKPDCore.OPENPKPD_VERSION)
    println("Event semantics: " * OpenPKPDCore.EVENT_SEMANTICS_VERSION)
    println("Solver semantics: " * OpenPKPDCore.SOLVER_SEMANTICS_VERSION)
    println("Artifact schema: " * OpenPKPDCore.ARTIFACT_SCHEMA_VERSION)
end

# ============================================================================
# Replay Command
# ============================================================================

const REPLAY_HELP = """
Replay a simulation from a saved artifact JSON file.

This command re-executes a simulation using the exact parameters stored in an
artifact file, validating reproducibility. Supports single, population, and
sensitivity artifacts.

Arguments:
  --artifact PATH    Path to the artifact JSON file (required)
  --out PATH         Optional output path for new artifact with replayed results

Examples:
  openpkpd replay --artifact validation/golden/pk_iv_bolus.json
  openpkpd replay --artifact my_sim.json --out replayed.json
"""

function cmd_replay(args)
    path = args["artifact"]
    out = get(args, "out", nothing)

    artifact = OpenPKPDCore.read_execution_json(path)
    atype = get(artifact, "artifact_type", "single")

    if atype == "population"
        res = OpenPKPDCore.replay_population_execution(artifact)
        _info("Replayed population simulation with $(length(res.individuals)) individuals")
        if out !== nothing
            parsed = OpenPKPDCore.deserialize_population_execution(artifact)
            OpenPKPDCore.write_population_json(
                out;
                population_spec = parsed.population_spec,
                grid = parsed.grid,
                solver = parsed.solver,
                result = res,
            )
            _info("Written to: $out")
        end
        return
    end

    if atype == "sensitivity_single"
        res = OpenPKPDCore.replay_sensitivity_execution(artifact)
        _info("Replayed single sensitivity analysis")
        if out !== nothing
            parsed = OpenPKPDCore.deserialize_sensitivity_execution(artifact)
            OpenPKPDCore.write_sensitivity_json(
                out;
                model_spec = parsed.model_spec,
                grid = parsed.grid,
                solver = parsed.solver,
                result = res,
            )
            _info("Written to: $out")
        end
        return
    end

    if atype == "sensitivity_population"
        res = OpenPKPDCore.replay_population_sensitivity_execution(artifact)
        _info("Replayed population sensitivity analysis")
        if out !== nothing
            parsed = OpenPKPDCore.deserialize_population_sensitivity_execution(artifact)
            OpenPKPDCore.write_population_sensitivity_json(
                out;
                population_spec = parsed.population_spec,
                grid = parsed.grid,
                solver = parsed.solver,
                result = res,
            )
            _info("Written to: $out")
        end
        return
    end

    # default: single execution
    res = OpenPKPDCore.replay_execution(artifact)
    _info("Replayed single simulation")
    if out !== nothing
        parsed = OpenPKPDCore.deserialize_execution(artifact)
        OpenPKPDCore.write_execution_json(
            out;
            model_spec = parsed.model_spec,
            grid = parsed.grid,
            solver = parsed.solver,
            result = res,
            pd_spec = parsed.pd_spec,
        )
        _info("Written to: $out")
    end
end

# ============================================================================
# Validate Golden Command
# ============================================================================

const VALIDATE_GOLDEN_HELP = """
Run golden validation tests.

Validates all golden artifacts in the validation/golden directory to ensure
reproducibility across versions and platforms.

Example:
  openpkpd validate-golden
"""

function cmd_validate_golden()
    cli_src = @__DIR__
    repo_root = normpath(joinpath(cli_src, "..", "..", ".."))

    runner = joinpath(repo_root, "validation", "scripts", "run_golden_validation.jl")

    if !isfile(runner)
        _die("Golden validation runner not found: " * runner)
    end

    cmd = `julia $runner`
    run(Cmd(cmd; dir=repo_root))
end

# ============================================================================
# Simulate Command
# ============================================================================

const SIMULATE_HELP = """
Run a PK or PKPD simulation from a JSON specification file.

SUPPORTED PK MODELS:
  OneCompIVBolus          - One-compartment IV bolus
    params: {CL, V}
  OneCompOralFirstOrder   - One-compartment oral first-order absorption
    params: {Ka, CL, V}
  TwoCompIVBolus          - Two-compartment IV bolus
    params: {CL, V1, Q, V2}
  TwoCompOral             - Two-compartment oral first-order absorption
    params: {Ka, CL, V1, Q, V2}
  ThreeCompIVBolus        - Three-compartment IV bolus (mammillary)
    params: {CL, V1, Q2, V2, Q3, V3}
  TransitAbsorption       - Transit compartment chain model
    params: {N, Ktr, Ka, CL, V}
  MichaelisMentenElimination - Saturable (nonlinear) elimination
    params: {Vmax, Km, V}

SUPPORTED PD MODELS:
  DirectEmax              - Direct Emax (hyperbolic)
    params: {E0, Emax, EC50}
  SigmoidEmax             - Sigmoid Emax (Hill equation)
    params: {E0, Emax, EC50, gamma}
  IndirectResponseTurnover - Indirect response turnover
    params: {Kin, Kout, R0, Imax, IC50}
  BiophaseEquilibration   - Effect compartment model
    params: {ke0, E0, Emax, EC50}

Input JSON format:
{
  "model": {
    "kind": "<model_kind>",
    "name": "my_simulation",
    "params": {...},
    "doses": [{"time": 0.0, "amount": 100.0}]
  },
  "grid": {"t0": 0.0, "t1": 24.0, "saveat": [0.0, 1.0, 2.0, ...]},
  "solver": {"alg": "Tsit5", "reltol": 1e-10, "abstol": 1e-12, "maxiters": 10000000},
  "pd": {  # optional
    "kind": "<pd_kind>",
    "name": "pd_model",
    "params": {...},
    "input_observation": "conc",
    "output_observation": "effect"
  }
}

Arguments:
  --spec PATH       Path to simulation specification JSON (required)
  --out PATH        Output path for results JSON (required)
  --format FORMAT   Output format: "artifact" or "simple" (default: artifact)

Examples:
  openpkpd simulate --spec pk_spec.json --out result.json
  openpkpd simulate --spec pkpd_spec.json --out result.json --format simple

  # Two-compartment IV bolus example spec:
  # {"model": {"kind": "TwoCompIVBolus", "params": {"CL": 10, "V1": 50, "Q": 5, "V2": 100}, "doses": [{"time": 0, "amount": 500}]}, ...}

  # Sigmoid Emax PD example:
  # {"pd": {"kind": "SigmoidEmax", "params": {"E0": 0, "Emax": 100, "EC50": 5, "gamma": 2}, ...}}
"""

function _parse_model_spec(spec::Dict)
    kind_str = spec["kind"]
    name = get(spec, "name", "simulation")
    params_dict = spec["params"]
    doses_raw = get(spec, "doses", [])

    doses = [OpenPKPDCore.DoseEvent(Float64(d["time"]), Float64(d["amount"])) for d in doses_raw]

    if kind_str == "OneCompIVBolus"
        params = OpenPKPDCore.OneCompIVBolusParams(
            Float64(params_dict["CL"]),
            Float64(params_dict["V"]),
        )
        return OpenPKPDCore.ModelSpec(OpenPKPDCore.OneCompIVBolus(), name, params, doses)

    elseif kind_str == "OneCompOralFirstOrder"
        params = OpenPKPDCore.OneCompOralFirstOrderParams(
            Float64(params_dict["Ka"]),
            Float64(params_dict["CL"]),
            Float64(params_dict["V"]),
        )
        return OpenPKPDCore.ModelSpec(OpenPKPDCore.OneCompOralFirstOrder(), name, params, doses)

    elseif kind_str == "TwoCompIVBolus"
        params = OpenPKPDCore.TwoCompIVBolusParams(
            Float64(params_dict["CL"]),
            Float64(params_dict["V1"]),
            Float64(params_dict["Q"]),
            Float64(params_dict["V2"]),
        )
        return OpenPKPDCore.ModelSpec(OpenPKPDCore.TwoCompIVBolus(), name, params, doses)

    elseif kind_str == "TwoCompOral"
        params = OpenPKPDCore.TwoCompOralParams(
            Float64(params_dict["Ka"]),
            Float64(params_dict["CL"]),
            Float64(params_dict["V1"]),
            Float64(params_dict["Q"]),
            Float64(params_dict["V2"]),
        )
        return OpenPKPDCore.ModelSpec(OpenPKPDCore.TwoCompOral(), name, params, doses)

    elseif kind_str == "ThreeCompIVBolus"
        params = OpenPKPDCore.ThreeCompIVBolusParams(
            Float64(params_dict["CL"]),
            Float64(params_dict["V1"]),
            Float64(params_dict["Q2"]),
            Float64(params_dict["V2"]),
            Float64(params_dict["Q3"]),
            Float64(params_dict["V3"]),
        )
        return OpenPKPDCore.ModelSpec(OpenPKPDCore.ThreeCompIVBolus(), name, params, doses)

    elseif kind_str == "TransitAbsorption"
        params = OpenPKPDCore.TransitAbsorptionParams(
            Int(params_dict["N"]),
            Float64(params_dict["Ktr"]),
            Float64(params_dict["Ka"]),
            Float64(params_dict["CL"]),
            Float64(params_dict["V"]),
        )
        return OpenPKPDCore.ModelSpec(OpenPKPDCore.TransitAbsorption(), name, params, doses)

    elseif kind_str == "MichaelisMentenElimination"
        params = OpenPKPDCore.MichaelisMentenEliminationParams(
            Float64(params_dict["Vmax"]),
            Float64(params_dict["Km"]),
            Float64(params_dict["V"]),
        )
        return OpenPKPDCore.ModelSpec(OpenPKPDCore.MichaelisMentenElimination(), name, params, doses)

    else
        _die("Unsupported model kind: $kind_str. Supported: OneCompIVBolus, OneCompOralFirstOrder, TwoCompIVBolus, TwoCompOral, ThreeCompIVBolus, TransitAbsorption, MichaelisMentenElimination")
    end
end

function _parse_pd_spec(spec::Dict)
    kind_str = spec["kind"]
    name = get(spec, "name", "pd_model")
    params_dict = spec["params"]
    input_obs = Symbol(get(spec, "input_observation", "conc"))
    output_obs = Symbol(get(spec, "output_observation", "effect"))

    if kind_str == "DirectEmax"
        params = OpenPKPDCore.DirectEmaxParams(
            Float64(params_dict["E0"]),
            Float64(params_dict["Emax"]),
            Float64(params_dict["EC50"]),
        )
        return OpenPKPDCore.PDSpec(OpenPKPDCore.DirectEmax(), name, params, input_obs, output_obs)

    elseif kind_str == "SigmoidEmax"
        params = OpenPKPDCore.SigmoidEmaxParams(
            Float64(params_dict["E0"]),
            Float64(params_dict["Emax"]),
            Float64(params_dict["EC50"]),
            Float64(params_dict["gamma"]),
        )
        return OpenPKPDCore.PDSpec(OpenPKPDCore.SigmoidEmax(), name, params, input_obs, output_obs)

    elseif kind_str == "IndirectResponseTurnover"
        params = OpenPKPDCore.IndirectResponseTurnoverParams(
            Float64(params_dict["Kin"]),
            Float64(params_dict["Kout"]),
            Float64(params_dict["R0"]),
            Float64(params_dict["Imax"]),
            Float64(params_dict["IC50"]),
        )
        return OpenPKPDCore.PDSpec(OpenPKPDCore.IndirectResponseTurnover(), name, params, input_obs, output_obs)

    elseif kind_str == "BiophaseEquilibration"
        params = OpenPKPDCore.BiophaseEquilibrationParams(
            Float64(params_dict["ke0"]),
            Float64(params_dict["E0"]),
            Float64(params_dict["Emax"]),
            Float64(params_dict["EC50"]),
        )
        return OpenPKPDCore.PDSpec(OpenPKPDCore.BiophaseEquilibration(), name, params, input_obs, output_obs)

    else
        _die("Unsupported PD kind: $kind_str. Supported: DirectEmax, SigmoidEmax, IndirectResponseTurnover, BiophaseEquilibration")
    end
end

function _parse_grid(spec::Dict)
    return OpenPKPDCore.SimGrid(
        Float64(spec["t0"]),
        Float64(spec["t1"]),
        [Float64(x) for x in spec["saveat"]],
    )
end

function _parse_solver(spec::Dict)
    return OpenPKPDCore.SolverSpec(
        Symbol(get(spec, "alg", "Tsit5")),
        Float64(get(spec, "reltol", 1e-10)),
        Float64(get(spec, "abstol", 1e-12)),
        Int(get(spec, "maxiters", 10^7)),
    )
end

function _simresult_to_dict(res::OpenPKPDCore.SimResult)
    return Dict(
        "t" => res.t,
        "states" => Dict(string(k) => v for (k, v) in res.states),
        "observations" => Dict(string(k) => v for (k, v) in res.observations),
        "metadata" => res.metadata,
    )
end

function cmd_simulate(args)
    spec_path = args["spec"]
    out_path = args["out"]
    format = get(args, "format", "artifact")

    if !isfile(spec_path)
        _die("Specification file not found: $spec_path")
    end

    spec = JSON.parsefile(spec_path; dicttype=Dict{String, Any})

    model_spec = _parse_model_spec(spec["model"])
    grid = _parse_grid(spec["grid"])
    solver = haskey(spec, "solver") ? _parse_solver(spec["solver"]) : OpenPKPDCore.SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd_spec = haskey(spec, "pd") ? _parse_pd_spec(spec["pd"]) : nothing

    if pd_spec !== nothing
        pd_kind = pd_spec.kind
        # Determine if this is a direct (pk_then_pd) or coupled PD model
        if pd_kind isa OpenPKPDCore.DirectEmax || pd_kind isa OpenPKPDCore.SigmoidEmax || pd_kind isa OpenPKPDCore.BiophaseEquilibration
            # Direct PD models - use simulate_pkpd (pk_then_pd mode)
            result = OpenPKPDCore.simulate_pkpd(model_spec, pd_spec, grid, solver)
            _info("Simulated PK then PD: $(length(result.t)) time points")
        elseif pd_kind isa OpenPKPDCore.IndirectResponseTurnover
            # Coupled ODE PD models - use simulate_pkpd_coupled
            result = OpenPKPDCore.simulate_pkpd_coupled(model_spec, pd_spec, grid, solver)
            _info("Simulated coupled PKPD: $(length(result.t)) time points")
        else
            _die("Unsupported PD model type for simulation")
        end
    else
        # PK only
        result = OpenPKPDCore.simulate(model_spec, grid, solver)
        _info("Simulated PK: $(length(result.t)) time points")
    end

    if format == "simple"
        _write_json(out_path, _simresult_to_dict(result))
    else
        OpenPKPDCore.write_execution_json(
            out_path;
            model_spec = model_spec,
            grid = grid,
            solver = solver,
            result = result,
            pd_spec = pd_spec,
        )
    end

    _info("Written to: $out_path")
end

# ============================================================================
# Population Command
# ============================================================================

const POPULATION_HELP = """
Run a population PK/PD simulation with inter-individual variability.

Input JSON format:
{
  "model": { ... },  # Same as simulate command
  "grid": { ... },
  "solver": { ... },
  "iiv": {
    "kind": "LogNormalIIV",
    "omegas": {"CL": 0.3, "V": 0.2},
    "seed": 12345,
    "n": 100
  },
  "iov": {  # optional
    "kind": "LogNormalIIV",
    "pis": {"CL": 0.1},
    "seed": 54321,
    "occasion_def": "dose_times"
  },
  "covariate_model": {  # optional
    "name": "wt_model",
    "effects": [
      {"kind": "PowerCovariate", "param": "CL", "covariate": "WT", "beta": 0.75, "ref": 70.0}
    ]
  },
  "covariates": [  # optional, one per individual
    {"values": {"WT": 70.0}},
    ...
  ]
}

Arguments:
  --spec PATH       Path to population specification JSON (required)
  --out PATH        Output path for results JSON (required)
  --format FORMAT   Output format: "artifact" or "simple" (default: artifact)

Examples:
  openpkpd population --spec pop_spec.json --out pop_result.json
"""

function _parse_iiv_spec(spec::Dict)
    omegas = Dict(Symbol(k) => Float64(v) for (k, v) in spec["omegas"])
    return OpenPKPDCore.IIVSpec(
        OpenPKPDCore.LogNormalIIV(),
        omegas,
        UInt64(spec["seed"]),
        Int(spec["n"]),
    )
end

function _parse_iov_spec(spec::Dict)
    pis = Dict(Symbol(k) => Float64(v) for (k, v) in spec["pis"])
    occasion_mode = Symbol(get(spec, "occasion_def", "dose_times"))
    return OpenPKPDCore.IOVSpec(
        OpenPKPDCore.LogNormalIIV(),
        pis,
        UInt64(spec["seed"]),
        OpenPKPDCore.OccasionDefinition(occasion_mode),
    )
end

function _parse_covariate_effect(spec::Dict)
    kind_str = spec["kind"]
    if kind_str == "LinearCovariate"
        kind = OpenPKPDCore.LinearCovariate()
    elseif kind_str == "PowerCovariate"
        kind = OpenPKPDCore.PowerCovariate()
    elseif kind_str == "ExpCovariate"
        kind = OpenPKPDCore.ExpCovariate()
    else
        _die("Unsupported covariate effect kind: $kind_str")
    end
    return OpenPKPDCore.CovariateEffect(
        kind,
        Symbol(spec["param"]),
        Symbol(spec["covariate"]),
        Float64(spec["beta"]),
        Float64(spec["ref"]),
    )
end

function _parse_covariate_model(spec::Dict)
    effects = [_parse_covariate_effect(e) for e in spec["effects"]]
    return OpenPKPDCore.CovariateModel(spec["name"], effects)
end

function _parse_individual_covariates(specs::Vector)
    return [
        OpenPKPDCore.IndividualCovariates(
            Dict(Symbol(k) => Float64(v) for (k, v) in get(s, "values", Dict())),
            nothing,  # time_varying not supported via CLI yet
        )
        for s in specs
    ]
end

function _popresult_to_dict(res)
    individuals = [_simresult_to_dict(r) for r in res.individuals]
    params = [Dict(string(k) => v for (k, v) in d) for d in res.params]

    summaries = Dict()
    for (k, s) in res.summaries
        summaries[string(k)] = Dict(
            "observation" => string(s.observation),
            "probs" => s.probs,
            "mean" => s.mean,
            "median" => s.median,
            "quantiles" => Dict(string(p) => v for (p, v) in s.quantiles),
        )
    end

    return Dict(
        "individuals" => individuals,
        "params" => params,
        "summaries" => summaries,
        "metadata" => res.metadata,
    )
end

function cmd_population(args)
    spec_path = args["spec"]
    out_path = args["out"]
    format = get(args, "format", "artifact")

    if !isfile(spec_path)
        _die("Specification file not found: $spec_path")
    end

    spec = JSON.parsefile(spec_path; dicttype=Dict{String, Any})

    model_spec = _parse_model_spec(spec["model"])
    grid = _parse_grid(spec["grid"])
    solver = haskey(spec, "solver") ? _parse_solver(spec["solver"]) : OpenPKPDCore.SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    iiv = haskey(spec, "iiv") ? _parse_iiv_spec(spec["iiv"]) : nothing
    iov = haskey(spec, "iov") ? _parse_iov_spec(spec["iov"]) : nothing
    cov_model = haskey(spec, "covariate_model") ? _parse_covariate_model(spec["covariate_model"]) : nothing
    covariates = haskey(spec, "covariates") ? _parse_individual_covariates(spec["covariates"]) : OpenPKPDCore.IndividualCovariates[]

    pop_spec = OpenPKPDCore.PopulationSpec(model_spec, iiv, iov, cov_model, covariates)

    result = OpenPKPDCore.simulate_population(pop_spec, grid, solver)
    n = length(result.individuals)
    _info("Simulated population: $n individuals")

    if format == "simple"
        _write_json(out_path, _popresult_to_dict(result))
    else
        OpenPKPDCore.write_population_json(
            out_path;
            population_spec = pop_spec,
            grid = grid,
            solver = solver,
            result = result,
        )
    end

    _info("Written to: $out_path")
end

# ============================================================================
# Sensitivity Command
# ============================================================================

const SENSITIVITY_HELP = """
Run sensitivity analysis by perturbing parameters.

Input JSON format:
{
  "model": { ... },  # Same as simulate command
  "grid": { ... },
  "solver": { ... },
  "perturbation": {
    "name": "cl_sensitivity",
    "param": "CL",
    "delta": 0.01
  },
  "observation": "conc"
}

Arguments:
  --spec PATH       Path to sensitivity specification JSON (required)
  --out PATH        Output path for results JSON (required)

Examples:
  openpkpd sensitivity --spec sens_spec.json --out sens_result.json
"""

function cmd_sensitivity(args)
    spec_path = args["spec"]
    out_path = args["out"]

    if !isfile(spec_path)
        _die("Specification file not found: $spec_path")
    end

    spec = JSON.parsefile(spec_path; dicttype=Dict{String, Any})

    model_spec = _parse_model_spec(spec["model"])
    grid = _parse_grid(spec["grid"])
    solver = haskey(spec, "solver") ? _parse_solver(spec["solver"]) : OpenPKPDCore.SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pert_spec = spec["perturbation"]
    plan = OpenPKPDCore.PerturbationPlan(
        pert_spec["name"],
        Symbol(pert_spec["param"]),
        Float64(pert_spec["delta"]),
    )

    observation = Symbol(get(spec, "observation", "conc"))

    result = OpenPKPDCore.run_sensitivity(model_spec, grid, solver, plan, observation)

    _info("Sensitivity analysis complete")
    _info("  Max absolute delta: $(result.metrics.max_abs_delta)")
    _info("  Max relative delta: $(result.metrics.max_rel_delta)")
    _info("  L2 norm delta: $(result.metrics.l2_norm_delta)")

    OpenPKPDCore.write_sensitivity_json(
        out_path;
        model_spec = model_spec,
        grid = grid,
        solver = solver,
        result = result,
    )

    _info("Written to: $out_path")
end

# ============================================================================
# Metrics Command
# ============================================================================

const METRICS_HELP = """
Compute PK/PD metrics from a simulation result or artifact.

Available metrics:
  - cmax: Maximum concentration
  - auc: Area under the curve (trapezoidal)
  - emin: Minimum response value
  - time_below: Time spent below a threshold
  - auc_above_baseline: AUC above a baseline value

Arguments:
  --artifact PATH      Path to artifact or result JSON (required)
  --observation NAME   Observation to analyze (default: conc)
  --metric NAME        Metric to compute: cmax, auc, emin, time_below, auc_above_baseline
  --threshold VALUE    Threshold for time_below or auc_above_baseline

Examples:
  openpkpd metrics --artifact result.json --metric cmax
  openpkpd metrics --artifact result.json --observation effect --metric emin
  openpkpd metrics --artifact result.json --metric time_below --threshold 10.0
"""

function cmd_metrics(args)
    artifact_path = args["artifact"]
    observation = Symbol(get(args, "observation", "conc"))
    metric_name = get(args, "metric", nothing)
    threshold = haskey(args, "threshold") && args["threshold"] !== nothing ? parse(Float64, args["threshold"]) : nothing

    if !isfile(artifact_path)
        _die("Artifact file not found: $artifact_path")
    end

    if metric_name === nothing
        _die("Metric name required. Use --metric cmax|auc|emin|time_below|auc_above_baseline")
    end

    artifact = JSON.parsefile(artifact_path)

    # Extract time and observation data
    t = nothing
    y = nothing

    if haskey(artifact, "result") && haskey(artifact["result"], "t")
        # Artifact format
        t = Float64.(artifact["result"]["t"])
        obs_data = artifact["result"]["observations"]
        obs_key = string(observation)
        if !haskey(obs_data, obs_key)
            _die("Observation '$obs_key' not found. Available: $(join(keys(obs_data), ", "))")
        end
        y = Float64.(obs_data[obs_key])
    elseif haskey(artifact, "t")
        # Simple format
        t = Float64.(artifact["t"])
        obs_data = artifact["observations"]
        obs_key = string(observation)
        if !haskey(obs_data, obs_key)
            _die("Observation '$obs_key' not found. Available: $(join(keys(obs_data), ", "))")
        end
        y = Float64.(obs_data[obs_key])
    else
        _die("Unrecognized file format")
    end

    result = nothing

    if metric_name == "cmax"
        result = OpenPKPDCore.cmax(t, y)
        println("Cmax: $result")
    elseif metric_name == "auc"
        result = OpenPKPDCore.auc_trapezoid(t, y)
        println("AUC: $result")
    elseif metric_name == "emin"
        result = OpenPKPDCore.emin(t, y)
        println("Emin: $result")
    elseif metric_name == "time_below"
        if threshold === nothing
            _die("--threshold required for time_below metric")
        end
        result = OpenPKPDCore.time_below(t, y, threshold)
        println("Time below $threshold: $result")
    elseif metric_name == "auc_above_baseline"
        if threshold === nothing
            _die("--threshold required for auc_above_baseline metric")
        end
        result = OpenPKPDCore.auc_above_baseline(t, y, threshold)
        println("AUC above baseline $threshold: $result")
    else
        _die("Unknown metric: $metric_name. Supported: cmax, auc, emin, time_below, auc_above_baseline")
    end
end

# ============================================================================
# Help Command
# ============================================================================

const MAIN_HELP = """
OpenPKPD - Professional PK/PD Simulation Platform

USAGE:
    openpkpd <command> [options]

COMMANDS:
    version          Display version information
    simulate         Run a PK or PKPD simulation from a JSON spec
    population       Run a population simulation with IIV/IOV
    sensitivity      Run parameter sensitivity analysis
    metrics          Compute PK/PD metrics from simulation results
    replay           Replay a simulation from an artifact file
    validate-golden  Run golden validation tests
    help             Show this help message

SUPPORTED PK MODELS:
    OneCompIVBolus, OneCompOralFirstOrder, TwoCompIVBolus, TwoCompOral,
    ThreeCompIVBolus, TransitAbsorption, MichaelisMentenElimination

SUPPORTED PD MODELS:
    DirectEmax, SigmoidEmax, IndirectResponseTurnover, BiophaseEquilibration

Use 'openpkpd help <command>' for detailed help on each command.

EXAMPLES:
    openpkpd version
    openpkpd simulate --spec pk_spec.json --out result.json
    openpkpd population --spec pop_spec.json --out pop_result.json
    openpkpd metrics --artifact result.json --metric cmax
    openpkpd replay --artifact validation/golden/pk_iv_bolus.json

For more information, visit: https://github.com/openpkpd/openpkpd
"""

function cmd_help(command::Union{String,Nothing})
    if command === nothing
        println(MAIN_HELP)
    elseif command == "version"
        println(VERSION_HELP)
    elseif command == "simulate"
        println(SIMULATE_HELP)
    elseif command == "population"
        println(POPULATION_HELP)
    elseif command == "sensitivity"
        println(SENSITIVITY_HELP)
    elseif command == "metrics"
        println(METRICS_HELP)
    elseif command == "replay"
        println(REPLAY_HELP)
    elseif command == "validate-golden"
        println(VALIDATE_GOLDEN_HELP)
    else
        println("Unknown command: $command")
        println()
        println(MAIN_HELP)
    end
end

# ============================================================================
# Main Entry Point
# ============================================================================

function main()
    if length(ARGS) < 1
        println(MAIN_HELP)
        return
    end

    cmd = ARGS[1]
    rest = ARGS[2:end]

    if cmd == "version"
        cmd_version()

    elseif cmd == "help"
        if length(rest) > 0
            cmd_help(rest[1])
        else
            cmd_help(nothing)
        end

    elseif cmd == "replay"
        rs = ArgParseSettings()
        @add_arg_table rs begin
            "--artifact"
                required = true
                help = "Path to the artifact JSON file"
            "--out"
                required = false
                help = "Output path for replayed artifact"
        end
        args = parse_args(rest, rs)
        cmd_replay(args)

    elseif cmd == "validate-golden"
        cmd_validate_golden()

    elseif cmd == "simulate"
        rs = ArgParseSettings()
        @add_arg_table rs begin
            "--spec"
                required = true
                help = "Path to simulation specification JSON"
            "--out"
                required = true
                help = "Output path for results"
            "--format"
                required = false
                default = "artifact"
                help = "Output format: artifact or simple"
        end
        args = parse_args(rest, rs)
        cmd_simulate(args)

    elseif cmd == "population"
        rs = ArgParseSettings()
        @add_arg_table rs begin
            "--spec"
                required = true
                help = "Path to population specification JSON"
            "--out"
                required = true
                help = "Output path for results"
            "--format"
                required = false
                default = "artifact"
                help = "Output format: artifact or simple"
        end
        args = parse_args(rest, rs)
        cmd_population(args)

    elseif cmd == "sensitivity"
        rs = ArgParseSettings()
        @add_arg_table rs begin
            "--spec"
                required = true
                help = "Path to sensitivity specification JSON"
            "--out"
                required = true
                help = "Output path for results"
        end
        args = parse_args(rest, rs)
        cmd_sensitivity(args)

    elseif cmd == "metrics"
        rs = ArgParseSettings()
        @add_arg_table rs begin
            "--artifact"
                required = true
                help = "Path to artifact or result JSON"
            "--observation"
                required = false
                default = "conc"
                help = "Observation to analyze"
            "--metric"
                required = true
                help = "Metric: cmax, auc, emin, time_below, auc_above_baseline"
            "--threshold"
                required = false
                help = "Threshold for time_below or auc_above_baseline"
        end
        args = parse_args(rest, rs)
        cmd_metrics(args)

    else
        _die("Unknown command: $cmd. Use 'openpkpd help' for usage.")
    end
end

# Run when invoked as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end # module
