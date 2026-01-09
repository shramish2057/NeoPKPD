<p align="center">
  <h1 align="center">NeoPKPD</h1>
  <p align="center">
    <strong>Transparent, validated pharmacokinetics and pharmacodynamics modeling infrastructure</strong>
  </p>
  <p align="center">
    <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="License: MIT"></a>
    <a href="https://julialang.org/"><img src="https://img.shields.io/badge/Julia-1.10+-purple.svg" alt="Julia"></a>
    <a href="https://python.org/"><img src="https://img.shields.io/badge/Python-3.10+-blue.svg" alt="Python"></a>
  </p>
  <p align="center">
    <a href="#features">Features</a> •
    <a href="#installation">Installation</a> •
    <a href="#quick-start">Quick Start</a> •
    <a href="#documentation">Documentation</a> •
    <a href="#contributing">Contributing</a>
  </p>
</p>

---

## Overview

**NeoPKPD** is a reference-grade PK/PD simulation platform designed for research, method development, and reproducible scientific computation. It emphasizes deterministic execution, transparent numerical semantics, and complete artifact serialization.

## Features

| Category | Features |
|----------|----------|
| **PK Models** | One/Two/Three-compartment IV & oral, transit absorption, Michaelis-Menten |
| **PD Models** | Direct Emax, sigmoid Emax, biophase equilibration, indirect response |
| **IV Infusion** | Zero-order infusion with duration support, overlapping infusions |
| **Population** | IIV, IOV, static & time-varying covariates |
| **Parameter Estimation** | NLME with FOCE-I, SAEM, and Laplacian methods |
| **NCA** | FDA/EMA-compliant non-compartmental analysis |
| **Trial Simulation** | Parallel, crossover, dose-escalation, bioequivalence designs |
| **Sensitivity** | Single-subject and population-level analysis |
| **VPC** | Visual Predictive Checks with pcVPC and stratification |
| **Model Import** | NONMEM (.ctl) and Monolix (.mlxtran) model parsing |
| **Data Import** | CDISC/SDTM format support (PC, EX, DM domains) |
| **Residual Error** | Additive, proportional, combined, exponential models |
| **Visualization** | Matplotlib/Plotly backends with estimation diagnostics |
| **Interfaces** | Julia API, Python bindings, CLI |
| **Reproducibility** | Versioned artifacts with deterministic replay |

## Installation

### Julia (Core)

```bash
git clone https://github.com/neopkpd/neopkpd.git
cd neopkpd

# Install dependencies
julia --project=packages/core -e 'using Pkg; Pkg.instantiate()'

# Verify installation
julia --project=packages/core -e 'using NeoPKPDCore; println("v", NEOPKPD_VERSION)'
```

### Python (Optional)

```bash
cd packages/python
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

### CLI

```bash
./packages/cli/bin/neopkpd version
```

## Quick Start

### Julia

```julia
using NeoPKPDCore

# Define model
spec = ModelSpec(
    OneCompIVBolus(),
    "example",
    OneCompIVBolusParams(5.0, 50.0),  # CL=5 L/h, V=50 L
    [DoseEvent(0.0, 100.0)]            # 100 mg at t=0
)

# Configure simulation
grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# Run and access results
result = simulate(spec, grid, solver)
println(result.observations[:conc])
```

### Python

```python
import neopkpd

neopkpd.init_julia()
result = neopkpd.simulate_pk_iv_bolus(
    cl=5.0, v=50.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0, t1=24.0,
    saveat=[float(t) for t in range(25)]
)
print(result["observations"]["conc"])
```

### CLI

```bash
# Check version
./packages/cli/bin/neopkpd version

# Run simulation from JSON spec
./packages/cli/bin/neopkpd simulate --spec simulation.json --out result.json

# Run NCA analysis
./packages/cli/bin/neopkpd nca --spec nca_spec.json --out nca_result.json

# Run clinical trial simulation
./packages/cli/bin/neopkpd trial --spec trial_spec.json --out trial_result.json

# Compute VPC
./packages/cli/bin/neopkpd vpc --spec vpc_spec.json --out vpc_result.json

# Import NONMEM model
./packages/cli/bin/neopkpd import --input model.ctl --format nonmem --out model.json

# Replay an artifact
./packages/cli/bin/neopkpd replay --artifact validation/golden/pk_iv_bolus.json

# Validate all golden artifacts
./packages/cli/bin/neopkpd validate-golden
```

## Repository Structure

```
neopkpd/
├── packages/
│   ├── core/                 # Julia simulation engine (NeoPKPDCore)
│   │   ├── src/
│   │   │   ├── models/       # PK/PD model definitions
│   │   │   ├── engine/       # ODE solving, population, infusion
│   │   │   ├── specs/        # Type specifications (ModelSpec, etc.)
│   │   │   ├── estimation/   # NLME: FOCE-I, SAEM, Laplacian
│   │   │   ├── import/       # NONMEM/Monolix parsers
│   │   │   ├── data/         # CDISC data handling
│   │   │   ├── analysis/     # VPC, NCA, sensitivity
│   │   │   └── serialization/# JSON artifact I/O
│   │   └── test/             # Comprehensive test suite
│   ├── python/               # Python bindings (neopkpd)
│   │   └── neopkpd/
│   │       ├── simulations/  # PK/PD simulation wrappers
│   │       ├── nca/          # Non-compartmental analysis
│   │       ├── trial/        # Clinical trial simulation
│   │       ├── estimation/   # NLME Python interface
│   │       ├── import_/      # Model import utilities
│   │       ├── data/         # CDISC data utilities
│   │       └── viz/          # Visualization (matplotlib/plotly)
│   └── cli/                  # Command-line interface (NeoPKPDCLI)
│       ├── src/              # CLI commands
│       └── bin/              # Entry point script
├── validation/               # Golden artifacts & validation
│   ├── golden/               # Reference outputs (deterministic)
│   └── scripts/              # Validation runners
├── docs/                     # Documentation (MkDocs)
│   ├── examples/             # Executable examples
│   └── *.md                  # API references
└── scripts/                  # Development tools
```

## Testing

```bash
# Julia unit tests
julia --project=packages/core -e 'using Pkg; Pkg.test()'

# Golden artifact validation
./packages/cli/bin/neopkpd validate-golden

# Python tests
cd packages/python && source .venv/bin/activate && pytest tests/

# Documentation build
mkdocs build --strict
```

## Semantic Versioning

NeoPKPD uses independent version numbers for numerical behavior:

| Version | Current | Scope |
|---------|---------|-------|
| Event Semantics | 1.0.0 | Dose handling |
| Solver Semantics | 1.0.0 | ODE solver behavior |
| Artifact Schema | 1.0.0 | JSON format |

Any change to numerical output requires a version bump.

## Documentation

Full documentation: [shramish2057.github.io/NeoPKPD](https://shramish2057.github.io/NeoPKPD/)

| Guide | Description |
|-------|-------------|
| [Models](docs/models.md) | PK/PD model reference |
| [Estimation](docs/estimation.md) | NLME parameter estimation (FOCE-I, SAEM) |
| [Population](docs/population.md) | IIV, IOV, covariates |
| [NCA](docs/nca.md) | Non-compartmental analysis |
| [VPC](docs/vpc.md) | Visual Predictive Checks |
| [Trial Simulation](docs/trial.md) | Clinical trial design |
| [Data Import](docs/data.md) | CDISC/SDTM format |
| [Model Import](docs/import.md) | NONMEM/Monolix parsing |
| [CLI Reference](docs/cli.md) | Command-line interface |
| [Python API](docs/python.md) | Python bindings |
| [Visualization](docs/visualization.md) | Plotting functions |

Build locally:
```bash
pip install -r docs/requirements.txt
mkdocs serve
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

MIT License - see [LICENSE](LICENSE) for details.

## Citation

```bibtex
@software{neopkpd,
  title = {NeoPKPD: Transparent PK/PD Modeling Infrastructure},
  url = {https://github.com/neopkpd/neopkpd},
  version = {0.1.0}
}
```

---

<p align="center">
  <sub>Built for reproducible pharmacometric research</sub>
</p>
