# OpenPKPD Documentation Restructure Plan

## Executive Summary

Complete restructure of OpenPKPD documentation to provide comprehensive, industry-standard, modularized documentation with separate Julia and Python sections, complete API coverage, and professional visual design.

---

## Current State Analysis

### Documentation Structure
```
docs/
├── index.md          # Home page with quick start
├── models.md         # Combined PK/PD models reference (762 lines)
├── nca.md            # NCA reference
├── estimation.md     # Parameter estimation
├── population.md     # Population simulation
├── vpc.md            # VPC documentation
├── trial.md          # Clinical trial simulation
├── visualization.md  # Visualization module
├── python.md         # Python bindings (1155 lines)
├── cli.md            # CLI reference
├── data.md           # Data import (CDISC)
├── import.md         # Model import (NONMEM/Monolix)
├── architecture.md   # System architecture
├── semantics.md      # Versioning semantics
├── reproducibility.md # Artifact reproducibility
└── examples/         # 80+ example directories
```

### Coverage Gaps Identified

| Category | Julia Exports | Documented | Gap |
|----------|--------------|------------|-----|
| PK Models | 8 models | 7 models | TMDD missing |
| PD Models | 15+ models | 4 models | IRM types, Disease progression missing |
| Population | 50+ items | ~20 items | Covariates, IOV details missing |
| Estimation | 40+ items | ~15 items | FOCE variants, diagnostics missing |
| NCA | 30+ items | ~15 items | BE analysis incomplete |
| VPC | 25+ items | ~10 items | pcVPC, stratification missing |
| Import | 20+ items | ~8 items | Monolix incomplete |
| Trial | 30+ items | ~15 items | Adaptive designs missing |
| **TOTAL Julia** | **~450** | **~150** | **~300 items undocumented** |

| Category | Python Exports | Documented | Gap |
|----------|---------------|------------|-----|
| Core Simulation | 20+ functions | 14 | 6 missing |
| NCA | 50+ items | ~20 | 30 missing |
| Trial | 70+ items | ~30 | 40 missing |
| Visualization | 55+ functions | ~20 | 35 missing |
| Data Import | 30+ items | ~10 | 20 missing |
| **TOTAL Python** | **~388** | **~120** | **~268 items undocumented** |

---

## Target Structure

### New Navigation Hierarchy

```
Home
├── Introduction
│   ├── What is OpenPKPD?
│   ├── Key Features
│   ├── Architecture Overview
│   └── Getting Started
│
├── Docs Julia
│   ├── Tutorial: Getting Started with Julia
│   ├── Models
│   │   ├── Pharmacokinetic Models
│   │   │   ├── One-Compartment IV Bolus
│   │   │   ├── One-Compartment IV Infusion
│   │   │   ├── One-Compartment Oral
│   │   │   ├── Two-Compartment IV
│   │   │   ├── Two-Compartment Oral
│   │   │   ├── Three-Compartment IV
│   │   │   ├── Transit Absorption
│   │   │   ├── Michaelis-Menten Elimination
│   │   │   └── TMDD Models
│   │   └── Pharmacodynamic Models
│   │       ├── Direct Emax
│   │       ├── Sigmoid Emax (Hill)
│   │       ├── Effect Compartment (Biophase)
│   │       ├── Indirect Response Type I
│   │       ├── Indirect Response Type II
│   │       ├── Indirect Response Type III
│   │       ├── Indirect Response Type IV
│   │       ├── Turnover Models
│   │       └── Disease Progression
│   ├── Population Modeling
│   │   ├── Inter-Individual Variability (IIV)
│   │   ├── Inter-Occasion Variability (IOV)
│   │   ├── Covariate Models
│   │   ├── Residual Error Models
│   │   └── Population Simulation
│   ├── Non-Compartmental Analysis
│   │   ├── Exposure Metrics (Cmax, Tmax, AUC)
│   │   ├── Terminal Phase Analysis
│   │   ├── Bioequivalence Analysis
│   │   └── Population NCA
│   ├── Parameter Estimation
│   │   ├── FOCE-I Method
│   │   ├── SAEM Algorithm
│   │   ├── Laplacian Estimation
│   │   ├── Estimation Diagnostics
│   │   └── Model Comparison
│   ├── Visual Predictive Check
│   │   ├── Standard VPC
│   │   ├── Prediction-Corrected VPC
│   │   ├── Stratified VPC
│   │   └── VPC Diagnostics
│   ├── Model Import
│   │   ├── NONMEM Control Stream
│   │   ├── Monolix Project Files
│   │   └── CDISC Data Formats
│   ├── Clinical Trials
│   │   ├── Parallel Design
│   │   ├── Crossover Design
│   │   ├── Dose Escalation (3+3)
│   │   ├── Bioequivalence Studies
│   │   └── Power Analysis
│   └── API Reference (Julia)
│
├── Docs Python
│   ├── Tutorial: Getting Started with Python
│   ├── Models
│   │   ├── Pharmacokinetic Models
│   │   │   ├── simulate_pk_iv_bolus
│   │   │   ├── simulate_pk_oral_first_order
│   │   │   ├── simulate_pk_twocomp_iv_bolus
│   │   │   ├── simulate_pk_twocomp_oral
│   │   │   ├── simulate_pk_threecomp_iv_bolus
│   │   │   ├── simulate_pk_transit_absorption
│   │   │   └── simulate_pk_michaelis_menten
│   │   └── Pharmacodynamic Models
│   │       ├── simulate_pkpd_direct_emax
│   │       ├── simulate_pkpd_sigmoid_emax
│   │       ├── simulate_pkpd_biophase_equilibration
│   │       └── simulate_pkpd_indirect_response
│   ├── Population Modeling
│   │   ├── simulate_population_iv_bolus
│   │   ├── simulate_population_oral
│   │   ├── Covariate Specification
│   │   └── Population Summaries
│   ├── Non-Compartmental Analysis
│   │   ├── run_nca
│   │   ├── NCA Configuration
│   │   ├── Individual NCA Functions
│   │   ├── Population NCA
│   │   └── Bioequivalence Functions
│   ├── Parameter Estimation
│   │   ├── run_estimation
│   │   ├── Estimation Methods
│   │   ├── Results Interpretation
│   │   └── Diagnostics
│   ├── Clinical Trials
│   │   ├── Study Designs
│   │   ├── Dosing Regimens
│   │   ├── Virtual Population
│   │   ├── Trial Simulation
│   │   └── Power Analysis
│   ├── Visualization (55+ functions)
│   │   ├── Backends & Themes
│   │   ├── PK Plots
│   │   │   ├── plot_conc_time
│   │   │   ├── plot_multi_conc_time
│   │   │   ├── plot_spaghetti
│   │   │   ├── plot_mean_ribbon
│   │   │   └── plot_individual_fits
│   │   ├── NCA Plots
│   │   │   ├── plot_lambda_z_fit
│   │   │   ├── plot_auc_visualization
│   │   │   └── plot_dose_proportionality
│   │   ├── PKPD Plots
│   │   │   ├── plot_effect_conc
│   │   │   ├── plot_hysteresis
│   │   │   └── plot_dose_response
│   │   ├── VPC Plots
│   │   │   ├── plot_vpc_detailed
│   │   │   ├── plot_pcvpc
│   │   │   ├── plot_stratified_vpc
│   │   │   ├── plot_vpc_with_blq
│   │   │   └── plot_vpc_ci
│   │   ├── Estimation Diagnostics
│   │   │   ├── plot_parameter_estimates
│   │   │   ├── plot_omega_matrix
│   │   │   ├── plot_convergence
│   │   │   ├── plot_parameter_convergence
│   │   │   ├── plot_shrinkage
│   │   │   ├── plot_eta_distributions
│   │   │   ├── plot_individual_parameters
│   │   │   ├── plot_ofv_comparison
│   │   │   ├── plot_correlation_matrix
│   │   │   └── plot_sigma_residuals
│   │   ├── Bootstrap Plots
│   │   │   ├── plot_bootstrap_distributions
│   │   │   ├── plot_bootstrap_ci
│   │   │   ├── plot_bootstrap_stability
│   │   │   └── plot_bootstrap_correlation
│   │   ├── Sensitivity Plots
│   │   │   ├── plot_tornado
│   │   │   ├── plot_spider
│   │   │   ├── plot_sensitivity_heatmap
│   │   │   ├── plot_waterfall
│   │   │   └── plot_one_at_a_time
│   │   ├── Population Plots
│   │   │   ├── plot_vpc
│   │   │   ├── plot_parameter_distributions
│   │   │   ├── plot_forest
│   │   │   ├── plot_boxplot
│   │   │   ├── plot_goodness_of_fit
│   │   │   ├── plot_estimation_summary
│   │   │   ├── plot_sensitivity
│   │   │   └── plot_sensitivity_tornado
│   │   └── Trial Plots
│   │       ├── plot_power_curve
│   │       ├── plot_trial_tornado
│   │       ├── plot_kaplan_meier
│   │       └── plot_endpoint_distribution
│   ├── Data Import
│   │   ├── CDISC Formats
│   │   ├── CSV/XPT Import
│   │   └── Data Preparation
│   └── API Reference (Python)
│
├── Examples
│   ├── Quickstart Examples
│   ├── Model Examples
│   ├── Population Examples
│   ├── NCA Examples
│   ├── Estimation Examples
│   ├── VPC Examples
│   ├── Trial Examples
│   └── Real-World Validation
│
├── Concepts
│   ├── Architecture
│   ├── Semantics & Versioning
│   ├── Reproducibility
│   └── Numerical Methods
│
└── Reference
    ├── CLI Reference
    ├── Glossary
    └── Changelog
```

---

## New Files to Create

### Introduction Section (4 new files)
| File | Description | Lines Est. |
|------|-------------|-----------|
| `docs/intro/index.md` | What is OpenPKPD? | 150 |
| `docs/intro/features.md` | Key features overview | 200 |
| `docs/intro/architecture.md` | System architecture visual | 250 |
| `docs/intro/getting-started.md` | Installation & first steps | 300 |

### Julia Documentation (35+ new files)
| Category | Files | Lines Est. |
|----------|-------|-----------|
| Tutorial | `docs/julia/tutorial.md` | 500 |
| PK Models | 9 files in `docs/julia/models/pk/` | 2,700 |
| PD Models | 9 files in `docs/julia/models/pd/` | 2,700 |
| Population | 5 files in `docs/julia/population/` | 1,500 |
| NCA | 4 files in `docs/julia/nca/` | 1,200 |
| Estimation | 5 files in `docs/julia/estimation/` | 1,500 |
| VPC | 4 files in `docs/julia/vpc/` | 1,200 |
| Import | 3 files in `docs/julia/import/` | 900 |
| Trial | 5 files in `docs/julia/trial/` | 1,500 |
| API Reference | `docs/julia/api-reference.md` | 800 |
| **Subtotal** | **~45 files** | **~14,500 lines** |

### Python Documentation (40+ new files)
| Category | Files | Lines Est. |
|----------|-------|-----------|
| Tutorial | `docs/python/tutorial.md` | 500 |
| PK Models | 7 files in `docs/python/models/pk/` | 2,100 |
| PD Models | 4 files in `docs/python/models/pd/` | 1,200 |
| Population | 4 files in `docs/python/population/` | 1,200 |
| NCA | 5 files in `docs/python/nca/` | 1,500 |
| Estimation | 4 files in `docs/python/estimation/` | 1,200 |
| Trial | 5 files in `docs/python/trial/` | 1,500 |
| Visualization | 10 files in `docs/python/viz/` | 4,000 |
| Data Import | 3 files in `docs/python/data/` | 900 |
| API Reference | `docs/python/api-reference.md` | 1,000 |
| **Subtotal** | **~47 files** | **~15,100 lines** |

### Total New Content
- **~96 new documentation files**
- **~30,000 lines of documentation**

---

## Proposed mkdocs.yml Structure

```yaml
site_name: OpenPKPD
site_description: Transparent, validated PK/PD modeling infrastructure
site_author: OpenPKPD contributors

repo_name: openpkpd/OpenPKPD
repo_url: https://github.com/openpkpd/OpenPKPD
edit_uri: edit/main/docs/

theme:
  name: material
  logo: assets/logo.svg
  favicon: assets/favicon.png
  features:
    - navigation.tabs
    - navigation.tabs.sticky
    - navigation.sections
    - navigation.expand
    - navigation.path
    - navigation.top
    - navigation.indexes
    - toc.follow
    - toc.integrate
    - search.suggest
    - search.highlight
    - content.tabs.link
    - content.code.copy
    - content.code.annotate
    - content.code.select
  palette:
    - scheme: default
      primary: indigo
      accent: indigo
      toggle:
        icon: material/brightness-7
        name: Switch to dark mode
    - scheme: slate
      primary: indigo
      accent: indigo
      toggle:
        icon: material/brightness-4
        name: Switch to light mode
  font:
    text: Inter
    code: JetBrains Mono

nav:
  - Home: index.md

  - Introduction:
    - intro/index.md
    - Key Features: intro/features.md
    - Architecture: intro/architecture.md
    - Getting Started: intro/getting-started.md

  - Docs Julia:
    - julia/index.md
    - Tutorial: julia/tutorial.md
    - Models:
      - julia/models/index.md
      - Pharmacokinetic Models:
        - One-Compartment IV Bolus: julia/models/pk/onecomp-iv-bolus.md
        - One-Compartment IV Infusion: julia/models/pk/onecomp-iv-infusion.md
        - One-Compartment Oral: julia/models/pk/onecomp-oral.md
        - Two-Compartment IV: julia/models/pk/twocomp-iv.md
        - Two-Compartment Oral: julia/models/pk/twocomp-oral.md
        - Three-Compartment IV: julia/models/pk/threecomp-iv.md
        - Transit Absorption: julia/models/pk/transit-absorption.md
        - Michaelis-Menten: julia/models/pk/michaelis-menten.md
        - TMDD Models: julia/models/pk/tmdd.md
      - Pharmacodynamic Models:
        - Direct Emax: julia/models/pd/direct-emax.md
        - Sigmoid Emax: julia/models/pd/sigmoid-emax.md
        - Effect Compartment: julia/models/pd/effect-compartment.md
        - Indirect Response Type I: julia/models/pd/indirect-response-1.md
        - Indirect Response Type II: julia/models/pd/indirect-response-2.md
        - Indirect Response Type III: julia/models/pd/indirect-response-3.md
        - Indirect Response Type IV: julia/models/pd/indirect-response-4.md
        - Turnover Models: julia/models/pd/turnover.md
        - Disease Progression: julia/models/pd/disease-progression.md
    - Population Modeling:
      - julia/population/index.md
      - Inter-Individual Variability: julia/population/iiv.md
      - Inter-Occasion Variability: julia/population/iov.md
      - Covariate Models: julia/population/covariates.md
      - Residual Error Models: julia/population/residual-error.md
      - Population Simulation: julia/population/simulation.md
    - Non-Compartmental Analysis:
      - julia/nca/index.md
      - Exposure Metrics: julia/nca/exposure-metrics.md
      - Terminal Phase Analysis: julia/nca/terminal-phase.md
      - Bioequivalence: julia/nca/bioequivalence.md
      - Population NCA: julia/nca/population-nca.md
    - Parameter Estimation:
      - julia/estimation/index.md
      - FOCE-I Method: julia/estimation/foce.md
      - SAEM Algorithm: julia/estimation/saem.md
      - Laplacian: julia/estimation/laplacian.md
      - Diagnostics: julia/estimation/diagnostics.md
      - Model Comparison: julia/estimation/comparison.md
    - Visual Predictive Check:
      - julia/vpc/index.md
      - Standard VPC: julia/vpc/standard.md
      - Prediction-Corrected VPC: julia/vpc/pcvpc.md
      - Stratified VPC: julia/vpc/stratified.md
      - VPC Diagnostics: julia/vpc/diagnostics.md
    - Model Import:
      - julia/import/index.md
      - NONMEM: julia/import/nonmem.md
      - Monolix: julia/import/monolix.md
      - CDISC Data: julia/import/cdisc.md
    - Clinical Trials:
      - julia/trial/index.md
      - Parallel Design: julia/trial/parallel.md
      - Crossover Design: julia/trial/crossover.md
      - Dose Escalation: julia/trial/dose-escalation.md
      - Bioequivalence Studies: julia/trial/bioequivalence.md
      - Power Analysis: julia/trial/power.md
    - API Reference: julia/api-reference.md

  - Docs Python:
    - python/index.md
    - Tutorial: python/tutorial.md
    - Models:
      - python/models/index.md
      - Pharmacokinetic Models:
        - IV Bolus: python/models/pk/iv-bolus.md
        - Oral First-Order: python/models/pk/oral.md
        - Two-Compartment IV: python/models/pk/twocomp-iv.md
        - Two-Compartment Oral: python/models/pk/twocomp-oral.md
        - Three-Compartment: python/models/pk/threecomp.md
        - Transit Absorption: python/models/pk/transit.md
        - Michaelis-Menten: python/models/pk/michaelis-menten.md
      - Pharmacodynamic Models:
        - Direct Emax: python/models/pd/direct-emax.md
        - Sigmoid Emax: python/models/pd/sigmoid-emax.md
        - Effect Compartment: python/models/pd/effect-compartment.md
        - Indirect Response: python/models/pd/indirect-response.md
    - Population Modeling:
      - python/population/index.md
      - Population IV Bolus: python/population/iv-bolus.md
      - Population Oral: python/population/oral.md
      - Covariates: python/population/covariates.md
      - Summaries: python/population/summaries.md
    - Non-Compartmental Analysis:
      - python/nca/index.md
      - run_nca Function: python/nca/run-nca.md
      - NCA Configuration: python/nca/config.md
      - Individual Functions: python/nca/functions.md
      - Population NCA: python/nca/population.md
      - Bioequivalence: python/nca/bioequivalence.md
    - Parameter Estimation:
      - python/estimation/index.md
      - Estimation Methods: python/estimation/methods.md
      - Results: python/estimation/results.md
      - Diagnostics: python/estimation/diagnostics.md
    - Clinical Trials:
      - python/trial/index.md
      - Study Designs: python/trial/designs.md
      - Dosing Regimens: python/trial/dosing.md
      - Virtual Population: python/trial/population.md
      - Trial Simulation: python/trial/simulation.md
      - Power Analysis: python/trial/power.md
    - Visualization:
      - python/viz/index.md
      - Backends & Themes: python/viz/backends.md
      - PK Plots: python/viz/pk.md
      - NCA Plots: python/viz/nca.md
      - PKPD Plots: python/viz/pkpd.md
      - VPC Plots: python/viz/vpc.md
      - Estimation Diagnostics: python/viz/estimation.md
      - Bootstrap Plots: python/viz/bootstrap.md
      - Sensitivity Plots: python/viz/sensitivity.md
      - Population Plots: python/viz/population.md
      - Trial Plots: python/viz/trial.md
    - Data Import:
      - python/data/index.md
      - CDISC Formats: python/data/cdisc.md
      - Data Preparation: python/data/preparation.md
    - API Reference: python/api-reference.md

  - Examples:
    - examples/index.md
    - Quickstart: examples/quickstart/README.md
    - Model Examples:
      - examples/models/README.md
      - PK Models:
        - One-Compartment IV: examples/models/pk/01_onecomp_iv_bolus/README.md
        - One-Compartment Infusion: examples/models/pk/02_onecomp_iv_infusion/README.md
        - One-Compartment Oral: examples/models/pk/03_onecomp_oral/README.md
        - Two-Compartment IV: examples/models/pk/04_twocomp_iv/README.md
        - Two-Compartment Oral: examples/models/pk/05_twocomp_oral/README.md
        - Three-Compartment: examples/models/pk/06_threecomp_iv/README.md
        - Transit Absorption: examples/models/pk/07_transit_absorption/README.md
        - Michaelis-Menten: examples/models/pk/08_michaelis_menten/README.md
      - PKPD Models:
        - Direct Emax: examples/models/pkpd/01_direct_emax/README.md
        - Sigmoid Emax: examples/models/pkpd/02_sigmoid_emax/README.md
        - Biophase: examples/models/pkpd/03_biophase_equilibration/README.md
        - Indirect Response: examples/models/pkpd/04_indirect_response/README.md
    - Population Examples:
      - examples/population/README.md
    - NCA Examples:
      - examples/nca/README.md
    - Estimation Examples:
      - examples/estimation/README.md
    - VPC Examples:
      - examples/vpc/README.md
    - Trial Examples:
      - examples/trial/README.md
    - Visualization Examples:
      - examples/visualization/README.md
    - Sensitivity Examples:
      - examples/sensitivity/README.md
    - Import Examples:
      - examples/import/README.md
    - Real-World Validation:
      - examples/real_world_validation/SOURCE.md
      - Theophylline SD: examples/real_world_validation/studies/theophylline_theo_sd/README.md
      - Theophylline MD: examples/real_world_validation/studies/theophylline_theo_md/README.md
      - Warfarin PKPD: examples/real_world_validation/studies/warfarin_pkpd/README.md
    - Use Cases:
      - Theophylline Analysis: examples/use_cases/real_world_theophylline/README.md
      - Bioequivalence Study: examples/use_cases/bioequivalence_study/README.md
      - Population PKPD: examples/use_cases/population_pkpd_analysis/README.md
      - NONMEM Migration: examples/use_cases/nonmem_migration/README.md

  - Concepts:
    - concepts/index.md
    - Architecture: architecture.md
    - Semantics: semantics.md
    - Reproducibility: reproducibility.md

  - Reference:
    - CLI Reference: cli.md
    - Glossary: reference/glossary.md
    - Changelog: reference/changelog.md

plugins:
  - search
  - mkdocstrings:
      handlers:
        python:
          options:
            show_source: true
            show_root_heading: true
  - git-revision-date-localized:
      enable_creation_date: true
  - minify:
      minify_html: true

markdown_extensions:
  - admonition
  - pymdownx.details
  - pymdownx.superfences:
      custom_fences:
        - name: mermaid
          class: mermaid
          format: !!python/name:pymdownx.superfences.fence_code_format
  - pymdownx.highlight:
      anchor_linenums: true
      line_spans: __span
      pygments_lang_class: true
  - pymdownx.inlinehilite
  - pymdownx.snippets
  - pymdownx.tabbed:
      alternate_style: true
  - pymdownx.arithmatex:
      generic: true
  - attr_list
  - md_in_html
  - tables
  - footnotes
  - toc:
      permalink: true
      toc_depth: 3

extra_javascript:
  - https://polyfill.io/v3/polyfill.min.js?features=es6
  - https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js

extra_css:
  - stylesheets/extra.css

extra:
  version:
    provider: mike
  analytics:
    provider: google
    property: G-XXXXXXXXXX
  social:
    - icon: fontawesome/brands/github
      link: https://github.com/openpkpd/OpenPKPD
    - icon: fontawesome/brands/python
      link: https://pypi.org/project/openpkpd/
```

---

## Page Templates

### Model Page Template (Julia)
```markdown
# [Model Name]

## Overview
Brief description of the model and its clinical applications.

## Model Structure

### Diagram
[Mermaid diagram of compartments]

### Parameters
| Parameter | Symbol | Units | Description | Constraints |
|-----------|--------|-------|-------------|-------------|
| CL | CL | L/h | Clearance | CL > 0 |

### State Variables
| State | Symbol | Description |
|-------|--------|-------------|
| A_central | A | Amount in central |

### Differential Equations
$$\frac{dA}{dt} = -\frac{CL}{V} \cdot A$$

## Julia API

### Type Definition
```julia
struct ModelNameParams
    cl::Float64
    v::Float64
end
```

### Simulation Example
```julia
using OpenPKPDCore

params = ModelNameParams(5.0, 50.0)
# ... complete example
```

## Clinical Applications
- Use case 1
- Use case 2

## Related Models
- [Link to related model 1]
- [Link to related model 2]

## References
1. Academic reference
```

### Function Page Template (Python)
```markdown
# function_name

## Signature
```python
def function_name(
    param1: Type,
    param2: Type,
    *,
    optional1: Type = default,
) -> ReturnType:
```

## Description
Detailed description of what the function does.

## Parameters
| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `param1` | `float` | Yes | Description |

## Returns
| Field | Type | Description |
|-------|------|-------------|
| `result` | `dict` | Description |

## Examples

### Basic Usage
```python
import openpkpd

result = openpkpd.function_name(...)
```

### Advanced Usage
```python
# More complex example
```

## Notes
- Important note 1
- Important note 2

## See Also
- [`related_function`](related.md)

## Source
[View source on GitHub](link)
```

---

## Implementation Phases

### Phase 1: Structure Setup (Week 1)
- [ ] Create new directory structure under `docs/`
- [ ] Update `mkdocs.yml` with new navigation
- [ ] Create index pages for all sections
- [ ] Add custom CSS for professional styling

### Phase 2: Introduction & Getting Started (Week 1)
- [ ] Write `intro/index.md` - What is OpenPKPD
- [ ] Write `intro/features.md` - Feature overview
- [ ] Write `intro/architecture.md` - Architecture diagram
- [ ] Write `intro/getting-started.md` - Installation guide

### Phase 3: Julia Tutorial & Models (Week 2)
- [ ] Write `julia/tutorial.md` - Complete tutorial
- [ ] Create all 9 PK model pages
- [ ] Create all 9 PD model pages
- [ ] Add model diagrams using Mermaid

### Phase 4: Julia Advanced Topics (Week 3)
- [ ] Write 5 population modeling pages
- [ ] Write 4 NCA pages
- [ ] Write 5 estimation pages
- [ ] Write 4 VPC pages
- [ ] Write 3 import pages
- [ ] Write 5 trial pages

### Phase 5: Python Tutorial & Models (Week 4)
- [ ] Write `python/tutorial.md` - Complete tutorial
- [ ] Create all 7 PK function pages
- [ ] Create all 4 PD function pages
- [ ] Link to Julia concepts where appropriate

### Phase 6: Python Advanced Topics (Week 5)
- [ ] Write 4 population pages
- [ ] Write 5 NCA pages
- [ ] Write 4 estimation pages
- [ ] Write 5 trial pages
- [ ] Write 3 data import pages

### Phase 7: Visualization Documentation (Week 6)
- [ ] Write viz backends page
- [ ] Write 5 viz section pages (PK, NCA, PKPD, VPC, etc.)
- [ ] Document all 55+ visualization functions
- [ ] Add example plots with outputs

### Phase 8: Polish & Review (Week 7)
- [ ] API reference generation
- [ ] Cross-linking between pages
- [ ] Add glossary
- [ ] Changelog
- [ ] Final review and testing

---

## Visual Design Guidelines

### Color Palette
- **Primary**: Indigo (#3F51B5)
- **Accent**: Deep Purple (#673AB7)
- **Julia**: Purple (#9558B2)
- **Python**: Blue (#3776AB)
- **Success**: Green (#4CAF50)
- **Warning**: Orange (#FF9800)
- **Error**: Red (#F44336)

### Typography
- **Headings**: Inter (sans-serif)
- **Body**: Inter (sans-serif)
- **Code**: JetBrains Mono

### Code Blocks
- Julia: Purple background tint
- Python: Blue background tint
- Shell: Gray background

### Admonitions
- **Note**: Blue info box
- **Tip**: Green success box
- **Warning**: Orange warning box
- **Danger**: Red error box
- **Example**: Purple example box

---

## Success Criteria

### Coverage
- [ ] 100% of Julia core exports documented
- [ ] 100% of Python exports documented
- [ ] All 55+ visualization functions documented
- [ ] All 8 PK models with dedicated pages
- [ ] All 9+ PD models with dedicated pages

### Quality
- [ ] Every function has signature, description, parameters, returns, examples
- [ ] All models have differential equations and diagrams
- [ ] Cross-references between related topics
- [ ] No broken links

### User Experience
- [ ] New user can start in <10 minutes
- [ ] API reference is searchable
- [ ] Examples are copy-pasteable
- [ ] Mobile-responsive

---

## Next Steps

1. **User Approval**: Review and approve this plan
2. **Phase 1 Start**: Create directory structure and updated mkdocs.yml
3. **Content Creation**: Begin writing documentation files in phases
4. **Continuous Validation**: Test documentation builds after each phase
