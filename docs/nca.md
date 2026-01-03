# Non-Compartmental Analysis (NCA) Reference

OpenPKPD provides comprehensive Non-Compartmental Analysis (NCA) functionality for pharmacokinetic analysis. The NCA module is FDA/EMA compliant and supports both single-dose and multiple-dose studies.

## Overview

Non-Compartmental Analysis is a model-independent approach for calculating pharmacokinetic parameters from concentration-time data. OpenPKPD's NCA module provides:

- **Exposure Metrics**: Cmax, Tmax, Cmin, Clast, Tlast, Cavg
- **AUC Calculations**: AUC0-t, AUC0-inf, AUC0-tau (multiple dose)
- **Lambda-z Estimation**: Terminal slope with R² quality assessment
- **Secondary Parameters**: t½, MRT, CL/F, Vz/F, Vss
- **Multiple Dose Metrics**: Accumulation index, PTF, Swing
- **Bioequivalence**: 90% CI, TOST, GMR analysis

---

## Python API

### Quick Start

```python
import openpkpd
from openpkpd.nca import run_nca, NCAConfig

# Initialize Julia
openpkpd.init_julia()

# Sample concentration-time data
times = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0]
conc = [0.0, 0.82, 1.44, 1.62, 1.28, 0.94, 0.68, 0.36, 0.08]
dose = 100.0

# Run full NCA
result = run_nca(times, conc, dose)

print(f"Cmax: {result.cmax:.2f}")
print(f"Tmax: {result.tmax:.1f} h")
print(f"AUC0-t: {result.auc_0_t:.2f} h*mg/L")
print(f"t½: {result.t_half:.2f} h")
```

---

## Exposure Metrics

### Primary Metrics

#### Cmax - Maximum Concentration

```python
from openpkpd.nca import nca_cmax

cmax = nca_cmax(conc)
# Returns: 1.62 (maximum value in concentration vector)
```

#### Tmax - Time to Maximum Concentration

```python
from openpkpd.nca import nca_tmax

tmax = nca_tmax(times, conc)
# Returns: 2.0 (time at which Cmax occurs)
```

#### Cmin - Minimum Concentration

Used for multiple-dose studies at steady state.

```python
from openpkpd.nca import nca_cmin

# Steady-state profile
ss_conc = [2.5, 4.2, 3.8, 3.0, 2.8, 2.6, 2.5, 2.5]
cmin = nca_cmin(ss_conc)
# Returns: 2.5
```

#### Clast and Tlast - Last Quantifiable Concentration

```python
from openpkpd.nca import nca_clast, nca_tlast

clast = nca_clast(times, conc)
tlast = nca_tlast(times, conc)

# With LLOQ
clast = nca_clast(times, conc, lloq=0.1)  # Excludes values below 0.1
tlast = nca_tlast(times, conc, lloq=0.1)
```

---

## AUC Calculations

### AUC0-t (AUC from Zero to Last Observation)

```python
from openpkpd.nca import auc_0_t

# Linear trapezoidal method
auc = auc_0_t(times, conc, method="linear")

# Log-linear trapezoidal method
auc = auc_0_t(times, conc, method="log_linear")

# Lin-log mixed (linear ascending, log-linear descending)
auc = auc_0_t(times, conc, method="lin_log_mixed")
```

**Method Selection Guidelines**:

| Method | When to Use | FDA Guidance |
|--------|-------------|--------------|
| `linear` | Simple estimation, IV bolus | Acceptable |
| `log_linear` | Monoexponential decline | Preferred for elimination phase |
| `lin_log_mixed` | Standard oral PK | Recommended for absorption + elimination |

### AUC Extrapolation to Infinity

Automatically calculated when lambda-z can be estimated:

```python
result = run_nca(times, conc, dose)

print(f"AUC0-inf: {result.auc_0_inf}")
print(f"Extrapolation %: {result.auc_extra_pct}%")

# FDA requires extrapolation < 20% for reliable AUC0-inf
if result.auc_extra_pct < 20:
    print("Acceptable extrapolation")
```

### AUC0-tau (Multiple Dose)

For steady-state studies with defined dosing interval:

```python
ss_times = [0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]
ss_conc = [2.5, 4.2, 3.8, 3.0, 2.8, 2.6, 2.5, 2.5]
tau = 12.0  # Dosing interval

result = run_nca(ss_times, ss_conc, dose, dosing_type="multiple", tau=tau)
print(f"AUC0-tau: {result.auc_0_tau}")
print(f"Cavg: {result.cavg}")  # AUC0-tau / tau
```

---

## Terminal Phase Analysis

### Lambda-z Estimation

The terminal elimination rate constant is estimated by log-linear regression of the terminal phase.

```python
from openpkpd.nca import estimate_lambda_z

result = estimate_lambda_z(times, conc)

print(f"Lambda-z: {result['lambda_z']:.4f} /h")
print(f"t½: {result['t_half']:.2f} h")
print(f"R²: {result['r_squared']:.4f}")
print(f"Points used: {result['n_points']}")
print(f"Quality: {result['quality_flag']}")
```

**Quality Flags**:

| Flag | Meaning | R² Threshold |
|------|---------|--------------|
| `good` | Reliable estimate | R² ≥ 0.95 |
| `warning` | Use with caution | 0.90 ≤ R² < 0.95 |
| `insufficient` | Not reliable | R² < 0.90 or < 3 points |

### Half-Life

```python
from openpkpd.nca import nca_half_life
import math

lambda_z = 0.1  # /h
t_half = nca_half_life(lambda_z)
# t½ = ln(2) / λz = 0.693 / 0.1 = 6.93 h
```

---

## Secondary PK Parameters

### Mean Residence Time (MRT)

```python
from openpkpd.nca import nca_mrt

# For extravascular administration
mrt = nca_mrt(aumc, auc)  # AUMC / AUC

# For IV infusion (adjusted for infusion time)
mrt = nca_mrt(aumc, auc, route="iv_infusion", t_inf=1.0)
```

### Clearance and Volume

```python
from openpkpd.nca import nca_cl_f, nca_vz_f

dose = 100.0
auc = 20.0
lambda_z = 0.1

# Apparent clearance: CL/F = Dose / AUC
cl_f = nca_cl_f(dose, auc)
# Returns: 5.0 L/h

# Apparent volume: Vz/F = CL/F / λz
vz_f = nca_vz_f(dose, lambda_z, auc)
# Returns: 50.0 L
```

---

## Multiple Dose Metrics

### Accumulation Index

```python
from openpkpd.nca import nca_accumulation_index

auc_ss = 25.0   # Steady-state AUC
auc_sd = 20.0   # Single-dose AUC

rac = nca_accumulation_index(auc_ss, auc_sd)
# Returns: 1.25 (25% accumulation)
```

### Peak-Trough Fluctuation (PTF)

```python
from openpkpd.nca import nca_ptf

cmax = 4.2
cmin = 2.5
cavg = 3.0

ptf = nca_ptf(cmax, cmin, cavg)
# Returns: 100 * (Cmax - Cmin) / Cavg = 56.7%
```

### Swing

```python
from openpkpd.nca import nca_swing

swing = nca_swing(cmax, cmin)
# Returns: 100 * (Cmax - Cmin) / Cmin = 68.0%
```

---

## Bioequivalence Analysis

OpenPKPD provides comprehensive bioequivalence analysis tools for regulatory submissions.

### Geometric Mean Ratio (GMR)

```python
from openpkpd.nca import geometric_mean_ratio, geometric_mean

test_values = [20.0, 22.0, 18.0, 25.0, 21.0, 19.0, 23.0, 20.0]
ref_values = [19.0, 21.0, 17.0, 24.0, 20.0, 18.0, 22.0, 19.0]

gmr = geometric_mean_ratio(test_values, ref_values)
print(f"GMR: {gmr:.4f}")  # ~1.05
```

### 90% Confidence Interval

```python
from openpkpd.nca import bioequivalence_90ci

result = bioequivalence_90ci(test_values, ref_values)

print(f"GMR: {result['gmr']:.4f}")
print(f"90% CI: [{result['ci_lower']:.4f}, {result['ci_upper']:.4f}]")
print(f"Intra-subject CV: {result['cv_intra']:.1f}%")
```

### TOST Analysis

Two One-Sided Tests for bioequivalence:

```python
from openpkpd.nca import tost_analysis

result = tost_analysis(test_values, ref_values)

print(f"Conclusion: {result['conclusion']}")
# Returns: "bioequivalent" or "not_bioequivalent"
```

### BE Conclusion Helper

```python
from openpkpd.nca import be_conclusion

# Standard limits (80-125%)
conclusion = be_conclusion(ci_lower=0.85, ci_upper=1.20)
# Returns: "bioequivalent"

# Highly variable drugs (EMA)
conclusion = be_conclusion(ci_lower=0.72, ci_upper=1.35,
                          lower_limit=0.6984, upper_limit=1.4319)
```

---

## NCA Configuration

Customize NCA calculations with `NCAConfig`:

```python
from openpkpd.nca import run_nca, NCAConfig

config = NCAConfig(
    method="lin_log_mixed",       # AUC calculation method
    lambda_z_min_points=3,        # Minimum points for λz
    lambda_z_r2_threshold=0.90,   # R² threshold for quality
    extrapolation_max_pct=20.0,   # Maximum acceptable extrapolation
    blq_handling="zero",          # Handle below LLOQ: "zero", "missing", "lloq_half"
)

result = run_nca(times, conc, dose, config=config)
```

---

## Full NCA Workflow

### Single-Dose Study

```python
import openpkpd
from openpkpd.nca import run_nca

openpkpd.init_julia()

# PK data from oral administration
times = [0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24]
conc = [0, 0.5, 1.2, 2.1, 2.4, 2.5, 2.2, 1.9, 1.4, 1.0, 0.5, 0.25, 0.06]
dose = 100.0  # mg

result = run_nca(times, conc, dose, route="oral")

# Print NCA report
print("=" * 50)
print("NON-COMPARTMENTAL ANALYSIS REPORT")
print("=" * 50)
print(f"\n{'Parameter':<20} {'Value':<15} {'Unit':<10}")
print("-" * 50)
print(f"{'Cmax':<20} {result.cmax:<15.3f} {'mg/L':<10}")
print(f"{'Tmax':<20} {result.tmax:<15.1f} {'h':<10}")
print(f"{'Clast':<20} {result.clast:<15.3f} {'mg/L':<10}")
print(f"{'Tlast':<20} {result.tlast:<15.1f} {'h':<10}")
print(f"{'AUC0-t':<20} {result.auc_0_t:<15.3f} {'h*mg/L':<10}")
if result.auc_0_inf:
    print(f"{'AUC0-inf':<20} {result.auc_0_inf:<15.3f} {'h*mg/L':<10}")
    print(f"{'AUC extrap %':<20} {result.auc_extra_pct:<15.1f} {'%':<10}")
if result.t_half:
    print(f"{'t½':<20} {result.t_half:<15.2f} {'h':<10}")
if result.cl_f:
    print(f"{'CL/F':<20} {result.cl_f:<15.2f} {'L/h':<10}")
if result.vz_f:
    print(f"{'Vz/F':<20} {result.vz_f:<15.1f} {'L':<10}")

# Quality assessment
print(f"\nQuality Flags: {result.quality_flags}")
if result.warnings:
    print(f"Warnings: {result.warnings}")
```

### Multiple-Dose Steady-State Study

```python
# Steady-state data (tau = 12 hours)
ss_times = [0, 0.5, 1, 2, 4, 6, 8, 10, 12]
ss_conc = [2.0, 3.5, 4.0, 3.6, 3.0, 2.6, 2.3, 2.1, 2.0]
dose = 50.0  # mg
tau = 12.0   # hours

result = run_nca(
    ss_times, ss_conc, dose,
    dosing_type="multiple",
    tau=tau,
    route="oral"
)

print(f"Cmax,ss: {result.cmax}")
print(f"Cmin,ss: {result.cmin}")
print(f"Cavg,ss: {result.cavg}")
print(f"AUC0-tau: {result.auc_0_tau}")
print(f"PTF: {result.ptf}%")
print(f"Swing: {result.swing}%")
```

---

## Julia API

### Basic NCA

```julia
using OpenPKPDCore

# Run full NCA
t = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0]
c = [0.0, 0.82, 1.44, 1.62, 1.28, 0.94, 0.68, 0.36, 0.08]
dose = 100.0

result = run_nca(t, c, dose)

println("Cmax: ", result.cmax)
println("AUC0-t: ", result.auc_0_t)
```

### NCA with Configuration

```julia
config = NCAConfig(
    method = :lin_log_mixed,
    lambda_z_min_points = 4,
    lambda_z_r2_threshold = 0.90
)

result = run_nca(t, c, dose; config = config)
```

---

## Regulatory Compliance

### FDA Guidance

OpenPKPD's NCA module follows FDA guidance for:

- **Bioequivalence Studies**: Log-transformation, 90% CI calculation
- **Lambda-z Estimation**: Minimum 3 points in terminal phase
- **AUC Extrapolation**: Warning if > 20% of AUC0-inf
- **Method Selection**: Lin-log mixed for oral absorption profiles

### EMA Guidance

Support for EMA-specific requirements:

- **Highly Variable Drugs**: Scaled average BE approach
- **Narrow Therapeutic Index**: Tightened limits (90.00-111.11%)
- **Reference-Scaled Average Bioequivalence**: When CVintra > 30%

---

## Quality Control

### Result Validation

```python
result = run_nca(times, conc, dose)

# Check quality flags
if 'insufficient_terminal_phase' in result.quality_flags:
    print("Warning: Lambda-z estimation may be unreliable")

if result.auc_extra_pct and result.auc_extra_pct > 20:
    print("Warning: High AUC extrapolation - interpret AUC0-inf with caution")

# Check warnings
for warning in result.warnings:
    print(f"Warning: {warning}")
```

### Export to Dictionary

```python
result_dict = result.to_dict()

# Contains all parameters with metadata
import json
print(json.dumps(result_dict, indent=2))
```

---

## Common Issues and Solutions

| Issue | Cause | Solution |
|-------|-------|----------|
| Lambda-z not calculated | < 3 points in terminal phase | Extend sampling time |
| Poor R² for terminal phase | Multi-exponential decline | Use only terminal points |
| High AUC extrapolation | Sampling stopped too early | Sample to 4-5 half-lives |
| Negative lambda-z | Rising concentrations | Check for enterohepatic recirculation |

---

## See Also

- [Models Reference](models.md) - PK model documentation
- [Visualization](visualization.md) - NCA plots
- [Trial Simulation](trial.md) - Using NCA in trial simulations
