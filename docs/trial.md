# Clinical Trial Simulation Reference

OpenPKPD provides comprehensive clinical trial simulation capabilities for pharmacometric studies, including study design, dosing regimens, virtual population generation, and statistical analysis.

## Overview

The trial simulation module supports:

- **Study Designs**: Parallel, crossover, dose-escalation, adaptive, bioequivalence
- **Dosing Regimens**: QD, BID, TID, QID, titration, custom schedules
- **Virtual Populations**: Demographic modeling, disease states, covariates
- **Trial Events**: Enrollment, dropout, compliance
- **Statistical Analysis**: Power analysis, sample size, alpha spending, BE analysis

---

## Quick Start

```python
from openpkpd import trial

# Create a parallel design study
design = trial.parallel_design(2)

# Create dosing regimens
placebo = trial.dosing_qd(0.0, 28)
active = trial.dosing_qd(100.0, 28)

# Create treatment arms
arms = [
    trial.TreatmentArm("Placebo", placebo, 50, placebo=True),
    trial.TreatmentArm("Active", active, 50),
]

# Create trial specification
spec = trial.TrialSpec(
    name="Phase 2 Study",
    design=design,
    arms=arms,
    duration_days=28,
    seed=42
)

# Run simulation
result = trial.simulate_trial(spec)

print(f"Completion rate: {result.overall_completion_rate:.1%}")
```

---

## Study Designs

### Parallel Group Design

Standard parallel group design with independent treatment arms.

```python
from openpkpd.trial import parallel_design, ParallelDesign

# Two-arm parallel design (1:1 randomization)
design = parallel_design(2)

# Three-arm with custom randomization (2:1:1)
design = parallel_design(3, randomization_ratio=[0.5, 0.25, 0.25])

# With stratification factors
design = parallel_design(
    2,
    stratification_factors=["age_group", "disease_severity"]
)
```

**Attributes**:

| Attribute | Type | Description |
|-----------|------|-------------|
| `n_arms` | int | Number of treatment arms |
| `randomization_ratio` | List[float] | Allocation ratio (sums to 1.0) |
| `stratification_factors` | List[str] | Stratification variables |

---

### Crossover Designs

#### 2x2 Crossover (AB, BA)

```python
from openpkpd.trial import crossover_2x2

design = crossover_2x2(washout_duration=14.0)
# Sequences: [[1, 2], [2, 1]]
```

#### 3x3 Latin Square

```python
from openpkpd.trial import crossover_3x3

design = crossover_3x3(washout_duration=21.0)
# Sequences: [[1, 2, 3], [2, 3, 1], [3, 1, 2]]
```

#### Williams Design

Balanced for first-order carryover effects.

```python
from openpkpd.trial import williams_design

# 2-treatment Williams (same as 2x2)
design = williams_design(2, washout_duration=7.0)

# 3-treatment Williams (6 sequences)
design = williams_design(3, washout_duration=14.0)

# 4-treatment Williams (8 sequences)
design = williams_design(4, washout_duration=21.0)
```

**Crossover Attributes**:

| Attribute | Type | Description |
|-----------|------|-------------|
| `n_periods` | int | Number of treatment periods |
| `n_sequences` | int | Number of treatment sequences |
| `washout_duration` | float | Washout period in days |
| `sequence_assignments` | List[List[int]] | Treatment order per sequence |

---

### Dose Escalation Designs

For Phase I first-in-human studies.

#### 3+3 Design

Traditional rule-based escalation.

```python
from openpkpd.trial import dose_escalation_3plus3

dose_levels = [10.0, 25.0, 50.0, 100.0, 200.0]
design = dose_escalation_3plus3(
    dose_levels,
    starting_dose=10.0,
    max_dlt_rate=0.33,
    cohort_size=3,
    max_subjects=30
)
```

#### mTPI (Modified Toxicity Probability Interval)

```python
from openpkpd.trial import dose_escalation_mtpi

design = dose_escalation_mtpi(
    dose_levels=[10.0, 25.0, 50.0, 100.0],
    target_dlt_rate=0.25,
    cohort_size=3
)
```

#### CRM (Continual Reassessment Method)

Model-based escalation.

```python
from openpkpd.trial import dose_escalation_crm

design = dose_escalation_crm(
    dose_levels=[10.0, 25.0, 50.0, 100.0, 200.0],
    target_dlt_rate=0.25,
    cohort_size=1  # CRM typically uses single-patient cohorts
)
```

**Escalation Attributes**:

| Attribute | Type | Description |
|-----------|------|-------------|
| `dose_levels` | List[float] | Available dose levels |
| `starting_dose` | float | Initial dose (default: first level) |
| `escalation_rule` | str | "3+3", "mTPI", or "CRM" |
| `cohort_size` | int | Subjects per cohort |
| `target_dlt_rate` | float | Target DLT probability |

---

### Adaptive Designs

Designs with pre-planned interim analyses.

```python
from openpkpd.trial import adaptive_design, parallel_design

base = parallel_design(2)

design = adaptive_design(
    base,
    interim_analyses=[0.5, 0.75],  # At 50% and 75% information
    alpha_spending="obrien_fleming",
    futility_boundary=0.10
)
```

**Alpha Spending Functions**:

| Function | Description | Early Stopping |
|----------|-------------|----------------|
| `obrien_fleming` | Conservative early, liberal late | Rare early stops |
| `pocock` | Equal spending across analyses | More early stops |
| `haybittle_peto` | Fixed boundary (p<0.001) for interim | Very conservative |

---

### Bioequivalence Designs

```python
from openpkpd.trial import bioequivalence_design

# FDA standard BE study
design = bioequivalence_design(
    n_periods=2,
    n_sequences=2,
    washout_duration=7.0,
    bioequivalence_limits=(0.80, 1.25),
    parameters=["cmax", "auc_0_inf"],
    regulatory_guidance="fda"
)

# EMA highly variable drugs
design = bioequivalence_design(
    bioequivalence_limits=(0.6984, 1.4319),
    regulatory_guidance="ema"
)
```

---

## Dosing Regimens

### Standard Regimens

```python
from openpkpd.trial import dosing_qd, dosing_bid, dosing_tid, dosing_qid

# Once daily - 100 mg for 28 days
regimen = dosing_qd(100.0, 28)

# Twice daily - 50 mg for 14 days
regimen = dosing_bid(50.0, 14)

# Three times daily - 25 mg for 7 days
regimen = dosing_tid(25.0, 7)

# Four times daily - 20 mg for 5 days
regimen = dosing_qid(20.0, 5)

# With loading dose
regimen = dosing_qd(100.0, 28, loading_dose=200.0)
```

**Default Dose Times** (hours after midnight):

| Frequency | Dose Times |
|-----------|------------|
| QD | [8:00] |
| BID | [8:00, 20:00] |
| TID | [8:00, 14:00, 20:00] |
| QID | [8:00, 12:00, 16:00, 20:00] |

### Custom Regimens

```python
from openpkpd.trial import dosing_custom

# Custom timing: 6 AM and 10 PM
regimen = dosing_custom(75.0, 21, dose_times=[6.0, 22.0])
```

### Titration Regimens

```python
from openpkpd.trial import titration_regimen

# Titrate from 25 to 100 mg in 4 steps, 7 days each
regimen = titration_regimen(
    start_dose=25.0,
    target_dose=100.0,
    n_steps=4,
    days_per_step=7
)

# With maintenance phase
regimen = titration_regimen(
    start_dose=25.0,
    target_dose=100.0,
    n_steps=4,
    days_per_step=7,
    maintenance_days=28  # 28 days at target dose
)
```

### Dose Event Calculations

```python
from openpkpd.trial import dose_event_times, total_regimen_duration

regimen = dosing_bid(50.0, 3)

# Get all dose times
times = dose_event_times(regimen)
# [8.0, 20.0, 32.0, 44.0, 56.0, 68.0]

# Get total duration
duration = total_regimen_duration(regimen)
# 3 (days)
```

---

## Virtual Population Generation

### Basic Population

```python
from openpkpd.trial import generate_virtual_population

# Generate 100 subjects with default demographics
population = generate_virtual_population(100, seed=42)

for subject in population[:3]:
    print(f"ID: {subject.id}, Age: {subject.age:.1f}, "
          f"Weight: {subject.weight:.1f} kg, Sex: {subject.sex}")
```

### Demographic Specifications

```python
from openpkpd.trial import DemographicSpec, generate_virtual_population

# Custom demographics
demo = DemographicSpec(
    age_mean=50.0,
    age_sd=12.0,
    age_range=(30.0, 70.0),
    weight_mean=80.0,
    weight_sd=15.0,
    weight_range=(50.0, 120.0),
    female_proportion=0.45,
    race_distribution={
        "caucasian": 0.65,
        "asian": 0.15,
        "black": 0.15,
        "hispanic": 0.05
    }
)

population = generate_virtual_population(100, demographics=demo, seed=42)
```

### Pre-defined Populations

```python
from openpkpd.trial import (
    default_demographic_spec,
    healthy_volunteer_spec,
    patient_population_spec
)

# Healthy volunteers (Phase I)
demo = healthy_volunteer_spec()
# Age 30±8, 18-45 years, 72±10 kg

# Disease-specific populations
demo, disease = patient_population_spec("diabetes")
demo, disease = patient_population_spec("renal")
demo, disease = patient_population_spec("hepatic")
demo, disease = patient_population_spec("oncology")
```

### Population with Disease State

```python
demo, disease = patient_population_spec("diabetes")

population = generate_virtual_population(
    100,
    demographics=demo,
    disease=disease,
    seed=42
)

for subject in population[:3]:
    print(f"ID: {subject.id}, Severity: {subject.disease_severity}, "
          f"Biomarker: {subject.baseline_biomarker:.1f}")
```

### Population Summary

```python
from openpkpd.trial import summarize_population

summary = summarize_population(population)

print(f"N: {summary['n']}")
print(f"Age: {summary['age_mean']:.1f} ± {summary['age_sd']:.1f}")
print(f"Weight: {summary['weight_mean']:.1f} ± {summary['weight_sd']:.1f}")
print(f"Female: {summary['female_proportion']:.1%}")
print(f"Race: {summary['race_distribution']}")
```

---

## Trial Simulation

### Trial Specification

```python
from openpkpd.trial import (
    TrialSpec, TreatmentArm, parallel_design, dosing_qd,
    VirtualPopulationSpec, DemographicSpec,
    DropoutSpec, ComplianceSpec
)

spec = TrialSpec(
    name="Phase 2 Study",
    design=parallel_design(2),
    arms=[
        TreatmentArm("Placebo", dosing_qd(0.0, 28), 50, placebo=True),
        TreatmentArm("Active", dosing_qd(100.0, 28), 50),
    ],
    population_spec=VirtualPopulationSpec(
        demographics=DemographicSpec(age_mean=55.0)
    ),
    duration_days=28,
    enrollment_rate=5.0,  # Subjects per day
    dropout=DropoutSpec(random_rate_per_day=0.005),
    compliance=ComplianceSpec(mean_compliance=0.90, pattern="decay"),
    pk_sampling_times=[0, 1, 2, 4, 8, 12, 24],
    endpoints=["pk_exposure"],
    n_replicates=1,
    seed=42
)
```

### Running Simulations

```python
from openpkpd.trial import simulate_trial, simulate_trial_replicates

# Single simulation
result = simulate_trial(spec)

print(f"Trial: {result.trial_name}")
print(f"Completion: {result.overall_completion_rate:.1%}")
print(f"Compliance: {result.overall_compliance:.1%}")

# Access arm results
for arm_name, arm_result in result.arms.items():
    print(f"\n{arm_name}:")
    print(f"  Enrolled: {arm_result.n_enrolled}")
    print(f"  Completed: {arm_result.n_completed}")
    print(f"  Compliance: {arm_result.mean_compliance:.1%}")

# Multiple replicates
results = simulate_trial_replicates(spec, n_replicates=100)
```

### Dropout Modeling

```python
from openpkpd.trial import simulate_dropout, DropoutSpec

# Random dropout at 1% per day
spec = DropoutSpec(random_rate_per_day=0.01)

dropout_days = simulate_dropout(100, 28.0, spec, seed=42)
n_dropouts = sum(1 for d in dropout_days if d is not None)
print(f"Dropouts: {n_dropouts}/100")
```

### Compliance Modeling

```python
from openpkpd.trial import apply_compliance, ComplianceSpec

times = list(range(0, 168, 24))  # 7 days
amounts = [100.0] * 7

# Random non-compliance
spec = ComplianceSpec(mean_compliance=0.80, pattern="random")

# Patterns:
# - "random": Random misses
# - "weekend_miss": Higher miss rate on weekends
# - "decay": Compliance decreases over time
# - "early_good": Better compliance early in study

actual = apply_compliance(times, amounts, spec, seed=42)
taken = sum(1 for d in actual if d > 0)
print(f"Doses taken: {taken}/7")
```

---

## Statistical Analysis

### Power Analysis

```python
from openpkpd.trial import estimate_power_analytical, PowerResult

# Calculate power for given sample size and effect size
result = estimate_power_analytical(
    n_per_arm=50,
    effect_size=0.5,  # Cohen's d
    sd=1.0,
    alpha=0.05,
    alternative="two-sided"
)

print(f"Power: {result.power:.1%}")
```

### Sample Size Estimation

```python
from openpkpd.trial import estimate_sample_size

result = estimate_sample_size(
    target_power=0.80,
    effect_size=0.5,
    alpha=0.05
)

print(f"N per arm: {result.n_per_arm}")
print(f"Total N: {result.total_n}")
print(f"Achieved power: {result.achieved_power:.1%}")
```

### Effect Size and Sample Size Table

| Cohen's d | N per arm (80% power) | N per arm (90% power) |
|-----------|----------------------|----------------------|
| 0.2 (small) | ~394 | ~526 |
| 0.5 (medium) | ~64 | ~86 |
| 0.8 (large) | ~26 | ~34 |

### Alpha Spending Functions

```python
from openpkpd.trial import alpha_spending_function, incremental_alpha

# At 50% information
alpha_50 = alpha_spending_function(0.5, 0.05, "obrien_fleming")
print(f"Cumulative alpha at 50%: {alpha_50:.4f}")

# Incremental alpha for multiple analyses
fractions = [0.5, 0.75, 1.0]
alphas = incremental_alpha(fractions, 0.05, "obrien_fleming")
print(f"Alpha per analysis: {alphas}")
# Approximately: [0.0031, 0.0106, 0.0363]
```

### Arm Comparison

```python
from openpkpd.trial import compare_arms

active_values = [12.5, 14.2, 11.8, 15.1, 13.0, 12.8, 14.5]
placebo_values = [10.2, 9.8, 11.1, 10.5, 9.5, 10.8, 10.0]

result = compare_arms(active_values, placebo_values, "Active", "Placebo")

print(f"Difference: {result.difference:.2f}")
print(f"95% CI: [{result.ci_lower:.2f}, {result.ci_upper:.2f}]")
print(f"p-value: {result.p_value:.4f}")
print(f"Significant: {result.significant}")
```

### Responder Analysis

```python
from openpkpd.trial import responder_analysis

values = [10.2, 15.5, 8.1, 12.3, 9.5, 14.2, 11.8, 13.0, 7.5, 16.0]

result = responder_analysis(
    values,
    threshold=10.0,
    direction="greater",  # Response if value > threshold
    confidence=0.95
)

print(f"Response rate: {result.response_rate:.1%}")
print(f"Responders: {result.n_responders}/{result.n_total}")
print(f"95% CI: [{result.ci_lower:.1%}, {result.ci_upper:.1%}]")
```

### Bioequivalence Assessment

```python
from openpkpd.trial import bioequivalence_90ci, assess_bioequivalence

test = [95.2, 102.1, 98.5, 105.3, 97.8, 100.2]
ref = [100.0, 98.5, 101.2, 99.8, 100.5, 97.2]

# 90% CI
lower, upper = bioequivalence_90ci(test, ref)
print(f"90% CI: [{lower:.4f}, {upper:.4f}]")

# Full assessment
result = assess_bioequivalence(test, ref)
print(f"GMR: {result['point_estimate']:.4f}")
print(f"Bioequivalent: {result['bioequivalent']}")
```

---

## Complete Workflow Examples

### Phase I Dose Escalation

```python
from openpkpd import trial

# Design 3+3 escalation
dose_levels = [10, 25, 50, 100, 200]
design = trial.dose_escalation_3plus3(dose_levels)

print(f"Starting at: {design.starting_dose} mg")
print(f"Cohort size: {design.cohort_size}")
print(f"Rule: {design.escalation_rule}")

# Create regimens for each dose level
for dose in dose_levels:
    regimen = trial.dosing_qd(dose, 7)
    print(f"Dose {dose} mg: {trial.total_regimen_duration(regimen)} days")
```

### Phase II Parallel Study with Power Analysis

```python
from openpkpd import trial

# 1. Determine sample size
sample_size = trial.estimate_sample_size(
    target_power=0.80,
    effect_size=0.5
)
n_per_arm = sample_size.n_per_arm
print(f"Need {n_per_arm} per arm for 80% power")

# 2. Set up study
spec = trial.TrialSpec(
    name="Phase 2 Efficacy",
    design=trial.parallel_design(2),
    arms=[
        trial.TreatmentArm("Placebo", trial.dosing_qd(0.0, 56), n_per_arm, placebo=True),
        trial.TreatmentArm("Active", trial.dosing_qd(100.0, 56), n_per_arm),
    ],
    duration_days=56,
    dropout=trial.DropoutSpec(random_rate_per_day=0.002),
    seed=42
)

# 3. Simulate multiple trials
results = trial.simulate_trial_replicates(spec, n_replicates=100)

# 4. Analyze completion rates
completion_rates = [r.overall_completion_rate for r in results]
print(f"Mean completion: {sum(completion_rates)/len(completion_rates):.1%}")
```

### Bioequivalence Study Simulation

```python
from openpkpd import trial
import random

# Design
design = trial.bioequivalence_design(regulatory_guidance="fda")

# Generate virtual BE data (crossover)
random.seed(42)
n_subjects = 24

test_auc = [random.gauss(100, 15) for _ in range(n_subjects)]
ref_auc = [random.gauss(100, 15) for _ in range(n_subjects)]

test_cmax = [random.gauss(20, 4) for _ in range(n_subjects)]
ref_cmax = [random.gauss(20, 4) for _ in range(n_subjects)]

# Assess BE for both parameters
for param, test_vals, ref_vals in [
    ("AUC", test_auc, ref_auc),
    ("Cmax", test_cmax, ref_cmax)
]:
    result = trial.assess_bioequivalence(test_vals, ref_vals)
    print(f"\n{param}:")
    print(f"  GMR: {result['point_estimate']:.4f}")
    print(f"  90% CI: [{result['ci_90_lower']:.4f}, {result['ci_90_upper']:.4f}]")
    print(f"  BE: {'Yes' if result['bioequivalent'] else 'No'}")
```

---

## See Also

- [Models Reference](models.md) - PK/PD model documentation
- [NCA Reference](nca.md) - Non-compartmental analysis
- [Visualization](visualization.md) - Trial result plots
- [Population Simulation](population.md) - IIV and covariate modeling
