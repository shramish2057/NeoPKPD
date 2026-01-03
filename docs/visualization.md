# Visualization Reference

OpenPKPD provides comprehensive visualization capabilities for pharmacometric data using both Matplotlib (static) and Plotly (interactive) backends.

## Overview

The visualization module includes:

- **PK Plots**: Concentration-time, spaghetti, mean with ribbon
- **PKPD Plots**: Effect-concentration, hysteresis, dose-response
- **Population Plots**: VPC, parameter distributions, forest plots, boxplots
- **NCA Plots**: Lambda-z fit, AUC visualization, dose proportionality
- **Trial Plots**: Power curves, tornado diagrams, Kaplan-Meier

---

## Quick Start

```python
from openpkpd import viz

# Set backend
viz.set_backend("matplotlib")  # or "plotly"

# Basic concentration-time plot
result = openpkpd.simulate_pk_iv_bolus(cl=5.0, v=50.0, ...)
viz.plot_conc_time(result)
```

---

## Backend Selection

### Matplotlib (Static)

Best for publication-quality figures and reports.

```python
from openpkpd import viz

viz.set_backend("matplotlib")

# All plots return matplotlib Figure/Axes
fig = viz.plot_conc_time(result)
fig.savefig("pk_profile.png", dpi=300)
```

### Plotly (Interactive)

Best for data exploration and web applications.

```python
viz.set_backend("plotly")

# All plots return Plotly Figure objects
fig = viz.plot_conc_time(result)
fig.show()  # Opens in browser
fig.write_html("pk_profile.html")
```

### Backend Detection

```python
backend = viz.get_backend()
print(f"Current backend: {backend}")

if viz.is_plotly_available():
    viz.set_backend("plotly")
```

---

## Themes

### OpenPKPD Theme

Professional theme optimized for pharmacometric visualizations.

```python
from openpkpd.viz import apply_openpkpd_theme, OPENPKPD_COLORS

apply_openpkpd_theme()  # Apply to all subsequent plots

# Access color palette
print(OPENPKPD_COLORS)
# {'primary': '#1f77b4', 'secondary': '#ff7f0e', 'tertiary': '#2ca02c', ...}
```

### Color Palette

| Name | Hex | Use Case |
|------|-----|----------|
| primary | #1f77b4 | Main data series |
| secondary | #ff7f0e | Comparison/reference |
| tertiary | #2ca02c | Third series |
| error | #d62728 | Errors/warnings |
| highlight | #9467bd | Emphasis |
| muted | #7f7f7f | Background/secondary |

---

## PK Plots

### Concentration-Time Plot

```python
from openpkpd import viz

fig = viz.plot_conc_time(
    result,
    log_scale=False,           # Linear or log y-axis
    show_doses=True,           # Show dose markers
    title="PK Profile",
    xlabel="Time (h)",
    ylabel="Concentration (mg/L)"
)
```

**Options**:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `log_scale` | bool | False | Log-scale y-axis |
| `show_doses` | bool | True | Mark dose times |
| `time_unit` | str | "h" | Time unit for label |
| `conc_unit` | str | "mg/L" | Concentration unit |

### Spaghetti Plot (Population)

Individual PK profiles overlaid.

```python
fig = viz.plot_spaghetti(
    pop_result,
    n_subjects=None,     # None = all subjects
    alpha=0.3,           # Line transparency
    highlight_mean=True, # Show population mean
    color_by=None        # Optional: 'sex', 'weight_group'
)
```

### Mean with Confidence Ribbon

```python
fig = viz.plot_mean_ribbon(
    pop_result,
    ci_levels=[0.05, 0.95],    # 90% CI
    show_median=True,
    show_individual=False,      # Overlay individual profiles
    n_bootstrap=1000            # Bootstrap samples for CI
)
```

### Individual Fits

Grid of individual PK profiles for model validation.

```python
fig = viz.plot_individual_fits(
    pop_result,
    observed_data=obs_df,    # DataFrame with observed data
    n_cols=4,                # Number of columns in grid
    n_subjects=16,           # Subjects to show
    show_residuals=False     # Add residual subplot
)
```

---

## PKPD Plots

### Effect-Concentration Plot

```python
fig = viz.plot_effect_conc(
    pkpd_result,
    color_by_time=True,      # Color gradient by time
    show_regression=True,    # Fit Emax curve
    ec50_line=True          # Show EC50 reference
)
```

### Hysteresis Plot

Visualize PK-PD temporal delay.

```python
fig = viz.plot_hysteresis(
    pkpd_result,
    arrow_frequency=5,      # Arrows every 5 points
    show_time_labels=True,  # Time annotations
    absorption_color="blue",
    elimination_color="red"
)
```

**Interpretation**:
- Clockwise loop: Effect lags concentration (indirect response)
- Counter-clockwise loop: Effect precedes concentration (tolerance)

### Dose-Response Curve

```python
# From multiple dose simulations
dose_results = [simulate_at_dose(d) for d in [10, 25, 50, 100, 200]]
doses = [10, 25, 50, 100, 200]

fig = viz.plot_dose_response(
    dose_results,
    doses,
    endpoint="effect",
    fit_model="emax",       # "emax", "emax_hill", "linear", "log_linear"
    show_ci=True,
    ed50_line=True
)
```

---

## Population Plots

### Visual Predictive Check (VPC)

Standard validation plot for population models.

```python
fig = viz.plot_vpc(
    sim_results,                    # List of simulation results
    observed_data=obs_df,          # Observed data DataFrame
    prediction_intervals=[0.05, 0.50, 0.95],
    bin_method="equal_n",          # "equal_n", "equal_width", "optimal"
    n_bins=10,
    show_ci=True,
    log_scale=False
)
```

**VPC Elements**:

| Element | Description |
|---------|-------------|
| Shaded regions | Predicted PI (5th, 50th, 95th) |
| Lines | Observed PI |
| Points | Individual observations |

### Parameter Distributions

```python
# Histogram of population parameters
fig = viz.plot_parameter_distributions(
    pop_result,
    parameters=["CL", "V", "Ka"],
    plot_type="histogram",    # "histogram", "kde", "violin"
    n_cols=3,
    show_typical=True,        # Mark typical values
    log_scale=False
)
```

### Forest Plot

For meta-analysis or subgroup effects.

```python
effects = [
    {"group": "Overall", "estimate": 1.0, "ci_lower": 0.85, "ci_upper": 1.18},
    {"group": "Age <65", "estimate": 0.92, "ci_lower": 0.75, "ci_upper": 1.12},
    {"group": "Age >=65", "estimate": 1.15, "ci_lower": 0.90, "ci_upper": 1.47},
    {"group": "Male", "estimate": 0.98, "ci_lower": 0.80, "ci_upper": 1.20},
    {"group": "Female", "estimate": 1.05, "ci_lower": 0.85, "ci_upper": 1.30},
]

fig = viz.plot_forest(
    effects,
    reference_line=1.0,
    show_weights=True,
    diamond_for_overall=True
)
```

### Boxplot / Violin Plot

Compare distributions across groups.

```python
fig = viz.plot_boxplot(
    pop_result,
    groups=["treatment", "sex"],
    metric="cmax",              # "cmax", "auc", "cl"
    violin=False,               # True for violin plot
    show_points=True,           # Overlay individual points
    show_means=True
)
```

---

## NCA Plots

### Lambda-z Fit Visualization

```python
from openpkpd import viz

fig = viz.plot_lambda_z_fit(
    nca_result,
    times,
    conc,
    show_excluded=True,      # Show points not used in fit
    show_equation=True,      # Display regression equation
    show_r2=True            # Show R² value
)
```

### AUC Visualization

```python
fig = viz.plot_auc_visualization(
    times,
    conc,
    nca_result,
    show_extrapolation=True,   # Show AUC extrapolated portion
    fill_alpha=0.3,
    method="lin_log_mixed"
)
```

### Dose Proportionality

```python
# NCA results at multiple doses
nca_results = [run_nca(t, c, d) for d in doses]
doses = [10, 25, 50, 100]

fig = viz.plot_dose_proportionality(
    nca_results,
    doses,
    metric="auc",           # "auc", "cmax", "cl_f"
    power_model=True,       # Fit power model
    log_scale=True,
    show_reference=True     # Proportionality reference line
)
```

**Dose Proportionality Assessment**:

| Slope (log-log) | Interpretation |
|-----------------|----------------|
| ~1.0 | Dose proportional |
| <1.0 | Less than proportional |
| >1.0 | Greater than proportional |

---

## Trial Plots

### Power Curve

```python
# Power at different sample sizes
sample_sizes = [20, 30, 40, 50, 60, 80, 100]
powers = [estimate_power(n=n, effect_size=0.5).power for n in sample_sizes]

fig = viz.plot_power_curve(
    sample_sizes,
    powers,
    target_power=0.80,
    target_n=None,           # Highlight specific N
    show_80_line=True,
    show_90_line=True
)
```

### Tornado Diagram (Sensitivity)

```python
sensitivity_results = [
    {"parameter": "CL", "low": 0.85, "high": 1.20},
    {"parameter": "V", "low": 0.90, "high": 1.15},
    {"parameter": "Ka", "low": 0.95, "high": 1.08},
    {"parameter": "F", "low": 0.75, "high": 1.30},
]

fig = viz.plot_tornado(
    sensitivity_results,
    baseline_value=1.0,
    sort_by_range=True,
    show_values=True,
    symmetric_axis=True
)
```

### Kaplan-Meier Survival Curve

```python
# Time-to-event data
time_to_event = [10, 15, 22, 30, 45, 60, 90, 120, 150, 180]
event_occurred = [1, 1, 1, 0, 1, 0, 1, 0, 0, 0]  # 1=event, 0=censored
groups = ["A", "A", "A", "A", "A", "B", "B", "B", "B", "B"]

fig = viz.plot_kaplan_meier(
    time_to_event,
    event_occurred,
    groups=groups,
    ci=0.95,
    show_censored=True,
    show_median=True,
    show_at_risk=True        # Number at risk table
)
```

### Endpoint Distribution by Arm

```python
fig = viz.plot_endpoint_distribution(
    trial_result,
    endpoint="pk_exposure",
    by_arm=True,
    plot_type="violin",      # "histogram", "kde", "violin", "box"
    show_stats=True
)
```

---

## Plot Customization

### Figure Size and DPI

```python
# Matplotlib
fig = viz.plot_conc_time(result, figsize=(10, 6), dpi=150)

# Plotly
fig = viz.plot_conc_time(result, width=800, height=600)
```

### Adding Annotations

```python
import matplotlib.pyplot as plt

fig = viz.plot_conc_time(result)
ax = fig.axes[0]

# Add annotations
ax.axhline(y=1.0, color='red', linestyle='--', label='Target')
ax.annotate('Cmax', xy=(2, 1.62), xytext=(3, 1.8),
            arrowprops=dict(arrowstyle='->'))
ax.legend()

plt.tight_layout()
```

### Multi-Panel Figures

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot to each axis
viz.plot_conc_time(result1, ax=axes[0, 0], title="Subject 1")
viz.plot_conc_time(result2, ax=axes[0, 1], title="Subject 2")
viz.plot_effect_conc(result1, ax=axes[1, 0], title="PKPD")
viz.plot_hysteresis(result1, ax=axes[1, 1], title="Hysteresis")

plt.tight_layout()
plt.savefig("multi_panel.png", dpi=300)
```

---

## Saving Plots

### Matplotlib

```python
fig = viz.plot_conc_time(result)

# PNG (raster)
fig.savefig("plot.png", dpi=300, bbox_inches='tight')

# SVG (vector)
fig.savefig("plot.svg", format='svg', bbox_inches='tight')

# PDF (publication)
fig.savefig("plot.pdf", format='pdf', bbox_inches='tight')
```

### Plotly

```python
fig = viz.plot_conc_time(result)

# Interactive HTML
fig.write_html("plot.html")

# Static image (requires kaleido)
fig.write_image("plot.png", scale=2)
fig.write_image("plot.svg")
fig.write_image("plot.pdf")

# JSON (for web apps)
fig.write_json("plot.json")
```

---

## Common Workflows

### Publication-Ready PK Figure

```python
from openpkpd import viz
import matplotlib.pyplot as plt

viz.set_backend("matplotlib")
viz.apply_openpkpd_theme()

# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Linear scale
viz.plot_conc_time(result, ax=ax1, log_scale=False)
ax1.set_title("A) Linear Scale")

# Log scale
viz.plot_conc_time(result, ax=ax2, log_scale=True)
ax2.set_title("B) Semi-log Scale")

plt.tight_layout()
fig.savefig("pk_figure.pdf", format='pdf', dpi=300)
```

### Interactive Dashboard Plot

```python
viz.set_backend("plotly")

fig = viz.plot_spaghetti(pop_result, alpha=0.5)

# Add dropdown for log scale
fig.update_layout(
    updatemenus=[
        dict(
            type="buttons",
            buttons=[
                dict(label="Linear", method="relayout",
                     args=[{"yaxis.type": "linear"}]),
                dict(label="Log", method="relayout",
                     args=[{"yaxis.type": "log"}]),
            ],
        )
    ]
)

fig.write_html("interactive_pk.html")
```

### Trial Report Figures

```python
from openpkpd import viz, trial
import matplotlib.pyplot as plt

# Run trial simulation
spec = trial.TrialSpec(...)
result = trial.simulate_trial(spec)

# Create report figure
fig = plt.figure(figsize=(14, 10))

# Sample size vs power
ax1 = fig.add_subplot(2, 2, 1)
viz.plot_power_curve(sample_sizes, powers, ax=ax1)
ax1.set_title("Power Analysis")

# Completion by arm
ax2 = fig.add_subplot(2, 2, 2)
arms = list(result.arms.keys())
completion = [result.arms[a].completion_rate for a in arms]
ax2.bar(arms, completion)
ax2.set_ylabel("Completion Rate")
ax2.set_title("Study Completion")

# Endpoint distribution
ax3 = fig.add_subplot(2, 2, 3)
viz.plot_endpoint_distribution(result, ax=ax3)
ax3.set_title("Endpoint Distribution")

# Forest plot
ax4 = fig.add_subplot(2, 2, 4)
viz.plot_forest(effects, ax=ax4)
ax4.set_title("Subgroup Analysis")

plt.tight_layout()
fig.savefig("trial_report.pdf", dpi=300)
```

---

## Troubleshooting

### Common Issues

| Issue | Solution |
|-------|----------|
| Plot not showing | Call `plt.show()` or `fig.show()` |
| Fonts not embedding in PDF | Use `fig.savefig(..., fonttype=42)` |
| Plotly images not exporting | Install `kaleido`: `pip install kaleido` |
| Legend overlapping | Use `loc='best'` or `bbox_to_anchor` |
| Axis labels cut off | Use `bbox_inches='tight'` when saving |

### Performance Tips

```python
# For large datasets, use rasterization
ax.set_rasterization_zorder(0)

# Reduce Plotly file size
fig.write_html("plot.html", include_plotlyjs='cdn')

# Sample data for quick previews
viz.plot_spaghetti(pop_result, n_subjects=50)
```

---

## See Also

- [Models Reference](models.md) - PK/PD models
- [NCA Reference](nca.md) - NCA calculations
- [Trial Reference](trial.md) - Trial simulation
- [Population](population.md) - Population variability
