"""
OpenPKPD Population Visualization

Population-level visualization functions including:
- Visual Predictive Check (VPC)
- Parameter distributions
- Forest plots
- Boxplots and violin plots
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union

from .backends import _get_plotter
from .themes import get_theme_config, get_color, get_colors


def plot_vpc(
    population_result: Dict[str, Any],
    observed_data: Optional[Dict[str, List[float]]] = None,
    observation: str = "conc",
    prediction_intervals: List[float] = [0.05, 0.50, 0.95],
    log_scale: bool = False,
    title: Optional[str] = "Visual Predictive Check",
    xlabel: str = "Time",
    ylabel: str = "Concentration",
    figsize: Tuple[float, float] = (12, 8),
    show_observed: bool = True,
    n_bins: int = 10,
) -> Any:
    """
    Create Visual Predictive Check (VPC) plot.

    VPC compares simulated prediction intervals with observed data to
    assess model adequacy.

    Args:
        population_result: Population simulation result dictionary
        observed_data: Optional dict with 't' and 'c' lists for observed data
        observation: Observation key to plot
        prediction_intervals: Quantiles for prediction intervals
        log_scale: Use log scale for y-axis
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size
        show_observed: Show observed data points
        n_bins: Number of time bins for prediction intervals

    Returns:
        Figure object

    Example:
        >>> pop_result = openpkpd.simulate_population_iv_bolus(...)
        >>> observed = {"t": [0, 1, 2, 4], "c": [10, 8, 6, 3]}
        >>> fig = plot_vpc(pop_result, observed_data=observed)
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    import numpy as np

    individuals = population_result["individuals"]
    t = np.array(individuals[0]["t"])

    # Collect all concentration data
    all_conc = []
    for ind in individuals:
        all_conc.append(ind["observations"][observation])
    all_conc = np.array(all_conc)

    # Calculate prediction intervals
    pi_lower = np.percentile(all_conc, prediction_intervals[0] * 100, axis=0)
    pi_median = np.percentile(all_conc, prediction_intervals[1] * 100, axis=0)
    pi_upper = np.percentile(all_conc, prediction_intervals[2] * 100, axis=0)

    fig = plotter.create_figure(figsize=figsize, title=title)

    # Plot outer prediction interval
    pi_label = f"{int((prediction_intervals[2] - prediction_intervals[0]) * 100)}% PI"
    plotter.fill_between(fig, list(t), list(pi_lower), list(pi_upper),
                         color=get_color(0), alpha=theme["alpha_ribbon"],
                         label=pi_label)

    # Plot median line
    plotter.line_plot(fig, list(t), list(pi_median), label="Median",
                      color=get_color(0), linewidth=theme["line_width"] * 1.5)

    # Plot observed data if provided
    if show_observed and observed_data is not None:
        t_obs = observed_data.get("t", [])
        c_obs = observed_data.get("c", observed_data.get("conc", []))
        if t_obs and c_obs:
            plotter.scatter_plot(fig, t_obs, c_obs, color=get_color(1),
                                 size=theme["marker_size"] * 8,
                                 label="Observed")

    plotter.set_labels(fig, xlabel=xlabel, ylabel=ylabel, title=title)

    if log_scale:
        plotter.set_log_scale(fig, y=True)

    plotter.add_legend(fig)

    return plotter.finalize(fig)


def plot_parameter_distributions(
    population_result: Dict[str, Any],
    parameters: Optional[List[str]] = None,
    plot_type: str = "histogram",
    n_cols: int = 3,
    subplot_size: Tuple[float, float] = (4, 3),
    show_stats: bool = True,
    log_scale: bool = False,
) -> Any:
    """
    Plot distributions of individual PK parameters.

    Args:
        population_result: Population simulation result dictionary
        parameters: List of parameter names to plot (None = all)
        plot_type: "histogram" or "kde"
        n_cols: Number of columns in subplot grid
        subplot_size: Size of each subplot
        show_stats: Show mean/median statistics
        log_scale: Use log scale for x-axis

    Returns:
        Figure object

    Example:
        >>> pop_result = openpkpd.simulate_population_iv_bolus(...)
        >>> fig = plot_parameter_distributions(pop_result, parameters=["CL", "V"])
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    import numpy as np

    individuals = population_result["individuals"]

    # Extract parameters
    if parameters is None:
        # Try to infer from first individual
        if "parameters" in individuals[0]:
            parameters = list(individuals[0]["parameters"].keys())
        else:
            parameters = []

    if not parameters:
        # Fallback: use summary statistics if available
        if "summaries" in population_result:
            parameters = list(population_result["summaries"].keys())

    if not parameters:
        raise ValueError("No parameters found in population result")

    # Collect parameter values
    param_values = {p: [] for p in parameters}
    for ind in individuals:
        if "parameters" in ind:
            for p in parameters:
                if p in ind["parameters"]:
                    param_values[p].append(ind["parameters"][p])

    # Filter out empty parameters
    parameters = [p for p in parameters if param_values[p]]

    n = len(parameters)
    n_rows = (n + n_cols - 1) // n_cols
    figsize = (subplot_size[0] * n_cols, subplot_size[1] * n_rows)

    backend = plotter.__class__.__name__

    if "Matplotlib" in backend:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
        if n > 1:
            axes = axes.flatten()
        else:
            axes = [axes]

        for i, param in enumerate(parameters):
            ax = axes[i]
            values = np.array(param_values[param])

            if plot_type == "histogram":
                ax.hist(values, bins=30, color=get_color(0), alpha=0.7,
                        edgecolor="white")
            else:  # kde
                from scipy import stats
                try:
                    kde = stats.gaussian_kde(values)
                    x_range = np.linspace(values.min(), values.max(), 100)
                    ax.fill_between(x_range, kde(x_range), alpha=0.7,
                                    color=get_color(0))
                except Exception:
                    ax.hist(values, bins=30, color=get_color(0), alpha=0.7)

            if show_stats:
                mean_val = np.mean(values)
                median_val = np.median(values)
                ax.axvline(mean_val, color=get_color(1), linestyle="--",
                           linewidth=1.5, label=f"Mean: {mean_val:.3g}")
                ax.axvline(median_val, color=get_color(2), linestyle=":",
                           linewidth=1.5, label=f"Median: {median_val:.3g}")
                ax.legend(fontsize=8)

            ax.set_xlabel(param)
            ax.set_ylabel("Frequency" if plot_type == "histogram" else "Density")
            ax.set_title(param)

            if log_scale:
                ax.set_xscale("log")

        # Hide empty subplots
        for i in range(n, len(axes)):
            axes[i].set_visible(False)

        fig.tight_layout()
        return fig

    else:  # Plotly
        from plotly.subplots import make_subplots
        import plotly.graph_objects as go

        fig = make_subplots(rows=n_rows, cols=n_cols,
                            subplot_titles=parameters)

        for i, param in enumerate(parameters):
            row = i // n_cols + 1
            col = i % n_cols + 1
            values = param_values[param]

            fig.add_trace(
                go.Histogram(x=values, name=param, marker_color=get_color(0),
                             showlegend=False),
                row=row, col=col
            )

            if show_stats:
                mean_val = np.mean(values)
                fig.add_vline(x=mean_val, line_dash="dash",
                              line_color=get_color(1), row=row, col=col)

        fig.update_layout(
            height=figsize[1] * 100,
            width=figsize[0] * 100,
            title="Parameter Distributions"
        )

        return fig


def plot_forest(
    effects: List[Dict[str, Any]],
    reference_line: float = 1.0,
    log_scale: bool = True,
    title: Optional[str] = "Forest Plot",
    xlabel: str = "Effect Size",
    figsize: Tuple[float, float] = (10, 8),
    show_summary: bool = True,
) -> Any:
    """
    Create forest plot for effect estimates with confidence intervals.

    Args:
        effects: List of dicts with keys 'name', 'estimate', 'lower', 'upper'
        reference_line: Reference line value (usually 1.0 for ratios)
        log_scale: Use log scale for x-axis
        title: Plot title
        xlabel: X-axis label
        figsize: Figure size
        show_summary: Show overall summary diamond

    Returns:
        Figure object

    Example:
        >>> effects = [
        ...     {"name": "Study 1", "estimate": 1.2, "lower": 0.9, "upper": 1.6},
        ...     {"name": "Study 2", "estimate": 0.95, "lower": 0.7, "upper": 1.3},
        ... ]
        >>> fig = plot_forest(effects)
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    import numpy as np

    n = len(effects)
    y_positions = list(range(n, 0, -1))

    backend = plotter.__class__.__name__

    if "Matplotlib" in backend:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize)

        for i, (effect, y) in enumerate(zip(effects, y_positions)):
            estimate = effect["estimate"]
            lower = effect["lower"]
            upper = effect["upper"]
            name = effect.get("name", f"Study {i+1}")

            # Draw confidence interval line
            ax.plot([lower, upper], [y, y], color=get_color(0),
                    linewidth=theme["line_width"], solid_capstyle="butt")

            # Draw point estimate
            ax.scatter([estimate], [y], color=get_color(0),
                       s=theme["marker_size"] * 15, zorder=3)

            # Add study name
            ax.text(ax.get_xlim()[0] if ax.get_xlim()[0] > 0 else 0.1, y,
                    f"  {name}", va="center", ha="left", fontsize=10)

        # Reference line
        ax.axvline(reference_line, color=get_color(2), linestyle="--",
                   linewidth=1, alpha=0.7, label=f"Reference ({reference_line})")

        # Summary diamond
        if show_summary and len(effects) > 1:
            all_estimates = [e["estimate"] for e in effects]
            all_weights = [1.0 / ((e["upper"] - e["lower"]) ** 2) for e in effects]
            total_weight = sum(all_weights)
            summary_estimate = sum(e * w for e, w in zip(all_estimates, all_weights)) / total_weight

            # Simplified summary CI (weighted average of CIs)
            summary_lower = sum(e["lower"] * w for e, w in zip(effects, all_weights)) / total_weight
            summary_upper = sum(e["upper"] * w for e, w in zip(effects, all_weights)) / total_weight

            diamond_y = 0
            diamond_height = 0.4
            diamond = plt.Polygon([
                [summary_lower, diamond_y],
                [summary_estimate, diamond_y + diamond_height],
                [summary_upper, diamond_y],
                [summary_estimate, diamond_y - diamond_height]
            ], color=get_color(1), alpha=0.8)
            ax.add_patch(diamond)
            ax.text(ax.get_xlim()[0] if ax.get_xlim()[0] > 0 else 0.1, diamond_y,
                    "  Summary", va="center", ha="left", fontsize=10, fontweight="bold")

        ax.set_yticks(y_positions + ([0] if show_summary else []))
        ax.set_yticklabels([e.get("name", f"Study {i+1}") for i, e in enumerate(effects)] +
                           (["Summary"] if show_summary else []))
        ax.set_xlabel(xlabel)
        ax.set_title(title)

        if log_scale:
            ax.set_xscale("log")

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()
        return fig

    else:  # Plotly
        import plotly.graph_objects as go

        fig = go.Figure()

        names = [e.get("name", f"Study {i+1}") for i, e in enumerate(effects)]
        estimates = [e["estimate"] for e in effects]
        lowers = [e["lower"] for e in effects]
        uppers = [e["upper"] for e in effects]

        # Error bars
        fig.add_trace(go.Scatter(
            x=estimates, y=y_positions,
            mode="markers",
            marker=dict(size=10, color=get_color(0)),
            error_x=dict(
                type="data",
                symmetric=False,
                array=[u - e for u, e in zip(uppers, estimates)],
                arrayminus=[e - l for e, l in zip(estimates, lowers)]
            ),
            name="Studies"
        ))

        # Reference line
        fig.add_vline(x=reference_line, line_dash="dash",
                      line_color=get_color(2), annotation_text="Reference")

        fig.update_layout(
            title=title,
            xaxis_title=xlabel,
            yaxis=dict(
                tickvals=y_positions,
                ticktext=names
            ),
            height=figsize[1] * 100,
            width=figsize[0] * 100
        )

        if log_scale:
            fig.update_xaxes(type="log")

        return fig


def plot_boxplot(
    population_result: Dict[str, Any],
    groups: Optional[List[str]] = None,
    metric: str = "cmax",
    violin: bool = False,
    show_points: bool = True,
    title: Optional[str] = None,
    xlabel: str = "Group",
    ylabel: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 6),
) -> Any:
    """
    Create boxplot or violin plot comparing groups.

    Args:
        population_result: Population simulation result dictionary
        groups: Group labels (one per individual, or derived from result)
        metric: Metric to plot (e.g., "cmax", "auc")
        violin: Use violin plot instead of boxplot
        show_points: Overlay individual data points
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figsize: Figure size

    Returns:
        Figure object

    Example:
        >>> pop_result = openpkpd.simulate_population_iv_bolus(...)
        >>> groups = ["Low Dose"] * 50 + ["High Dose"] * 50
        >>> fig = plot_boxplot(pop_result, groups=groups, metric="cmax")
        >>> fig.show()
    """
    plotter = _get_plotter()
    theme = get_theme_config()

    import numpy as np

    individuals = population_result["individuals"]

    # Extract metric values
    values = []
    for ind in individuals:
        obs = ind["observations"]
        if metric == "cmax":
            values.append(max(obs.get("conc", obs.get(list(obs.keys())[0]))))
        elif metric == "tmax":
            t = ind["t"]
            c = obs.get("conc", obs.get(list(obs.keys())[0]))
            values.append(t[np.argmax(c)])
        elif metric == "auc":
            t = np.array(ind["t"])
            c = np.array(obs.get("conc", obs.get(list(obs.keys())[0])))
            values.append(np.trapz(c, t))
        elif metric in obs:
            values.append(max(obs[metric]))
        else:
            values.append(0.0)

    values = np.array(values)

    # Handle groups
    if groups is None:
        groups = ["All"] * len(values)

    unique_groups = list(dict.fromkeys(groups))  # Preserve order
    n_groups = len(unique_groups)

    if ylabel is None:
        ylabel = metric.upper()
    if title is None:
        title = f"{ylabel} by Group"

    backend = plotter.__class__.__name__

    if "Matplotlib" in backend:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize)

        # Prepare data by group
        group_data = [values[[i for i, g in enumerate(groups) if g == ug]]
                      for ug in unique_groups]

        positions = list(range(1, n_groups + 1))

        if violin:
            parts = ax.violinplot(group_data, positions=positions,
                                  showmeans=True, showmedians=True)
            for pc in parts["bodies"]:
                pc.set_facecolor(get_color(0))
                pc.set_alpha(0.7)
        else:
            bp = ax.boxplot(group_data, positions=positions, patch_artist=True)
            for patch in bp["boxes"]:
                patch.set_facecolor(get_color(0))
                patch.set_alpha(0.7)

        if show_points:
            for i, (ug, data) in enumerate(zip(unique_groups, group_data)):
                jitter = np.random.uniform(-0.15, 0.15, len(data))
                ax.scatter([i + 1 + j for j in jitter], data,
                           color=get_color(1), alpha=0.5, s=20, zorder=3)

        ax.set_xticks(positions)
        ax.set_xticklabels(unique_groups)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

        fig.tight_layout()
        return fig

    else:  # Plotly
        import plotly.graph_objects as go

        fig = go.Figure()

        for i, ug in enumerate(unique_groups):
            group_values = values[[j for j, g in enumerate(groups) if g == ug]]

            if violin:
                fig.add_trace(go.Violin(
                    y=group_values, name=ug,
                    box_visible=True, meanline_visible=True,
                    fillcolor=get_color(i), line_color=get_color(i)
                ))
            else:
                fig.add_trace(go.Box(
                    y=group_values, name=ug,
                    marker_color=get_color(i),
                    boxpoints="all" if show_points else False
                ))

        fig.update_layout(
            title=title,
            xaxis_title=xlabel,
            yaxis_title=ylabel,
            height=figsize[1] * 100,
            width=figsize[0] * 100
        )

        return fig
