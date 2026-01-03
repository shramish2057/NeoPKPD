"""
OpenPKPD Trial Designs

Study design definitions for clinical trials including:
- Parallel group designs
- Crossover designs (2x2, 3x3, Williams)
- Dose escalation (3+3, mTPI, CRM)
- Adaptive designs
- Bioequivalence designs
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, Union
from enum import Enum

class EscalationRuleType(str, Enum):
    """Dose escalation rule types."""
    THREE_PLUS_THREE = "3+3"
    MTPI = "mTPI"
    CRM = "CRM"


@dataclass
class ParallelDesign:
    """
    Parallel group study design.

    Attributes:
        n_arms: Number of treatment arms
        randomization_ratio: Randomization ratio for each arm
        stratification_factors: Factors for stratified randomization

    Example:
        >>> design = ParallelDesign(n_arms=2)
        >>> design = ParallelDesign(n_arms=3, randomization_ratio=[0.5, 0.25, 0.25])
    """
    n_arms: int
    randomization_ratio: List[float] = field(default_factory=list)
    stratification_factors: List[str] = field(default_factory=list)

    def __post_init__(self):
        if not self.randomization_ratio:
            self.randomization_ratio = [1.0 / self.n_arms] * self.n_arms


@dataclass
class CrossoverDesign:
    """
    Crossover study design.

    Attributes:
        n_periods: Number of treatment periods
        n_sequences: Number of treatment sequences
        washout_duration: Washout period duration in days
        sequence_assignments: Treatment assignments per sequence

    Example:
        >>> design = CrossoverDesign(n_periods=2, n_sequences=2)  # 2x2 crossover
    """
    n_periods: int
    n_sequences: int
    washout_duration: float = 7.0
    sequence_assignments: List[List[int]] = field(default_factory=list)


@dataclass
class DoseEscalationDesign:
    """
    Dose escalation study design for Phase I trials.

    Attributes:
        dose_levels: Available dose levels
        starting_dose: Starting dose level
        escalation_rule: Rule type ('3+3', 'mTPI', 'CRM')
        cohort_size: Number of subjects per cohort
        max_subjects: Maximum total subjects
        target_dlt_rate: Target DLT rate (for mTPI/CRM)

    Example:
        >>> design = DoseEscalationDesign(
        ...     dose_levels=[10.0, 25.0, 50.0, 100.0, 200.0],
        ...     escalation_rule='3+3'
        ... )
    """
    dose_levels: List[float]
    starting_dose: Optional[float] = None
    escalation_rule: str = "3+3"
    cohort_size: int = 3
    max_subjects: int = 30
    target_dlt_rate: float = 0.25

    def __post_init__(self):
        if self.starting_dose is None:
            self.starting_dose = self.dose_levels[0]


@dataclass
class AdaptiveDesign:
    """
    Adaptive trial design with interim analyses.

    Attributes:
        base_design: Underlying study design
        interim_analyses: Information fractions for interim analyses
        alpha_spending: Alpha spending function type
        futility_boundary: Futility boundary threshold

    Example:
        >>> base = ParallelDesign(n_arms=2)
        >>> design = AdaptiveDesign(base_design=base, interim_analyses=[0.5])
    """
    base_design: Union[ParallelDesign, CrossoverDesign]
    interim_analyses: List[float] = field(default_factory=lambda: [0.5])
    alpha_spending: str = "obrien_fleming"
    futility_boundary: float = 0.10


@dataclass
class BioequivalenceDesign:
    """
    Bioequivalence study design.

    Attributes:
        n_periods: Number of periods (default: 2)
        n_sequences: Number of sequences (default: 2)
        washout_duration: Washout period in days
        bioequivalence_limits: BE acceptance limits
        parameters: PK parameters for BE assessment
        regulatory_guidance: Regulatory guidance ('fda' or 'ema')

    Example:
        >>> design = BioequivalenceDesign(regulatory_guidance='fda')
    """
    n_periods: int = 2
    n_sequences: int = 2
    washout_duration: float = 7.0
    bioequivalence_limits: Tuple[float, float] = (0.80, 1.25)
    parameters: List[str] = field(default_factory=lambda: ["cmax", "auc_0_inf"])
    regulatory_guidance: str = "fda"


def parallel_design(
    n_arms: int,
    randomization_ratio: Optional[List[float]] = None,
    stratification_factors: Optional[List[str]] = None,
) -> ParallelDesign:
    """
    Create a parallel group study design.

    Args:
        n_arms: Number of treatment arms
        randomization_ratio: Randomization ratio (default: equal)
        stratification_factors: Stratification factors

    Returns:
        ParallelDesign object

    Example:
        >>> design = parallel_design(2)
        >>> design = parallel_design(3, randomization_ratio=[0.5, 0.25, 0.25])
    """
    return ParallelDesign(
        n_arms=n_arms,
        randomization_ratio=randomization_ratio or [],
        stratification_factors=stratification_factors or [],
    )


def crossover_2x2(washout_duration: float = 7.0) -> CrossoverDesign:
    """
    Create a 2-period, 2-sequence crossover design (AB, BA).

    Args:
        washout_duration: Washout period in days

    Returns:
        CrossoverDesign object

    Example:
        >>> design = crossover_2x2(washout_duration=14.0)
    """
    return CrossoverDesign(
        n_periods=2,
        n_sequences=2,
        washout_duration=washout_duration,
        sequence_assignments=[[1, 2], [2, 1]],
    )


def crossover_3x3(washout_duration: float = 7.0) -> CrossoverDesign:
    """
    Create a 3-period, 3-sequence crossover design (Latin square).

    Args:
        washout_duration: Washout period in days

    Returns:
        CrossoverDesign object

    Example:
        >>> design = crossover_3x3(washout_duration=21.0)
    """
    return CrossoverDesign(
        n_periods=3,
        n_sequences=3,
        washout_duration=washout_duration,
        sequence_assignments=[[1, 2, 3], [2, 3, 1], [3, 1, 2]],
    )


def williams_design(n_treatments: int, washout_duration: float = 7.0) -> CrossoverDesign:
    """
    Create a Williams design for crossover studies.

    Williams designs are balanced for first-order carryover effects.

    Args:
        n_treatments: Number of treatments (2, 3, or 4)
        washout_duration: Washout period in days

    Returns:
        CrossoverDesign object

    Example:
        >>> design = williams_design(4, washout_duration=14.0)
    """
    if n_treatments == 2:
        sequences = [[1, 2], [2, 1]]
    elif n_treatments == 3:
        sequences = [[1, 2, 3], [2, 3, 1], [3, 1, 2],
                     [1, 3, 2], [3, 2, 1], [2, 1, 3]]
    elif n_treatments == 4:
        sequences = [[1, 2, 4, 3], [2, 3, 1, 4], [3, 4, 2, 1], [4, 1, 3, 2],
                     [1, 4, 2, 3], [2, 1, 3, 4], [3, 2, 4, 1], [4, 3, 1, 2]]
    else:
        raise ValueError("Williams design only supports 2, 3, or 4 treatments")

    return CrossoverDesign(
        n_periods=n_treatments,
        n_sequences=len(sequences),
        washout_duration=washout_duration,
        sequence_assignments=sequences,
    )


def dose_escalation_3plus3(
    dose_levels: List[float],
    starting_dose: Optional[float] = None,
    max_dlt_rate: float = 0.33,
    cohort_size: int = 3,
    max_subjects: int = 30,
) -> DoseEscalationDesign:
    """
    Create a 3+3 dose escalation design.

    Args:
        dose_levels: Available dose levels
        starting_dose: Starting dose (default: first level)
        max_dlt_rate: Maximum acceptable DLT rate
        cohort_size: Cohort size (default: 3)
        max_subjects: Maximum total subjects

    Returns:
        DoseEscalationDesign object

    Example:
        >>> design = dose_escalation_3plus3([10.0, 25.0, 50.0, 100.0, 200.0])
    """
    return DoseEscalationDesign(
        dose_levels=dose_levels,
        starting_dose=starting_dose,
        escalation_rule="3+3",
        cohort_size=cohort_size,
        max_subjects=max_subjects,
        target_dlt_rate=max_dlt_rate,
    )


def dose_escalation_mtpi(
    dose_levels: List[float],
    starting_dose: Optional[float] = None,
    target_dlt_rate: float = 0.25,
    cohort_size: int = 3,
    max_subjects: int = 30,
) -> DoseEscalationDesign:
    """
    Create an mTPI dose escalation design.

    Args:
        dose_levels: Available dose levels
        starting_dose: Starting dose (default: first level)
        target_dlt_rate: Target DLT rate
        cohort_size: Cohort size
        max_subjects: Maximum total subjects

    Returns:
        DoseEscalationDesign object

    Example:
        >>> design = dose_escalation_mtpi([10.0, 25.0, 50.0, 100.0], target_dlt_rate=0.30)
    """
    return DoseEscalationDesign(
        dose_levels=dose_levels,
        starting_dose=starting_dose,
        escalation_rule="mTPI",
        cohort_size=cohort_size,
        max_subjects=max_subjects,
        target_dlt_rate=target_dlt_rate,
    )


def dose_escalation_crm(
    dose_levels: List[float],
    starting_dose: Optional[float] = None,
    target_dlt_rate: float = 0.25,
    cohort_size: int = 1,
    max_subjects: int = 30,
) -> DoseEscalationDesign:
    """
    Create a CRM (Continual Reassessment Method) dose escalation design.

    Args:
        dose_levels: Available dose levels
        starting_dose: Starting dose (default: first level)
        target_dlt_rate: Target DLT rate
        cohort_size: Cohort size (usually 1 for CRM)
        max_subjects: Maximum total subjects

    Returns:
        DoseEscalationDesign object

    Example:
        >>> design = dose_escalation_crm([10.0, 25.0, 50.0, 100.0, 200.0])
    """
    return DoseEscalationDesign(
        dose_levels=dose_levels,
        starting_dose=starting_dose,
        escalation_rule="CRM",
        cohort_size=cohort_size,
        max_subjects=max_subjects,
        target_dlt_rate=target_dlt_rate,
    )


def adaptive_design(
    base_design: Union[ParallelDesign, CrossoverDesign],
    interim_analyses: List[float] = [0.5],
    alpha_spending: str = "obrien_fleming",
    futility_boundary: float = 0.10,
) -> AdaptiveDesign:
    """
    Create an adaptive trial design.

    Args:
        base_design: Underlying study design
        interim_analyses: Information fractions for interim analyses
        alpha_spending: Alpha spending function ('obrien_fleming', 'pocock', 'haybittle_peto')
        futility_boundary: Futility boundary (conditional power threshold)

    Returns:
        AdaptiveDesign object

    Example:
        >>> base = parallel_design(2)
        >>> design = adaptive_design(base, interim_analyses=[0.5])
    """
    return AdaptiveDesign(
        base_design=base_design,
        interim_analyses=interim_analyses,
        alpha_spending=alpha_spending,
        futility_boundary=futility_boundary,
    )


def bioequivalence_design(
    n_periods: int = 2,
    n_sequences: int = 2,
    washout_duration: float = 7.0,
    bioequivalence_limits: Tuple[float, float] = (0.80, 1.25),
    parameters: List[str] = ["cmax", "auc_0_inf"],
    regulatory_guidance: str = "fda",
) -> BioequivalenceDesign:
    """
    Create a bioequivalence study design.

    Args:
        n_periods: Number of periods
        n_sequences: Number of sequences
        washout_duration: Washout period in days
        bioequivalence_limits: BE acceptance limits
        parameters: PK parameters to assess
        regulatory_guidance: Regulatory guidance ('fda' or 'ema')

    Returns:
        BioequivalenceDesign object

    Example:
        >>> design = bioequivalence_design(regulatory_guidance='fda')
        >>> # For highly variable drugs (EMA)
        >>> design = bioequivalence_design(
        ...     bioequivalence_limits=(0.6984, 1.4319),
        ...     regulatory_guidance='ema'
        ... )
    """
    return BioequivalenceDesign(
        n_periods=n_periods,
        n_sequences=n_sequences,
        washout_duration=washout_duration,
        bioequivalence_limits=bioequivalence_limits,
        parameters=parameters,
        regulatory_guidance=regulatory_guidance,
    )


def get_design_description(design: Any) -> str:
    """
    Get a human-readable description of a study design.

    Args:
        design: Study design object

    Returns:
        Design description string

    Example:
        >>> design = parallel_design(2)
        >>> print(get_design_description(design))
        '2-arm parallel group design'
    """
    if isinstance(design, ParallelDesign):
        return f"{design.n_arms}-arm parallel group design"
    elif isinstance(design, CrossoverDesign):
        return f"{design.n_periods}-period {design.n_sequences}-sequence crossover design"
    elif isinstance(design, DoseEscalationDesign):
        return f"Dose escalation ({design.escalation_rule}) with {len(design.dose_levels)} dose levels"
    elif isinstance(design, AdaptiveDesign):
        base_desc = get_design_description(design.base_design)
        n_interim = len(design.interim_analyses)
        return f"Adaptive {base_desc} with {n_interim} interim analysis(es)"
    elif isinstance(design, BioequivalenceDesign):
        return f"Bioequivalence {design.n_periods}x{design.n_sequences} crossover ({design.regulatory_guidance})"
    else:
        return str(type(design).__name__)
