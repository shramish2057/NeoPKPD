# PKPD biomarker turnover

This use case models exposure driven biomarker suppression using a coupled PKPD indirect response turnover model.

Purpose:
- Validate coupled PKPD, regimen comparison, and IOV behavior under realistic workflows.
- Provide a reproducible, replayable contract for response metrics used in decision making.

PK model:
- One-compartment IV bolus

PD model:
- Indirect response turnover
- Baseline R0 = Kin / Kout
- Drug effect driven by concentration with an inhibitory Emax form

Scenarios:
- QD dosing: 100 mg at 0, 24, 48
- BID dosing: 50 mg at 0, 12, 24, 36, 48, 60
- Each scenario is simulated in two modes:
  - No IOV
  - IOV on CL across occasions

Outputs:
- Mean and quantile response summaries across population
- Decision metrics on mean response:
  - Emin
  - Time below threshold (80 percent of baseline)
  - Suppression AUC (baseline - response)

Validation:
- Generated artifacts are compared against committed expected artifacts.
- Metrics are compared with strict tolerance.
