#!/usr/bin/env bash
set -euo pipefail

julia use_cases/fih_dose_exploration/run.jl
julia use_cases/fih_dose_exploration/validate.jl
