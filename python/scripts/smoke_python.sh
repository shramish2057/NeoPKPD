#!/usr/bin/env bash
set -euo pipefail

# Use existing venv if present, otherwise create one
if [ -d "python/.venv" ]; then
    source python/.venv/bin/activate
else
    python3 -m venv python/.venv
    source python/.venv/bin/activate
fi

pip install -e python
pip install pytest

pytest -q python/tests

python python/examples/write_pk_iv_bolus.py
