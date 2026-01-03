#!/usr/bin/env bash
set -euo pipefail

./bin/openpkpd version
./bin/openpkpd replay --artifact validation/golden/pk_iv_bolus.json
./bin/openpkpd validate-golden
