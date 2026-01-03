#!/usr/bin/env bash
set -euo pipefail

# First instantiate the core package
julia -e 'using Pkg; Pkg.activate("core/OpenPKPDCore"); Pkg.instantiate()'

# Develop the local OpenPKPDCore package and instantiate CLI dependencies
julia -e '
using Pkg
Pkg.activate("cli/OpenPKPDCLI")
Pkg.develop(path="core/OpenPKPDCore")
Pkg.instantiate()
'

./bin/openpkpd version
./bin/openpkpd replay --artifact validation/golden/pk_iv_bolus.json
./bin/openpkpd validate-golden
