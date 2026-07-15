#!/usr/bin/env bash
set -euo pipefail

export PYTHONPATH="${PYTHONPATH:-}:$(pwd)/src"

python -m fast_glycan_masking.glycan_rotamer_generator \
  --pdb examples/reference_glycoprotein.pdb \
  --chain A \
  --resid 200 \
  --config configs/man5_sampling.yaml \
  --n-conformers 10000 \
  --out-prefix Man5_library
