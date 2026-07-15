#!/usr/bin/env bash
set -euo pipefail

export PYTHONPATH="${PYTHONPATH:-}:$(pwd)/src"

python -m fast_glycan_masking.glycan_library_placer \
  --library Man5_library.npz \
  --protein examples/designed_antigen.pdb \
  --site A:87 \
  --site A:200 \
  --n-models 100 \
  --clash-mode vdw \
  --vdw-scale 0.70 \
  --out-dir outputs/A87_A200 \
  --prefix Glyc_Des_A87_A200
