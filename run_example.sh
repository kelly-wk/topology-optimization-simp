#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
PYTHON_BIN="${PYTHON_BIN:-python3}"
"$PYTHON_BIN" topology_simp.py \
  --nelx 18 \
  --nely 6 \
  --volfrac 0.40 \
  --max-iter 80 \
  --output outputs
mkdir -p docs/assets
cp outputs/final_topology.svg docs/assets/final_topology.svg
cp outputs/summary.json docs/summary.json
"$PYTHON_BIN" -m unittest discover -s tests -v
