#!/usr/bin/env bash
# Run the whole pipeline. Steps 1-4, 6 and 7 are CPU work; step 5 (MD) wants a GPU.
# Every step skips work that is already done, so re-running after a failure is safe.
#   ./run_pipeline.sh           full run (3 x 10 ns per ligand)
#   ./run_pipeline.sh --quick   smoke test (100 ps per replica, separate output folders)
set -euo pipefail
cd "$(dirname "$0")"

QUICK=()
[[ "${1:-}" == "--quick" ]] && QUICK=(--quick)

python scripts/01_prepare_receptor.py
python scripts/02_parameterise_ligands.py
python scripts/03_dock.py
python scripts/04_build_systems.py
python scripts/05_run_md.py "${QUICK[@]}"
python scripts/06_mmgbsa.py "${QUICK[@]}"
python scripts/07_analyse.py "${QUICK[@]}"
