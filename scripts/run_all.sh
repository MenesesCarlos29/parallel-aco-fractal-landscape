#!/usr/bin/env bash
set -euo pipefail

RACINE_PROJET="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

cd "${RACINE_PROJET}"

echo "[GLOBAL] Execution de Q0..."
make -C Q0 run

echo "[GLOBAL] Execution de Q1..."
make -C Q1 run

echo "[GLOBAL] Termine. Resultats disponibles dans results/."
