#!/usr/bin/env bash
set -euo pipefail

RACINE_PROJET="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

cd "${RACINE_PROJET}"

echo "[Q1] Compilation de Q1..."
make -C Q1 all

echo "[Q1] Creation du dossier results/..."
mkdir -p results

echo "[Q1] Lancement du benchmark officiel Q1..."
./Q1/ant_simu_q1.exe --gui 0 --iter 3000 --rep 5

echo "[Q1] Termine. Resultat principal: results/Q1_timings_breakdown.csv"
