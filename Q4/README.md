# Q1 - Instrumentation et timings par etape

Ce dossier contient la version instrumentee de la simulation ACO.

## Compilation

Depuis la racine du projet:

```bash
cd Q1
make all
```

## Execution benchmark Q1

```bash
cd Q1
make run
```

La cible `run` execute:

- `./Q1/ant_simu_q1.exe --gui 0 --iter 3000 --rep 5`

Le CSV est genere dans:

- `../results/Q1_timings_breakdown.csv`

## Execution manuelle

```bash
cd ..
./Q1/ant_simu_q1.exe --gui 0 --iter 3000 --rep 5
```
