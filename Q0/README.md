# Q0 - Baseline monocoeur

Ce dossier contient la version baseline de la simulation ACO (mesure globale `T_total`).

## Compilation

Depuis la racine du projet:

```bash
cd Q0
make all
```

## Execution benchmark Q0

```bash
cd Q0
make run
```

La cible `run` execute:

- `./Q0/ant_simu_q0.exe --gui 0 --iter 3000 --rep 5`

Le CSV est genere dans:

- `../results/Q0_baseline.csv`

## Execution manuelle

```bash
cd ..
./Q0/ant_simu_q0.exe --gui 0 --iter 3000 --rep 5
```
