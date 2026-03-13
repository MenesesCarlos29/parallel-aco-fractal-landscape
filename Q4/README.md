# Q4 - MPI (Environnement replique)

## Resume de l'implementation

`Q4` implemente un schema MPI pur (sans OpenMP) avec environnement replique:

1. Chaque rang MPI construit la meme carte (terrain + pheromones).
2. Les fourmis sont reparties par sous-ensembles (`[start_ant, end_ant)` par rang).
3. Chaque rang simule uniquement ses fourmis locales.
4. Le compteur de nourriture est agrege globalement.
5. La carte de pheromones est synchronisee a chaque iteration par:
   `MPI_Allreduce(..., MPI_MAX, ...)`.
6. Les timings sont mesures avec `MPI_Wtime()`.
7. Seul `rank 0` ecrit les CSV de sortie.

Ce choix est volontairement celui de "Q4 - approche 1 (environnement replique)".

## Compilation et execution

Depuis la racine du projet:

```bash
cd Q4
make all
```

Execution simple:

```bash
cd Q4
make run P=4 ITER=4000 ANTS=5000
```

Benchmark officiel (sans hwthreads, avec binding coeur):

```bash
cd Q4
make run_scaling P_LIST="1 2 4" ITER=4000 REP=1 ANTS=5000
```

Benchmark ponctuel sur un `P` fixe:

```bash
cd Q4
make benchmark P_BENCH=4 REP_BENCH=5 ITER=4000 ANTS=5000
```

Smoke test non officiel (si `np=8` n'est pas possible en coeurs reels):

```bash
cd Q4
make run_scaling_smoke P_LIST_SMOKE="8" ITER=4000 REP=1 ANTS=5000
```

Resume automatique des CSV de scaling:

```bash
cd Q4
make summarize_results
```

Sorties CSV:

1. `Q4/results/Q4_timings_breakdown.csv` (sortie immediate d'un run).
2. `results/Q4_scaling_np*.csv` (copie pour comparaison inter-processus).
3. `results/Q4_scaling_summary.csv` (resume calcule par `summarize_results`).

## Resultats observes

Table issue de:

1. `results/Q4_scaling_np1.csv`
2. `results/Q4_scaling_np2.csv`
3. `results/Q4_scaling_np4.csv`
4. `results/Q4_scaling_np8.csv`

Baseline `T1 = 5106.248294 ms` (`np=1`).

`results/Q4_timings_breakdown.csv` reste utile pour inspecter une execution ponctuelle,
mais la comparaison inter-processus doit se faire avec `Q4_scaling_np*.csv`.

| P | T_compute (ms) | T_comm (ms) | T_total (ms) | Speedup vs np=1 | Efficiency | Comm % |
|---|---------------:|------------:|-------------:|----------------:|-----------:|-------:|
| 1 | 5093.754156 | 10.251623 | 5106.248294 | 1.0000 | 1.0000 | 0.20 |
| 2 | 5155.694814 | 4927.521773 | 10071.528715 | 0.5070 | 0.2535 | 48.93 |
| 4 | 6398.531615 | 13341.296475 | 19556.918251 | 0.2611 | 0.0653 | 68.22 |
| 8 | 10943.162797 | 38330.687776 | 48693.749867 | 0.1049 | 0.0131 | 78.72 |

Lecture rapide:

1. Le temps total augmente quand `P` augmente.
2. `T_comm` devient dominant tres vite.
3. Le gain sur `move ants` ne compense pas les couts de synchronisation.

## Analyse technique: pourquoi ca ne scale pas

1. La reduction MPI des pheromones porte sur toute la grille a chaque iteration.
2. Le volume reduit par iteration est important (`~4.05 MiB` pour `grid=513`).
3. L'evaporation reste calculee sur la grille complete sur chaque rang: c'est du calcul replique.
4. Le travail qui se parallelise vraiment (`advance_one`) baisse avec `P`, mais la communication globale augmente plus vite.
5. Le cas `np=8` lance avec `--use-hwthread-cpus` ne doit pas etre interprete comme benchmark officiel "1 core/process".
6. Le comportement observe est coherent avec l'approche "environnement replique", pas necessairement un bug de mesure.

## Limitations actuelles

1. L'approche repliquee a une limite naturelle de scalabilite.
2. `T_comm` finit par dominer le temps total.
3. `np=1` reste la reference la plus rapide dans ce contexte.
4. `np>1` est surtout utile ici pour caracteriser le cout de communication.
5. Le `Checksum` peut differer pour `np>1`; ce point est connu et hors objectif de cette passe.

