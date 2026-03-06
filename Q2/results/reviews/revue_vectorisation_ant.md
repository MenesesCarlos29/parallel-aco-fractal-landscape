# Revue vectorisation Q2

- Fichier analyse: `ant.cpp`
- Signaux detectes: **25**
- Signaux `optimized`: **0**
- Signaux `missed`: **25**
- Signaux utiles (loops/causes actionnables): **9**
- Utiles `optimized`: **0**
- Utiles `missed`: **9**

## Table prioritaire (a corriger en premier)

| Ligne | Colonne | Statut | Message |
|---:|---:|:---:|---|
| 61 | 31 | missed | couldn't vectorize loop |
| 97 | 39 | missed | couldn't vectorize loop |
| 97 | 39 | missed | not vectorized: loop nest containing two or more consecutive inner loops cannot be vectorized |
| 126 | 54 | missed | couldn't vectorize loop |
| 126 | 54 | missed | not vectorized: number of iterations cannot be computed. |
| 126 | 54 | missed | couldn't vectorize loop |
| 126 | 54 | missed | not vectorized: number of iterations cannot be computed. |
| 79 | 33 | missed | couldn't vectorize loop |
| 70 | 6 | missed | not vectorized: unsupported use in stmt. |

## Causes principales (prioritaires)

| Occurrences | Cause |
|---:|---|
| 5 | couldn't vectorize loop |
| 2 | not vectorized: number of iterations cannot be computed. |
| 1 | not vectorized: unsupported use in stmt. |
| 1 | not vectorized: loop nest containing two or more consecutive inner loops cannot be vectorized |

## Table complete (tous signaux)

| Ligne | Colonne | Statut | Message |
|---:|---:|:---:|---|
| 61 | 31 | missed | couldn't vectorize loop |
| 61 | 31 | missed | not vectorized: latch block not empty. |
| 97 | 39 | missed | couldn't vectorize loop |
| 97 | 39 | missed | not vectorized: loop nest containing two or more consecutive inner loops cannot be vectorized |
| 126 | 54 | missed | couldn't vectorize loop |
| 126 | 54 | missed | not vectorized: number of iterations cannot be computed. |
| 126 | 54 | missed | couldn't vectorize loop |
| 126 | 54 | missed | not vectorized: number of iterations cannot be computed. |
| 79 | 33 | missed | couldn't vectorize loop |
| 70 | 6 | missed | not vectorized: unsupported use in stmt. |
| 10 | 11 | missed | not vectorized: statement can throw an exception: AntSwarm::resize (this_2(D) |
| 8 | 1 | missed | not vectorized: statement can throw an exception: resx 15 |
| 8 | 1 | missed | not vectorized: statement can throw an exception: resx 14 |
| 8 | 1 | missed | not vectorized: statement can throw an exception: resx 13 |
| 8 | 1 | missed | not vectorized: statement can throw an exception: resx 12 |
| 8 | 1 | missed | not vectorized: statement can throw an exception: resx 11 |
| 8 | 1 | missed | not vectorized: statement can throw an exception: resx 10 |
| 8 | 1 | missed | not vectorized: statement can throw an exception: resx 9 |
| 8 | 1 | missed | not vectorized: statement can throw an exception: resx 8 |
| 8 | 1 | missed | not vectorized: statement can throw an exception: resx 7 |
| 8 | 1 | missed | not vectorized: statement can throw an exception: resx 6 |
| 8 | 1 | missed | not vectorized: statement can throw an exception: resx 5 |
| 8 | 1 | missed | not vectorized: statement can throw an exception: resx 4 |
| 8 | 1 | missed | not vectorized: statement can throw an exception: resx 3 |
| 8 | 1 | missed | not vectorized: statement can throw an exception: resx 2 |

## Causes principales (complet)

| Occurrences | Cause |
|---:|---|
| 5 | couldn't vectorize loop |
| 2 | not vectorized: number of iterations cannot be computed. |
| 1 | not vectorized: unsupported use in stmt. |
| 1 | not vectorized: statement can throw an exception: resx 9 |
| 1 | not vectorized: statement can throw an exception: resx 8 |
| 1 | not vectorized: statement can throw an exception: resx 7 |
| 1 | not vectorized: statement can throw an exception: resx 6 |
| 1 | not vectorized: statement can throw an exception: resx 5 |
| 1 | not vectorized: statement can throw an exception: resx 4 |
| 1 | not vectorized: statement can throw an exception: resx 3 |
