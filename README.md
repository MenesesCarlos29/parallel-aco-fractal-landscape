# parallel-aco-fractal-landscape
Parallel Ant Colony Optimization (ACO) on a fractal landscape, with OpenMP/MPI acceleration and performance analysis.
# Projet Final — ACO (Ant Colony Optimisation) sur paysage fractal

## 1) Objectif (résumé)
Optimiser et paralléliser la simulation ACO :
1) Mesurer le temps par itération (timers + breakdown)
2) Vectorisation (SoA)
3) Parallélisation mémoire partagée (OpenMP)
4) Parallélisation distribuée (MPI) : Méthode 1 (implémentation)
5) MPI Méthode 2 : stratégie (bonus : implémentation)

---

## 2) Base fournie (fichiers principaux)
- `ant_simu.cpp` : main (init, boucle temporelle, render)
- `ant.hpp/.cpp` : logique de mouvement + état (chargée/non chargée)
- `pheronome.hpp` : V1/V2, propagation (alpha), évaporation (beta), bords
- `fractal_land.hpp/.cpp` : génération du terrain fractal (plasma)
- `rand_generator.hpp` : RNG déterministe (seed + indices)
- `window.hpp/.cpp`, `renderer.hpp/.cpp` : SDL (affichage)

---

## 3) Organisation du dépôt (1 dossier par question)
Chaque dossier est une version autonome, compilable et testable :

- `Q0_baseline/` : code fourni (référence)
- `Q1_timers/` : instrumentation + mode headless
- `Q2_vector/` : vectorisation (SoA)
- `Q3_openmp/` : OpenMP + mesures speedup/efficacité
- `Q4_mpi_method1/` : MPI méthode 1 (Allreduce MAX) + mesures
- `Q5_mpi_method2_bonus/` : STRATEGY (obligatoire) + bonus code si possible

Dans chaque `Qx/` :
- `Makefile`
- `README.md` local (build/run + paramètres)
- `results/` (CSV + figures)

---

## 4) Répartition du travail (binôme) + rapport croisé

### Personne A — Implémentation + mesures
**À faire**
- `Q1_timers/`
  - Ajouter mode headless (sans SDL) pour benchmark
  - Ajouter timers (total + phases)
  - Export CSV (moyenne + écart-type sur N runs)
- `Q3_openmp/`
  - OpenMP sur les boucles adaptées (évaporation + dépôts thread-local si besoin)
  - Courbes temps/speedup/efficacité vs #threads

**Rapport**
- Personne B rédige les sections du rapport correspondant à `Q1` et `Q3`.

---

### Personne B — Implémentation + mesures
**À faire**
- `Q2_vector/`
  - Vectoriser (SoA) : positions/états/graines en tableaux
  - Comparer vs baseline
- `Q4_mpi_method1/`
  - MPI méthode 1 : carte complète par processus + fourmis réparties
  - Synchro phéromones avec `MPI_Allreduce(..., MPI_MAX, ...)`
- `Q5_mpi_method2_bonus/`
  - Écrire une stratégie claire (obligatoire)
  - Bonus : code (ghost cells + migration fourmis) si temps

**Rapport**
- Personne A rédige les sections du rapport correspondant à `Q2`, `Q4`, `Q5`.

---

## 5) Règles pour éviter les blocages (important)
- Chacun ne modifie QUE ses dossiers `Qx/`.
- Fichiers partagés autorisés :
  - `README.md` (racine) : organisation globale uniquement
  - `report/` : assemblage final (à faire à la fin, via PR)
- Résultats :
  - CSV : `results/*.csv`
  - Figures : `results/fig_*.png`
- Paramètres du benchmark : doivent être IDENTIQUES entre dossiers (voir contrat).

---

## 6) Contrat de paramètres (BENCHMARK PROTOCOL) — à compléter

> Objectif : garantir des mesures comparables entre versions (baseline, vector, OMP, MPI).
> Beaucoup d’infos sont probablement déjà dans le code séquentiel : on doit les extraire et figer ces valeurs.

### 6.1 Mesure : définition exacte
- Mode benchmark : **HEADLESS** (sans SDL/render) : **[TBD: oui/non]**
- Est-ce qu’on exclut init + génération fractal du timing ? **[TBD]**
- Une “itération” = (update phéromones + move fourmis + évaporation + …) : **[TBD]**

**Où regarder dans le code (baseline)**
- Boucle principale : `ant_simu.cpp` (rechercher `while`, `for` de simulation, ou compteur de pas)
- Render : `renderer.cpp` / `window.cpp` (si inclus dans le loop)

---

### 6.2 Paramètres de l’environnement (terrain fractal)
- Taille grille (W×H) : **W = [TBD]**, **H = [TBD]**
- Paramètres “plasma” / fractal :
  - `n` (nb sous-grilles) : **[TBD]**
  - `ns = 2^k + 1` : **k = [TBD]**, donc **ns = [TBD]**
  - déviation `d` : **[TBD]**
- Normalisation terrain dans [0,1] : **[TBD: confirmé]**

**Où regarder**
- `fractal_land.cpp` : rechercher `ns`, `k`, `deviation`, `seed`, `normalize`, `subdiv`
- `fractal_land.hpp` : constructeur + paramètres

---

### 6.3 Paramètres ACO (modèle)
- Nombre de fourmis `m` : **[TBD]**
- Init positions :
  - ( ) toutes à la fourmilière
  - ( ) uniformes sur la grille
  - Choix : **[TBD]**
- α (bruit / propagation) : **[TBD]**
- β (évaporation, proche de 1) : **[TBD]**
- ε (taux d’exploration) : **[TBD]**
- Nombre de sources de nourriture : **[TBD]**
- Placement obstacles (-1) : **[TBD]** (fixe / aléatoire / %)

**Où regarder**
- `pheronome.hpp` : alpha/beta + update V1/V2
- `ant.hpp/.cpp` : epsilon + règle de mouvement + état chargée/non chargée
- `ant_simu.cpp` : création fourmilière/nourriture/obstacles + `m`

---

### 6.4 RNG et reproductibilité
- Seed globale : **[TBD]**
- Seed dépendante index (fourmi / cellule) : **[TBD: oui/non]**
- Nombre de seeds testées :
  - ( ) 1 seed fixe
  - ( ) plusieurs seeds (ex: 3) + moyenne
  - Choix : **[TBD]**

**Où regarder**
- `rand_generator.hpp` : fonctions `rand_int32`, `rand_double`, usage de seed
- appels RNG dans `fractal_land.cpp` et `ant.cpp`

---

### 6.5 Protocole d’expériences
- Nb d’itérations (pas de temps) par run : **[TBD]**
- Nb de runs par configuration (pour stats) : **[TBD]**
- Stats à reporter : moyenne + écart-type : **[TBD: oui/non]**

---

### 6.6 OpenMP (Q3)
- Liste #threads testés : **[TBD]** (ex: 1,2,4,8,…)
- Scheduling (static/dynamic) : **[TBD]**
- Limite : jusqu’à ( ) cœurs physiques ( ) threads logiques : **[TBD]**

---

### 6.7 MPI (Q4)
- Méthode 1 : carte complète sur chaque processus : **OK**
- #processus testés : **[TBD]** (ex: 1,2,4,8,…)
- Contrat CPU : 1 cœur / processus : **[TBD: oui/non]**
- Synchronisation phéromones : `MAX` global : **OK**
- Fréquence synchro :
  - ( ) chaque itération
  - ( ) toutes les K itérations (K = [TBD])
- Données synchronisées :
  - ( ) V1 seulement
  - ( ) V2 seulement
  - ( ) V1 + V2
  - Choix : **[TBD]**

---

## 7) Checklist “Done” (pour chaque Qx)
- `make` fonctionne
- `README.md` local indique build + run + paramètres
- `results/` contient au minimum :
  - 1 CSV (timing ou speedup)
  - 1 figure OU tableau
- Paramètres utilisés = contrat (section 6)

---

## 8) Notes pour la soutenance
Chaque personne doit savoir expliquer l’autre partie :
- Personne A : comprend Q2/Q4/Q5
- Personne B : comprend Q1/Q3
Objectif : défense fluide + cross-check des résultats.