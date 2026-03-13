# Q5 - MPI Domain Decomposition

`Q5` is the second MPI implementation. Unlike `Q4`, it does not replicate the whole domain on every rank. Instead, it decomposes the grid by rows and exchanges only halo data plus migrating ants.

## Execution model

- MPI only, no OpenMP
- 1D row-wise domain decomposition
- each rank owns a contiguous range of global rows
- local pheromone storage includes two halo rows
- halo exchange is performed only with neighbor ranks
- ants migrate when they cross a subdomain boundary
- only rank `0` writes CSV files
- timings are measured with `MPI_Wtime()`

## What the CSV reports

- `First_Iteration`
- `Food_KPI`
- `Checksum`
- `Total_Ants_Global`
- compute time split:
  - `T_compute_move_ms`
  - `T_compute_evap_ms`
  - `T_compute_update_ms`
  - `T_compute_ms`
- communication time split:
  - `T_comm_halo_ms`
  - `T_comm_migration_ms`
  - `T_comm_food_ms`
  - `T_comm_ms`
- load-balance metrics:
  - `Ants_Min`
  - `Ants_Avg`
  - `Ants_Max`
  - `Ants_Stddev`
  - `Imbalance_Ratio`

## Requirements

- `make`
- `mpicxx` and `mpirun`
- SDL2 if you want GUI mode, although MPI benchmarks should use `--gui 0`

## Build

```bash
make -C Q5 all
```

Executable produced:

- `Q5/ant_simu.exe`

## Main commands

### One MPI run

```bash
make -C Q5 run P=2 ITER=4000 REP=1 ANTS=5000
```

Generated outputs:

- local CSV: `Q5/results/Q5_timings_breakdown.csv`
- copied CSV: `results/Q5_run_np2.csv`

### Fixed-process benchmark

```bash
make -C Q5 benchmark P_BENCH=4 ITER=4000 REP_BENCH=5 ANTS=5000
```

Generated output:

- `results/Q5_benchmark_p4.csv`

### Scaling campaign

```bash
make -C Q5 run_scaling P_LIST="1 2 4" ITER=4000 REP=1 ANTS=5000
```

Generated outputs:

- `results/Q5_scaling_np1.csv`
- `results/Q5_scaling_np2.csv`
- `results/Q5_scaling_np4.csv`

If you want explicit core binding, pass it through `MPI_BIND_FLAGS`:

```bash
make -C Q5 run_scaling P_LIST="1 2 4" ITER=4000 REP=1 ANTS=5000 MPI_BIND_FLAGS="--bind-to core --map-by core"
```

### Summarize scaling CSVs

```bash
make -C Q5 summarize_results
```

This requires the `results/Q5_scaling_np*.csv` files to exist first.

To summarize the archived datasets stored in `final_results/` instead of the runtime `results/` directory:

```bash
make -C Q5 summarize_results RESULTS_DIR=../final_results
```

## Manual execution

From the repository root:

```bash
mpirun -np 4 Q5/ant_simu.exe --gui 0 --iter 4000 --rep 1 --ants 5000
```

With explicit binding:

```bash
mpirun --bind-to core --map-by core -np 4 Q5/ant_simu.exe --gui 0 --iter 4000 --rep 1 --ants 5000
```

## Validation notes

`Q5` uses `Q2/results/Q2_timings_breakdown.csv` as its baseline reference when that file is available. The validation logic reports the comparison status explicitly instead of silently assuming a match.

## Cleaning

```bash
make -C Q5 clean
```
