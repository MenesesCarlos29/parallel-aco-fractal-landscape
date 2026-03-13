# Q4 - MPI Replicated Environment

`Q4` is the first MPI implementation. It follows the replicated-environment strategy: every MPI rank owns the full terrain and pheromone map, but only simulates a deterministic subset of the ants.

## Execution model

- MPI only, no OpenMP
- each rank builds the same global environment
- ant ownership is partitioned by index range
- pheromones are synchronized globally with `MPI_Allreduce(..., MPI_MAX, ...)`
- `Food_KPI` is aggregated across ranks
- only rank `0` writes the CSV files
- timings are measured with `MPI_Wtime()`

## Requirements

- `make`
- `mpicxx` and `mpirun`
- SDL2 if you want GUI mode, although benchmarking should use `--gui 0`

## Build

```bash
make -C Q4 all
```

Executable produced:

- `Q4/ant_simu.exe`

## Main commands

### One MPI run

```bash
make -C Q4 run P=4 ITER=4000 ANTS=5000
```

This uses the `run` target with:

- `P`: number of MPI processes
- `ITER`: iteration count
- `ANTS`: number of ants

The target runs with core binding by default:

```bash
mpirun --bind-to core --map-by core -np 4 ./ant_simu.exe --gui 0 --ants 5000 --iter 4000 --rep 1
```

Generated outputs:

- local CSV: `Q4/results/Q4_timings_breakdown.csv`
- copied CSV: `results/Q4_run_np4.csv`

### Official scaling campaign

```bash
make -C Q4 run_scaling P_LIST="1 2 4" ITER=4000 REP=1 ANTS=5000
```

Generated outputs:

- `results/Q4_scaling_np1.csv`
- `results/Q4_scaling_np2.csv`
- `results/Q4_scaling_np4.csv`

### Fixed-process benchmark

```bash
make -C Q4 benchmark P_BENCH=4 REP_BENCH=5 ITER=4000 ANTS=5000
```

Generated output:

- `results/Q4_mpi_timings_p4.csv`

### Non-official smoke test with hardware threads

```bash
make -C Q4 run_scaling_smoke P_LIST_SMOKE="8" ITER=4000 REP=1 ANTS=5000
```

Generated output:

- `results/Q4_scaling_smoke_np8.csv`

### Summarize archived scaling CSVs

```bash
make -C Q4 summarize_results
```

This requires the `results/Q4_scaling_np*.csv` files to exist first.

## Manual execution

From the repository root:

```bash
mpirun --bind-to core --map-by core -np 4 Q4/ant_simu.exe --gui 0 --iter 4000 --rep 1 --ants 5000
```

## Validation notes

`Q4` validates against the `Q2` baseline for `np=1` using `First_Iteration`, `Food_KPI`, and the current checksum logic implemented in `Q4`.

## Cleaning

```bash
make -C Q4 clean
```
