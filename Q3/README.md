# Q3 - OpenMP Shared-Memory Parallelization

`Q3` is the OpenMP version of the simulation. It keeps the `Q2`-style timing breakdown and adds shared-memory scaling experiments over multiple thread counts.

## What this stage adds

- OpenMP parallel regions in the hot path
- per-stage timing CSV
- scaling workflow for `1`, `2`, `4`, and `8` threads
- consolidated OpenMP scaling CSV

## Requirements

- `make`
- `g++` with C++17 support
- OpenMP support (`-fopenmp`)
- SDL2 if you want GUI mode

## Build

```bash
make -C Q3 all
```

Executable produced:

- `Q3/ant_simu.exe`

## Main commands

### Single benchmark run

```bash
make -C Q3 run
```

This writes:

- local CSV: `Q3/results/Q3_timings_breakdown.csv`
- copied CSV: `results/Q3_timings_breakdown.csv`

### Full OpenMP scaling campaign

```bash
make -C Q3 run_q3
```

This runs the benchmark for `OMP_NUM_THREADS=1 2 4 8` and updates:

- `results/Q3_openmp_scaling.csv`

### One thread count only

```bash
make -C Q3 run_q3_thread THREADS=4
```

This also stores:

- `results/Q3_t4_timings_breakdown.csv`

When `THREADS=1`, the target also refreshes:

- `results/Q3_threads1_reference.csv`

## Manual execution

```bash
cd Q3
OMP_NUM_THREADS=4 ./ant_simu.exe --gui 0 --iter 4000 --rep 5
```

Manual execution writes the local file:

- `Q3/results/Q3_timings_breakdown.csv`

## Useful CLI options

```bash
cd Q3
./ant_simu.exe --help
```

Main options:

- `--grid-size N`
- `--ants N`
- `--alpha X`
- `--beta X`
- `--eps X`
- `--iter N`
- `--rep N`
- `--seed N`
- `--gui 0|1`

## Cleaning

```bash
make -C Q3 clean
```
