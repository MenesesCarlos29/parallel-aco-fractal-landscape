# Parallel ACO on a Fractal Landscape

This repository contains the same ant-colony simulation implemented across six progressively more advanced stages:

- `Q0`: single-core baseline
- `Q1`: single-core instrumentation with per-stage timings
- `Q2`: data-layout and hot-loop optimization (the main sequential baseline for correctness)
- `Q3`: OpenMP shared-memory parallelization
- `Q4`: MPI replicated-environment parallelization
- `Q5`: MPI domain decomposition with halos and ant migration

Each `Qx/` folder has its own `README.md` with exact build and run instructions for that stage.

## Repository layout

- `Q0/` to `Q5/`: source code and stage-specific Makefiles
- `final_results/`: archived CSV files and figures used by the final report
- `report/`: LaTeX report (`rapport.tex`, `rapport.pdf`)
- `scripts/`: helper scripts for batch runs and Q3 plotting
- `results/`: top-level scratch/output directory populated by some `make run` targets

## Requirements

The exact toolchain depends on the stage you want to run.

### Common requirements

- Linux or WSL with `bash`
- `make`
- a C++17 compiler (`g++` for `Q0` to `Q3`, `mpicxx` for `Q4` and `Q5`)
- SDL2 for GUI mode (`--gui 1`)

### Additional requirements by stage

- `Q3`: OpenMP support (`g++` with `-fopenmp`)
- `Q4`, `Q5`: OpenMPI or another MPI implementation providing `mpicxx` and `mpirun`
- `scripts/plot_q3_results.py`: Python 3 with `numpy`, `pandas`, and `matplotlib`

## Quick start

Build and run any stage from the repository root.

```bash
make -C Q0 run
make -C Q1 run
make -C Q2 run
make -C Q3 run_q3
make -C Q4 run_scaling P_LIST="1 2 4" ITER=4000 REP=1 ANTS=5000
make -C Q5 run_scaling P_LIST="1 2 4" ITER=4000 REP=1 ANTS=5000
```

For a list of targets in any stage:

```bash
make -C Q4 help
```

## Stage summary

| Stage | Main idea | Build | Typical run |
|---|---|---|---|
| `Q0` | single-core baseline | `make -C Q0 all` | `make -C Q0 run` |
| `Q1` | stage-level timing instrumentation | `make -C Q1 all` | `make -C Q1 run` |
| `Q2` | optimized sequential/vector-friendly layout | `make -C Q2 all` | `make -C Q2 run` |
| `Q3` | OpenMP scaling | `make -C Q3 all` | `make -C Q3 run_q3` |
| `Q4` | MPI replicated environment | `make -C Q4 all` | `make -C Q4 run_scaling P_LIST="1 2 4" ITER=4000 REP=1 ANTS=5000` |
| `Q5` | MPI domain decomposition | `make -C Q5 all` | `make -C Q5 run_scaling P_LIST="1 2 4" ITER=4000 REP=1 ANTS=5000` |

## Common CLI options

All executables accept the same main command-line interface, although the default values can differ slightly between stages:

- `--grid-size N`
- `--ants N`
- `--alpha X`
- `--beta X`
- `--eps X`
- `--iter N`
- `--rep N`
- `--seed N`
- `--gui 0|1`
- `--help`

Headless benchmarking should always use `--gui 0`.

## Notes about outputs

Output paths are not fully uniform across stages:

- `Q0` and `Q1` write their benchmark CSV to the top-level `results/` directory when you use `make run`.
- `Q2` to `Q5` generate a local CSV inside `Qx/results/`, then some Makefile targets also copy selected files to the top-level `results/` directory.
- `final_results/` stores archived datasets and figures used in the report; it is not the default runtime output directory.

## Report

Compile the final report with:

```bash
cd report
pdflatex -interaction=nonstopmode -halt-on-error rapport.tex
```

The resulting PDF is `report/rapport.pdf`.

## Recommended reading order

If you are reviewing the project for the first time, use this order:

1. Read this top-level README.
2. Read the README inside the stage you want to run.
3. Use `make -C Qx help` to inspect the available targets.
4. Open `report/rapport.pdf` for the methodological and performance analysis.
