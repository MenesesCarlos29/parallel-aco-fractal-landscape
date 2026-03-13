# Q0 - Single-Core Baseline

`Q0` is the plain sequential reference version of the ant-colony simulation. It measures only the total runtime and exports the main functional indicators used throughout the project.

## What this stage contains

- single-process, single-thread execution
- optional SDL visualization
- baseline CSV with:
  - runtime
  - `First_Iteration`
  - `Food_KPI`
  - `Checksum`

## Requirements

- `make`
- `g++` with C++17 support
- SDL2 if you want `--gui 1`

## Build

From the repository root:

```bash
make -C Q0 all
```

Executable produced:

- `Q0/ant_simu_q0.exe`

## Official benchmark command

```bash
make -C Q0 run
```

This target runs:

```bash
./Q0/ant_simu_q0.exe --gui 0 --iter 3000 --rep 5
```

When launched through `make -C Q0 run`, the program runs from the repository root, so the CSV is written to:

- `results/Q0_baseline.csv`

## Manual execution

From inside `Q0/`:

```bash
cd Q0
./ant_simu_q0.exe --gui 0 --iter 3000 --rep 5
```

In that case, the output path becomes local to the folder:

- `Q0/results/Q0_baseline.csv`

## Useful CLI options

```bash
./ant_simu_q0.exe --help
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
make -C Q0 clean
```
