# Q1 - Sequential Instrumentation

`Q1` keeps the simulation sequential but adds per-stage timing instrumentation. It is the first version that breaks the execution time into separate phases.

## What this stage measures

The CSV reports the runtime split into:

- `T_move_ants`
- `T_evap`
- `T_pher_update`
- `T_render`
- `T_total`

It also keeps the functional indicators:

- `First_Iteration`
- `Food_KPI`
- `Checksum`

## Requirements

- `make`
- `g++` with C++17 support
- SDL2 if you want GUI mode

## Build

```bash
make -C Q1 all
```

Executable produced:

- `Q1/ant_simu_q1.exe`

## Official benchmark command

```bash
make -C Q1 run
```

This target runs:

```bash
./Q1/ant_simu_q1.exe --gui 0 --iter 3000 --rep 5
```

When launched through `make -C Q1 run`, the CSV is written to:

- `results/Q1_timings_breakdown.csv`

## Manual execution

```bash
cd Q1
./ant_simu_q1.exe --gui 0 --iter 3000 --rep 5
```

Manual execution from inside `Q1/` writes to:

- `Q1/results/Q1_timings_breakdown.csv`

## Useful CLI options

```bash
./Q1/ant_simu_q1.exe --help
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
make -C Q1 clean
```
