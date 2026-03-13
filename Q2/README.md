# Q2 - Optimized Sequential Baseline

`Q2` is the optimized sequential implementation used as the main baseline for later correctness checks. It keeps the same observable behavior as the sequential simulation while restructuring the hot loop and the data layout.

## Why Q2 matters

- it is the main sequential reference for `Q4` and `Q5`
- it reports per-stage timings
- it exports the reference `First_Iteration`, `Food_KPI`, and `Checksum`

## Requirements

- `make`
- `g++` with C++17 support
- SDL2 if you want GUI mode

## Build

```bash
make -C Q2 all
```

Executable produced:

- `Q2/ant_simu.exe`

## Official benchmark command

```bash
make -C Q2 run
```

This target:

1. runs `Q2/ant_simu.exe` in headless mode,
2. writes a local CSV to `Q2/results/Q2_timings_breakdown.csv`,
3. copies that CSV to `results/Q2_timings_breakdown.csv`.

Default benchmark command used by the Makefile:

```bash
./ant_simu.exe --gui 0 --iter 3000 --rep 5
```

## Manual execution

```bash
cd Q2
./ant_simu.exe --gui 0 --iter 3000 --rep 5
```

Local CSV output:

- `Q2/results/Q2_timings_breakdown.csv`

## Optional review script

`Q2/scripts/` contains helper material related to the vectorization review.

## Useful CLI options

```bash
cd Q2
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
make -C Q2 clean
```
