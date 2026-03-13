# Helper Scripts

This folder contains convenience scripts used during development and analysis. They are not required to build the core executables.

## Available scripts

### `run_all.sh`

Runs the default `Q0` and `Q1` benchmark targets from the repository root.

```bash
./scripts/run_all.sh
```

### `run_q1.sh`

Builds and runs the default `Q1` benchmark.

```bash
./scripts/run_q1.sh
```

### `plot_q3_results.py`

Generates Q3 scaling figures from CSV files stored in `final_results/` or `results/`.

Requirements:

- Python 3
- `numpy`
- `pandas`
- `matplotlib`

Usage:

```bash
python3 scripts/plot_q3_results.py
```

Output figures are written to the top-level `results/` directory.
