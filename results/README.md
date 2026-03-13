# Runtime Results Directory

This folder is the default scratch/output directory used by several `make run` targets.

## What goes here

Depending on the stage, this directory may receive:

- copied benchmark CSV files from `Q0` to `Q5`
- consolidated scaling CSV files
- generated plots from helper scripts

## Important distinction

- `results/` is the working output directory used during local runs.
- `final_results/` contains the archived datasets and figures referenced by the final report.

If you rerun benchmarks, expect files in this directory to be replaced.
