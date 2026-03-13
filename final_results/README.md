# Archived Final Results

This folder stores the frozen CSV files and figures used by the final report.

## Typical contents

- final timing breakdown CSV files
- scaling summaries for Q3, Q4, and Q5
- PNG figures included in the report

## Intended use

Use this folder when you want to inspect the exact datasets cited in `report/rapport.tex` and `report/rapport.pdf`.

Standard `make run` targets usually write to the top-level `results/` directory first. Some summary commands can also be pointed here explicitly, for example:

```bash
make -C Q5 summarize_results RESULTS_DIR=../final_results
```
