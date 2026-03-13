# Report Folder

This folder contains the final written report for the project.

## Main files

- `rapport.tex`: LaTeX source of the report
- `rapport.pdf`: compiled PDF
- `imgs/`: figures and static assets used by the report

## Build the PDF

```bash
cd report
pdflatex -interaction=nonstopmode -halt-on-error rapport.tex
```

Run the command a second time if LaTeX asks for another pass to refresh the table of contents or references.

## Notes

- The report is written in French because it was prepared for the course submission.
- The numerical results cited in the report come from `../final_results/`.
