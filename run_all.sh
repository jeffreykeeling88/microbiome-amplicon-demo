#!/usr/bin/env bash
# Run the full pipeline: data -> QC -> DADA2 -> mock check -> diversity -> report
set -euo pipefail
cd "$(dirname "$0")"

bash scripts/00_get_data.sh
Rscript scripts/01_qc.R
Rscript scripts/02_dada2.R
Rscript scripts/03_mock_validation.R
Rscript scripts/04_diversity.R
Rscript -e 'rmarkdown::render("report/report.Rmd", output_dir = "results", quiet = TRUE)'
echo "All done. Open results/report.html"
