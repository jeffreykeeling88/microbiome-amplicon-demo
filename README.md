# 16S rRNA Amplicon Analysis: From Raw Reads to Publication Figures

A complete, reproducible microbiome analysis pipeline built on **DADA2**, **phyloseq**, and **vegan**. One command takes paired-end Illumina FASTQ files to an ASV table, taxonomy, statistics, figures, and a client-ready HTML report.

**This is a demo of the analysis I deliver for clients.** If you have 16S, ITS, or 18S amplicon data that needs processing, [contact me on Upwork](#contact).

---

## What you get

| Deliverable | Description |
|---|---|
| **HTML report** | Plain-language summary, every figure, statistical results, and a methods paragraph you can paste into a manuscript |
| **ASV table + taxonomy** | Counts and SILVA v138.2 classifications for every amplicon sequence variant |
| **phyloseq object** | Ready to load in R for any further analysis |
| **Figures** | Quality profiles, read tracking, alpha diversity, PCoA, taxonomic composition, differential abundance (200 dpi PNG) |
| **Statistics tables** | Alpha diversity tests, PERMANOVA, dispersion tests, differential abundance with FDR correction |
| **Code** | Commented scripts so you or a reviewer can reproduce every number |

## Demo dataset

The pipeline runs on the mothur MiSeq SOP dataset (Kozich et al. 2013): 16S V4 amplicons from mouse fecal samples collected early (days 0 to 9) and late (days 141 to 150) after weaning, plus a 20-strain **mock community** used to verify accuracy. The mock community check confirms the pipeline recovers the known strains without inflating diversity with sequencing errors.

## Pipeline

```
Raw FASTQ
  │  01_qc.R              quality profiles → choose truncation lengths
  ▼
  │  02_dada2.R           filter → learn errors → denoise → merge → remove chimeras → assign taxonomy
  ▼
ASV table + taxonomy
  │  03_mock_validation.R exact-match check against known mock strains
  │  04_diversity.R       alpha · beta (PERMANOVA, betadisper) · composition · differential abundance
  ▼
report/report.Rmd → results/report.html
```

## Methods at a glance

- **Denoising:** DADA2 error modeling and exact sequence variant inference (no OTU clustering)
- **Taxonomy:** naive Bayesian classifier, SILVA v138.2, minimum bootstrap 50, exact-match species
- **Alpha diversity:** observed ASVs and Shannon index on rarefied data; Wilcoxon or Kruskal-Wallis
- **Beta diversity:** Bray-Curtis PCoA; PERMANOVA and homogeneity of dispersions (9,999 permutations)
- **Differential abundance:** CLR-transformed genus counts, Wilcoxon test, Benjamini-Hochberg FDR
- **Quality control:** per-step read tracking, length filtering, mock community validation

## Run it yourself

Requires R (≥ 4.3) with `dada2`, `phyloseq`, `vegan`, `ggplot2`, `dplyr`, `tidyr`, `patchwork`, `rmarkdown`.

```bash
bash run_all.sh          # downloads data + SILVA, runs everything, renders results/report.html
```

For a new dataset, put FASTQ files in `data/raw/`, write `data/metadata.tsv` (first column = sample ID), and edit `config.R`: truncation lengths, primer trimming, expected amplicon length, and the grouping column. No other code changes are needed.

## Repository layout

```
config.R                 all project settings in one place
run_all.sh               one-command pipeline
scripts/                 numbered analysis steps
report/report.Rmd        report template (numbers fill in automatically)
results/                 created on run: report, tables, figures, phyloseq object
```

## Contact

**Jeff Keeling** · PhD candidate, Biological Sciences, University of Texas at El Paso
Dissertation research: bacterial, archaeal, and fungal communities of Arctic plants using amplicon sequencing and metatranscriptomics.

Available for freelance amplicon analysis, figure preparation, and methods writing. Fixed-price quotes based on sample number and scope.
