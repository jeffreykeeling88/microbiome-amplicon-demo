# 01_qc.R: read quality profiles, used to choose truncation lengths
suppressPackageStartupMessages({ library(dada2); library(ggplot2) })
source("config.R")

fnFs <- sort(list.files(cfg$fastq_dir, pattern = cfg$fwd_pattern, full.names = TRUE))
fnRs <- sort(list.files(cfg$fastq_dir, pattern = cfg$rev_pattern, full.names = TRUE))
stopifnot(length(fnFs) > 0, length(fnFs) == length(fnRs))
message(length(fnFs), " paired samples found")

# Aggregate profile across all samples (one panel per read direction)
pF <- plotQualityProfile(fnFs, aggregate = TRUE) + ggtitle("Forward reads (all samples)")
pR <- plotQualityProfile(fnRs, aggregate = TRUE) + ggtitle("Reverse reads (all samples)")

ggsave(file.path(cfg$out_dir, "figures", "qc_quality_forward.png"), pF, width = 7, height = 4.5, dpi = 200)
ggsave(file.path(cfg$out_dir, "figures", "qc_quality_reverse.png"), pR, width = 7, height = 4.5, dpi = 200)
message("Quality profiles written to ", file.path(cfg$out_dir, "figures"))
