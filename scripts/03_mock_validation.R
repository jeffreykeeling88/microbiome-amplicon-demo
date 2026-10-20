# 03_mock_validation.R: check pipeline accuracy against a known mock community
# Every ASV found in the mock sample is compared to the known reference sequences.
suppressPackageStartupMessages({ library(phyloseq); library(Biostrings) })
source("config.R")

ps <- readRDS(file.path(cfg$out_dir, "phyloseq.rds"))
mocks <- intersect(cfg$mock_samples, sample_names(ps))
if (length(mocks) == 0 || is.na(cfg$mock_reference) || !file.exists(cfg$mock_reference)) {
  message("No mock community configured; skipping.")
  quit(save = "no")
}

ref  <- as.character(readDNAStringSet(cfg$mock_reference))
seqs <- as.character(refseq(ps))
res <- lapply(mocks, function(m) {
  counts <- as.numeric(otu_table(ps)[m, ])
  names(counts) <- taxa_names(ps)
  present <- counts[counts > 0]
  hit <- sapply(seqs[names(present)], function(s) any(grepl(s, ref, fixed = TRUE)))
  data.frame(
    sample              = m,
    asvs_detected       = length(present),
    asvs_exact_match    = sum(hit),
    asvs_unexpected     = sum(!hit),
    reference_sequences = length(ref),
    pct_reads_exact     = round(100 * sum(present[hit]) / sum(present), 2)
  )
})
res <- do.call(rbind, res)
print(res)
write.table(res, file.path(cfg$out_dir, "tables", "mock_validation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
