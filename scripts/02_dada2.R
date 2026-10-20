# 02_dada2.R: filter, denoise, merge, remove chimeras, assign taxonomy
# Output: ASV table, taxonomy, ASV sequences, read-tracking table, phyloseq object
suppressPackageStartupMessages({ library(dada2); library(phyloseq); library(Biostrings) })
source("config.R")

fnFs <- sort(list.files(cfg$fastq_dir, pattern = cfg$fwd_pattern, full.names = TRUE))
fnRs <- sort(list.files(cfg$fastq_dir, pattern = cfg$rev_pattern, full.names = TRUE))
samples <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
names(fnFs) <- names(fnRs) <- samples

# ---- 1. Quality filtering ---------------------------------------------------
filt_dir <- file.path(cfg$out_dir, "filtered")
filtFs <- file.path(filt_dir, paste0(samples, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(samples, "_R_filt.fastq.gz"))
names(filtFs) <- names(filtRs) <- samples

filt <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                      trimLeft = cfg$trim_left, truncLen = cfg$trunc_len,
                      maxN = 0, maxEE = cfg$max_ee, truncQ = 2, rm.phix = TRUE,
                      compress = TRUE, multithread = cfg$threads)
rownames(filt) <- samples
keep <- filt[, "reads.out"] > 0
if (any(!keep)) message("Dropping samples with no reads after filtering: ", paste(samples[!keep], collapse = ", "))
filtFs <- filtFs[keep]; filtRs <- filtRs[keep]

# ---- 2. Learn error models and denoise -------------------------------------
errF <- learnErrors(filtFs, multithread = cfg$threads, randomize = TRUE)
errR <- learnErrors(filtRs, multithread = cfg$threads, randomize = TRUE)
ggplot2::ggsave(file.path(cfg$out_dir, "figures", "qc_error_model_forward.png"),
                plotErrors(errF, nominalQ = TRUE), width = 8, height = 6, dpi = 150)

dadaFs <- dada(filtFs, err = errF, multithread = cfg$threads)
dadaRs <- dada(filtRs, err = errR, multithread = cfg$threads)

# ---- 3. Merge pairs, build table, remove chimeras ---------------------------
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = cfg$min_overlap)
seqtab  <- makeSequenceTable(mergers)
len_dist <- table(nchar(getSequences(seqtab)))

in_range <- nchar(colnames(seqtab)) >= cfg$expected_len[1] & nchar(colnames(seqtab)) <= cfg$expected_len[2]
seqtab <- seqtab[, in_range, drop = FALSE]
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = cfg$threads)

message(sprintf("%d ASVs after chimera removal; %.1f%% of merged reads retained",
                ncol(seqtab_nochim), 100 * sum(seqtab_nochim) / sum(seqtab)))

# ---- 4. Read tracking -------------------------------------------------------
getN <- function(x) sum(getUniques(x))
s <- names(filtFs)
track <- data.frame(
  sample     = samples,
  input      = filt[samples, "reads.in"],
  filtered   = filt[samples, "reads.out"],
  denoisedF  = NA_integer_, denoisedR = NA_integer_, merged = NA_integer_,
  length_ok  = NA_integer_, nonchim = NA_integer_, row.names = samples
)
track[s, "denoisedF"] <- sapply(dadaFs[s], getN)
track[s, "denoisedR"] <- sapply(dadaRs[s], getN)
track[s, "merged"]    <- sapply(mergers[s], getN)
track[s, "length_ok"] <- rowSums(seqtab)[s]
track[s, "nonchim"]   <- rowSums(seqtab_nochim)[s]
track$pct_retained <- round(100 * track$nonchim / track$input, 1)
write.table(track, file.path(cfg$out_dir, "tables", "read_tracking.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(length = names(len_dist), n_asvs = as.integer(len_dist)),
            file.path(cfg$out_dir, "tables", "asv_length_distribution.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 5. Taxonomy ------------------------------------------------------------
taxa <- assignTaxonomy(seqtab_nochim, cfg$silva_train, multithread = cfg$threads, minBoot = 50)
if (!is.na(cfg$silva_species) && file.exists(cfg$silva_species)) {
  taxa <- addSpecies(taxa, cfg$silva_species)
}

# ---- 6. Export --------------------------------------------------------------
asv_ids <- sprintf("ASV%04d", seq_len(ncol(seqtab_nochim)))
seqs    <- colnames(seqtab_nochim)
names(seqs) <- asv_ids

otu <- seqtab_nochim; colnames(otu) <- asv_ids
tax <- taxa;          rownames(tax) <- asv_ids

write.table(data.frame(ASV = asv_ids, t(otu), check.names = FALSE),
            file.path(cfg$out_dir, "tables", "asv_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(ASV = asv_ids, tax, sequence = seqs, check.names = FALSE),
            file.path(cfg$out_dir, "tables", "asv_taxonomy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeXStringSet(DNAStringSet(seqs), file.path(cfg$out_dir, "asv_sequences.fasta"))

meta <- read.delim(cfg$metadata, row.names = 1, check.names = FALSE)
meta <- meta[rownames(otu), , drop = FALSE]
ps <- phyloseq(otu_table(otu, taxa_are_rows = FALSE), tax_table(tax),
               sample_data(meta), refseq(DNAStringSet(seqs)))
saveRDS(ps, file.path(cfg$out_dir, "phyloseq.rds"))
message("Done: ", nsamples(ps), " samples, ", ntaxa(ps), " ASVs saved to ", cfg$out_dir)
