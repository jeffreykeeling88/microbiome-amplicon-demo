# 04_diversity.R: alpha diversity, beta diversity, composition, differential abundance
suppressPackageStartupMessages({
  library(phyloseq); library(vegan); library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})
source("config.R")

# ---- Style ------------------------------------------------------------------
pal_cat <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
col_other <- "#b9b8b2"
theme_set(theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title.position = "plot",
        panel.grid.major = element_line(colour = "#ecebe7", linewidth = 0.4),
        axis.text  = element_text(colour = "#52514e"),
        axis.title = element_text(colour = "#0b0b0b"),
        plot.title = element_text(face = "bold", colour = "#0b0b0b"),
        plot.subtitle = element_text(colour = "#52514e"),
        legend.position = "top", legend.justification = "left",
        strip.text = element_text(face = "bold")))
save_fig <- function(p, name, w = 7, h = 4.5) {
  ggsave(file.path(cfg$out_dir, "figures", name), p, width = w, height = h, dpi = 200, bg = "white")
}
fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))

# ---- Data -------------------------------------------------------------------
ps_all <- readRDS(file.path(cfg$out_dir, "phyloseq.rds"))
ps <- prune_samples(!sample_names(ps_all) %in% cfg$mock_samples, ps_all)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
# Remove non-bacterial host signal: chloroplast and mitochondrial 16S, unassigned kingdom
tt <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
offtarget <- is.na(tt$Kingdom) | tt$Kingdom == "Eukaryota" |
  (!is.na(tt$Order) & tt$Order == "Chloroplast") | (!is.na(tt$Family) & tt$Family == "Mitochondria")
n_off_reads <- sum(taxa_sums(ps)[offtarget])
ps <- prune_taxa(!offtarget, ps)
message(sprintf("Removed %d off-target ASVs (chloroplast, mitochondria, unassigned): %.2f%% of reads",
                sum(offtarget), 100 * n_off_reads / (n_off_reads + sum(taxa_sums(ps)))))
offtarget_summary <- data.frame(asvs_removed = sum(offtarget),
                                pct_reads_removed = round(100 * n_off_reads / (n_off_reads + sum(taxa_sums(ps))), 2))
grp <- cfg$group_col
groups <- sort(unique(as.character(sample_data(ps)[[grp]])))
sample_data(ps)[[grp]] <- factor(sample_data(ps)[[grp]], levels = groups)
group_cols <- setNames(pal_cat[seq_along(groups)], groups)
two_groups <- length(groups) == 2

depth <- if (is.na(cfg$rarefy_depth)) min(sample_sums(ps)) else cfg$rarefy_depth
ps_rare <- rarefy_even_depth(ps, sample.size = depth, rngseed = cfg$seed, verbose = FALSE)
message(sprintf("%d samples, %d ASVs; rarefied to %d reads/sample for alpha diversity", nsamples(ps), ntaxa(ps), depth))

# ---- Alpha diversity --------------------------------------------------------
alpha <- estimate_richness(ps_rare, measures = c("Observed", "Shannon", "InvSimpson"))
alpha$sample <- rownames(alpha)
alpha[[grp]] <- sample_data(ps_rare)[[grp]]
write.table(alpha, file.path(cfg$out_dir, "tables", "alpha_diversity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

alpha_long <- alpha %>%
  pivot_longer(c(Observed, Shannon), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, Observed = "Observed ASVs", Shannon = "Shannon index"))
alpha_tests <- alpha_long %>% group_by(metric) %>% summarise(
  p = if (two_groups) wilcox.test(value ~ .data[[grp]], exact = FALSE)$p.value
      else kruskal.test(value ~ .data[[grp]])$p.value,
  test = if (two_groups) "Wilcoxon rank-sum" else "Kruskal-Wallis", .groups = "drop")
write.table(alpha_tests, file.path(cfg$out_dir, "tables", "alpha_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

p_alpha <- ggplot(alpha_long, aes(x = .data[[grp]], y = value)) +
  geom_boxplot(aes(colour = .data[[grp]]), width = 0.5, outlier.shape = NA, linewidth = 0.6, fill = NA) +
  geom_jitter(aes(fill = .data[[grp]]), width = 0.12, height = 0, shape = 21, size = 2.6,
              colour = "white", stroke = 0.6) +
  geom_text(data = alpha_tests, aes(x = -Inf, y = Inf, label = paste(test, fmt_p(p))),
            hjust = -0.05, vjust = 1.5, size = 3.2, colour = "#52514e", inherit.aes = FALSE) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_colour_manual(values = group_cols, guide = "none") +
  scale_fill_manual(values = group_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(title = "Alpha diversity", subtitle = sprintf("Rarefied to %s reads per sample", format(depth, big.mark = ",")),
       x = NULL, y = NULL)
save_fig(p_alpha, "alpha_diversity.png", w = 7.5, h = 4.2)

# ---- Beta diversity: Bray-Curtis PCoA + PERMANOVA ---------------------------
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
bc <- phyloseq::distance(ps_rel, method = "bray")
meta <- data.frame(sample_data(ps_rel))
perm <- adonis2(bc ~ meta[[grp]], permutations = 9999)
disp <- permutest(betadisper(bc, meta[[grp]]), permutations = 9999)
beta_tests <- data.frame(
  test = c("PERMANOVA (adonis2)", "Dispersion (betadisper)"),
  R2   = c(round(perm$R2[1], 3), NA),
  F    = c(round(perm$F[1], 2), round(disp$tab$F[1], 2)),
  p    = c(perm$`Pr(>F)`[1], disp$tab$`Pr(>F)`[1]))
write.table(beta_tests, file.path(cfg$out_dir, "tables", "beta_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

ord <- ordinate(ps_rel, method = "PCoA", distance = bc)
ev  <- ord$values$Relative_eig[1:2] * 100
pcoa <- data.frame(ord$vectors[, 1:2], meta)
p_beta <- ggplot(pcoa, aes(Axis.1, Axis.2)) +
  stat_ellipse(aes(colour = .data[[grp]]), type = "t", level = 0.95, linewidth = 0.6, linetype = "dashed") +
  geom_point(aes(fill = .data[[grp]]), shape = 21, size = 3.2, colour = "white", stroke = 0.7) +
  scale_colour_manual(values = group_cols, name = NULL) +
  scale_fill_manual(values = group_cols, name = NULL) +
  labs(title = "Community composition (Bray-Curtis PCoA)",
       subtitle = sprintf("PERMANOVA R2 = %.2f, %s; dispersion test %s",
                          perm$R2[1], fmt_p(perm$`Pr(>F)`[1]), fmt_p(disp$tab$`Pr(>F)`[1])),
       x = sprintf("PCoA 1 (%.1f%%)", ev[1]), y = sprintf("PCoA 2 (%.1f%%)", ev[2])) +
  coord_equal()
save_fig(p_beta, "beta_pcoa_bray.png", w = 6.5, h = 5.5)

# ---- Composition: top genera ------------------------------------------------
glom_rank <- function(physeq, rank) {
  tt <- as.data.frame(tax_table(physeq), stringsAsFactors = FALSE)
  lab <- tt[[rank]]
  # Label unclassified ASVs by their deepest known rank
  ranks <- colnames(tt)[seq_len(which(colnames(tt) == rank))]
  for (i in which(is.na(lab))) {
    known <- unlist(tt[i, ranks]); known <- known[!is.na(known)]
    lab[i] <- if (length(known)) paste0("Unclassified ", tail(known, 1)) else "Unclassified"
  }
  otu <- as(otu_table(physeq), "matrix"); if (taxa_are_rows(physeq)) otu <- t(otu)
  t(rowsum(t(otu), lab))   # samples x taxa
}
gen <- glom_rank(ps_rel, "Genus")
top_n <- 7
top <- names(sort(colMeans(gen), decreasing = TRUE))[seq_len(min(top_n, ncol(gen)))]
comp <- data.frame(sample = rownames(gen), gen[, top, drop = FALSE],
                   Other = 1 - rowSums(gen[, top, drop = FALSE]), check.names = FALSE) %>%
  pivot_longer(-sample, names_to = "Genus", values_to = "abundance") %>%
  left_join(data.frame(sample = rownames(meta), group = meta[[grp]]), by = "sample") %>%
  mutate(Genus = factor(Genus, levels = c(top, "Other")))
ord_samples <- comp %>% filter(Genus == top[1]) %>% arrange(group, desc(abundance)) %>% pull(sample)
comp$sample <- factor(comp$sample, levels = ord_samples)
genus_cols <- c(setNames(pal_cat[seq_along(top)], top), Other = col_other)

p_comp <- ggplot(comp, aes(sample, abundance, fill = Genus)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.25) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = genus_cols, name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(title = "Most abundant genera", subtitle = sprintf("Top %d genera by mean relative abundance; remainder pooled as Other", length(top)),
       x = NULL, y = "Relative abundance") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
        panel.grid.major.x = element_blank(), legend.position = "right")
save_fig(p_comp, "taxa_barplot_genus.png", w = 9, h = 4.8)

# ---- Differential abundance: CLR + Wilcoxon, BH-adjusted --------------------
# Genera present in >= 20% of samples; CLR transform with pseudocount.
gen_counts <- glom_rank(ps, "Genus")
prev <- colMeans(gen_counts > 0)
gc <- gen_counts[, prev >= 0.2, drop = FALSE]
clr <- t(apply(gc + 0.5, 1, function(x) log(x) - mean(log(x))))
g <- meta[rownames(clr), grp]
da <- lapply(colnames(clr), function(tx) {
  x <- clr[, tx]
  p <- if (two_groups) wilcox.test(x ~ g, exact = FALSE)$p.value else kruskal.test(x ~ g)$p.value
  diff <- if (two_groups) median(x[g == groups[2]]) - median(x[g == groups[1]]) else NA
  data.frame(genus = tx, prevalence = round(prev[tx], 2), clr_median_diff = round(diff, 2), p = p)
})
da <- do.call(rbind, da)
da$q <- p.adjust(da$p, method = "BH")
da <- da[order(da$q, da$p), ]
write.table(da, file.path(cfg$out_dir, "tables", "differential_abundance_genus.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

if (two_groups) {
  sig <- head(da[da$q < 0.05, ], 15)
  if (nrow(sig) > 0) {
    sig$direction <- ifelse(sig$clr_median_diff > 0, paste("Higher in", groups[2]), paste("Higher in", groups[1]))
    dir_cols <- setNames(group_cols[c(2, 1)], paste("Higher in", groups[c(2, 1)]))
    p_da <- ggplot(sig, aes(clr_median_diff, reorder(genus, clr_median_diff))) +
      geom_vline(xintercept = 0, colour = "#52514e", linewidth = 0.4) +
      geom_segment(aes(x = 0, xend = clr_median_diff, yend = reorder(genus, clr_median_diff)),
                   colour = "#c3c2b7", linewidth = 0.6) +
      geom_point(aes(fill = direction), shape = 21, size = 3.2, colour = "white", stroke = 0.7) +
      scale_fill_manual(values = dir_cols, name = NULL) +
      labs(title = "Differentially abundant genera",
           subtitle = sprintf("%d of %d genera significant (Wilcoxon on CLR counts, BH q < 0.05)%s",
                              sum(da$q < 0.05), nrow(da), if (sum(da$q < 0.05) > 15) "; top 15 shown" else ""),
           x = sprintf("Median CLR difference (%s minus %s)", groups[2], groups[1]), y = NULL)
    save_fig(p_da, "differential_abundance.png", w = 7.5, h = 0.6 + 0.32 * nrow(sig) + 1.4)
  }
}

# ---- Read tracking figure ---------------------------------------------------
track <- read.delim(file.path(cfg$out_dir, "tables", "read_tracking.tsv"))
steps <- c("input", "filtered", "merged", "nonchim")
tl <- track %>% select(sample, all_of(steps)) %>% pivot_longer(-sample, names_to = "step", values_to = "reads") %>%
  mutate(step = factor(step, levels = steps, labels = c("Raw", "Quality-filtered", "Merged", "Chimera-free")))
p_track <- ggplot(tl, aes(step, reads, group = sample)) +
  geom_line(colour = "#c3c2b7", linewidth = 0.5) +
  geom_point(colour = pal_cat[1], size = 1.6) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Reads retained at each step", subtitle = sprintf("Median %.0f%% of raw reads kept through chimera removal",
       median(track$pct_retained, na.rm = TRUE)), x = NULL, y = "Reads per sample")
save_fig(p_track, "read_tracking.png", w = 6.5, h = 4)

saveRDS(list(offtarget = offtarget_summary, alpha_tests = alpha_tests, beta_tests = beta_tests, da = da, depth = depth,
             groups = groups, top = top), file.path(cfg$out_dir, "summary.rds"))
message("Diversity analysis complete.")
