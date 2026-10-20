# ---------------------------------------------------------------
# Project configuration
# Edit this file for a new dataset. Every script reads it.
# ---------------------------------------------------------------

cfg <- list(
  # Input
  fastq_dir      = "data/raw",                 # folder of paired FASTQ files
  fwd_pattern    = "_R1_001.fastq",            # forward-read filename suffix
  rev_pattern    = "_R2_001.fastq",            # reverse-read filename suffix
  metadata       = "data/metadata.tsv",        # sample metadata, first column = sample ID

  # Reference database (DADA2-formatted SILVA training set)
  silva_train    = "ref/silva_nr99_v138.2_toGenus_trainset.fa.gz",
  silva_species  = "ref/silva_v138.2_assignSpecies.fa.gz",   # optional, set NA to skip

  # Read trimming / filtering (set from the quality profiles in 01_qc)
  trim_left      = c(0, 0),     # bases removed from start (primers, if still present)
  trunc_len      = c(240, 160), # forward, reverse truncation lengths
  max_ee         = c(2, 2),     # max expected errors per read
  min_overlap    = 12,          # minimum overlap when merging pairs
  expected_len   = c(250, 256), # expected amplicon length range (V4 = ~253 bp)

  # Analysis
  group_col      = "Group",     # metadata column used for comparisons
  mock_samples   = c("Mock"),   # mock community sample IDs (excluded from ecology stats)
  mock_reference = "data/HMP_MOCK.v35.fasta",  # known mock sequences, NA if none
  rarefy_depth   = NA,          # NA = use minimum sample depth
  seed           = 100,

  # Output
  out_dir        = "results",
  threads        = TRUE
)

dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(cfg$out_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(cfg$out_dir, "tables"),  showWarnings = FALSE)
set.seed(cfg$seed)
