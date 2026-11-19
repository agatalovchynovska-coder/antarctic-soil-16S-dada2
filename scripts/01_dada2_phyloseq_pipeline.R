# ==============================================================================
# 16S rRNA Amplicon Sequencing Pipeline (DADA2 -> Phyloseq)
# Project: Culturable and unculturable microbiota of Thala Hills Oasis, East Antarctica
# ==============================================================================

# 1. Load required libraries ---------------------------------------------------
library(dada2)
library(phyloseq)
library(decontam)
library(DECIPHER)
library(phangorn)
library(tidyverse)

# Set paths (relative paths for reproducibility)
raw_reads_dir <- "data/clean_reads"
tax_db_path   <- "data/silva_nr99_v138.2_toGenus_trainset.fa.gz"
results_dir   <- "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Match forward and reverse FASTQ files -------------------------------------
fnFs <- sort(list.files(raw_reads_dir, pattern = "_R1_001.16S.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(raw_reads_dir, pattern = "_R2_001.16S.fastq.gz", full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# 3. Quality filtering and trimming --------------------------------------------
# Note: Primers were already removed prior to this step (no trimLeft applied)
filtFs <- file.path(results_dir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(results_dir, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  truncLen = c(240, 200),
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

# 4. Error rate learning & ASV inference ---------------------------------------
set.seed(42)
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

# 5. Merge paired reads & construct sequence table -----------------------------
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
seqtab  <- makeSequenceTable(mergers)

# Filter amplicon lengths (V3-V4 target window: 395-450 bp)
seqtab.len <- seqtab[, nchar(colnames(seqtab)) %in% 395:450]

# 6. Chimera removal -----------------------------------------------------------
seqtab.nochim <- removeBimeraDenovo(
  seqtab.len,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)

# Track read retention across pipeline steps
getN <- function(x) sum(getUniques(x))
track <- cbind(
  out,
  sapply(dadaFs, getN),
  sapply(dadaRs, getN),
  sapply(mergers, getN),
  rowSums(seqtab.nochim)
)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.csv(track, file.path(results_dir, "read_tracking_summary.csv"))

# 7. Taxonomic assignment (SILVA v138.2) ---------------------------------------
taxa <- assignTaxonomy(
  seqtab.nochim,
  tax_db_path,
  minBoot = 80,
  multithread = TRUE
)

# 8. Build initial Phyloseq object & Decontamination ---------------------------
# Define sample metadata indicating negative controls
sample_df <- data.frame(
  SampleID = rownames(seqtab.nochim),
  Sample_or_Control = ifelse(grepl("^K", rownames(seqtab.nochim)), "Control", "True Sample"),
  row.names = rownames(seqtab.nochim)
)

ps_raw <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  sample_data(sample_df),
  tax_table(taxa)
)

# Identify and remove contaminants using prevalence method (decontam)
sample_data(ps_raw)$is.neg <- sample_data(ps_raw)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(ps_raw, method = "prevalence", neg = "is.neg", threshold = 0.5)
ps_noncontam  <- prune_taxa(!contamdf.prev$contaminant, ps_raw)

# Remove negative and positive control samples from downstream analysis
controls_to_remove <- c("K1", "K2_1", "K2_2", "Lsp1654_Psp1767")
ps_samples <- prune_samples(!(sample_names(ps_noncontam) %in% controls_to_remove), ps_noncontam)

# 9. Taxonomic filtering (remove non-bacterial & organelle sequences) ----------
ps_clean <- ps_samples %>%
  subset_taxa(
    Kingdom == "Bacteria" &
    !is.na(Phylum) &
    (is.na(Order)  | Order != "Chloroplast") &
    (is.na(Family) | Family != "Mitochondria")
  )

# Remove empty taxa and zero-read samples
ps_clean <- prune_taxa(taxa_sums(ps_clean) > 0, ps_clean)
ps_clean <- prune_samples(sample_sums(ps_clean) > 0, ps_clean)

# 10. Phylogenetic tree construction (Neighbor-Joining) ------------------------
seqs <- getSequences(otu_table(ps_clean))
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA, verbose = FALSE)

phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm          <- dist.ml(phang.align)
treeNJ      <- NJ(dm)

# Merge tree into final Phyloseq object and save
ps_final <- merge_phyloseq(ps_clean, phy_tree(treeNJ))
saveRDS(ps_final, file.path(results_dir, "phyloseq_final_soils.rds"))
