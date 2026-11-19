# antarctic-soil-16S-dada2
16S rRNA (V3-V4) amplicon sequencing analysis pipeline (DADA2 in R) for 51 East Antarctic soil samples (Thala Hills oasis)
# 16S rRNA Amplicon Sequencing Pipeline: East Antarctic Soil Microbiota

Bioinformatic workflow for processing paired-end Illumina 16S rRNA (V3–V4) amplicon sequencing data from 51 soil samples collected in the **Thala Hills Oasis (East Antarctica)**.

## Pipeline Overview
1. **Quality Control & Filtering (`dada2`)**: Filtering primer-trimmed paired-end reads (`maxEE = c(2, 2)`), error rate modeling, and ASV inference.
2. **Merging & Length Filtering**: Merging forward and reverse reads and selecting target V3–V4 amplicon lengths (`395–450 bp`).
3. **Chimera Removal**: Consensus-based chimera filtering (`removeBimeraDenovo`).
4. **Taxonomic Assignment**: Classification against the **SILVA v138.2** reference database with an 80% bootstrap confidence threshold (`minBoot = 80`).
5. **Decontamination (`decontam`)**: Prevalence-based contaminant identification using extraction and PCR negative controls, followed by removal of control samples.
6. **Taxonomic Curation & Phylogeny (`phyloseq`, `DECIPHER`, `phangorn`)**: Removal of chloroplast, mitochondrial, and unassigned phylum sequences, followed by multiple sequence alignment and Neighbor-Joining phylogenetic tree reconstruction.

## Dependencies
* R (>= 4.3.0)
* `dada2`, `phyloseq`, `decontam`, `DECIPHER`, `phangorn`, `tidyverse`

## Output
* `results/read_tracking_summary.csv` — read counts retained at each step of the DADA2 workflow.
* `results/phyloseq_final_soils.rds` — curated `phyloseq` object (ASV table, taxonomy, metadata, and phylogenetic tree) ready for alpha/beta diversity and statistical analyses.
