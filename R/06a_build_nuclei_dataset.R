# 06a_build_nuclei_dataset.R — single-nucleus version of the network table (experiment B: no metacells)
#
# Same 50 variables as the metacell table (CPS_Local, total_UMIs and the 48 genes from step 04),
# but each row is ONE nucleus with its raw UMI counts instead of a sum over 30 nuclei.
# To keep donors balanced, we sample the same number of nuclei from every donor (~12 -> ~1000 rows).
#
# Needs: the .h5ad file, data/processed/mtg_microglia_meta.rds (step 02), results/04_selected_genes.csv (step 04)
# Output: data/processed/network_dataset_nuclei.csv / .tsv — table for ARACNE, PC, MIIC (rows = nuclei)
#         data/processed/nucleus_info.csv                — donor, supertype and covariates of each nucleus
#         results/06a_build_nuclei_dataset_output.txt   — report

library(Matrix)
source("R/read_h5ad.R")
set.seed(2026)

path        <- "data/SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_RNAseq_final-nuclei.2026-06-22.h5ad"
N_PER_DONOR <- 12      # nuclei sampled per donor (83 donors x 12 ~ 1000 rows)

sink("results/06a_build_nuclei_dataset_output.txt", split = TRUE)

md     <- read_h5ad_df(path, "obs")
genes  <- read_h5ad_df(path, "var")
counts <- read_h5ad_matrix(path, "UMIs", n_cells = nrow(md), n_genes = nrow(genes))

# Same nuclei pool as the metacells: MTG microglia after QC, without Micro-PVM_1 (perivascular macrophages)
meta <- readRDS("data/processed/mtg_microglia_meta.rds")
meta <- meta[meta$Supertype != "Micro-PVM_1", ]
meta$donor <- droplevels(meta[["Donor ID"]])
sel_genes <- read.csv("results/04_selected_genes.csv")$gene      # the same 48 genes as the metacell table
stopifnot(length(sel_genes) == 48, all(sel_genes %in% rownames(genes)))

# Balanced sample: N_PER_DONOR random nuclei from every donor (all nuclei if a donor has fewer)
pick <- unlist(lapply(split(seq_len(nrow(meta)), meta$donor), function(r)
  r[sample.int(length(r), min(N_PER_DONOR, length(r)))]))
nuc <- meta[pick, ]
nuc <- nuc[order(nuc$cell_index), ]
cat("Nuclei sampled:", nrow(nuc), "| donors:", length(unique(nuc$donor)), "\n")
cat("\nNuclei per supertype:\n"); print(table(droplevels(nuc$Supertype)))

X <- as(counts[, nuc$cell_index], "dgCMatrix")                   # genes x nuclei, raw counts
rownames(X) <- rownames(genes)
cat("\nAll values are integers:", all(X@x == round(X@x)), "\n")

# total_UMIs is defined exactly as for metacells: sum over ALL genes of the row
dataset <- data.frame(CPS_Local  = nuc$CPS_Local,
                      total_UMIs = colSums(X),
                      as.matrix(t(X[sel_genes, ])), check.names = FALSE)

cat("\nTable size (nuclei x variables):", dim(dataset), "\n")
cat("\nUMIs per nucleus:\n"); print(summary(dataset$total_UMIs))
cat("\nFraction of nuclei where each gene is non-zero (metacells had >= 0.8 for all):\n")
print(round(sort(colMeans(dataset[, sel_genes] > 0)), 2))
cat("\nDistinct CPS values (= donors):", length(unique(dataset$CPS_Local)), "\n")

write.csv(dataset, "data/processed/network_dataset_nuclei.csv", row.names = FALSE)
write.table(dataset, "data/processed/network_dataset_nuclei.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
info_cols <- c("cell_index", "donor", "Supertype", "CPS_Local", "CPS_Local_ABeta", "CPS_Local_pTau",
               "Age at Death", "Sex", "APOE Genotype", "Overall AD neuropathological Change", "Braak",
               "Cognitive Status", "Severely Affected Donor")
write.csv(nuc[, info_cols], "data/processed/nucleus_info.csv", row.names = FALSE)
cat("\nSaved: data/processed/network_dataset_nuclei.csv/.tsv and nucleus_info.csv\n")
sink()
