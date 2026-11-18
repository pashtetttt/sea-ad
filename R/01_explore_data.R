# 01_explore_data.R — first look at the SEA-AD data (immune nuclei, 10 brain regions)
# Run: open alzheimers-project.Rproj in RStudio, open this file, click Source.
# All output is also written to results/01_explore_output.txt.

library(Matrix)
source("R/read_h5ad.R")

path <- "data/SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_RNAseq_final-nuclei.2026-06-22.h5ad"

dir.create("results", showWarnings = FALSE)
sink("results/01_explore_output.txt", split = TRUE)   # output to both console and file

cat("\n######## 0. FILE CONTENTS (top levels) ########\n")
print(h5ls(path, recursive = 2, datasetinfo = FALSE)[, c("group", "name", "otype")])

cat("\n######## 1. NUCLEUS (obs) AND GENE (var) METADATA ########\n")
md    <- read_h5ad_df(path, "obs")
genes <- read_h5ad_df(path, "var")
cat("Nuclei:", nrow(md), "| genes:", nrow(genes), "\n")
cat("\nFirst gene names:\n"); print(head(rownames(genes), 10))
cat("\nGene annotation columns:\n"); print(colnames(genes))
cat("\nNucleus metadata columns:\n"); print(colnames(md))

cat("\nColumns that look like CPS / pathology / donor / region / cell type:\n")
print(grep("pseudo|cps|braak|thal|adnc|apoe|age|sex|donor|region|class|type",
           colnames(md), ignore.case = TRUE, value = TRUE))

# For each column: a frequency table if it has few distinct values, otherwise a short summary
for (col in colnames(md)) {
  x <- md[[col]]
  n <- length(unique(x))
  cat("\n====", col, "| type:", class(x)[1], "| distinct values:", n, "\n")
  if (n <= 40) print(table(x, useNA = "ifany")) else if (is.numeric(x)) print(summary(x)) else print(head(unique(x), 5))
}

cat("\n######## 2. RAW COUNTS IN 5000 RANDOM NUCLEI ########\n")
counts <- read_h5ad_matrix(path, "UMIs", n_cells = nrow(md), n_genes = nrow(genes))
cat("Matrix size (genes x nuclei):", dim(counts), "\n")

set.seed(42)
cells <- sort(sample(ncol(counts), 5000))
umi <- as(counts[, cells], "dgCMatrix")
rownames(umi) <- rownames(genes)

cat("Example values (should be integers if these are raw counts):\n")
print(as.matrix(umi[1:5, 1:5]))

cat("\nGenes detected per nucleus:\n");     print(summary(colSums(umi > 0)))
cat("\nUMIs (molecules) per nucleus:\n");    print(summary(colSums(umi)))

detect <- rowMeans(umi > 0)   # fraction of nuclei in which the gene is seen at least once
cat("\nNumber of genes detected in at least this fraction of nuclei:\n")
print(sapply(c(any = 0, `10%` = 0.1, `30%` = 0.3, `50%` = 0.5), function(t) sum(detect > t)))

cat("\n30 most frequently detected genes:\n")
print(round(head(sort(detect, decreasing = TRUE), 30), 2))

cand <- c("P2RY12", "CX3CR1", "TMEM119", "CSF1R",                    # homeostatic microglia
          "TREM2", "TYROBP", "APOE", "SPP1", "GPNMB", "LPL", "CD9", "CD74",  # disease response
          "PLCG2", "INPP5D", "MS4A6A",                                # AD risk genes
          "SPI1", "IRF8")                                             # transcription factors
cat("\nFraction of nuclei in which candidate genes are detected:\n")
print(round(detect[intersect(cand, names(detect))], 2))
cat("\nCandidates not found among gene names (maybe different IDs):\n")
print(setdiff(cand, names(detect)))

sink()
cat("\nDone. Output: results/01_explore_output.txt\n")
