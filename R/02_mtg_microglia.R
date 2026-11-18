# 02_mtg_microglia.R — select MTG microglia, quality control, metacell check, CPS
# Before the first run: install.packages("readxl", type = "binary"); renv::snapshot()
# Output: results/02_mtg_microglia_output.txt and data/processed/mtg_microglia_meta.rds

library(Matrix)
library(readxl)
source("R/read_h5ad.R")

path <- "data/SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_RNAseq_final-nuclei.2026-06-22.h5ad"
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
sink("results/02_mtg_microglia_output.txt", split = TRUE)

md     <- read_h5ad_df(path, "obs")
genes  <- read_h5ad_df(path, "var")
counts <- read_h5ad_matrix(path, "UMIs", n_cells = nrow(md), n_genes = nrow(genes))

# ---------------------------------------------------------------- 1. Filtering
cat("\n######## 1. FILTERING ########\n")
f_region <- md$`Brain Region` == "MTG"
f_type   <- !md$Supertype %in% c("Monocyte", "Lymphocyte")   # blood-derived cells, not microglia
f_qc     <- md$experiment_component_failed != "TRUE" & md$rna_amplification_pass_fail == "Pass"
cat("All nuclei:                      ", nrow(md), "\n")
cat("  in MTG:                        ", sum(f_region), "\n")
cat("  + without monocytes/lymphocytes:", sum(f_region & f_type), "\n")
keep <- f_region & f_type & f_qc
cat("  + passed quality control:      ", sum(keep), "\n")

cat("\nSupertypes after filtering:\n"); print(table(droplevels(md$Supertype[keep])))
n_per_donor <- table(droplevels(md$`Donor ID`[keep]))
cat("\nDonors:", length(n_per_donor), "| nuclei per donor:\n"); print(summary(as.vector(n_per_donor)))

# ---------------------------------------------------------------- 2. What the supertypes are
cat("\n######## 2. MARKERS BY SUPERTYPE (fraction of nuclei in which the gene is detected) ########\n")
markers <- c("MRC1", "LYVE1", "F13A1", "CD163",      # perivascular macrophages (PVM)
             "P2RY12", "CX3CR1", "TMEM119",          # homeostatic microglia
             "APOE", "SPP1", "GPNMB", "TREM2")       # disease response
set.seed(1)
by_type <- split(which(keep), droplevels(md$Supertype[keep]))
samp <- sort(unlist(lapply(by_type, function(i) i[sample.int(length(i), min(300, length(i)))])))
sub <- as(counts[match(markers, rownames(genes)), samp], "dgCMatrix")
rownames(sub) <- markers
det_type <- sapply(split(seq_along(samp), droplevels(md$Supertype[samp])),
                   function(j) rowMeans(sub[, j, drop = FALSE] > 0))
print(round(det_type, 2))

# ---------------------------------------------------------------- 3. Raw counts and sparsity
cat("\n######## 3. 5000 RANDOM MTG MICROGLIA NUCLEI ########\n")
s5  <- sort(sample(which(keep), 5000))
umi <- as(counts[, s5], "dgCMatrix")
rownames(umi) <- rownames(genes)
cat("All values are integers (i.e. raw counts):", all(umi@x == round(umi@x)), "\n")

cand <- c("P2RY12", "CX3CR1", "TMEM119", "CSF1R", "TREM2", "TYROBP", "APOE", "SPP1",
          "GPNMB", "LPL", "CD9", "CD74", "PLCG2", "INPP5D", "MS4A6A", "SPI1", "IRF8")
det_nuc <- rowMeans(umi[cand, ] > 0)

# Trial metacells: within each donor, split nuclei into groups of 30 and sum their counts
don <- droplevels(md$`Donor ID`[s5])
grp <- ave(seq_along(s5), don, FUN = function(i) ceiling(seq_along(i) / 30))
key <- factor(paste(don, grp))
full <- levels(key)[table(key) == 30]                  # complete groups only
mc  <- umi[cand, ] %*% t(fac2sparse(key))
mc  <- mc[, full, drop = FALSE]
det_mc <- rowMeans(mc > 0)
cat("\nFraction with the gene detected: single nuclei vs trial metacells (", length(full), "metacells of 30 nuclei)\n")
print(round(cbind(nuclei = det_nuc, metacells = det_mc), 2))

# ---------------------------------------------------------------- 4. CPS and metadata
cat("\n######## 4. CPS (Travaglini 2026, local to MTG) ########\n")
cps <- read_excel("data/sea-ad_global_and_local_cps_050126.xlsx", sheet = "Travaglini 2026")
cps_mtg <- cps[cps$`Brain Region` == "MTG", ]

meta_cols <- c("Donor ID", "Supertype", "Number of UMIs", "Genes detected", "Fraction mitochondrial UMIs",
               "Doublet score", "Age at Death", "Sex", "APOE Genotype", "Overall AD neuropathological Change",
               "Braak", "Thal", "Cognitive Status", "Last MMSE Score", "Severely Affected Donor")
meta <- md[keep, meta_cols]
meta$cell_index <- which(keep)                          # nucleus index in the original file
m <- match(as.character(meta$`Donor ID`), cps_mtg$`Donor ID`)
for (cc in c("CPS_Local", "CPS_Local_ABeta", "CPS_Local_pTau", "CPS_Global")) meta[[cc]] <- cps_mtg[[cc]][m]
cat("Nuclei without CPS:", sum(is.na(meta$CPS_Local)), "\n")

# Donor-level sanity check: CPS should increase with ADNC and Braak stage
donors <- unique(meta[, c("Donor ID", "CPS_Local", "Overall AD neuropathological Change", "Braak",
                          "Cognitive Status", "Severely Affected Donor")])
cat("\nCPS_Local by ADNC group (donor level):\n")
print(aggregate(CPS_Local ~ `Overall AD neuropathological Change`, donors, function(x) round(c(n = length(x), median = median(x)), 2)))
cat("\nCPS_Local by Braak stage:\n")
print(aggregate(CPS_Local ~ Braak, donors, function(x) round(c(n = length(x), median = median(x)), 2)))
cat("\nCPS_Local in severely affected donors (Severely Affected Donor):\n")
print(aggregate(CPS_Local ~ `Severely Affected Donor`, donors, function(x) round(c(n = length(x), median = median(x)), 2)))

saveRDS(meta, "data/processed/mtg_microglia_meta.rds")
cat("\nSaved: data/processed/mtg_microglia_meta.rds (", nrow(meta), "nuclei )\n")
sink()
