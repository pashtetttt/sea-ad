# 03_build_dataset.R — MTG microglia metacells and the final table for the networks (~1000 x ~50)
#
# What it does:
#   1. Takes the MTG microglia nuclei from step 02 (without Micro-PVM_1, which are perivascular macrophages).
#   2. Within each "donor x supertype" pair, randomly splits nuclei into groups of 30
#      and SUMS their raw counts — these are the metacells (a remainder of < 30 nuclei is dropped).
#   3. Caps the number of metacells per donor so that the total is ~1000 and donors with many
#      nuclei do not dominate.
#   4. Selects 48 genes: some from biology, some from the data (most variable).
#      NOTE: this gene selection is superseded by 04_select_genes.R (the final gene set).
#   5. Adds the CPS_Local and total_UMIs nodes (technical control for sequencing depth).
#
# Output:
#   data/processed/network_dataset.csv / .tsv    — table for ARACNE, PC, MIIC (rows = metacells)
#   data/processed/metacell_info.csv             — donor, supertype and covariates of each metacell
#   data/processed/metacell_counts_all_genes.rds — all genes, to change the gene set without rebuilding
#   results/03_build_dataset_output.txt          — report

library(Matrix)
source("R/read_h5ad.R")
set.seed(2026)

path         <- "data/SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_RNAseq_final-nuclei.2026-06-22.h5ad"
MC_SIZE      <- 30     # nuclei per metacell
TARGET_N     <- 1000   # desired number of metacells
N_GENES      <- 48     # genes in the network (+ CPS_Local + total_UMIs = 50 variables)
MIN_DETECT   <- 0.9    # a gene must be non-zero in at least 90% of metacells

sink("results/03_build_dataset_output.txt", split = TRUE)

md     <- read_h5ad_df(path, "obs")
genes  <- read_h5ad_df(path, "var")
counts <- read_h5ad_matrix(path, "UMIs", n_cells = nrow(md), n_genes = nrow(genes))

meta <- readRDS("data/processed/mtg_microglia_meta.rds")
meta <- meta[meta$Supertype != "Micro-PVM_1", ]
meta$Supertype  <- droplevels(meta$Supertype)
meta$donor      <- droplevels(meta[["Donor ID"]])
cat("MTG microglia nuclei (without Micro-PVM_1):", nrow(meta), "| donors:", nlevels(meta$donor), "\n")

# ---------------------------------------------------------------- 1. Split into metacells
cat("\n######## 1. METACELLS ########\n")
meta$group <- NA_character_
cells_by_ds <- split(seq_len(nrow(meta)), interaction(meta$donor, meta$Supertype, drop = TRUE))
for (nm in names(cells_by_ds)) {
  rows <- cells_by_ds[[nm]]
  k <- floor(length(rows) / MC_SIZE)
  if (k < 1) next
  rows <- rows[sample.int(length(rows))][seq_len(k * MC_SIZE)]   # shuffle and take k*30 nuclei
  meta$group[rows] <- paste0(nm, "_mc", rep(seq_len(k), each = MC_SIZE))
}
mc_info <- unique(meta[!is.na(meta$group), c("group", "donor", "Supertype")])
cat("Metacells available before capping:", nrow(mc_info), "\n")

# Per-donor cap: the smallest cap for which the total is >= TARGET_N
avail <- table(mc_info$donor)
caps  <- seq_len(max(avail))
tot   <- sapply(caps, function(cp) sum(pmin(avail, cp)))
cap   <- if (any(tot >= TARGET_N)) caps[which(tot >= TARGET_N)[1]] else max(avail)
keep_mc <- unlist(lapply(split(mc_info$group, mc_info$donor), function(g)
  if (length(g) > cap) g[sample.int(length(g), cap)] else g))
mc_info <- mc_info[mc_info$group %in% keep_mc, ]
cat("Cap: at most", cap, "metacells per donor -> total", nrow(mc_info), "metacells\n")
cat("Donors with at least one metacell:", length(unique(mc_info$donor)), "\n")
cat("\nMetacells per donor:\n"); print(summary(as.vector(table(droplevels(mc_info$donor)))))
cat("\nMetacells per supertype:\n"); print(table(droplevels(mc_info$Supertype)))

# ---------------------------------------------------------------- 2. Summing counts
cat("\n######## 2. SUMMING COUNTS (donor by donor) ########\n")
mc_list <- list()
donors  <- unique(as.character(mc_info$donor))
for (i in seq_along(donors)) {
  rows <- which(meta$group %in% mc_info$group & meta$donor == donors[i])
  o    <- order(meta$cell_index[rows])
  X    <- as(counts[, meta$cell_index[rows][o]], "dgCMatrix")       # genes x nuclei of this donor
  g    <- factor(meta$group[rows][o])
  mc_list[[i]] <- X %*% t(fac2sparse(g))                            # genes x metacells (sums)
  if (i %% 10 == 0) cat("  donors processed:", i, "of", length(donors), "\n")
}
MC <- do.call(cbind, mc_list)
rownames(MC) <- rownames(genes)
mc_info <- mc_info[match(colnames(MC), mc_info$group), ]
mc_info$total_UMIs <- colSums(MC)
cat("Metacell matrix (genes x metacells):", dim(MC), "\n")
cat("All sums are integers:", all(MC@x == round(MC@x)), "\n")
cat("\nUMIs (molecules) per metacell:\n"); print(summary(mc_info$total_UMIs))

# Donor metadata for each metacell
dmeta <- unique(meta[, c("donor", "CPS_Local", "CPS_Local_ABeta", "CPS_Local_pTau", "Age at Death", "Sex",
                         "APOE Genotype", "Overall AD neuropathological Change", "Braak",
                         "Cognitive Status", "Severely Affected Donor")])
mc_info <- cbind(mc_info, dmeta[match(mc_info$donor, dmeta$donor), -1])

# ---------------------------------------------------------------- 3. Gene selection
cat("\n######## 3. GENE SELECTION ########\n")
det <- rowMeans(MC > 0)

bio <- c("P2RY12", "CX3CR1", "TMEM119", "CSF1R",                     # homeostatic microglia
         "TREM2", "TYROBP", "APOE", "SPP1", "GPNMB", "CD74", "HLA-DRA", # disease response
         "C1QA", "C1QB", "CD163",                                     # complement, macrophage-like profile
         "PLCG2", "INPP5D", "MS4A6A", "PICALM", "SORL1", "BIN1", "CD33", "ABI3",  # AD risk genes
         "SPI1", "IRF8", "MEF2C")                                     # transcription factors
cat("Biology candidates missing from the data:", setdiff(bio, rownames(MC)), "\n")
bio <- intersect(bio, rownames(MC))
cat("Dropped as too sparse (<", MIN_DETECT, "of metacells):\n")
print(round(det[bio][det[bio] < MIN_DETECT], 2))
bio_ok <- bio[det[bio] >= MIN_DETECT]

# Data-driven: most variable genes among well-detected ones.
# Normalization (counts per million + log) is used ONLY to rank genes;
# the network gets raw summed counts.
excl <- c("MALAT1", "NEAT1", "XIST", "TSIX", "RPS4Y1", "DDX3Y", "UTY", "KDM5D", "USP9Y",
          "EIF1AY", "ZFY", "NLGN4Y", "TXLNGY")                       # nuclear lncRNAs and sex genes
pool <- names(det)[det >= MIN_DETECT & !grepl("^MT-|^RPL|^RPS|^HB[AB]", names(det)) &
                     !names(det) %in% c(excl, bio_ok)]
lc  <- log1p(as.matrix(MC[pool, ] %*% Diagonal(x = 1e6 / mc_info$total_UMIs)))
mu  <- rowMeans(lc); va <- apply(lc, 1, var)
res <- va - predict(loess(va ~ mu))                                  # variability above what is expected
n_data     <- N_GENES - length(bio_ok)
data_genes <- names(sort(res, decreasing = TRUE))[seq_len(n_data)]
cat("\nPool for data-driven selection:", length(pool), "genes; taking the", n_data, "most variable\n")

sel <- c(bio_ok, data_genes)

# ---------------------------------------------------------------- 4. Final table
dataset <- data.frame(CPS_Local = mc_info$CPS_Local,
                      total_UMIs = mc_info$total_UMIs,
                      as.matrix(t(MC[sel, ])), check.names = FALSE)
rownames(dataset) <- mc_info$group

cat("\n######## 4. FINAL TABLE ########\n")
cat("Size (metacells x variables):", dim(dataset), "\n")

rep_tab <- data.frame(
  gene        = sel,
  source      = ifelse(sel %in% bio_ok, "biology", "data"),
  detect      = round(det[sel], 3),
  median_umi  = apply(as.matrix(MC[sel, ]), 1, median),
  rho_CPS     = round(sapply(sel, function(g) cor(dataset[[g]], dataset$CPS_Local, method = "spearman")), 2),
  rho_depth   = round(sapply(sel, function(g) cor(dataset[[g]], dataset$total_UMIs, method = "spearman")), 2),
  row.names = NULL)
cat("\nSelected genes (rho = Spearman rank correlation, for inspection only):\n")
print(rep_tab, row.names = FALSE)

cat("\nCPS_Local across metacells:\n"); print(summary(dataset$CPS_Local))
cat("Distinct CPS values (= donors):", length(unique(dataset$CPS_Local)), "\n")

write.csv(dataset, "data/processed/network_dataset.csv", row.names = FALSE)
write.table(dataset, "data/processed/network_dataset.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(mc_info, "data/processed/metacell_info.csv", row.names = FALSE)
saveRDS(MC, "data/processed/metacell_counts_all_genes.rds")
write.csv(rep_tab, "results/03_selected_genes.csv", row.names = FALSE)
cat("\nSaved to data/processed/ and results/\n")
sink()
