# 08_gene_sets.R — two alternative gene sets for the network, each on metacells and on single nuclei
#
#   Set "lit"  (literature only): microglial genes reported in Gabitto et al. 2024 (Nat Neurosci, Fig. 6c-e)
#              and Adeoye et al. 2025 (Aging Dis, Fig. 8 microglial modules). No data-driven choice.
#   Set "data" (data only):       genes with the largest variability BETWEEN METACELLS OF THE SAME DONOR.
#              Using within-donor variability targets cell-state programmes and avoids genes that only
#              differ between people (genotype, post-mortem handling) — the problem we saw in step 03.
#
# Both sets: genes must be non-zero in >= MIN_DETECT of metacells; CPS_Local and total_UMIs are added.
# Normalization is used only to rank genes; the tables contain raw summed (or single-nucleus) counts.
#
# Output: data/processed/network_dataset_{lit,data}.csv/.tsv          (metacells)
#         data/processed/network_dataset_nuclei_{lit,data}.csv/.tsv   (12 nuclei per donor, same sample as 06b rep 1)
#         results/08_gene_sets_output.txt, results/08_gene_sets.csv
#         results/networks_{lit,data}[_nuclei]/                       (if RUN_NETWORKS)
# Run time: ~10 min to read the nuclei from the .h5ad (cached), plus ~10 min per network run.

library(Matrix)
set.seed(2026)

N_DATA_GENES <- 60      # size of the data-only set (try 48 / 60 / 80)
MIN_DETECT   <- 0.8
N_PER_DONOR  <- 12
RUN_NETWORKS <- TRUE    # run 05_networks.R on all four tables at the end

sink("results/08_gene_sets_output.txt", split = TRUE)
MC      <- readRDS("data/processed/metacell_counts_all_genes.rds")
mc_info <- read.csv("data/processed/metacell_info.csv", check.names = FALSE)
stopifnot(identical(colnames(MC), mc_info$group))
det <- rowMeans(MC > 0)

# ======================================================================= Set "lit"
lit <- rbind(
  data.frame(source = "Gabitto Fig6c: DAM (Micro-PVM_3) markers",
             gene = c("PTPRG","TRHDE","CADM1","TPRG1","WIPF3","PARVG","RALGAPA2","KCNJ5","MGLL","FPR3","ZNF804A",
                      "DSCAM","SAMD4A","LIPA","NHSL1","RASGRP3","TBC1D16","TANC2","ATG7","DOP1B","SNX29","GLDN",
                      "LINC02698","AC131956.2","SLC47A1","AL138756.1")),
  data.frame(source = "Gabitto Fig6d: up early along CPS",
             gene = c("IL1B","CSF1R","STAB1","NINJ1","JAK3","IRF1","IRF7","IFI16","FCGR1A","FCGR1B","FCGR2A","FCGR3B",
                      "CD74","HLA-DRB5","C1QA","C1QB","CTSC","LY86","FCGR3A")),
  data.frame(source = "Gabitto Fig6d: up late along CPS", gene = c("CTSD","CTSS","LYZ","APOE")),
  data.frame(source = "Gabitto Fig6e: microglial TFs up early", gene = c("RUNX1","IKZF1","NFATC2","MAF")),
  data.frame(source = "Adeoye Fig8A: MTG MIC-M2 (up in AD)", gene = c("APOE","B2M","CST3","CD81")),
  data.frame(source = "Adeoye Fig8A: MTG MIC-M5 (lipid/immune)", gene = c("GPNMB","LRRK2","MITF","ABCA1","STARD13","ZBTB16","PRKAG2")),
  data.frame(source = "Adeoye Fig8A: MTG MIC-M4 (activation)", gene = c("PICALM","DOCK2","TAB2")),
  data.frame(source = "Adeoye Fig8A: MTG MIC-M1 (down in AD)", gene = c("NLGN1","RYR2","PLCB1")),
  data.frame(source = "Adeoye: hubs shared across regions", gene = c("CHST11","FKBP5","GNAQ","ITPR2","LINGO1","MEF2A",
                                                                     "MEF2C","NAV3","TBC1D5","MAML3","MTHFD1L","TGFBR1")))
# (Adeoye MIC-M2 ribosomal hubs RPS27A/RPS15/RPS19/RPS28 are left out: ribosomal genes track sequencing depth)
lit <- lit[!duplicated(lit$gene), ]
lit$in_data <- lit$gene %in% rownames(MC)
lit$detect  <- round(ifelse(lit$in_data, det[lit$gene], NA), 3)
lit_genes   <- lit$gene[lit$in_data & lit$detect >= MIN_DETECT]
cat("######## SET lit ########\n")
cat("Genes listed in the papers:", nrow(lit), "| in our data:", sum(lit$in_data),
    "| pass the detection filter:", length(lit_genes), "\n")
cat("Dropped as too sparse:", lit$gene[lit$in_data & lit$detect < MIN_DETECT], "\n")

# ======================================================================= Set "data"
noncoding <- grepl("^(AC|AL|AP|AF|BX|CR|FP|Z)[0-9]+\\.[0-9]+$|^LINC|-AS[0-9]*$|^MIR|-DT$", rownames(MC))
technical <- grepl("^MT-|^RPL|^RPS|^HB[AB]|^HSPA|^HSPH|^HSP90|^DNAJB1$|^FOS|^JUN|^EGR[0-9]$|^IER[0-9]|^DUSP1$|^ATF3$",
                   rownames(MC))
sex_misc  <- rownames(MC) %in% c("MALAT1", "NEAT1", "XIST", "TSIX", "XACT", "RPS4Y1", "DDX3Y", "UTY", "KDM5D",
                                 "USP9Y", "EIF1AY", "ZFY", "NLGN4Y", "TXLNGY", "ARL17A", "ARL17B")
pool <- rownames(MC)[det >= MIN_DETECT & !noncoding & !technical & !sex_misc]
lc <- log1p(as.matrix(MC[pool, ] %*% Diagonal(x = 1e6 / colSums(MC))))   # ranking only

f     <- factor(mc_info$donor)
means <- rowsum(t(lc), f) / as.vector(table(f))          # donor x gene means
resid <- t(lc) - means[as.integer(f), ]                   # metacell x gene, donor mean removed
mu    <- rowMeans(lc)
w_var <- colSums(resid^2) / (ncol(lc) - nlevels(f))       # within-donor variance
t_var <- apply(lc, 1, var)                                # total variance
score <- w_var - predict(loess(w_var ~ mu))               # within-donor variability above the mean trend
data_genes <- names(sort(score, decreasing = TRUE))[seq_len(N_DATA_GENES)]
cat("\n######## SET data ########\n")
cat("Pool:", length(pool), "genes | taking the", N_DATA_GENES, "with the highest within-donor variability\n")
dtab <- data.frame(gene = data_genes, detect = round(det[data_genes], 3),
                   share_within_donor = round(w_var[data_genes] / t_var[data_genes], 2),
                   in_lit_set = data_genes %in% lit_genes, row.names = NULL)
print(dtab, row.names = FALSE)
cat("\nShare of variance within donors: median", round(median(dtab$share_within_donor), 2),
    "(1 = all variability is between metacells of the same person, 0 = only between people)\n")
cat("Overlap with the lit set:", sum(dtab$in_lit_set), "genes\n")

write.csv(rbind(data.frame(set = "lit", gene = lit$gene, source = lit$source, detect = lit$detect,
                           selected = lit$gene %in% lit_genes),
                data.frame(set = "data", gene = data_genes, source = "within-donor variability",
                           detect = dtab$detect, selected = TRUE)),
          "results/08_gene_sets.csv", row.names = FALSE)

# ======================================================================= Metacell tables
save_table <- function(df, name) {
  write.csv(df, paste0("data/processed/", name, ".csv"), row.names = FALSE)
  write.table(df, paste0("data/processed/", name, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Saved", name, ":", nrow(df), "rows x", ncol(df), "variables\n")
}
mc_table <- function(g) data.frame(CPS_Local = mc_info$CPS_Local, total_UMIs = colSums(MC),
                                   as.matrix(t(MC[g, ])), check.names = FALSE)
cat("\n")
save_table(mc_table(lit_genes),  "network_dataset_lit")
save_table(mc_table(data_genes), "network_dataset_data")

# ======================================================================= Single-nucleus tables
cache <- "data/processed/nuclei_all_genesets.rds"
need  <- unique(c(lit_genes, data_genes))
if (!file.exists(cache) || !all(need %in% colnames(readRDS(cache)$G))) {
  cat("\nReading nuclei from the .h5ad for", length(need), "genes (one time only)...\n")
  source("R/read_h5ad.R")
  path   <- "data/SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_RNAseq_final-nuclei.2026-06-22.h5ad"
  md     <- read_h5ad_df(path, "obs")
  genes  <- read_h5ad_df(path, "var")
  counts <- read_h5ad_matrix(path, "UMIs", n_cells = nrow(md), n_genes = nrow(genes))
  meta <- readRDS("data/processed/mtg_microglia_meta.rds")
  meta <- meta[meta$Supertype != "Micro-PVM_1", ]
  meta$donor <- droplevels(meta[["Donor ID"]])
  meta <- meta[order(meta$cell_index), ]
  gi <- match(need, rownames(genes))
  blocks <- split(seq_len(nrow(meta)), ceiling(seq_len(nrow(meta)) / 2000))
  parts <- lapply(seq_along(blocks), function(b) {
    X <- as(counts[, meta$cell_index[blocks[[b]]]], "dgCMatrix")
    if (b %% 5 == 0) cat("  blocks read:", b, "of", length(blocks), "\n")
    list(G = as.matrix(t(X[gi, ])), total = colSums(X))
  })
  G <- do.call(rbind, lapply(parts, `[[`, "G")); colnames(G) <- need
  saveRDS(list(info = meta[, c("cell_index", "donor", "Supertype", "CPS_Local")], G = G,
               total = unlist(lapply(parts, `[[`, "total"))), cache)
}
nuc <- readRDS(cache)
set.seed(1001)                                            # same sampling as 06b rep 1
pick <- unlist(lapply(split(seq_len(nrow(nuc$G)), nuc$info$donor), function(r)
  r[sample.int(length(r), min(N_PER_DONOR, length(r)))]))
nuc_table <- function(g) data.frame(CPS_Local = nuc$info$CPS_Local[pick], total_UMIs = nuc$total[pick],
                                    nuc$G[pick, g, drop = FALSE], check.names = FALSE)
save_table(nuc_table(lit_genes),  "network_dataset_nuclei_lit")
save_table(nuc_table(data_genes), "network_dataset_nuclei_data")
sink()

# ======================================================================= Networks
if (RUN_NETWORKS) {
  for (nm in c("lit", "data", "nuclei_lit", "nuclei_data")) {
    INPUT   <- paste0("data/processed/network_dataset_", nm, ".csv")
    OUT_DIR <- paste0("results/networks_", sub("^nuclei_(.*)$", "\\1_nuclei", nm))
    REPORT  <- file.path(OUT_DIR, "networks_report.txt")
    message("Running networks on ", INPUT)
    source("R/05_networks.R")
  }
  rm(INPUT, OUT_DIR, REPORT)
}
