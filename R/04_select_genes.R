# 04_select_genes.R — rebuild the gene set (fast: no .h5ad reading, uses the saved metacells)
#
# Why: in step 03 the "most variable" genes turned out to be mostly
#   - unnamed non-coding RNAs (AC..., AL..., LINC..., ...-AS1),
#   - tissue-handling stress genes (HSPA1A, HSPA1B, DUSP1),
#   - genes of other cell types (FLT1 — endothelial cells) and donor-genotype-driven genes (ARL17B).
# They vary between donors for technical/genetic reasons and are hard to interpret.
# New data-driven criterion: genes that best separate microglial STATES
#   (homeostatic Micro-PVM_2 vs reactive Micro-PVM_3/4) — the transition we are interested in.
# Normalization is used only for ranking; the network gets raw summed counts.

library(Matrix)
set.seed(2026)
N_GENES    <- 48
MIN_DETECT <- 0.8   # was 0.9: this way GPNMB (0.85), a key disease-response gene, enters the network

sink("results/04_select_genes_output.txt", split = TRUE)
MC      <- readRDS("data/processed/metacell_counts_all_genes.rds")
mc_info <- read.csv("data/processed/metacell_info.csv", check.names = FALSE)
stopifnot(identical(colnames(MC), mc_info$group))
det <- rowMeans(MC > 0)

bio <- c("P2RY12", "CX3CR1", "TMEM119", "CSF1R",
         "TREM2", "TYROBP", "APOE", "SPP1", "GPNMB", "CD74", "HLA-DRA",
         "C1QA", "C1QB", "CD163",
         "PLCG2", "INPP5D", "MS4A6A", "PICALM", "SORL1", "BIN1", "CD33", "ABI3",
         "SPI1", "IRF8", "MEF2C",
         # known human homeostatic microglia genes (used instead of data-driven selection for this side:
         # the genes "higher in homeostatic" included other cell types' genes — MMRN1, TNFRSF13C, CLEC9A)
         "P2RY13", "SALL1", "GPR34", "OLFML3", "HEXB", "SELPLG", "TGFBR1", "SLC2A5", "CST3")
bio_ok <- bio[det[bio] >= MIN_DETECT]
cat("Biology genes:", length(bio_ok), "| dropped:", setdiff(bio, bio_ok), "\n")

# Genes excluded from data-driven selection
noncoding <- grepl("^(AC|AL|AP|AF|BX|CR|FP|Z)[0-9]+\\.[0-9]+$|^LINC|-AS[0-9]*$|^MIR|-DT$", rownames(MC))
technical <- grepl("^MT-|^RPL|^RPS|^HB[AB]|^HSPA|^HSPH|^HSP90|^DNAJB1$|^FOS|^JUN|^EGR[0-9]$|^IER[0-9]|^DUSP1$|^ATF3$",
                   rownames(MC))
sex_misc  <- rownames(MC) %in% c("MALAT1", "NEAT1", "XIST", "TSIX", "XACT", "RPS4Y1", "DDX3Y", "UTY", "KDM5D",
                                 "USP9Y", "EIF1AY", "ZFY", "NLGN4Y", "TXLNGY", "ARL17A", "ARL17B")
pool <- rownames(MC)[det >= MIN_DETECT & !noncoding & !technical & !sex_misc & !rownames(MC) %in% bio_ok]
cat("Pool for data-driven selection:", length(pool), "genes\n")

# Counts per million + log — only to compare the states
lc <- log1p(as.matrix(MC[pool, ] %*% Diagonal(x = 1e6 / colSums(MC))))
homeo    <- mc_info$Supertype == "Micro-PVM_2"
reactive <- mc_info$Supertype %in% c("Micro-PVM_3-SEAAD", "Micro-PVM_4-SEAAD")
cat("Metacells: homeostatic", sum(homeo), "| reactive", sum(reactive), "\n")
diff <- rowMeans(lc[, reactive]) - rowMeans(lc[, homeo])   # log fold difference

n_data <- N_GENES - length(bio_ok)
up   <- names(sort(diff, decreasing = TRUE))[seq_len(n_data)]   # higher in reactive microglia
down <- character(0)
sel  <- c(bio_ok, up)

dataset <- data.frame(CPS_Local = mc_info$CPS_Local, total_UMIs = colSums(MC),
                      as.matrix(t(MC[sel, ])), check.names = FALSE)

rep_tab <- data.frame(
  gene      = sel,
  source    = c(rep("biology", length(bio_ok)), rep("up_in_reactive", length(up)), rep("up_in_homeostatic", length(down))),
  detect    = round(det[sel], 3),
  logFC_reactive_vs_homeo = round(c(rep(NA, length(bio_ok)), diff[up], diff[down]), 2),
  rho_CPS   = round(sapply(sel, function(g) cor(dataset[[g]], dataset$CPS_Local, method = "spearman")), 2),
  rho_depth = round(sapply(sel, function(g) cor(dataset[[g]], dataset$total_UMIs, method = "spearman")), 2),
  row.names = NULL)
cat("\nFinal genes:\n"); print(rep_tab, row.names = FALSE)
cat("\nTable size:", dim(dataset), "\n")

if (!file.exists("data/processed/network_dataset_v1_hvg.csv"))
  file.copy("data/processed/network_dataset.csv", "data/processed/network_dataset_v1_hvg.csv")
write.csv(dataset, "data/processed/network_dataset.csv", row.names = FALSE)
write.table(dataset, "data/processed/network_dataset.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(rep_tab, "results/04_selected_genes.csv", row.names = FALSE)
cat("Saved: data/processed/network_dataset.csv/.tsv (previous version: network_dataset_v1_hvg.csv)\n")
sink()
