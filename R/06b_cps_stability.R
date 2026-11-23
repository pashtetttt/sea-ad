# 06b_cps_stability.R — how stable are the CPS_Local neighbours, and how many appear "just from donors"?
#
# Three checks, all with the same method settings as 05_networks.R:
#   1. Stability:        repeat experiment B (12 random nuclei per donor) with N_SEEDS different samples.
#                        If PTPRG / TGFBI are CPS neighbours in every sample, the result is robust.
#   2. Shuffle, nuclei:  same samples, but CPS values are randomly reassigned BETWEEN DONORS
#                        (each donor still has one value, but someone else's). Any CPS neighbour found
#                        then comes from donor identity, not from disease severity.
#   3. Shuffle, metacells: the same shuffle on the metacell table — tests the 19 MIIC neighbours.
#
# Output: results/06b_cps_stability_output.txt — report
#         results/06b_cps_stability_runs.csv   — one row per run x method (CPS neighbours, edge counts)
# Cache:  data/processed/nuclei_all_48genes.rds — all MTG microglia nuclei, 48 genes + total UMIs
#         (built once from the .h5ad, ~10 min; later runs skip this step)
# Run time: roughly 30-45 min (PC on nuclei ~2.5 min and MIIC on metacells ~5 min per run).

library(Matrix)
library(minet)
library(pcalg)
library(miic)
set.seed(2026)

N_SEEDS            <- 5      # independent nucleus samples for experiment B
N_PER_DONOR        <- 12     # nuclei per donor in each sample
RUN_SHUFFLE_NUCLEI <- TRUE   # check 2
N_SHUFFLE_MC       <- 2      # check 3: number of shuffled metacell runs (0 to skip)
ALPHA_PC           <- 0.01   # same as 05_networks.R
NULL_Q             <- 0.999  # same as 05_networks.R

sink("results/06b_cps_stability_output.txt", split = TRUE)

# ---------------------------------------------------------------- Cache of all nuclei (48 genes)
cache <- "data/processed/nuclei_all_48genes.rds"
sel_genes <- read.csv("results/04_selected_genes.csv")$gene
if (!file.exists(cache)) {
  cat("Building nucleus cache from the .h5ad (one time only)...\n")
  source("R/read_h5ad.R")
  path   <- "data/SEAAD_Immune_HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_AnG_V1C_RNAseq_final-nuclei.2026-06-22.h5ad"
  md     <- read_h5ad_df(path, "obs")
  genes  <- read_h5ad_df(path, "var")
  counts <- read_h5ad_matrix(path, "UMIs", n_cells = nrow(md), n_genes = nrow(genes))
  meta <- readRDS("data/processed/mtg_microglia_meta.rds")
  meta <- meta[meta$Supertype != "Micro-PVM_1", ]
  meta$donor <- droplevels(meta[["Donor ID"]])
  meta <- meta[order(meta$cell_index), ]
  gi <- match(sel_genes, rownames(genes))
  blocks <- split(seq_len(nrow(meta)), ceiling(seq_len(nrow(meta)) / 2000))   # read 2000 nuclei at a time
  parts <- lapply(seq_along(blocks), function(b) {
    X <- as(counts[, meta$cell_index[blocks[[b]]]], "dgCMatrix")
    if (b %% 5 == 0) cat("  blocks read:", b, "of", length(blocks), "\n")
    list(G = as.matrix(t(X[gi, ])), total = colSums(X))
  })
  G <- do.call(rbind, lapply(parts, `[[`, "G")); colnames(G) <- sel_genes
  saveRDS(list(info = meta[, c("cell_index", "donor", "Supertype", "CPS_Local")],
               G = G, total = unlist(lapply(parts, `[[`, "total"))), cache)
}
nuc <- readRDS(cache)
cat("Nucleus pool:", nrow(nuc$G), "nuclei from", length(unique(nuc$info$donor)), "donors\n")

make_nuclei_table <- function(seed) {
  set.seed(seed)
  pick <- unlist(lapply(split(seq_len(nrow(nuc$G)), nuc$info$donor), function(r)
    r[sample.int(length(r), min(N_PER_DONOR, length(r)))]))
  list(df = data.frame(CPS_Local = nuc$info$CPS_Local[pick], total_UMIs = nuc$total[pick],
                       nuc$G[pick, ], check.names = FALSE),
       donor = as.character(nuc$info$donor[pick]))
}

# CPS values reassigned between donors: each donor keeps one value, but a random other donor's
shuffle_cps <- function(df, donor, seed) {
  set.seed(seed)
  per_donor <- tapply(df$CPS_Local, donor, `[`, 1)
  new_val   <- setNames(sample(per_donor), names(per_donor))
  df$CPS_Local <- unname(new_val[donor])
  df
}

# ---------------------------------------------------------------- The three methods (as in 05)
run_aracne <- function(df) {
  mim <- build.mim(df, estimator = "spearman")
  null_mi <- unlist(lapply(1:20, function(i) {
    m <- build.mim(as.data.frame(lapply(df, sample)), estimator = "spearman"); m[upper.tri(m)] }))
  mim[mim < quantile(null_mi, NULL_Q)] <- 0
  ar <- aracne(mim, eps = 0)
  idx <- which(ar > 0 & upper.tri(ar), arr.ind = TRUE)
  data.frame(from = colnames(df)[idx[, 1]], to = colnames(df)[idx[, 2]])
}
run_pc <- function(df) {
  C <- 2 * sin(pi / 6 * cor(df, method = "spearman"))
  fit <- pc(suffStat = list(C = C, n = nrow(df)), indepTest = gaussCItest, alpha = ALPHA_PC,
            labels = colnames(df), skel.method = "stable", u2pd = "relaxed", solve.confl = TRUE)
  A <- as(fit@graph, "matrix"); A <- (A + t(A)) > 0
  idx <- which(A & upper.tri(A), arr.ind = TRUE)
  data.frame(from = colnames(df)[idx[, 1]], to = colnames(df)[idx[, 2]])
}
run_miic <- function(df) {
  s <- miic(input_data = df, n_threads = max(1, parallel::detectCores() - 1))$summary
  if ("type" %in% names(s)) s <- s[s$type == "P", ]
  data.frame(from = s$x, to = s$y)
}
cps_nb <- function(e) sort(unique(c(e$to[e$from == "CPS_Local"], e$from[e$to == "CPS_Local"])))

runs <- list()
run_all <- function(df, variant, rep) {
  for (m in c("ARACNE", "PC", "MIIC")) {
    t0 <- Sys.time()
    e  <- switch(m, ARACNE = run_aracne(df), PC = run_pc(df), MIIC = run_miic(df))
    nb <- setdiff(cps_nb(e), "total_UMIs")
    runs[[length(runs) + 1]] <<- data.frame(variant = variant, rep = rep, method = m,
      n_edges = nrow(e), n_cps_neighbours = length(nb), cps_neighbours = paste(nb, collapse = ", "))
    cat(sprintf("  %-18s rep %d  %-6s edges %3d | CPS neighbours %2d: %s  (%.1f min)\n", variant, rep, m,
                nrow(e), length(nb), paste(nb, collapse = ", "),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}

# ---------------------------------------------------------------- 1 + 2. Nuclei: real and shuffled CPS
cat("\n######## 1-2. SINGLE NUCLEI:", N_SEEDS, "samples of", N_PER_DONOR, "nuclei per donor ########\n")
for (s in seq_len(N_SEEDS)) {
  tab <- make_nuclei_table(1000 + s)
  run_all(tab$df, "nuclei_real", s)
  if (RUN_SHUFFLE_NUCLEI) run_all(shuffle_cps(tab$df, tab$donor, 2000 + s), "nuclei_shuffled", s)
}

# ---------------------------------------------------------------- 3. Metacells: shuffled CPS
if (N_SHUFFLE_MC > 0) {
  cat("\n######## 3. METACELLS WITH SHUFFLED CPS ########\n")
  mc_df   <- read.csv("data/processed/network_dataset.csv", check.names = FALSE)
  mc_info <- read.csv("data/processed/metacell_info.csv", check.names = FALSE)
  stopifnot(nrow(mc_df) == nrow(mc_info), all(mc_df$CPS_Local == mc_info$CPS_Local))
  for (r in seq_len(N_SHUFFLE_MC)) run_all(shuffle_cps(mc_df, as.character(mc_info$donor), 3000 + r), "metacells_shuffled", r)
}

# ---------------------------------------------------------------- Summary
res <- do.call(rbind, runs)
write.csv(res, "results/06b_cps_stability_runs.csv", row.names = FALSE)

cat("\n######## SUMMARY: number of CPS neighbours (excluding total_UMIs) ########\n")
cat("Reference (05): metacells real CPS -> ARACNE 1, PC 0, MIIC 18 (+ total_UMIs)\n\n")
print(aggregate(n_cps_neighbours ~ variant + method, res, function(x)
  c(mean = round(mean(x), 1), min = min(x), max = max(x))))

cat("\nHow often each gene is a CPS neighbour across the", N_SEEDS, "real nucleus samples:\n")
real <- res[res$variant == "nuclei_real", ]
for (m in c("ARACNE", "PC", "MIIC")) {
  g <- unlist(strsplit(real$cps_neighbours[real$method == m], ", "))
  f <- sort(table(g[g != ""]), decreasing = TRUE)
  cat(sprintf("  %-6s %s\n", m, if (length(f)) paste0(names(f), " ", f, "/", N_SEEDS, collapse = ", ") else "none"))
}
cat("\nSaved: results/06b_cps_stability_runs.csv\n")
sink()
