# 05_networks.R — network reconstruction with the three course methods (ARACNE, PC, MIIC) and comparison
#
# Input:  data/processed/network_dataset.csv  (983 metacells x 50 variables: CPS_Local, total_UMIs, 48 genes)
# Output: results/05_networks_output.txt       — report
#         results/networks/edges_*.csv          — edge list of each method
#         results/networks/edges_comparison.csv — all edges and which methods found them
#         results/networks/*.png                — network plots (same node layout in all plots)
#
# The data are raw summed counts. ARACNE and PC work on RANKS (Spearman), so the scale of
# counts does not matter to them; MIIC discretizes continuous values itself.

library(minet)
library(pcalg)
library(miic)
library(igraph)
set.seed(2026)

ALPHA_PC  <- 0.01    # significance level of the conditional independence tests in PC
NULL_Q    <- 0.999   # ARACNE threshold: quantile of mutual information on shuffled data
RUN_MIIC  <- TRUE    # MIIC is the slowest (minutes); can be switched off temporarily

dir.create("results/networks", recursive = TRUE, showWarnings = FALSE)
sink("results/05_networks_output.txt", split = TRUE)

df   <- read.csv("data/processed/network_dataset.csv", check.names = FALSE)
vars <- colnames(df); n <- nrow(df); p <- ncol(df)
rho  <- cor(df, method = "spearman")
cat("Data:", n, "metacells x", p, "variables\n")

key  <- function(a, b) ifelse(a < b, paste(a, b, sep = " -- "), paste(b, a, sep = " -- "))
timed <- function(expr) { t0 <- Sys.time(); r <- expr; cat("  time:", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n"); r }

# ======================================================================= ARACNE
# 1) mutual information (MI) for all pairs, estimated from rank correlation: MI = -0.5*log(1 - rho^2)
# 2) significance threshold: MI on data with every column shuffled (no dependencies) -> 99.9% quantile
# 3) ARACNE: in every triangle, removes the weakest edge (likely indirect)
cat("\n######## ARACNE ########\n")
edges_aracne <- timed({
  mim <- build.mim(df, estimator = "spearman")
  null_mi <- unlist(lapply(1:20, function(i) {
    m <- build.mim(as.data.frame(lapply(df, sample)), estimator = "spearman"); m[upper.tri(m)] }))
  thr <- quantile(null_mi, NULL_Q)
  cat("  MI threshold:", signif(thr, 3), "(corresponds to |rho| ~", round(sqrt(1 - exp(-2 * thr)), 3), ")\n")
  mim_thr <- mim; mim_thr[mim_thr < thr] <- 0
  ar  <- aracne(mim_thr, eps = 0)
  idx <- which(ar > 0 & upper.tri(ar), arr.ind = TRUE)
  data.frame(from = vars[idx[, 1]], to = vars[idx[, 2]], type = "--", weight = ar[idx])
})
cat("  edges:", nrow(edges_aracne), "\n")

# ======================================================================= PC
# Starts from the complete graph; removes edge X–Y if X and Y are independent given other variables;
# then orients some of the edges. Test: partial correlations on ranks
# (Spearman -> 2*sin(pi/6*rho): this makes the Gaussian test valid for any monotone relationships).
cat("\n######## PC ########\n")
edges_pc <- timed({
  C <- 2 * sin(pi / 6 * rho)
  pc_fit <- pc(suffStat = list(C = C, n = n), indepTest = gaussCItest, alpha = ALPHA_PC,
               labels = vars, skel.method = "stable", u2pd = "relaxed", solve.confl = TRUE)
  A <- as(pc_fit@graph, "matrix")                 # A[i, j] = 1 means edge i -> j
  out <- list()
  for (i in 1:(p - 1)) for (j in (i + 1):p) {
    if (A[i, j] == 1 && A[j, i] == 1) out[[length(out) + 1]] <- c(vars[i], vars[j], "--")
    else if (A[i, j] == 1)            out[[length(out) + 1]] <- c(vars[i], vars[j], "->")
    else if (A[j, i] == 1)            out[[length(out) + 1]] <- c(vars[j], vars[i], "->")
  }
  e <- as.data.frame(do.call(rbind, out), stringsAsFactors = FALSE)
  colnames(e) <- c("from", "to", "type"); e
})
cat("  edges:", nrow(edges_pc), "| of which oriented:", sum(edges_pc$type == "->"), "\n")

# ======================================================================= MIIC
# Like PC, but uses mutual information with a complexity penalty instead of significance tests;
# orients edges and can flag a latent common cause (type "<->").
edges_miic <- NULL
if (RUN_MIIC) {
  cat("\n######## MIIC ########\n")
  edges_miic <- timed({
    res <- miic(input_data = df, n_threads = max(1, parallel::detectCores() - 1))
    s <- res$summary
    if ("type" %in% names(s)) s <- s[s$type == "P", ]
    o <- s$ort_inferred
    data.frame(from = ifelse(o == -2, s$y, s$x),
               to   = ifelse(o == -2, s$x, s$y),
               type = ifelse(abs(o) == 2, "->", ifelse(o == 6, "<->", "--")),
               stringsAsFactors = FALSE)
  })
  cat("  edges:", nrow(edges_miic), "| oriented:", sum(edges_miic$type == "->"),
      "| with a possible latent cause (<->):", sum(edges_miic$type == "<->"), "\n")
}

nets <- list(ARACNE = edges_aracne, PC = edges_pc, MIIC = edges_miic)
nets <- nets[!sapply(nets, is.null)]
for (m in names(nets)) {
  nets[[m]]$rho <- round(rho[cbind(nets[[m]]$from, nets[[m]]$to)], 3)
  write.csv(nets[[m]], paste0("results/networks/edges_", tolower(m), ".csv"), row.names = FALSE)
}

# ======================================================================= COMPARISON
cat("\n######## METHOD COMPARISON ########\n")
keys <- lapply(nets, function(e) key(e$from, e$to))
all_keys <- sort(unique(unlist(keys)))
comp <- data.frame(edge = all_keys, sapply(keys, function(k) all_keys %in% k), check.names = FALSE)
ab <- do.call(rbind, strsplit(all_keys, " -- "))
comp$rho <- round(rho[ab], 3)
comp$n_methods <- rowSums(comp[, names(nets), drop = FALSE])
write.csv(comp, "results/networks/edges_comparison.csv", row.names = FALSE)

cat("\nNumber of edges:\n"); print(sapply(keys, length))
cat("\nEdge overlap between methods (Jaccard index = shared / union):\n")
J <- outer(names(keys), names(keys), Vectorize(function(a, b)
  round(length(intersect(keys[[a]], keys[[b]])) / length(union(keys[[a]], keys[[b]])), 2)))
dimnames(J) <- list(names(keys), names(keys)); print(J)
cat("\nEdges found by all methods:", sum(comp$n_methods == length(nets)), "\n")
print(comp[comp$n_methods == length(nets), c("edge", "rho")], row.names = FALSE)

neigh <- function(e, v) sort(unique(c(e$to[e$from == v], e$from[e$to == v])))
for (v in c("CPS_Local", "total_UMIs")) {
  cat("\nNeighbours of", v, ":\n")
  for (m in names(nets)) {
    nb <- neigh(nets[[m]], v)
    cat(sprintf("  %-7s (%2d): %s\n", m, length(nb),
                paste0(nb, " (", round(rho[v, nb], 2), ")", collapse = ", ")))
  }
}

cat("\nHubs — 10 nodes with the most connections (excluding total_UMIs):\n")
for (m in names(nets)) {
  d <- table(c(nets[[m]]$from, nets[[m]]$to)); d <- d[names(d) != "total_UMIs"]
  d <- sort(d, decreasing = TRUE)[1:min(10, length(d))]
  cat(sprintf("  %-7s %s\n", m, paste0(names(d), "=", d, collapse = ", ")))
}

if (all(c("PC", "MIIC") %in% names(nets))) {
  dir_pc <- nets$PC[nets$PC$type == "->", ]; dir_mi <- nets$MIIC[nets$MIIC$type == "->", ]
  both <- intersect(key(dir_pc$from, dir_pc$to), key(dir_mi$from, dir_mi$to))
  same <- intersect(paste(dir_pc$from, dir_pc$to), paste(dir_mi$from, dir_mi$to))
  cat("\nOrientation: edges oriented in both PC and MIIC:", length(both),
      "| of which in the same direction:", length(same), "\n")
}

# ======================================================================= PLOTS
grp <- setNames(rep("reactive (data)", p), vars)
grp[c("P2RY12","CX3CR1","TMEM119","CSF1R","P2RY13","SALL1","GPR34","OLFML3","HEXB","SELPLG","TGFBR1","SLC2A5","CST3")] <- "homeostatic"
grp[c("TREM2","TYROBP","APOE","SPP1","GPNMB","CD74","HLA-DRA","C1QA","C1QB","CD163")] <- "disease response"
grp[c("PLCG2","INPP5D","MS4A6A","PICALM","SORL1","BIN1","CD33","ABI3")] <- "AD risk"
grp[c("SPI1","IRF8","MEF2C")] <- "TF"
grp["CPS_Local"] <- "CPS"; grp["total_UMIs"] <- "technical"
pal <- c(homeostatic = "#4C9BE8", "disease response" = "#E8833A", "AD risk" = "#5DBB63",
         TF = "#9B6BD0", "reactive (data)" = "#F2C14E", CPS = "#D64545", technical = "#9E9E9E")

g_all <- graph_from_data_frame(ab, vertices = data.frame(name = vars), directed = FALSE)
lay <- layout_with_fr(g_all)

plot_net <- function(e, title, file) {
  g <- graph_from_data_frame(e[, c("from", "to")], vertices = data.frame(name = vars), directed = TRUE)
  png(file, width = 2000, height = 2000, res = 220)
  plot(g, layout = lay, vertex.color = pal[grp[V(g)$name]], vertex.size = 7,
       vertex.label.cex = 0.55, vertex.label.color = "black", vertex.frame.color = NA,
       edge.arrow.mode = ifelse(e$type == "->", 2, ifelse(e$type == "<->", 3, 0)),
       edge.arrow.size = 0.35, edge.color = ifelse(e$rho > 0, "#555555", "#C0392B"), main = title)
  legend("bottomleft", legend = names(pal), fill = pal, bty = "n", cex = 0.7)
  dev.off()
}
for (m in names(nets)) plot_net(nets[[m]], paste(m, "—", nrow(nets[[m]]), "edges"),
                                paste0("results/networks/net_", tolower(m), ".png"))
cons <- comp[comp$n_methods == length(nets), ]
if (nrow(cons) > 0) {
  ce <- data.frame(do.call(rbind, strsplit(cons$edge, " -- ")), type = "--", rho = cons$rho)
  colnames(ce)[1:2] <- c("from", "to")
  plot_net(ce, paste("Consensus of all methods —", nrow(ce), "edges"), "results/networks/net_consensus.png")
}
cat("\nPlots: results/networks/*.png (grey edges = positive association, red = negative)\n")
sink()
