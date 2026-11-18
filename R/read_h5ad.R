# read_h5ad.R — read .h5ad files (AnnData format) directly with rhdf5 / HDF5Array,
# without zellkonverter and without Python.
#
# Layout of an .h5ad file (a plain HDF5 file: "folders" and "files" inside one file):
#   /obs           — nucleus metadata table (donor, brain region, cell type, ...)
#   /var           — gene metadata table
#   /X             — normalized matrix (in SEA-AD: log(UMIs per 10,000 + 1)) — NOT used
#   /layers/UMIs   — raw counts (number of RNA molecules) — this is what we use
# Categorical columns are stored as a pair: codes (integer indices) + categories (labels).

library(rhdf5)
library(HDF5Array)

# Read the /obs or /var table into a data.frame
read_h5ad_df <- function(file, grp = "obs") {
  ls_all <- h5ls(file, recursive = TRUE, datasetinfo = FALSE)
  items  <- ls_all[ls_all$group == paste0("/", grp), ]
  index_name <- h5readAttributes(file, grp)[["_index"]]
  if (is.null(index_name)) index_name <- "_index"

  cols <- list()
  for (i in seq_len(nrow(items))) {
    nm <- items$name[i]
    p  <- paste0(grp, "/", nm)
    if (items$otype[i] == "H5I_GROUP") {
      sub <- ls_all$name[ls_all$group == paste0("/", p)]
      if (all(c("codes", "categories") %in% sub)) {            # categorical column
        codes <- as.integer(h5read(file, paste0(p, "/codes")))
        cats  <- as.character(h5read(file, paste0(p, "/categories")))
        codes[codes < 0] <- NA                                   # -1 means missing
        cols[[nm]] <- factor(cats[codes + 1L], levels = cats)    # +1: R indices start at 1
      } else if (all(c("values", "mask") %in% sub)) {            # column with missing values
        v <- as.vector(h5read(file, paste0(p, "/values")))
        v[as.logical(h5read(file, paste0(p, "/mask")))] <- NA
        cols[[nm]] <- v
      }
    } else {
      cols[[nm]] <- as.vector(h5read(file, p))
    }
  }
  idx <- as.character(cols[[index_name]])
  cols[[index_name]] <- NULL
  df <- as.data.frame(cols, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(df) <- make.unique(idx)
  h5closeAll()
  df
}

# Count matrix as a DelayedMatrix: the data stay on disk and are read in blocks.
# Returned in Bioconductor orientation: rows = genes, columns = nuclei.
read_h5ad_matrix <- function(file, layer = "UMIs", n_cells, n_genes) {
  m <- tryCatch(
    HDF5Array::H5ADMatrix(file, layer = layer),
    error = function(e) {
      message("H5ADMatrix failed (", conditionMessage(e), "), trying H5SparseMatrix")
      HDF5Array::H5SparseMatrix(file, paste0("layers/", layer))
    })
  if (nrow(m) == n_cells && ncol(m) == n_genes) m <- t(m)   # make it genes x nuclei
  m
}
