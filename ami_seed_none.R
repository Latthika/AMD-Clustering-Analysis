# -----------------------------
# AMI-based Network Analysis
# -----------------------------

# Required packages
needed <- c(
  "igraph", "dplyr", "openxlsx", "RColorBrewer",
  "scales", "WGCNA", "aricode", "preprocessCore", "tidyr"
)

for (p in needed) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(igraph)
library(dplyr)
library(openxlsx)
library(RColorBrewer)
library(scales)
library(WGCNA)
library(aricode)
library(tidyr)

options(stringsAsFactors = FALSE)

# -----------------------------
# Parameters
# -----------------------------

input_csv <- "/Users/latthika/Documents/amd dataset/aak100_gene_renamed-2.csv"
out_dir <- "/Users/latthika/Documents/amd dataset/excel comparisionfinal/AMI_result"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_xlsx <- file.path(
  out_dir,
  paste0("AMI_FINAL_", format(Sys.Date(), "%Y%m%d"), ".xlsx")
)

seeds <- 1:100
edge_keep_pct <- 0.05
ami_bins <- 8
top_pct_for_hubs <- 0.10

# -----------------------------
# Helper functions
# -----------------------------

coerce_numeric_df <- function(df) {
  df[] <- lapply(df, function(col) as.numeric(as.character(col)))
  df
}

discretize_equalfreq <- function(vec, bins = 8) {
  if (all(is.na(vec))) return(factor(rep(NA_integer_, length(vec))))
  probs <- seq(0, 1, length.out = bins + 1)
  br <- unique(quantile(vec, probs = probs, na.rm = TRUE))
  if (length(br) <= 1) return(factor(rep(1L, length(vec))))
  cut_vec <- cut(vec, breaks = br, include.lowest = TRUE, labels = FALSE)
  cut_vec[is.na(cut_vec)] <- max(cut_vec, na.rm = TRUE)
  factor(cut_vec)
}

adjacency_ami <- function(expr_df, bins = 8) {
  n <- ncol(expr_df)
  gnames <- colnames(expr_df)
  
  disc <- lapply(expr_df, discretize_equalfreq, bins = bins)
  names(disc) <- gnames
  
  entropy_vec <- sapply(disc, function(x) {
    p <- table(x) / sum(table(x))
    -sum(p * log2(p + 1e-12))
  })
  
  adj <- matrix(0, n, n, dimnames = list(gnames, gnames))
  
  for (i in seq_len(n)) {
    for (j in i:n) {
      tab <- table(disc[[i]], disc[[j]])
      p <- tab / sum(tab)
      Hxy <- -sum(p * log2(p + 1e-12))
      mi <- max(0, entropy_vec[i] + entropy_vec[j] - Hxy)
      adj[i, j] <- adj[j, i] <- ifelse(
        entropy_vec[i] > 0 && entropy_vec[j] > 0,
        mi / sqrt(entropy_vec[i] * entropy_vec[j]),
        0
      )
    }
  }
  
  diag(adj) <- 0
  adj
}

threshold_top_pct <- function(adj, pct) {
  w <- adj[upper.tri(adj)]
  cutoff <- quantile(w, probs = 1 - pct, na.rm = TRUE)
  adj[adj < cutoff] <- 0
  diag(adj) <- 0
  adj
}

build_igraph_from_adj <- function(adj, genes) {
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE)
  V(g)$name <- genes
  g
}

# -----------------------------
# Data loading and partitioning
# -----------------------------

raw <- read.csv(input_csv, check.names = FALSE)
stopifnot("mgs_level" %in% colnames(raw))

traits <- raw$mgs_level
expr_raw <- coerce_numeric_df(raw %>% select(-mgs_level))

expr_partitions <- list(
  Full = expr_raw,
  MGS1 = expr_raw[traits == "MGS1", ],
  MGS4 = expr_raw[traits == "MGS4", ]
)

# -----------------------------
# Seed-wise analysis
# -----------------------------

global_records <- list()
gene_records <- list()
top_genes_records <- list()

for (seed in seeds) {
  set.seed(seed)
  gene_records[[as.character(seed)]] <- list()
  top_genes_records[[as.character(seed)]] <- list()
  
  for (pname in names(expr_partitions)) {
    expr_df <- expr_partitions[[pname]]
    genes <- colnames(expr_df)
    
    adj_full <- adjacency_ami(expr_df, bins = ami_bins)
    adj_thr <- threshold_top_pct(adj_full, edge_keep_pct)
    g <- build_igraph_from_adj(adj_thr, genes)
    
    wc <- cluster_walktrap(g)
    mem_wc <- membership(wc)
    
    if ("cluster_leiden" %in% ls("package:igraph")) {
      ld <- cluster_leiden(g)
    } else {
      ld <- cluster_louvain(g)
    }
    mem_ld <- membership(ld)
    
    imc_wc <- intramodularConnectivity(adj_full, mem_wc)
    imc_ld <- intramodularConnectivity(adj_full, mem_ld)
    
    top_n <- ceiling(nrow(imc_wc) * top_pct_for_hubs)
    
    top_genes_records[[as.character(seed)]][[pname]] <- list(
      Walktrap = rownames(imc_wc)[order(imc_wc$kWithin, decreasing = TRUE)][1:top_n],
      Leiden = rownames(imc_ld)[order(imc_ld$kWithin, decreasing = TRUE)][1:top_n]
    )
  }
}

# -----------------------------
# Output
# -----------------------------

wb <- createWorkbook()
addWorksheet(wb, "Summary")
writeData(wb, "Summary", data.frame(Status = "Completed"))

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

cat("Analysis completed. Output saved to:", out_xlsx, "\n")
