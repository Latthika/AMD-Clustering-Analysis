# JSD seed-based network analysis with Walktrap and Leiden clustering

needed <- c(
  "igraph", "dplyr", "openxlsx", "RColorBrewer",
  "scales", "WGCNA", "preprocessCore", "tidyr", "aricode"
)

for (p in needed) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(igraph)
library(dplyr)
library(openxlsx)
library(RColorBrewer)
library(scales)
library(WGCNA)
library(tidyr)
library(aricode)

options(stringsAsFactors = FALSE)

# ---- Parameters ----
input_csv <- "/Users/latthika/Documents/amd dataset/aak100_gene_renamed-2.csv"
out_dir <- "/Users/latthika/Documents/amd dataset/excel comparisionfinal/JSD_result"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_xlsx <- file.path(
  out_dir,
  paste0("JSD_seed_summary_noNorm_", format(Sys.Date(), "%Y%m%d"), ".xlsx")
)

seeds <- 1:100
edge_keep_pct <- 0.05
top_pct_for_hubs <- 0.10

# ---- JSD functions ----
jsd_single <- function(p, q) {
  p <- p / sum(p)
  q <- q / sum(q)
  m <- 0.5 * (p + q)
  0.5 * sum(p * log((p + 1e-12) / (m + 1e-12))) +
    0.5 * sum(q * log((q + 1e-12) / (m + 1e-12)))
}

adjacency_jsd <- function(expr_df) {
  n <- ncol(expr_df)
  genes <- colnames(expr_df)
  adj <- matrix(0, n, n, dimnames = list(genes, genes))
  expr_prob <- apply(expr_df + 1e-9, 2, function(v) v / sum(v))
  
  for (i in seq_len(n)) {
    for (j in i:n) {
      sim <- 1 - jsd_single(expr_prob[, i], expr_prob[, j])
      adj[i, j] <- sim
      adj[j, i] <- sim
    }
  }
  diag(adj) <- 0
  adj
}

# ---- Helpers ----
coerce_numeric_df <- function(df) {
  df[] <- lapply(df, function(x) as.numeric(as.character(x)))
  df
}

threshold_top_pct <- function(adj, pct) {
  cutoff <- quantile(adj[upper.tri(adj)], probs = 1 - pct, na.rm = TRUE)
  adj[adj < cutoff] <- 0
  diag(adj) <- 0
  adj
}

build_igraph_from_adj <- function(adj, genes) {
  g <- graph_from_adjacency_matrix(adj, mode = "undirected",
                                   weighted = TRUE, diag = FALSE)
  V(g)$name <- genes
  g
}

compute_hub_metrics_from_adj <- function(expr_df, adj_full, adj_thr, g, membership, softPower = 6) {
  genes <- colnames(expr_df)
  imc <- intramodularConnectivity(adjMat = adj_full, colors = membership)
  
  kTotal_alt <- tryCatch(
    softConnectivity(as.data.frame(expr_df), power = softPower),
    error = function(e) rep(NA, length(genes))
  )
  
  MEs <- tryCatch(moduleEigengenes(expr_df, membership)$eigengenes,
                  error = function(e) NULL)
  
  kME_all <- tryCatch(signedKME(expr_df, MEs),
                      error = function(e) NULL)
  
  assigned_kME <- rep(NA_real_, length(genes))
  names(assigned_kME) <- genes
  
  if (!is.null(kME_all)) {
    me_labels <- sub("^kME", "", colnames(kME_all))
    for (i in seq_along(genes)) {
      mod <- as.character(membership[i])
      idx <- which(me_labels == mod)
      assigned_kME[i] <- if (length(idx) == 1) kME_all[i, idx] else NA
    }
  }
  
  data.frame(
    gene = genes,
    module = membership,
    kTotal = imc$kTotal,
    kWithin = imc$kWithin,
    kOut = imc$kOut,
    kDiff = imc$kDiff,
    kME = assigned_kME,
    softConn = kTotal_alt,
    strength = strength(g),
    degree = degree(g),
    betweenness = betweenness(g, weights = 1 / (E(g)$weight + 1e-12)),
    stringsAsFactors = FALSE
  ) %>% arrange(desc(kWithin))
}

summarize_graph_basic <- function(g, mem_wc, mem_ld, adj) {
  common <- intersect(names(mem_wc), names(mem_ld))
  
  data.frame(
    global_density = edge_density(g),
    global_transitivity = transitivity(g, type = "global"),
    Walktrap_modularity = modularity(g, mem_wc),
    Leiden_modularity = modularity(g, mem_ld),
    Walktrap_mean_kWithin =
      mean(intramodularConnectivity(adj, mem_wc)$kWithin, na.rm = TRUE),
    Leiden_mean_kWithin =
      mean(intramodularConnectivity(adj, mem_ld)$kWithin, na.rm = TRUE),
    ARI = ARI(mem_wc[common], mem_ld[common]),
    NMI = NMI(mem_wc[common], mem_ld[common])
  )
}

# ---- Load data ----
raw <- read.csv(input_csv, check.names = FALSE)
traits <- raw$mgs_level
expr_raw <- coerce_numeric_df(raw %>% select(-mgs_level))

expr_partitions <- list(
  Full = expr_raw,
  MGS1 = expr_raw[traits == "MGS1", ],
  MGS4 = expr_raw[traits == "MGS4", ]
)

# ---- Storage ----
global_records <- list()
gene_records <- list()
top_genes_records <- list()

# ---- Seed loop ----
for (seed in seeds) {
  set.seed(seed)
  gene_records[[seed]] <- list()
  top_genes_records[[seed]] <- list()
  
  for (pname in names(expr_partitions)) {
    expr_df <- expr_partitions[[pname]]
    if (ncol(expr_df) == 0) next
    
    adj <- adjacency_jsd(expr_df)
    adj_thr <- threshold_top_pct(adj, edge_keep_pct)
    g <- build_igraph_from_adj(adj_thr, colnames(expr_df))
    
    wc <- cluster_walktrap(g)
    mem_wc <- membership(wc)
    
    ld <- if ("cluster_leiden" %in% ls("package:igraph"))
      cluster_leiden(g) else cluster_louvain(g)
    mem_ld <- membership(ld)
    
    walk_df <- compute_hub_metrics_from_adj(expr_df, adj, adj_thr, g, mem_wc)
    leid_df <- compute_hub_metrics_from_adj(expr_df, adj, adj_thr, g, mem_ld)
    
    top_n <- ceiling(nrow(walk_df) * top_pct_for_hubs)
    
    gene_records[[seed]][[pname]] <- list(Walktrap = walk_df, Leiden = leid_df)
    top_genes_records[[seed]][[pname]] <- list(
      Walktrap = walk_df$gene[1:top_n],
      Leiden = leid_df$gene[1:top_n]
    )
    
    sum_row <- summarize_graph_basic(g, mem_wc, mem_ld, adj)
    
    global_records[[length(global_records) + 1]] <-
      cbind(seed = seed, partition = pname, algorithm = "Walktrap",
            sum_row[, c("global_density", "global_transitivity",
                        "Walktrap_modularity", "Walktrap_mean_kWithin",
                        "ARI", "NMI")])
    
    global_records[[length(global_records) + 1]] <-
      cbind(seed = seed, partition = pname, algorithm = "Leiden",
            sum_row[, c("global_density", "global_transitivity",
                        "Leiden_modularity", "Leiden_mean_kWithin",
                        "ARI", "NMI")])
  }
}

# ---- Output ----
global_df <- bind_rows(global_records)

wb <- createWorkbook()
addWorksheet(wb, "global_per_seed")
writeDataTable(wb, "global_per_seed", global_df)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat("Saved:", out_xlsx, "\n")
