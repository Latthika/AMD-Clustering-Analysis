# Robust Gene Co-expression Network Analysis for AMD Using AMI & JSD

This repository contains the full workflow and datasets used to perform robustness-based gene network analysis for **Age-related Macular Degeneration (AMD)**. Two information-theoretic methods—**Adjusted Mutual Information (AMI)** and **Jensen–Shannon Divergence (JSD)**—are used to construct co-expression networks, followed by community detection across **100 stochastic seeds**.

The analysis is performed on:

- **81-gene curated AMD dataset**
- **181-gene dataset** containing **81 AMD genes + 100 real noise genes**

Both AMI and JSD scripts share the same structure; only the input dataset is changed when switching between 81-gene and 181-gene runs.

---

## Repository Contents

| File | Description |
|------|-------------|
| `aak100_gene_renamed-2.csv` | Preprocessed AMD dataset with **81 curated genes**. |
| `noise data.csv` | Dataset with **181 genes** (81 AMD genes + 100 real noise genes). |
| `ami_seed_none.R` | Complete AMI-based pipeline (100 seeds). |
| `jsd_seed_none.R` | Complete JSD-based pipeline (100 seeds). |

## Running the AMI Pipeline

Inside **ami_seed_none.R**, modify the input line:

```
input_csv <- "path/to/aak100_gene_renamed-2.csv"
```
input_csv <- "path/to/aak100_gene_renamed-2.csv"

To run the 181-gene noise dataset, change it to:

