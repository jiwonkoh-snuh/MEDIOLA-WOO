# PAM50 molecular subtyping (reference-population centering)
#
# PART 1 produces the reported PAM50 calls: each PAM50 gene is centered on the
# per-gene median of the full TCGA-BRCA primary-tumor cohort (rather than on the
# skewed study cohort), then samples are assigned to the genefu PAM50 centroids by
# Spearman nearest-centroid, so every tumor receives a call.
#
# PART 2 is the benchmark that justified this centering scheme. Using TCGA-BRCA
# tumors with established PAM50 labels as ground truth, two schemes are compared
# on basal-enriched test sets matching the study composition:
#   A = percentile subgroup-specific centering (Zhao et al.)
#   B = full-reference-median centering (adopted)
# Over 30 random partitions, full-reference-median centering was more accurate
# (96.4% vs 91.4%) and misassigned far fewer non-Basal tumors to Basal
# (4.5% vs 48.0%), so it is used in PART 1.
#
# Inputs  (place under data/):
#   rnaseq_tpm.txt   gene x sample TPM (RNA-seq; see GEO)
#   tcga_brca_reference.rds     TCGA-BRCA SummarizedExperiment (TCGAbiolinks STAR-Counts)
# Outputs:
#   PAM50_calls.csv                 reported per-sample PAM50 classifications
#   PAM50_centering_benchmark.csv   30-seed benchmark supporting the scheme choice

library(genefu)
library(org.Hs.eg.db)
library(SummarizedExperiment)
library(TCGAbiolinks)

data_dir <- "data"
out_dir  <- "."
exclude  <- c("S25_B", "S38_B", "S40_B", "S44_B", "S60_B")   # RNA-seq PCA outliers

# ---- shared setup: PAM50 centroids (Entrez -> symbol) and TCGA-BRCA reference ----
data(pam50)
cent <- pam50$centroids
rownames(cent) <- mapIds(org.Hs.eg.db,
                         keys = as.character(pam50$centroids.map$EntrezGene.ID),
                         column = "SYMBOL", keytype = "ENTREZID")
cent <- cent[!is.na(rownames(cent)) & !duplicated(rownames(cent)), ]

se       <- readRDS(file.path(data_dir, "tcga_brca_reference.rds"))
tcga_log <- log2(assay(se, "tpm_unstrand") + 1)
rownames(tcga_log) <- rowData(se)$gene_name
tcga_log <- tcga_log[!duplicated(rownames(tcga_log)) & rownames(tcga_log) != "", ]


# ===========================================================================
# PART 1 - reported PAM50 calls (full-reference-median centering)
# ===========================================================================
tpm     <- read.table(file.path(data_dir, "rnaseq_tpm.txt"), header = TRUE, row.names = 1)
tpm     <- tpm[, !colnames(tpm) %in% exclude]
log_tpm <- log2(as.matrix(tpm) + 1)

full_med <- apply(tcga_log, 1, median)                    # full-cohort reference median

g <- Reduce(intersect, list(rownames(cent), rownames(log_tpm), names(full_med)))
if (length(g) < 45)
  stop(sprintf("only %d/50 PAM50 genes matched - check gene identifiers", length(g)))

x       <- sweep(log_tpm[g, ], 1, full_med[g], "-")
cor_mat <- cor(x, cent[g, ], method = "spearman")
calls   <- data.frame(sample = colnames(x),
                      PAM50   = colnames(cor_mat)[max.col(cor_mat, ties.method = "first")],
                      max_cor = round(apply(cor_mat, 1, max), 3))
write.csv(calls, file.path(out_dir, "PAM50_calls.csv"), row.names = FALSE)

calls$timepoint <- sub(".*_(B|OT1|OT2)$", "\\1", calls$sample)
cat("PAM50 call proportions by timepoint:\n")
print(round(prop.table(table(calls$timepoint, calls$PAM50), 1), 3))


# ===========================================================================
# PART 2 - centering-scheme benchmark (30 seeds)
# ===========================================================================
sub <- TCGAquery_subtype("brca")
lab <- sub$BRCA_Subtype_PAM50[match(substr(colnames(se), 1, 12), sub$patient)]
ok  <- !is.na(lab) & lab %in% c("Basal", "Her2", "LumA", "LumB", "Normal")
ref <- tcga_log[, ok]; lab <- lab[ok]

run_seed <- function(seed) {
  set.seed(seed)
  train <- sample(ncol(ref), floor(ncol(ref) / 2))
  test  <- setdiff(seq_len(ncol(ref)), train)

  med  <- apply(ref[, train], 1, median)
  qpct <- rowMeans(ref[, train[lab[train] == "Basal"]] <= med)

  # basal-enriched test set (~9:1 basal:non-basal), matching the study composition
  basal_te <- test[lab[test] == "Basal"]
  mix   <- c(basal_te, sample(test[lab[test] != "Basal"], round(length(basal_te) / 9)))
  gg    <- Reduce(intersect, list(rownames(cent), rownames(ref), names(qpct)))
  M     <- ref[gg, mix]; truth <- lab[mix]; nb <- truth != "Basal"

  cA <- sapply(gg, function(gene) quantile(M[gene, ], probs = qpct[gene], names = FALSE))
  pA <- colnames(cent)[max.col(cor(sweep(M, 1, cA, "-"),      cent[gg, ], method = "spearman"))]
  pB <- colnames(cent)[max.col(cor(sweep(M, 1, med[gg], "-"), cent[gg, ], method = "spearman"))]

  data.frame(seed = seed,
             acc_percentile = mean(pA == truth), acc_refmedian = mean(pB == truth),
             nonBasal_to_Basal_percentile = mean(pA[nb] == "Basal"),
             nonBasal_to_Basal_refmedian  = mean(pB[nb] == "Basal"))
}

bench <- do.call(rbind, lapply(1:30, run_seed))
write.csv(bench, file.path(out_dir, "PAM50_centering_benchmark.csv"), row.names = FALSE)
cat("\nmean across 30 seeds (percentile vs reference-median):\n")
print(round(colMeans(bench[-1]), 3))
