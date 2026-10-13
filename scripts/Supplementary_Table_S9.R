# Supplementary Table S9
# Gene set enrichment analysis (GSEA) comparing groups by pCR status and by
# recurrence status within each timepoint (MEDIOLA WOO trial). Differential
# expression is fit with voom/limma per timepoint and contrast, ranked by log
# fold change, and tested against the Hallmark gene sets with fgsea. Two models
# are run: a tumor-area-adjusted model (primary) and an unadjusted model
# (sensitivity). Benjamini-Hochberg-adjusted P values are returned by fgsea.
#
# The long output below (one row per gene set x timepoint x contrast x model) is
# the source for Supplementary Table S9 (pivoted to place the adjusted and
# unadjusted models side by side) and for Figure 3C (adjusted model, adjusted
# P < 0.01).
#
# Inputs  (place under data/):
#   rnaseq_counts.txt   gene x sample raw counts (RNA-seq; see GEO)
#   tissue_composition.csv     AI tissue composition (Tumor area), pCR, Recur
# Output: gsea_hallmark_results.csv

library(edgeR)
library(limma)
library(fgsea)
library(msigdbr)
library(dplyr)
library(tidyr)

data_dir <- "data"
out_dir  <- "."
exclude  <- c("S25_B", "S38_B", "S40_B", "S44_B", "S60_B")

counts <- as.matrix(read.table(file.path(data_dir, "rnaseq_counts.txt"),
                               header = TRUE, row.names = 1, check.names = FALSE))
counts <- counts[, !colnames(counts) %in% exclude]

ld <- read.csv(file.path(data_dir, "tissue_composition.csv"))
ld$sample     <- paste(trimws(ld$Sample), trimws(ld$Class), sep = "_")
meta          <- ld[, c("sample", "Class", "Tumor", "pCR", "Recur")]
rownames(meta) <- meta$sample

hall     <- tryCatch(msigdbr(species = "Homo sapiens", category = "H"),
                     error = function(e) msigdbr(species = "Homo sapiens", collection = "H"))
pathways <- split(hall$gene_symbol, hall$gs_name)

contrasts <- list(
  pCR   = list(col = "pCR",   ref = "nonPCR", pos = "PCR"),
  Recur = list(col = "Recur", ref = "N",      pos = "Y")
)

run_tp <- function(tp, cn, adjust) {
  cd <- contrasts[[cn]]
  s  <- colnames(counts)[sub(".*_", "", colnames(counts)) == tp]
  s  <- s[s %in% rownames(meta)]
  m  <- meta[s, , drop = FALSE]
  m$grp <- factor(trimws(m[[cd$col]]), levels = c(cd$ref, cd$pos))
  m  <- m[!is.na(m$grp) & !is.na(m$Tumor), ]; s <- rownames(m)
  if (nlevels(droplevels(m$grp)) < 2 || min(table(m$grp)) < 2) return(NULL)

  y    <- DGEList(counts[, s])
  keep <- filterByExpr(y, group = m$grp)
  y    <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])

  design <- if (adjust) model.matrix(~ grp + Tumor, m) else model.matrix(~ grp, m)
  fit    <- eBayes(lmFit(voom(y, design), design))
  tt     <- topTable(fit, coef = paste0("grp", cd$pos), number = Inf, sort.by = "none")

  ranks <- sort(setNames(tt$logFC, rownames(tt)))
  set.seed(1)
  fg <- fgsea(pathways, ranks, minSize = 15, maxSize = 500, eps = 0)
  data.frame(contrast = cn, timepoint = tp,
             pathway = sub("^HALLMARK_", "", fg$pathway),
             NES = fg$NES, padj = fg$padj,
             model = ifelse(adjust, "purityadj", "unadj"))
}

grid <- expand.grid(tp = c("B", "OT1", "OT2"), cn = names(contrasts),
                    adjust = c(FALSE, TRUE), stringsAsFactors = FALSE)
res  <- do.call(rbind, Map(run_tp, grid$tp, grid$cn, grid$adjust))

write.csv(res, file.path(out_dir, "gsea_hallmark_results.csv"), row.names = FALSE)
cat("rows:", nrow(res), "| models:", paste(unique(res$model), collapse = ", "), "\n")
