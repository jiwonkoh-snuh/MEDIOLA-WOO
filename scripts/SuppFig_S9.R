# Supplementary Figure S9 — RNA-seq quality control and cohort assembly
# MEDIOLA WOO trial
#
# QC values (mapping rate, library size, PC1/PC2) are read from Supplementary
# Table S5. A library is excluded if PC1 is a robust outlier (MAD z-score |z| > 15)
# or the library size is below 20 million mapped reads.
#
# Inputs (place under data/):
#   Supplementary_Table_S5.csv   Sample, Timepoint, Mapping rate (%),
#                                Library size (mapped reads), PC1, PC2, QC status
#   rnaseq_counts.txt            gene x sample raw count matrix (panel E re-PCA)
# Output: FigureS9_RNAseq_QC.pdf

library(edgeR)
library(tidyverse)
library(ggrepel)
library(patchwork)

data_dir <- "data"
out_dir  <- "."

## --- read QC from the deposited Supplementary Table S5 -----------------------
s5 <- read.csv(file.path(data_dir, "Supplementary_Table_S5.csv"), check.names = FALSE)
qc <- data.frame(
  sample       = s5[["Sample"]],
  timepoint    = s5[["Timepoint"]],
  mapping_rate = as.numeric(s5[["Mapping rate (%)"]]),
  lib_size     = as.numeric(s5[["Library size (mapped reads)"]]),
  PC1          = as.numeric(s5[["PC1"]]),
  PC2          = as.numeric(s5[["PC2"]]),
  stringsAsFactors = FALSE
)

## exclusion rule
Z_THRESH   <- 15
DEPTH_MIN  <- 20e6
pc1_med <- median(qc$PC1)
pc1_mad <- mad(qc$PC1)                              # R mad() includes the 1.4826 constant
qc <- qc %>%
  mutate(robz_pc1 = (PC1 - pc1_med) / pc1_mad,
         excluded = abs(robz_pc1) > Z_THRESH | lib_size < DEPTH_MIN,
         status   = ifelse(excluded, "Excluded", "Retained"),
         Timepoint = factor(timepoint, levels = c("B","OT1","OT2")),
         lib_M = lib_size / 1e6)

# expected exclusions
target <- c("S25_B","S38_B","S40_B","S44_B","S60_B")
stopifnot(setequal(qc$sample[qc$excluded], target))
cat("Excluded (rule-derived):", paste(sort(qc$sample[qc$excluded]), collapse=", "),
    "\nRetained:", sum(!qc$excluded), "of", nrow(qc), "\n")

# PC1 cutoff lines (for the PCA panel) corresponding to |z| > 15
pc1_lo <- pc1_med - Z_THRESH * pc1_mad
pc1_hi <- pc1_med + Z_THRESH * pc1_mad

pal_time <- c("B"="#90728F","OT1"="#9D983D","OT2"="#E15759")
base_theme <- theme_classic(base_size = 11) +
  theme(panel.grid.major.y = element_line(color = "grey90"),
        panel.grid.minor   = element_blank(),
        plot.title = element_text(face = "bold", size = 11))

## (a) cohort-assembly flow diagram -------------------------------------------
flow <- data.frame(
  x = 1, y = c(3,2,1),
  lab = c("160 sequenced libraries passing initial QC\n(baseline 54 \u00b7 OT1 54 \u00b7 OT2 52)",
          "Exclude 5 baseline expression outliers\n(S25_B, S38_B, S40_B, S44_B, S60_B)",
          "155 analytic samples\n(baseline 49 \u00b7 OT1 54 \u00b7 OT2 52)"))
seg <- data.frame(x = 1, xend = 1, y = c(2.72,1.72), yend = c(2.28,1.28))
p_flow <- ggplot(flow, aes(x, y)) +
  geom_label(aes(label = lab), size = 3.0, label.size = 0.3,
             fill = c("grey96","#FBE4E4","#E7EEF7"), lineheight = 0.95) +
  geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.4) +
  xlim(0.4,1.6) + ylim(0.6,3.4) +
  labs(title = "A") + theme_void() +
  theme(plot.title = element_text(face = "bold", size = 12))

## (b) library size (dashed line = 20M-read depth threshold) -------------------
p_lib <- ggplot(qc, aes(Timepoint, lib_M)) +
  geom_violin(aes(fill = Timepoint), trim = TRUE, alpha = 0.6, color = NA, scale = "width") +
  geom_hline(yintercept = DEPTH_MIN/1e6, linetype = "dashed", color = "grey40") +
  geom_jitter(aes(color = status), width = 0.12, size = 1.6, alpha = 0.7) +
  geom_text_repel(data = subset(qc, excluded),
                  aes(label = sample), size = 2.6, color = "#C0392B",
                  max.overlaps = Inf, box.padding = 0.4) +
  scale_fill_manual(values = pal_time, guide = "none") +
  scale_color_manual(values = c("Excluded"="#C0392B","Retained"="grey35")) +
  labs(title = "B", x = "Timepoint", y = "Library size (million reads)", color = NULL) +
  base_theme

## (c) Salmon mapping rate -----------------------------------------------------
p_map <- ggplot(qc, aes(Timepoint, mapping_rate)) +
  geom_violin(aes(fill = Timepoint), trim = TRUE, alpha = 0.6, color = NA, scale = "width") +
  geom_jitter(aes(color = status), width = 0.12, size = 1.6, alpha = 0.7) +
  geom_text_repel(data = subset(qc, excluded),
                  aes(label = sample), size = 2.6, color = "#C0392B",
                  max.overlaps = Inf, box.padding = 0.4) +
  scale_fill_manual(values = pal_time, guide = "none") +
  scale_color_manual(values = c("Excluded"="#C0392B","Retained"="grey35")) +
  labs(title = "C", x = "Timepoint", y = "Salmon mapping rate (%)", color = NULL) +
  base_theme

## (d) PCA before exclusion (dashed lines = PC1 robust-z cutoffs) ---------------
p_pca <- ggplot(qc, aes(PC1, PC2)) +
  geom_vline(xintercept = c(pc1_lo, pc1_hi), linetype = "dashed", color = "grey60") +
  geom_point(aes(color = Timepoint, shape = status), size = 2, alpha = 0.8) +
  geom_text_repel(data = subset(qc, excluded),
                  aes(label = sample), size = 2.6, color = "#C0392B",
                  max.overlaps = Inf, box.padding = 0.5) +
  scale_color_manual(values = pal_time) +
  scale_shape_manual(values = c("Retained"=16, "Excluded"=4)) +
  labs(title = "D", x = "PC1", y = "PC2", color = "Timepoint", shape = NULL) +
  base_theme

## (e) PCA after exclusion (recomputed on the 155 retained libraries) -----------
# Rebuild log-CPM after dropping the excluded columns, so the PCA is
# recomputed on the retained libraries.
counts <- read.table(file.path(data_dir, "rnaseq_counts.txt"), header = TRUE,
                     row.names = 1, check.names = FALSE)
counts <- as.matrix(round(counts))
counts <- counts[, !colnames(counts) %in% target]      # 155 retained samples
stopifnot(ncol(counts) == 155)

dge5    <- DGEList(counts = counts)
tp5     <- factor(sub(".*_", "", colnames(counts)), levels = c("B","OT1","OT2"))
dge5    <- dge5[filterByExpr(dge5, group = tp5), , keep.lib.sizes = FALSE]
dge5    <- calcNormFactors(dge5)
logCPM5 <- cpm(dge5, log = TRUE)
pca5    <- prcomp(t(logCPM5), center = TRUE, scale. = TRUE)
var5    <- round(100 * summary(pca5)$importance[2, 1:2], 1)

pca5_df <- data.frame(PC1 = pca5$x[,1], PC2 = pca5$x[,2], Timepoint = tp5)
p_pca5 <- ggplot(pca5_df, aes(PC1, PC2, color = Timepoint)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_manual(values = pal_time) +
  labs(title = "E", x = paste0("PC1 (", var5[1], "%)"),
       y = paste0("PC2 (", var5[2], "%)"), color = "Timepoint") +
  base_theme

## assemble --------------------------------------------------------------------
design <- (p_flow | p_lib  | p_map) /
          (p_pca  | p_pca5 | plot_spacer()) +
          plot_layout(heights = c(1, 1.05))
pdf(file.path(out_dir, "FigureS9_RNAseq_QC.pdf"), width = 16, height = 8)   # landscape
print(design)
dev.off()
