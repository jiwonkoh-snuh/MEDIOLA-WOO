# Figure 3  (RNA-seq pathway analysis, purity-adjusted)
#   3A: GSVA heatmap of Hallmark gene sets across timepoints (whole cohort)
#   3B: GSVA violins for the highlighted gene sets, stratified by pCR, with
#       purity-adjusted mixed-model timepoint contrasts (non-significant brackets greyed)
#   3C: GSEA balloon plots (pCR and recurrence contrasts at each timepoint)
#
# Inputs  (place under data/):
#   rnaseq_tpm.txt             gene x sample TPM matrix (RNA-seq; see GEO)
#   sample_metadata.csv               per-patient metadata (pCR, PET2wk, ...)
#   tissue_composition.csv            AI tissue composition (Tumor area = purity covariate)
#   gsea_hallmark_results.csv precomputed GSEA results (for panel 3C)
# Outputs: GSVA_Fig3A.pdf, GSVA_Fig3B.pdf, GSEA_Fig3C.pdf + score/p-value CSVs.

library(GSVA)
library(msigdbr)
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(ggh4x)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(purrr)
library(data.table)
library(patchwork)
library(scales)
library(grid)
library(lmerTest)
library(emmeans)

data_dir <- "data"
out_dir  <- "."

# five baseline RNA-seq libraries flagged as global expression outliers
exclude    <- c("S25_B", "S38_B", "S40_B", "S44_B", "S60_B")
timepoints <- c("B", "OT1", "OT2")
tp_col     <- c(B = "#90728F", OT1 = "#9D983D", OT2 = "#E15759")
pcr_col    <- c("PCR" = "#8FAADC", "nonPCR" = "#F4707C")

tpm      <- read.table(file.path(data_dir, "rnaseq_tpm.txt"), header = TRUE, row.names = 1)
log_tpm  <- as.matrix(log2(tpm + 1))
metadata <- read.csv(file.path(data_dir, "sample_metadata.csv"), check.names = FALSE)

# Tumor-area proportion (AI segmentation) used as the purity covariate in 3B
ld <- read.csv(file.path(data_dir, "tissue_composition.csv"))
ld$sample  <- paste(trimws(ld$Sample), trimws(ld$Class), sep = "_")
purity_vec <- setNames(ld$Tumor, ld$sample)


# GSVA on Hallmark gene sets
hallmark      <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)

gsvaPar     <- gsvaParam(log_tpm, hallmark_list, kcdf = "Gaussian", maxDiff = TRUE)
gsva_result <- gsva(gsvaPar)
write.csv(gsva_result, file.path(out_dir, "gsva_hallmark_scores.csv"), row.names = TRUE)


# 3A : heatmap
gsva_hm <- gsva_result
rownames(gsva_hm) <- sub("^HALLMARK_", "", rownames(gsva_hm))
gsva_hm <- gsva_hm[, !colnames(gsva_hm) %in% exclude]

sample_info <- data.frame(ID2 = colnames(gsva_hm),
                          timepoint = sub(".*_", "", colnames(gsva_hm)),
                          stringsAsFactors = FALSE)

# keep gene sets with SD > 0.2 within at least one timepoint
keep_rows <- apply(gsva_hm, 1, function(row) {
  any(sapply(timepoints, function(tp)
    sd(row[sample_info$ID2[sample_info$timepoint == tp]], na.rm = TRUE)) > 0.2)
})
gsva_filtered <- gsva_hm[keep_rows, , drop = FALSE]

# order samples within each timepoint by Ward clustering, timepoints kept together
sample_clusters <- unlist(lapply(timepoints, function(tp) {
  cols <- sample_info$ID2[sample_info$timepoint == tp]
  hc   <- hclust(dist(t(gsva_filtered[, cols, drop = FALSE])), method = "ward.D2")
  cols[hc$order]
}))

annotation_col <- data.frame(Timepoint = sample_info$timepoint, row.names = sample_info$ID2)

pdf(file.path(out_dir, "GSVA_Fig3A.pdf"), width = 14, height = 8)
pheatmap(gsva_filtered[, sample_clusters],
         annotation_col = annotation_col[sample_clusters, , drop = FALSE],
         annotation_colors = list(Timepoint = tp_col),
         cluster_cols = FALSE, cluster_rows = TRUE, scale = "row",
         color = colorRampPalette(c("saddlebrown", "white", "darkcyan"))(50),
         fontsize_row = 6, show_colnames = FALSE, main = "GSVA scores across timepoints")
dev.off()


# 3B : purity-adjusted timepoint contrasts, stratified by pCR
gsva_long <- as.data.frame(gsva_result[, !colnames(gsva_result) %in% exclude]) %>%
  rownames_to_column("pathway") %>%
  pivot_longer(-pathway, names_to = "sample", values_to = "gsva_score") %>%
  mutate(timepoint  = str_extract(sample, "(B|OT1|OT2)$"),
         patient_id = str_replace(sample, "(_B|_OT1|_OT2)$", ""),
         timepoint  = factor(timepoint, levels = timepoints))

metadata_unique <- metadata %>% distinct(ID, .keep_all = TRUE)

df_all <- gsva_long %>%
  filter(pathway %in% paste0("HALLMARK_", rownames(gsva_filtered))) %>%
  left_join(metadata_unique %>% dplyr::select(any_of(c("ID", "pCR", "PET2wk"))),
            by = c("patient_id" = "ID")) %>%
  mutate(Group = timepoint, score = gsva_score, Sample = patient_id,
         pCR = factor(pCR, levels = c("PCR", "nonPCR")),
         pathway_label = sub("^HALLMARK_", "", pathway),
         purity = purity_vec[sample]) %>%
  filter(!is.na(pCR), !is.na(purity))

# per gene set x pCR group: mixed model with purity covariate, EMM timepoint contrasts
get_pvals <- function(data) {
  comps <- list(c("B", "OT1"), c("B", "OT2"), c("OT1", "OT2"))
  fit <- tryCatch(suppressMessages(lmer(score ~ Group + purity + (1 | Sample), data = data)),
                  error = function(e) NULL)
  pc  <- if (!is.null(fit)) as.data.frame(pairs(emmeans(fit, ~ Group), adjust = "none")) else NULL
  bind_rows(lapply(comps, function(cc) {
    p <- NA_real_
    if (!is.null(pc)) {
      hit <- pc[pc$contrast %in% c(paste(cc[1], "-", cc[2]), paste(cc[2], "-", cc[1])), ]
      if (nrow(hit)) p <- hit$p.value[1]
    }
    data.frame(group1 = cc[1], group2 = cc[2], p = p)
  }))
}

# BH across the full 228-test set (38 sets x 3 contrasts x 2 groups)
pvals_all <- df_all %>%
  group_by(pathway_label, pCR) %>%
  group_modify(~ get_pvals(.x)) %>%
  ungroup() %>%
  mutate(p_adj_global = p.adjust(p, method = "BH"))
write.csv(pvals_all, file.path(out_dir, "gsva_pcr_contrasts.csv"), row.names = FALSE)

fmt_p   <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                                         formatC(p, format = "f", digits = 3))
relabel <- function(x) recode_factor(x, "DNA_REPAIR" = "DNA repair",
                                        "E2F_TARGETS" = "E2F targets",
                                        "G2M_CHECKPOINT" = "G2M checkpoint")
wanted <- c("DNA_REPAIR", "E2F_TARGETS", "G2M_CHECKPOINT")

df_plot <- df_all %>% filter(pathway_label %in% wanted) %>% mutate(pathway_clean = relabel(pathway_label))

annot <- pvals_all %>%
  filter(pathway_label %in% wanted) %>%
  mutate(pathway_clean = relabel(pathway_label),
         group1 = as.character(group1), group2 = as.character(group2),
         label = fmt_p(p_adj_global),
         sig = !is.na(p_adj_global) & p_adj_global < 0.05,
         y.position = case_when(group1 == "B"   & group2 == "OT1" ~ 0.75,
                                group1 == "B"   & group2 == "OT2" ~ 0.95,
                                group1 == "OT1" & group2 == "OT2" ~ 1.15))

p_3b <- ggplot(df_plot, aes(Group, score, fill = Group)) +
  geom_violin(trim = TRUE, alpha = 0.7, color = NA, width = 0.65, scale = "width") +
  geom_line(aes(group = Sample), color = "grey75", linetype = "dotted") +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.25, size = 1.8, color = "black") +
  labs(x = "Timepoint", y = "GSVA score") +
  scale_fill_manual(values = tp_col) +
  facet_grid2(rows = vars(pathway_clean), cols = vars(pCR),
              strip = strip_themed(background_x = elem_list_rect(fill = unname(pcr_col), color = c(NA, NA)),
                                   text_x = elem_list_text(colour = c("black", "black"), face = c("bold", "bold")))) +
  theme_classic() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
        text = element_text(size = 14), legend.position = "none")

if (any(annot$sig))
  p_3b <- p_3b + stat_pvalue_manual(filter(annot, sig), label = "label",
                                    xmin = "group1", xmax = "group2", y.position = "y.position",
                                    tip.length = 0.01, size = 3, color = "black")
if (any(!annot$sig))
  p_3b <- p_3b + stat_pvalue_manual(filter(annot, !sig), label = "label",
                                    xmin = "group1", xmax = "group2", y.position = "y.position",
                                    tip.length = 0.01, size = 3, color = "grey60", bracket.size = 0.4)

pdf(file.path(out_dir, "GSVA_Fig3B.pdf"), width = 9, height = 9); print(p_3b); dev.off()


# 3C : GSEA balloon plots (purity-adjusted model; padj < 0.01 shown)
gsea <- fread(file.path(data_dir, "gsea_hallmark_results.csv"))
gsea <- gsea[gsea$model == "purityadj", ]

keys <- list(list(c = "pCR", t = "B"), list(c = "pCR", t = "OT1"), list(c = "pCR", t = "OT2"),
             list(c = "Recur", t = "B"), list(c = "Recur", t = "OT1"), list(c = "Recur", t = "OT2"))
tp_lab      <- c(B = "Baseline", OT1 = "OT1", OT2 = "OT2")
panel_title <- function(k) {
  if (k$c == "pCR") sprintf("PCR vs non-PCR at %s", tp_lab[k$t])
  else              sprintf("Recur vs non-recur at %s", tp_lab[k$t])
}
dir_labs <- function(k) {
  if (k$c == "pCR") list(left = "\u2190 non-PCR", right = "PCR \u2192")
  else              list(left = "\u2190 non-recur", right = "Recur \u2192")
}

theme_natcomm <- function(base = 11) {
  theme_classic(base_size = base) %+replace%
    theme(panel.grid.major.y = element_line(color = "#E6E6E6", linewidth = 0.35),
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          axis.line = element_line(color = "black", linewidth = 0.5),
          axis.ticks = element_line(color = "black", linewidth = 0.4), axis.ticks.length = unit(2, "pt"),
          plot.title = element_text(size = base + 2, hjust = 0.5), axis.title.y = element_blank(),
          axis.text = element_text(color = "black"), strip.background = element_blank())
}

all_sig <- gsea %>%
  filter(is.finite(padj), !is.na(NES), padj < 0.01) %>%
  filter(!str_detect(pathway, regex("UV|SPERMATOGENESIS|MYOGENESIS", ignore_case = TRUE))) %>%
  mutate(neglog10 = -log10(padj))
if (nrow(all_sig) == 0) stop("No significant pathways (padj < 0.01) in the selected model.")

size_lim <- range(all_sig$neglog10, na.rm = TRUE)

plots <- map(keys, function(k) {
  sig <- all_sig %>% filter(contrast == k$c, timepoint == k$t) %>%
    arrange(desc(NES), padj) %>% mutate(pathway = factor(pathway, levels = rev(unique(pathway))))
  if (nrow(sig) == 0) return(ggplot() + theme_void() + ggtitle(panel_title(k)))
  dl <- dir_labs(k)
  ggplot(sig, aes(NES, pathway, size = neglog10, color = NES)) +
    geom_point(alpha = 0.9, stroke = 0) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    scale_size_continuous(name = expression(-log[10]("adjP")), range = c(2.2, 10), limits = size_lim) +
    scale_color_gradient2(low = "saddlebrown", mid = "white", high = "darkcyan",
                          midpoint = 0, limits = c(-3, 3), name = "NES", oob = squish) +
    scale_y_discrete(labels = label_wrap(28)) +
    labs(x = "Normalized enrichment score (NES)", title = panel_title(k)) +
    xlim(c(-4, 4)) +
    annotate("text", x = -4, y = Inf, label = dl$left,  hjust = 0, vjust = -0.5, size = 3.2) +
    annotate("text", x =  4, y = Inf, label = dl$right, hjust = 1, vjust = -0.5, size = 3.2) +
    coord_cartesian(clip = "off") +
    theme_natcomm(11) +
    theme(axis.text.y = element_text(size = max(6, min(11, 18 - 0.05 * n_distinct(sig$pathway)))),
          plot.margin = margin(14, 6, 6, 6))
})

final_plot <- wrap_plots(plots, ncol = 3, nrow = 2, guides = "collect") & theme(legend.position = "right")
ggsave(file.path(out_dir, "GSEA_Fig3C.pdf"), final_plot, width = 18, height = 10, useDingbats = FALSE)
