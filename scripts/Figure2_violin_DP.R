# Figure 2
# Digital-pathology tissue composition (tumor / stromal / necrosis / inflammatory
# area) over serial biopsies (B / OT1 / OT2) for the MEDIOLA WOO trial.
#   Panel 2A: whole cohort, longitudinal change per compartment
#             (12-test set: 4 compartments x 3 timepoint contrasts, n = 54)
#   Panel 2B: same values stratified by pCR
#             (24-test set: 4 compartments x 3 contrasts x 2 response groups)
# Benjamini-Hochberg correction is applied once across each panel's set of paired
# Wilcoxon tests. The time-by-response interaction model supporting Supplementary
# Table S4 is in Supplementary_Table_S4.R.
#
# Input  (place under data/):
#   tissue_composition.csv
# Output: fig2A_*.pdf, fig2B_*.pdf + the p-value tables.

library(tidyverse)
library(ggpubr)
library(ggh4x)

data_dir <- "data"
out_dir  <- "."

df <- read.csv(file.path(data_dir, "tissue_composition.csv"))
df$Class <- factor(df$Class, levels = c("B", "OT1", "OT2"))
df$pCR_f <- factor(df$pCR, levels = c("PCR", "nonPCR"), labels = c("pCR", "non-pCR"))

comp  <- c("Tumor", "Stroma", "Necrosis", "Inflammatory")
lab   <- c(Tumor = "Tumor area", Stroma = "Stromal area",
           Necrosis = "Necrosis area", Inflammatory = "Inflammatory infiltrate area")
pairs <- list(c("B", "OT1"), c("B", "OT2"), c("OT1", "OT2"))

tp_col  <- c(B = "#90728F", OT1 = "#9D983D", OT2 = "#E15759")
pcr_col <- c("pCR" = "#8FAADC", "non-pCR" = "#F4707C")
ypos    <- c("B_OT1" = 100, "OT1_OT2" = 108, "B_OT2" = 116)

fmt_p <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                                       formatC(p, format = "f", digits = 3))

# raw paired-Wilcoxon p-values on complete pairs
compute_raw <- function(data, comps, pairs) {
  do.call(rbind, lapply(comps, function(v) {
    w <- data %>%
      transmute(Sample, Class, value = .data[[v]]) %>%
      pivot_wider(names_from = Class, values_from = value)
    do.call(rbind, lapply(pairs, function(pr) {
      ok <- complete.cases(w[[pr[1]]], w[[pr[2]]])
      p  <- wilcox.test(w[[pr[1]]][ok], w[[pr[2]]][ok], paired = TRUE, exact = FALSE)$p.value
      data.frame(Compartment = v, group1 = pr[1], group2 = pr[2], n = sum(ok), p = p)
    }))
  }))
}

# BH across every row bound together -> the set passed in defines the correction
finalize_stats <- function(res) {
  res$p.adj      <- p.adjust(res$p, method = "BH")
  res$y.position <- ypos[paste(res$group1, res$group2, sep = "_")]
  res$label      <- fmt_p(res$p.adj)
  res$sig        <- !is.na(res$p.adj) & res$p.adj < 0.05
  res
}
build_stats <- function(data, comps, pairs) finalize_stats(compute_raw(data, comps, pairs))

add_pvals <- function(p, st) {
  s  <- st[st$sig, ]; ns <- st[!st$sig, ]
  if (nrow(s) > 0)
    p <- p + stat_pvalue_manual(s, label = "label", tip.length = 0.01, size = 3.3, color = "black")
  if (nrow(ns) > 0)
    p <- p + stat_pvalue_manual(ns, label = "label", tip.length = 0.01, size = 3.3,
                                color = "grey60", bracket.size = 0.4)
  p
}

base_violin <- function(v) {
  ggplot(df, aes(Class, .data[[v]], fill = Class)) +
    geom_violin(trim = TRUE, alpha = 0.7, color = NA, width = 0.65, scale = "width") +
    geom_line(aes(group = Sample), color = "grey", linetype = "dotted") +
    geom_point(position = position_jitter(width = 0.1), alpha = 0.2) +
    scale_y_continuous(limits = c(0, 120), breaks = c(0, 25, 50, 75, 100),
                       labels = scales::percent_format(scale = 1)) +
    scale_fill_manual(values = tp_col) +
    labs(x = "Timepoint", y = lab[[v]]) +
    theme_classic() +
    theme(panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          text = element_text(size = 14), legend.position = "none")
}


# ---- Figure 2A : whole cohort (12-test set) ----
statsA <- build_stats(df, comp, pairs)
write.csv(statsA, file.path(out_dir, "fig2A_pvalues.csv"), row.names = FALSE)

for (v in comp) {
  p <- add_pvals(base_violin(v), statsA[statsA$Compartment == v, ])
  ggsave(file.path(out_dir, sprintf("fig2A_%s.pdf", v)), p, width = 3.3, height = 3.7)
}


# ---- Figure 2B : stratified by pCR (24-test set, pooled across both groups) ----
statsB <- do.call(rbind, lapply(levels(df$pCR_f), function(g) {
  s <- compute_raw(df[df$pCR_f == g, ], comp, pairs); s$pCR_f <- g; s
}))
statsB <- finalize_stats(statsB)
write.csv(statsB, file.path(out_dir, "fig2B_pvalues.csv"), row.names = FALSE)

for (v in comp) {
  st <- statsB[statsB$Compartment == v, ]
  p <- base_violin(v) +
    facet_wrap2(~ factor(pCR_f, levels = c("pCR", "non-pCR")),
                strip = strip_themed(background_x = elem_list_rect(fill = pcr_col[c("pCR", "non-pCR")]),
                                     text_x = elem_list_text(color = "black", face = "bold")))
  p <- add_pvals(p, st)
  ggsave(file.path(out_dir, sprintf("fig2B_%s.pdf", v)), p, width = 5.0, height = 3.7)
}

