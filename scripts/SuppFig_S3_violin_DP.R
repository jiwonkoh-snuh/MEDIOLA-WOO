# Supplementary Figure S3
# Sequential digital-pathology tissue composition (tumor / stromal / necrosis /
# inflammatory area) over serial biopsies (B / OT1 / OT2), stratified by metabolic
# response, for the MEDIOLA WOO trial. Two stratifications are shown:
#   - metabolic response at OT1 (PET, 2 weeks)
#   - metabolic response at OT2 (PET, 4 weeks)
#
# The plotted values are the same per-sample tissue-composition measurements as in
# Figure 2; here they are re-stratified by metabolic response. All comparisons in
# Supplementary Figure S3 form a single analysis: Benjamini-Hochberg is applied
# once across the complete set of 48 paired Wilcoxon tests (2 metabolic-response
# endpoints x 4 compartments x 2 response groups x 3 timepoint contrasts).
#
# Input  (place under data/):
#   tissue_composition.csv
# Output: figS3_mResOT1_*.pdf, figS3_mResOT2_*.pdf + the p-value tables.

library(tidyverse)
library(ggpubr)
library(ggh4x)

data_dir <- "data"
out_dir  <- "."

df <- read.csv(file.path(data_dir, "tissue_composition.csv"))
df$Class <- factor(df$Class, levels = c("B", "OT1", "OT2"))
# metabolic response at each on-treatment timepoint (Y = responder)
df$mRes_OT1 <- factor(df$PET2wk, levels = c("Y", "N"), labels = c("mRes", "non-mRes"))
df$mRes_OT2 <- factor(df$PET4wk, levels = c("Y", "N"), labels = c("mRes", "non-mRes"))

comp  <- c("Tumor", "Stroma", "Necrosis", "Inflammatory")
lab   <- c(Tumor = "Tumor area", Stroma = "Stromal area",
           Necrosis = "Necrosis area", Inflammatory = "Inflammatory infiltrate area")
pairs <- list(c("B", "OT1"), c("B", "OT2"), c("OT1", "OT2"))

tp_col   <- c(B = "#90728F", OT1 = "#9D983D", OT2 = "#E15759")
mres_col <- c("mRes" = "#8FAADC", "non-mRes" = "#F4707C")
ypos     <- c("B_OT1" = 100, "OT1_OT2" = 108, "B_OT2" = 116)

fmt_p <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                                       formatC(p, format = "f", digits = 3))

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
finalize_stats <- function(res) {
  res$p.adj      <- p.adjust(res$p, method = "BH")
  res$y.position <- ypos[paste(res$group1, res$group2, sep = "_")]
  res$label      <- fmt_p(res$p.adj)
  res$sig        <- !is.na(res$p.adj) & res$p.adj < 0.05
  res
}
add_pvals <- function(p, st) {
  s  <- st[st$sig, ]; ns <- st[!st$sig, ]
  if (nrow(s) > 0)
    p <- p + stat_pvalue_manual(s, label = "label", tip.length = 0.01, size = 3.3, color = "black")
  if (nrow(ns) > 0)
    p <- p + stat_pvalue_manual(ns, label = "label", tip.length = 0.01, size = 3.3,
                                color = "grey60", bracket.size = 0.4)
  p
}

# All S3 comparisons form a single analysis: compute raw paired-Wilcoxon p-values for
# both metabolic-response endpoints, then apply Benjamini-Hochberg once across the
# full set of 48 tests.
endpoints <- c(mRes_OT1 = "mResOT1", mRes_OT2 = "mResOT2")

raw_for_endpoint <- function(group_col) {
  d <- df[!is.na(df[[group_col]]), ]
  d$response_group <- factor(d[[group_col]], levels = c("mRes", "non-mRes"))
  do.call(rbind, lapply(c("mRes", "non-mRes"), function(g) {
    s <- compute_raw(d[d$response_group == g, ], comp, pairs)
    s$response_group <- g; s$endpoint <- group_col; s
  }))
}

stats_all <- do.call(rbind, lapply(names(endpoints), raw_for_endpoint))
stats_all$p.adj      <- p.adjust(stats_all$p, method = "BH")   # BH across the 48-test set
stats_all$y.position <- ypos[paste(stats_all$group1, stats_all$group2, sep = "_")]
stats_all$label      <- fmt_p(stats_all$p.adj)
stats_all$sig        <- !is.na(stats_all$p.adj) & stats_all$p.adj < 0.05
stats_all$response_group <- factor(stats_all$response_group, levels = c("mRes", "non-mRes"))

# Plot each endpoint using the pooled-adjusted p-values
plot_endpoint <- function(group_col) {
  tag <- endpoints[[group_col]]
  d <- df[!is.na(df[[group_col]]), ]
  d$response_group <- factor(d[[group_col]], levels = c("mRes", "non-mRes"))
  st_ep <- stats_all[stats_all$endpoint == group_col, ]
  write.csv(st_ep, file.path(out_dir, sprintf("figS3_%s_pvalues.csv", tag)), row.names = FALSE)

  for (v in comp) {
    st <- st_ep[st_ep$Compartment == v, ]
    p <- ggplot(d, aes(Class, .data[[v]], fill = Class)) +
      geom_violin(trim = TRUE, alpha = 0.7, color = NA, width = 0.65, scale = "width") +
      geom_line(aes(group = Sample), color = "grey", linetype = "dotted") +
      geom_point(position = position_jitter(width = 0.1), alpha = 0.2) +
      facet_wrap2(~ response_group,
                  strip = strip_themed(background_x = elem_list_rect(fill = mres_col[c("mRes", "non-mRes")]),
                                       text_x = elem_list_text(color = "black", face = "bold"))) +
      scale_y_continuous(limits = c(0, 120), breaks = c(0, 25, 50, 75, 100),
                         labels = scales::percent_format(scale = 1)) +
      scale_fill_manual(values = tp_col) +
      labs(x = "Timepoint", y = lab[[v]]) +
      theme_classic() +
      theme(panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
            panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
            text = element_text(size = 14), legend.position = "none")
    p <- add_pvals(p, st)
    ggsave(file.path(out_dir, sprintf("figS3_%s_%s.pdf", tag, v)), p, width = 5.0, height = 3.7)
  }
}

for (gc in names(endpoints)) plot_endpoint(gc)   # OT1 (PET 2 weeks), OT2 (PET 4 weeks)
