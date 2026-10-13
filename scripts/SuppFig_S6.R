# Supplementary Figure S6
# Temporal dynamics of CD8+ T-cell subsets in the tumor and stromal compartments
# over serial biopsies (B / OT1 / OT2), stratified by treatment response, for the
# MEDIOLA WOO trial.
#
# The figure is a single multiple-testing set: Benjamini-Hochberg is applied once
# across all 72 paired Wilcoxon tests (6 markers x 2 compartments x 2 response
# groups x 3 timepoint contrasts).
#
# Inputs  (place under data/):
#   mihc_density_P1.csv   CD8+PD1+ / CD8+CD244+ subsets
#   mihc_density_P2.csv   CD8 lineage / activation / memory subsets
# Output: SuppFig_S6.pdf

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

data_dir <- "data"
out_dir  <- "."
timepoints <- c("B", "OT1", "OT2")

read_density <- function(f) {
  d <- read.csv(file.path(data_dir, f))
  d <- d[d$Group %in% timepoints & d$Segment %in% c("Tumor", "Stroma"), ]
  d$Group   <- factor(d$Group,   levels = timepoints)
  d$Segment <- factor(d$Segment, levels = c("Tumor", "Stroma"))
  d
}
P1 <- read_density("mihc_density_P1.csv")
P2 <- read_density("mihc_density_P2.csv")

# marker key -> display label, and the source panel for each key
markers <- c(CD8_T = "CD8+ T", CD8_PD1 = "CD8+PD1+", CD8_CD244 = "CD8+CD244+",
             CD8_T_acti = "CD8+ T activated",
             CD8_T_acti_term = "Activated memory CD8+ T cell",
             CD8_T_memory = "CD8+ T memory")
source_of <- function(k) if (k %in% c("CD8_PD1", "CD8_CD244")) P1 else P2


# Paired Wilcoxon tests across the 72-test set
paired_p <- function(d, marker, pair) {
  w <- d %>%
    filter(Group %in% pair) %>%
    dplyr::select(ID, Group, all_of(marker)) %>%
    pivot_wider(names_from = Group, values_from = all_of(marker)) %>%
    drop_na()
  if (nrow(w) < 2) return(NA)
  wilcox.test(w[[pair[1]]], w[[pair[2]]], paired = TRUE, exact = FALSE)$p.value
}

contrasts <- list(c("B", "OT1"), c("B", "OT2"), c("OT1", "OT2"))
stats <- bind_rows(lapply(names(markers), function(k) {
  d0 <- source_of(k)
  bind_rows(lapply(c("Tumor", "Stroma"), function(comp) {
    bind_rows(lapply(c("pCR", "non-pCR"), function(resp) {
      d <- d0[d0$Segment == comp & d0$pCR == ifelse(resp == "pCR", 1, 0), ]
      data.frame(marker = k, Compartment = comp, Response = resp,
                 g1 = sapply(contrasts, `[`, 1), g2 = sapply(contrasts, `[`, 2),
                 p = sapply(contrasts, function(p) paired_p(d, k, p)),
                 stringsAsFactors = FALSE)
    }))
  }))
}))
stats$p_BH <- p.adjust(stats$p, method = "BH")

cat("Tests:", nrow(stats), "| significant (p_BH < 0.05):", sum(stats$p_BH < 0.05, na.rm = TRUE), "\n")


# One violin panel (a single marker x compartment x response combination)
fmt_p <- function(p) ifelse(is.na(p), "NA",
                     ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                                       formatC(p, format = "f", digits = 3)))
pal_time <- c(B = "#90728F", OT1 = "#9D983D", OT2 = "#E15759")

marker_panel <- function(resp, comp, key) {
  d <- source_of(key) %>%
    filter(Segment == comp, pCR == ifelse(resp == "pCR", 1, 0)) %>%
    dplyr::select(ID, Group, value = all_of(key))
  ymax <- suppressWarnings(max(d$value, na.rm = TRUE))
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1

  annot <- stats %>%
    filter(marker == key, Compartment == comp, Response == resp) %>%
    rename(group1 = g1, group2 = g2) %>%
    arrange(match(paste0(group1, group2), c("BOT1", "OT1OT2", "BOT2"))) %>%
    mutate(label = paste0("P=", fmt_p(p_BH)),
           lab_color = ifelse(!is.na(p_BH) & p_BH < 0.05, "black", "grey60"),
           y.position = ymax * c(1.05, 1.19, 1.33))

  ggplot(d, aes(Group, value, fill = Group)) +
    geom_violin(trim = TRUE, alpha = 0.7, color = NA, scale = "width", width = 0.7) +
    geom_line(aes(group = ID), color = "grey70", linetype = "dotted", linewidth = 0.25) +
    geom_point(position = position_jitter(width = 0.12, seed = 1), alpha = 0.35, size = 0.5) +
    scale_fill_manual(values = pal_time) +
    stat_pvalue_manual(annot, label = "label", color = "lab_color",
                       tip.length = 0.01, bracket.size = 0.3, size = 2.2) +
    scale_color_identity() +
    labs(title = paste(strwrap(markers[[key]], 16), collapse = "\n"), x = NULL, y = NULL) +
    coord_cartesian(ylim = c(0, ymax * 1.55), clip = "off") +
    theme_minimal(base_size = 9) +
    theme(legend.position = "none",
          plot.title = element_text(size = 8),
          axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 6),
          panel.grid.minor = element_blank(), plot.margin = margin(2, 4, 2, 2))
}

# a 6-marker (2 x 3) block for one response x compartment
marker_block <- function(resp, comp)
  wrap_plots(lapply(names(markers), function(k) marker_panel(resp, comp, k)),
             nrow = 2, ncol = 3)

txt <- function(t, angle = 0, size = 4.3)
  ggplot() + annotate("text", x = 0, y = 0, label = t, angle = angle,
                      fontface = "bold", size = size) +
  theme_void() + coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1))


# Assemble: rows = response group, columns = compartment
w_lab <- 0.045
header <- (txt("") | txt("Tumor area") | txt("Stromal area")) +
          plot_layout(widths = c(w_lab, 1, 1))
row_p  <- (txt("pCR", 90) | marker_block("pCR", "Tumor") | marker_block("pCR", "Stroma")) +
          plot_layout(widths = c(w_lab, 1, 1))
row_n  <- (txt("non-pCR", 90) | marker_block("non-pCR", "Tumor") | marker_block("non-pCR", "Stroma")) +
          plot_layout(widths = c(w_lab, 1, 1))

fig <- header / row_p / row_n +
  plot_layout(heights = c(0.05, 1, 1)) +
  plot_annotation(
    caption = paste("Paired Wilcoxon signed-rank test; P values Benjamini-Hochberg",
                    "adjusted across the 72-test set. Black = adjusted P < 0.05; grey = n.s."),
    theme = theme(plot.caption = element_text(size = 8, hjust = 0)))

ggsave(file.path(out_dir, "SuppFig_S6.pdf"), fig, width = 20, height = 12, device = cairo_pdf)
cat("Wrote SuppFig_S6.pdf\n")
