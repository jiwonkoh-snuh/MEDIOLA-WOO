# Supplementary Figure S7
# Baseline CXCL13+ cell densities (whole tissue) compared by treatment response
# and by baseline tertiary lymphoid structure (TLS) status, for the MEDIOLA WOO
# trial.
#
# Two comparisons (by pCR, by TLS) over two markers form a single 4-test set;
# Benjamini-Hochberg is applied once across all four Wilcoxon rank-sum tests.
#
# Inputs  (place under data/):
#   mihc_density_P4.csv   CXCL13+ subsets
#   sample_metadata.csv          baseline TLS status
# Output: SuppFig_S7_pCR.pdf, SuppFig_S7_TLS.pdf

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

data_dir <- "data"
out_dir  <- "."
markers  <- c("CD4_CXCL13", "CD8_CXCL13")
marker_labels <- c(CD4_CXCL13 = "CD4+CXCL13+", CD8_CXCL13 = "CD8+CXCL13+")

# baseline whole-tissue densities + baseline TLS status
P4  <- read.csv(file.path(data_dir, "mihc_density_P4.csv"))
met <- read.csv(file.path(data_dir, "sample_metadata.csv"), check.names = FALSE)
tls_baseline <- met %>% filter(timepoint == "B") %>% dplyr::select(ID, TLS, pCR)

long <- P4 %>%
  filter(Segment == "Total", Group == "B") %>%
  dplyr::select(ID, all_of(markers)) %>%
  pivot_longer(all_of(markers), names_to = "marker", values_to = "cell_density") %>%
  left_join(tls_baseline, by = "ID") %>%
  mutate(marker = factor(marker, levels = markers),
         pCR = factor(pCR, levels = c("PCR", "nonPCR")),
         TLS = factor(TLS, levels = c("Y", "N")))


# Wilcoxon rank-sum tests (4-test set)
rank_sum <- function(d, group_var, lv) {
  x <- d$cell_density[d[[group_var]] == lv[1]]
  y <- d$cell_density[d[[group_var]] == lv[2]]
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 1 || length(y) < 1) return(NA)
  wilcox.test(x, y, exact = FALSE)$p.value
}

stats <- bind_rows(lapply(markers, function(m) {
  d <- long %>% filter(marker == m)
  data.frame(
    comparison_set = c("pCR", "TLS"),
    marker = m,
    p = c(rank_sum(d, "pCR", c("PCR", "nonPCR")),
          rank_sum(d, "TLS", c("Y", "N"))),
    stringsAsFactors = FALSE)
}))
stats$p_BH <- p.adjust(stats$p, method = "BH")
print(stats)


# Plot: violin faceted by marker, adjusted P on the bracket
pal_pcr <- c(PCR = "#8FAADC", nonPCR = "#F4707C")
pal_tls <- c(Y = "lightgreen", N = "grey70")

make_fig <- function(set_name, group_var, group_levels, palette) {
  d <- long %>%
    transmute(marker, cell_density, grp = .data[[group_var]]) %>%
    drop_na()

  stat.test <- stats %>%
    filter(comparison_set == set_name) %>%
    transmute(marker = factor(marker, levels = markers),
              group1 = group_levels[1], group2 = group_levels[2],
              label = paste0("P = ", formatC(p_BH, format = "f", digits = 3)),
              col_lab = ifelse(p_BH < 0.05, "black", "grey50")) %>%
    left_join(d %>% group_by(marker) %>%
                summarise(y.position = max(cell_density, na.rm = TRUE) * 1.12, .groups = "drop"),
              by = "marker")

  ggplot(d, aes(grp, cell_density, fill = grp)) +
    geom_violin(trim = TRUE, alpha = 0.7, color = NA, width = 0.65, scale = "width") +
    geom_point(position = position_jitter(width = 0.12), alpha = 0.25, size = 1.8) +
    scale_fill_manual(values = palette) +
    facet_wrap(~ marker, ncol = 2, scales = "free_y",
               labeller = as_labeller(marker_labels)) +
    stat_pvalue_manual(stat.test, label = "label", color = "col_lab",
                       tip.length = 0.01, bracket.size = 0.5, size = 4.3) +
    scale_color_identity() +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
    labs(x = NULL, y = expression(paste("Density (cells/mm"^2, ")"))) +
    theme_classic(base_size = 14) +
    theme(panel.grid.major.y = element_line(color = "grey90"),
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold", size = 14),
          legend.position = "none")
}

p_pcr <- make_fig("pCR", "pCR", c("PCR", "nonPCR"), pal_pcr)
p_tls <- make_fig("TLS", "TLS", c("Y", "N"), pal_tls)

ggsave(file.path(out_dir, "SuppFig_S7_pCR.pdf"), p_pcr, width = 6.5, height = 3.8, useDingbats = FALSE)
ggsave(file.path(out_dir, "SuppFig_S7_TLS.pdf"), p_tls, width = 6.5, height = 3.8, useDingbats = FALSE)
cat("Wrote SuppFig_S7_pCR.pdf and SuppFig_S7_TLS.pdf\n")
