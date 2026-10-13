# Figure 4, Figure 5A-B and Supplementary Figure S5
# Multiplex IHC T-cell densities over serial biopsies (B / OT1 / OT2) and the
# TLS-vs-imaging-response analysis for the MEDIOLA WOO trial.
#
# Figure 4A, Figure 4B, Figure 5A (CXCL13+ subsets) and Supplementary Figure S5
# (CD8+TIM3+ / CD8+LAG3+) share a single multiple-testing correction:
# Benjamini-Hochberg is applied once across the whole set of 81 paired Wilcoxon
# tests. Figure 5B is independent and uses Fisher's exact test.
#
# Inputs  (place under data/):
#   mihc_density_P1.csv   CD8 exhaustion/co-signaling subsets
#   mihc_density_P2.csv   CD3/CD4/CD8 lineage and memory subsets
#   mihc_density_P4.csv   CXCL13+ subsets
#   sample_metadata.csv        universal per-sample metadata incl. TLS status (Figure 5B)
# Outputs: Figure4A, Figure4B, Figure5A, Figure5B, FigureS5 PDFs
#          + the adjusted p-value table and the Fisher results table.
#
# Output sizes: Figure 4A and 5B are saved at final print size (place at 100%).
# Figure 4B and 5A keep a large canvas with enlarged text and are scaled down to fit.

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

data_dir <- "data"
out_dir  <- "."
timepoints <- c("B", "OT1", "OT2")

# Print-size settings for panels saved at final size (4A, 5B)
mm  <- 1 / 25.4      # mm -> inch
fs4 <- 6             # base font size (pt) at print size
pv4 <- 5 / .pt       # p-value label size = 5 pt

read_density <- function(f) {
  d <- read.csv(file.path(data_dir, f))
  d <- d[d$Group %in% timepoints, ]
  d$Group   <- factor(d$Group,   levels = timepoints)
  d$Segment <- factor(d$Segment, levels = c("Total", "Tumor", "Stroma"))
  d
}

P1 <- read_density("mihc_density_P1.csv")
P2 <- read_density("mihc_density_P2.csv")
P4 <- read_density("mihc_density_P4.csv")
metadata <- read.csv(file.path(data_dir, "sample_metadata.csv"), check.names = FALSE)


# Paired Wilcoxon tests (Figure 4A, 4B, Figure 5A, Supplementary Figure S5)
# Each marker is compared across the three timepoint pairs on complete pairs.

paired_wilcox <- function(d, marker, pair) {
  w <- d %>%
    filter(Group %in% pair) %>%
    dplyr::select(ID, Group, all_of(marker)) %>%
    pivot_wider(names_from = Group, values_from = all_of(marker)) %>%
    drop_na()
  if (nrow(w) < 2) return(NA)
  wilcox.test(w[[pair[1]]], w[[pair[2]]], paired = TRUE, exact = FALSE)$p.value
}

three_contrasts <- function(d, marker) {
  pairs <- list(c("B", "OT1"), c("B", "OT2"), c("OT1", "OT2"))
  sapply(pairs, function(p) paired_wilcox(d, marker, p))
}

# helper that loops a marker set over one or both response groups
collect <- function(df, markers, segment, panel, by_pcr = TRUE) {
  d0 <- df[df$Segment == segment, ]
  groups <- if (by_pcr) c("pCR", "nonPCR") else "All"
  bind_rows(lapply(groups, function(sp) {
    d <- if (sp == "All") d0 else d0[d0$pCR == ifelse(sp == "pCR", 1, 0), ]
    bind_rows(lapply(markers, function(m) {
      data.frame(panel = panel, subpanel = sp, marker = m,
                 g1 = c("B", "B", "OT1"), g2 = c("OT1", "OT2", "OT2"),
                 p = three_contrasts(d, m), stringsAsFactors = FALSE)
    }))
  }))
}

# Figure 4A: whole-cohort lineage densities (Total segment, not split by response)
pvals_4A <- collect(P2, c("CD3", "CD4_T", "CD8_T"), "Total", "4A", by_pcr = FALSE)

# Figure 4B: response-stratified CD8/CD4 subsets. Most come from P2; the CD8+PD1+
# and CD8+CD244+ subsets come from P1.
pvals_4B <- bind_rows(
  collect(P2, c("CD8_T", "CD8_T_acti", "CD8_T_acti_term", "CD8_T_memory",
                "CD4_T", "CD4_T_memory"), "Total", "4B"),
  collect(P1, c("CD8_PD1", "CD8_CD244"), "Total", "4B")
)

# Figure 5A: CXCL13+ subsets, measured in the stromal segment
pvals_5A <- collect(P4, c("CD4_CXCL13", "CD8_CXCL13"), "Stroma", "5A")

# Supplementary Figure S5: CD8+TIM3+ / CD8+LAG3+ (Total segment), same design
# as Figure 4B and therefore part of the same correction set.
pvals_S5 <- collect(P1, c("CD8_TIM3", "CD8_LAG3"), "Total", "S5")

all_pvals <- bind_rows(pvals_4A, pvals_4B, pvals_5A, pvals_S5) %>%
  mutate(p_adj = p.adjust(p, method = "BH"))

cat("Total paired Wilcoxon tests:", nrow(all_pvals), "\n")   # 81
write.csv(all_pvals, file.path(out_dir, "Fig4_5A_S5_all_pvalues_BHcorrected.csv"),
          row.names = FALSE)


# Violin panel builder
# Adjusted p-values are printed on the brackets: black when < 0.05, grey when not.
# base / psize / ptsize / lw / ypos default to the original values (used for Figure S5).

fmt_p <- function(p) ifelse(is.na(p), "NA",
                            ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                                   formatC(p, format = "f", digits = 3)))

pal_time <- c(B = "#90728F", OT1 = "#9D983D", OT2 = "#E15759")

violin_panel <- function(d, marker, pvals, y_max, title = NULL,
                         base = 14, psize = 3.5, ptsize = 1.8, lw = 0.5,
                         ypos = c(0.80, 0.87, 0.94)) {
  annot <- pvals %>%
    rename(group1 = g1, group2 = g2) %>%
    arrange(match(paste0(group1, group2), c("BOT1", "BOT2", "OT1OT2"))) %>%
    mutate(label     = fmt_p(p_adj),
           lab_color = ifelse(!is.na(p_adj) & p_adj < 0.05, "black", "grey60"),
           y.position = y_max * ypos)
  
  ggplot(d, aes(Group, .data[[marker]], fill = Group)) +
    geom_violin(trim = TRUE, alpha = 0.7, color = NA, scale = "width", width = 0.65) +
    geom_line(aes(group = ID), color = "grey75", linetype = "dotted", linewidth = lw) +
    geom_point(position = position_jitter(width = 0.1), alpha = 0.2, size = ptsize) +
    labs(x = "Timepoint", y = expression(paste("Density (cells/mm"^2, ")")),
         title = title) +
    scale_fill_manual(values = pal_time) +
    scale_y_continuous(limits = c(0, y_max), breaks = seq(0, y_max, length.out = 5),
                       expand = expansion(mult = c(0, 0.02))) +
    stat_pvalue_manual(annot, label = "label", color = "lab_color",
                       tip.length = 0.01, size = psize, bracket.size = lw) +
    scale_color_identity() +
    theme_classic(base_size = base) +
    theme(legend.position = "none",
          panel.grid.major.y = element_line(color = "grey90"),
          panel.grid.major.x = element_blank(),
          plot.margin = margin(base * 0.43, base * 0.71, base * 0.43, base * 0.71)) +
    coord_cartesian(ylim = c(0, y_max), clip = "off")
}

# round a marker's maximum up to a clean axis limit
axis_max <- function(v) {
  v <- max(v, na.rm = TRUE)
  step <- 10 ^ floor(log10(v))
  ceiling(v / step) * step
}

save_grid <- function(plots, file, ncol, nrow, width, height) {
  pdf(file.path(out_dir, file), width = width, height = height)
  print(ggarrange(plotlist = plots, ncol = ncol, nrow = nrow, align = "hv"))
  dev.off()
}


# Figure 4A  (whole cohort) - saved at final print size
titles_A <- c(CD3 = "CD3+ T-cell", CD4_T = "CD4+ T-cell", CD8_T = "CD8+ T-cell")
ymax_A   <- c(CD3 = 8000, CD4_T = 4000, CD8_T = 8000)   # fixed y-axis limits

plots_A <- lapply(names(titles_A), function(m) {
  d  <- P2[P2$Segment == "Total", ]
  pv <- filter(all_pvals, panel == "4A", marker == m)
  violin_panel(d, m, pv, ymax_A[[m]], titles_A[[m]],
               base = fs4, psize = pv4, ptsize = 0.5, lw = 0.25)
})
save_grid(plots_A, "Figure4A_adjP.pdf", ncol = 3, nrow = 1,
          width = 180 * mm, height = 55 * mm)


# Figure 4B  (response-stratified, 8 markers in a 3 x 3 grid, last cell blank)
# Row 1: CD8+, CD8+PD1+, CD8+CD244+
# Row 2: CD8+ activated, CD8+ activated memory, CD8+ memory
# Row 3: CD4+, CD4+ memory, (blank)
markers_B <- c("CD8_T", "CD8_PD1", "CD8_CD244",
               "CD8_T_acti", "CD8_T_acti_term", "CD8_T_memory",
               "CD4_T", "CD4_T_memory")
titles_B <- c(CD8_T = "CD8+ T-cell", CD8_PD1 = "CD8+PD1+ T-cell",
              CD8_CD244 = "CD8+CD244+ T-cell",
              CD8_T_acti = "CD8+ activated\nT-cell",
              CD8_T_acti_term = "CD8+ activated\nmemory T-cell",
              CD8_T_memory = "CD8+ memory\nT-cell",
              CD4_T = "CD4+ T-cell", CD4_T_memory = "CD4+ memory\nT-cell")

density_source <- function(m) if (m %in% c("CD8_PD1", "CD8_CD244")) P1 else P2
ymax_B <- sapply(markers_B, function(m) {
  d <- density_source(m); axis_max(d[d$Segment == "Total", m])
})

for (sp in c("pCR", "nonPCR")) {
  plots <- lapply(markers_B, function(m) {
    d  <- density_source(m)
    d  <- d[d$Segment == "Total" & d$pCR == ifelse(sp == "pCR", 1, 0), ]
    pv <- filter(all_pvals, panel == "4B", subpanel == sp, marker == m)
    violin_panel(d, m, pv, ymax_B[[m]], titles_B[[m]],
                 base = 18, psize = 6, ptsize = 1.8, lw = 0.5,
                 ypos = c(0.72, 0.84, 0.96))
  })
  save_grid(c(plots, list(NULL)), paste0("Figure4B_", sp, "_adjP.pdf"),
            ncol = 3, nrow = 3, width = 12, height = 12)
}


# Figure 5A  (CXCL13+, stromal segment) - large canvas, enlarged text
markers_5A <- c("CD4_CXCL13", "CD8_CXCL13")
titles_5A  <- c(CD4_CXCL13 = "CD4+CXCL13+ cells", CD8_CXCL13 = "CD8+CXCL13+ cells")
ymax_5A    <- sapply(markers_5A, function(m) axis_max(P4[P4$Segment == "Stroma", m]))

for (sp in c("pCR", "nonPCR")) {
  plots <- lapply(markers_5A, function(m) {
    d  <- P4[P4$Segment == "Stroma" & P4$pCR == ifelse(sp == "pCR", 1, 0), ]
    pv <- filter(all_pvals, panel == "5A", subpanel == sp, marker == m)
    violin_panel(d, m, pv, ymax_5A[[m]], titles_5A[[m]],
                 base = 18, psize = 6, ptsize = 1.8, lw = 0.5,
                 ypos = c(0.72, 0.84, 0.96))
  })
  save_grid(plots, paste0("Figure5A_", sp, "_adjP.pdf"),
            ncol = 2, nrow = 1, width = 8, height = 4)
}


# Supplementary Figure S5  (CD8+TIM3+ / CD8+LAG3+, Total segment) - unchanged
markers_S5 <- c("CD8_TIM3", "CD8_LAG3")
titles_S5  <- c(CD8_TIM3 = "CD8+TIM3+ cells", CD8_LAG3 = "CD8+LAG3+ cells")
ymax_S5    <- sapply(markers_S5, function(m) axis_max(P1[P1$Segment == "Total", m]))

for (sp in c("pCR", "nonPCR")) {
  plots <- lapply(markers_S5, function(m) {
    d  <- P1[P1$Segment == "Total" & P1$pCR == ifelse(sp == "pCR", 1, 0), ]
    pv <- filter(all_pvals, panel == "S5", subpanel == sp, marker == m)
    violin_panel(d, m, pv, ymax_S5[[m]], titles_S5[[m]])
  })
  save_grid(plots, paste0("FigureS5_", sp, "_adjP.pdf"),
            ncol = 2, nrow = 1, width = 8, height = 4)
}


# Figure 5B  (TLS vs imaging response at OT1, Fisher's exact test) - final print size
df_ot1 <- metadata %>% filter(timepoint == "OT1")

fisher_pet <- fisher.test(table(df_ot1$TLS, df_ot1$PET4wk))
fisher_ct  <- fisher.test(table(df_ot1$TLS, df_ot1$CT4wk))

# BH correction across the two response readouts (mRes at OT2, RECIST at OT2)
p_adj_5B <- p.adjust(c(fisher_pet$p.value, fisher_ct$p.value), method = "BH")

fisher_results <- data.frame(
  modality   = c("PET4wk", "CT4wk"),
  p_value    = c(fisher_pet$p.value,  fisher_ct$p.value),
  p_adj      = p_adj_5B,
  odds_ratio = c(fisher_pet$estimate, fisher_ct$estimate)
)
write.csv(fisher_results, file.path(out_dir, "Fig5B_Fisher_results.csv"),
          row.names = FALSE)
print(fisher_results)

df_long <- df_ot1 %>%
  pivot_longer(c(PET4wk, CT4wk), names_to = "Modality", values_to = "Response") %>%
  mutate(TLS      = factor(TLS, levels = c("Y", "N")),
         Response = factor(ifelse(Response == "Y", "Response", "No Response"),
                           levels = c("Response", "No Response")),
         Modality = factor(Modality, levels = c("PET4wk", "CT4wk")))

pval_annot <- data.frame(
  Modality = factor(c("PET4wk", "CT4wk"), levels = c("PET4wk", "CT4wk")),
  label = paste0("P = ", signif(p_adj_5B, 3)),
  y = 1.05)

resp_colors <- c("Response" = "#8FAADC", "No Response" = "#F4707C")

p_5B <- ggplot(df_long, aes(TLS, fill = Response)) +
  geom_bar(position = "fill", width = 0.7) +
  facet_wrap(~ Modality, nrow = 1) +
  scale_fill_manual(values = resp_colors) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25),
                     labels = c("0", "25", "50", "75", "100"), expand = c(0, 0)) +
  geom_text(data = pval_annot, aes(x = 1.5, y = y, label = label),
            inherit.aes = FALSE, size = pv4) +
  labs(x = "TLS at OT1", y = "Proportion (%)", fill = "Response") +
  theme_classic(base_size = fs4) +
  theme(strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom",
        panel.grid.major.y = element_line(color = "grey90"),
        panel.grid.minor = element_blank())

pdf(file.path(out_dir, "Figure5B_TLS_response.pdf"), width = 55 * mm, height = 60 * mm)
print(p_5B)
dev.off()