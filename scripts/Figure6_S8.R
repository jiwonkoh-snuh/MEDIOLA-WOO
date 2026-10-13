# Figure 6 and Supplementary Figure S8
# Spatial proximity of immune cells to tumor cells (nearest-neighbor median
# distance) over serial biopsies for the MEDIOLA WOO trial.
#
# Figure 6 (panels A, B, D) and Supplementary Figure S8 share a single
# multiple-testing correction: Benjamini-Hochberg is applied once across the
# whole set of 62 Wilcoxon tests.
#   Panel A : P2 distances, four markers, response-stratified time course
#   Panel B : P1 CD8+CD244+ distance, stratified by response and by recurrence
#   Panel D : P2 memory T-cell distance at baseline, recurrence vs no recurrence
#   Fig S8  : P2 distances (same four markers as A), recurrence-stratified
#
# Inputs  (place under data/):
#   mihc_distance_P1.csv   CD8 co-signaling subset distances
#   mihc_distance_P2.csv   lineage / memory subset distances
# Outputs: Figure6A/B/D, FigureS8 PDFs + the adjusted p-value table.

library(data.table)
library(tidyverse)
library(ggpubr)

data_dir <- "data"
out_dir  <- "."
timepoints <- c("B", "OT1", "OT2")
pal_time  <- c(B = "#90728F", OT1 = "#9D983D", OT2 = "#E15759")
pal_recur <- c(N = "#8FAADC", Y = "#F4707C")

read_distance <- function(f) {
  d <- fread(file.path(data_dir, f))
  setnames(d, make.names(names(d), unique = TRUE))
  d %>% mutate(pCR   = factor(pCR),
               Time  = factor(Time, levels = timepoints),
               Recur = factor(Recur, levels = c(0, 1), labels = c("N", "Y")))
}

df1 <- read_distance("mihc_distance_P1.csv")
df2 <- read_distance("mihc_distance_P2.csv")


# Wilcoxon tests
# Timepoint contrasts within a stratum are paired (same patient across biopsies);
# recurrence contrasts (Panel D) are unpaired between two independent groups.

time_contrasts <- function(d) {
  pairs <- list(c("B", "OT1"), c("B", "OT2"), c("OT1", "OT2"))
  sapply(pairs, function(pair) {
    w <- d %>%
      filter(Time %in% pair) %>%
      dplyr::select(ID, Time, Median) %>%
      pivot_wider(names_from = Time, values_from = Median) %>%
      drop_na()
    if (nrow(w) < 2) return(NA)
    wilcox.test(w[[pair[1]]], w[[pair[2]]], paired = TRUE, exact = FALSE)$p.value
  })
}

markers_A <- c("CD3", "CD8_T", "CD8_T_acti", "CD8_T_acti_term")

# time-course helper for one stratifying variable (pCR or Recur)
time_course_tests <- function(df, markers, strat, panel) {
  levels_map <- if (strat == "pCR")
    list(c("1", "pCR"), c("0", "nonPCR")) else
    list(c("Y", "Recur"), c("N", "NoRecur"))
  bind_rows(lapply(levels_map, function(lv) {
    bind_rows(lapply(markers, function(m) {
      d <- df %>% filter(.data[[strat]] == lv[1], To == m, is.finite(Median))
      data.frame(panel = panel, subpanel = lv[2], marker = m,
                 g1 = c("B", "B", "OT1"), g2 = c("OT1", "OT2", "OT2"),
                 p = time_contrasts(d), stringsAsFactors = FALSE)
    }))
  }))
}

# Panel A: P2, 4 markers x 2 response groups x 3 contrasts = 24
pvals_A <- time_course_tests(df2, markers_A, "pCR", "A")

# Panel B: P1 CD8+CD244+, response- and recurrence-stratified = 12
pvals_B <- bind_rows(
  time_course_tests(df1, "CD8_CD244", "pCR",   "B"),
  time_course_tests(df1, "CD8_CD244", "Recur", "B")
)

# Panel D: P2 memory subsets at baseline, recurrence groups (unpaired) = 2
markers_D <- c("T_memory", "CD8_T_memory")
pvals_D <- bind_rows(lapply(markers_D, function(m) {
  d  <- df2 %>% filter(Time == "B", To == m, !is.na(Recur), is.finite(Median))
  g1 <- d$Median[d$Recur == "N"]; g2 <- d$Median[d$Recur == "Y"]
  p  <- if (length(g1) < 2 | length(g2) < 2) NA else
        wilcox.test(g1, g2, paired = FALSE, exact = FALSE)$p.value
  data.frame(panel = "D", subpanel = "B", marker = m,
             g1 = "N", g2 = "Y", p = p, stringsAsFactors = FALSE)
}))

# Supplementary Figure S8: P2, same 4 markers, recurrence-stratified = 24
pvals_S8 <- time_course_tests(df2, markers_A, "Recur", "S8")

all_pvals <- bind_rows(pvals_A, pvals_B, pvals_D, pvals_S8) %>%
  mutate(p_adj = p.adjust(p, method = "BH"))

cat("Total Wilcoxon tests:", nrow(all_pvals), "\n")   # 62
write.csv(all_pvals, file.path(out_dir, "Fig6_S8_all_pvalues_BHcorrected.csv"),
          row.names = FALSE)


# Plot builders
fmt_p <- function(p) ifelse(is.na(p), "NA",
                     ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                                       formatC(p, format = "f", digits = 3)))

nice_breaks <- function(ymax, k = 4) {
  if (!is.finite(ymax) || ymax <= 0) return(NULL)
  step <- ymax / k
  pow  <- floor(log10(step))
  base <- step / 10 ^ pow
  base <- if (base < 2) 2 else if (base < 5) 5 else 10
  step <- base * 10 ^ pow
  seq(0, ceiling(ymax / step) * step, by = step)
}

timecourse_panel <- function(d, pvals, y_max, title) {
  annot <- pvals %>%
    dplyr::rename(group1 = g1, group2 = g2) %>%
    arrange(match(paste0(group1, group2), c("BOT1", "BOT2", "OT1OT2"))) %>%
    mutate(label      = fmt_p(p_adj),
           lab_color  = ifelse(!is.na(p_adj) & p_adj < 0.05, "black", "grey60"),
           y.position = c(y_max * 0.80, y_max * 0.87, y_max * 0.94))

  ggplot(d, aes(Time, Median, fill = Time)) +
    geom_violin(trim = TRUE, alpha = 0.7, color = NA, width = 0.65, scale = "width") +
    geom_line(aes(group = ID), color = "grey75", linetype = "dotted") +
    geom_point(position = position_jitter(width = 0.1), alpha = 0.2, size = 1.8) +
    scale_fill_manual(values = pal_time) +
    scale_y_continuous(limits = c(0, y_max), breaks = nice_breaks(y_max),
                       expand = expansion(mult = c(0, 0.02))) +
    stat_pvalue_manual(annot, label = "label", color = "lab_color",
                       tip.length = 0.01, size = 2.8) +
    scale_color_identity() +
    labs(title = title, x = "Timepoint", y = "Median distance (\u00b5m)") +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          panel.grid.major.y = element_line(color = "grey90"),
          panel.grid.major.x = element_blank(),
          plot.margin = margin(6, 10, 6, 10)) +
    coord_cartesian(ylim = c(0, y_max), clip = "off")
}

recur_panel <- function(d, pvals, y_max, title) {
  annot <- pvals %>%
    dplyr::rename(group1 = g1, group2 = g2) %>%
    mutate(label = fmt_p(p_adj),
           lab_color = ifelse(!is.na(p_adj) & p_adj < 0.05, "black", "grey60"),
           y.position = y_max * 0.92)

  ggplot(d, aes(Recur, Median, fill = Recur)) +
    geom_violin(trim = TRUE, alpha = 0.7, color = NA, width = 0.65) +
    geom_jitter(width = 0.12, alpha = 0.6, size = 1.8, shape = 21, color = "black") +
    scale_fill_manual(values = pal_recur) +
    scale_y_continuous(limits = c(0, y_max), breaks = nice_breaks(y_max),
                       expand = expansion(mult = c(0, 0.02))) +
    stat_pvalue_manual(annot, label = "label", color = "lab_color",
                       tip.length = 0.01, size = 2.8) +
    scale_color_identity() +
    labs(title = title, x = "Recurrence", y = "Median distance (\u00b5m)") +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          panel.grid.major.y = element_line(color = "grey90"),
          panel.grid.major.x = element_blank(),
          plot.margin = margin(6, 10, 6, 10)) +
    coord_cartesian(ylim = c(0, y_max), clip = "off")
}

save_grid <- function(plots, file, ncol, nrow, width, height) {
  pdf(file.path(out_dir, file), width = width, height = height)
  print(ggarrange(plotlist = plots, ncol = ncol, nrow = nrow, align = "hv"))
  dev.off()
}

# shared per-marker y-limits so response, recurrence and S8 panels are comparable
finite2 <- df2 %>% filter(To %in% markers_A, Time %in% timepoints, is.finite(Median))
ymax_A  <- tapply(finite2$Median, finite2$To, max, na.rm = TRUE)

titles_A <- c(CD3 = "CD3+ T cell", CD8_T = "CD8+ T cell",
              CD8_T_acti = "CD8+ activated T cell",
              CD8_T_acti_term = "CD8+ activated memory T cell")


# Panel A  (response-stratified time course)
for (sp in c("pCR", "nonPCR")) {
  grp <- ifelse(sp == "pCR", "1", "0")
  plots <- lapply(markers_A, function(m) {
    d  <- df2 %>% filter(pCR == grp, To == m, is.finite(Median))
    pv <- filter(all_pvals, panel == "A", subpanel == sp, marker == m)
    timecourse_panel(d, pv, ymax_A[[m]], titles_A[[m]])
  })
  save_grid(plots, paste0("Figure6A_", sp, "_adjP.pdf"),
            ncol = 2, nrow = 2, width = 10, height = 8)
}


# Panel B  (CD8+CD244+, response and recurrence)
finite1  <- df1 %>% filter(To == "CD8_CD244", is.finite(Median))
ymax_B   <- max(finite1$Median, na.rm = TRUE)

plots_B_pcr <- lapply(c("pCR", "nonPCR"), function(sp) {
  grp <- ifelse(sp == "pCR", "1", "0")
  d   <- df1 %>% filter(pCR == grp, To == "CD8_CD244", is.finite(Median))
  pv  <- filter(all_pvals, panel == "B", subpanel == sp, marker == "CD8_CD244")
  timecourse_panel(d, pv, ymax_B, paste0("CD8+CD244+ cell (", sp, ")"))
})
save_grid(plots_B_pcr, "Figure6B_pCR_adjP.pdf",
          ncol = 2, nrow = 1, width = 10, height = 4)

plots_B_recur <- lapply(c("NoRecur", "Recur"), function(sp) {
  grp <- ifelse(sp == "Recur", "Y", "N")
  d   <- df1 %>% filter(Recur == grp, To == "CD8_CD244", is.finite(Median))
  pv  <- filter(all_pvals, panel == "B", subpanel == sp, marker == "CD8_CD244")
  timecourse_panel(d, pv, ymax_B, paste0("CD8+CD244+ cell (", sp, ")"))
})
save_grid(plots_B_recur, "Figure6B_recur_adjP.pdf",
          ncol = 2, nrow = 1, width = 10, height = 4)


# Panel D  (memory T cells at baseline, recurrence vs no recurrence)
finiteD <- df2 %>% filter(To %in% markers_D, Time == "B", is.finite(Median))
ymax_D  <- tapply(finiteD$Median, finiteD$To, max, na.rm = TRUE)
titles_D <- c(T_memory = "Memory T cell", CD8_T_memory = "CD8+ memory T cell")

plots_D <- list()
for (m in markers_D) {
  d  <- df2 %>% filter(Time == "B", To == m, !is.na(Recur), is.finite(Median))
  pv <- filter(all_pvals, panel == "D", subpanel == "B", marker == m)
  plots_D[[paste("B", m)]] <-
    recur_panel(d, pv, ymax_D[[m]], paste0(titles_D[[m]], " (B)"))
}
save_grid(plots_D, "Figure6D_memory_adjP.pdf",
          ncol = 2, nrow = 1, width = 10, height = 4)


# Supplementary Figure S8  (recurrence-stratified time course)
# The recurrence group is small (n = 5); its contrasts are shown for
# completeness, were uniformly non-significant, and are described as
# exploratory in the legend. Same per-marker scale as Panel A.
for (sp in c("NoRecur", "Recur")) {
  grp <- ifelse(sp == "Recur", "Y", "N")
  plots <- lapply(markers_A, function(m) {
    d  <- df2 %>% filter(Recur == grp, To == m, is.finite(Median))
    pv <- filter(all_pvals, panel == "S8", subpanel == sp, marker == m)
    timecourse_panel(d, pv, ymax_A[[m]], titles_A[[m]])
  })
  save_grid(plots, paste0("FigureS8_", sp, "_adjP.pdf"),
            ncol = 2, nrow = 2, width = 10, height = 8)
}
