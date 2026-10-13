# Supplementary Figure S4
# Molecular subtypes across treatment: alluvial plots of PAM50 (top) and Burstein
# TNBC (bottom) subtype classifications over serial biopsies (B / OT1 / OT2), each
# faceted by the whole cohort and by pCR status.
#
# Reads the deposited Source Data directly. The five RNA-seq outlier baselines are
# absent from the analytic cohort (B = 49, OT1 = 54, OT2 = 52); those patients render
# as an "NA" stratum at
# baseline and flow into their OT1/OT2 calls.
#
# Input  (place under data/):
#   Supplementary_Figure_S4.csv   columns: sample, patient, timepoint, pCR,
#                                 PAM50_subtype, Burstein_TNBC_subtype
# Output: SuppFig_S4_PAM50_alluvial.pdf, SuppFig_S4_Burstein_alluvial.pdf

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(ggh4x)

data_dir <- "data"
out_dir  <- "."

s4 <- read.csv(file.path(data_dir, "Supplementary_Figure_S4.csv"), check.names = FALSE)

facet_fill <- c("All population" = "#D9D9D9", "pCR" = "#8FAADC", "non-pCR" = "#F4707C")

# reshape to one row per patient x timepoint; missing timepoints become an "NA"
# stratum so the alluvium is continuous
alluvial_data <- function(subtype_col, subtype_levels) {
  wide <- s4 %>%
    dplyr::select(patient, timepoint, subtype = all_of(subtype_col), pCR) %>%
    pivot_wider(names_from = timepoint, values_from = subtype)

  long <- wide %>%
    pivot_longer(c("B", "OT1", "OT2"), names_to = "timepoint", values_to = "subtype") %>%
    mutate(timepoint = factor(timepoint, levels = c("B", "OT1", "OT2")),
           subtype   = factor(ifelse(is.na(subtype), "NA", subtype), levels = subtype_levels))

  bind_rows(mutate(long, pCR = "All population"), long) %>%
    mutate(pCR = factor(pCR, levels = c("All population", "pCR", "non-pCR")))
}

plot_alluvial <- function(data, colors, file) {
  p <- ggplot(data, aes(x = timepoint, stratum = subtype, alluvium = patient)) +
    geom_flow(stat = "alluvium", aes(fill = after_stat(stratum)),
              alpha = 0.7, lode.guidance = "frontback") +
    geom_stratum(aes(fill = after_stat(stratum))) +
    facet_wrap2(~ pCR, strip = strip_themed(
      background_x = elem_list_rect(fill = facet_fill, color = "black"),
      text_x = elem_list_text(color = "black", face = "bold"))) +
    scale_fill_manual(values = colors, drop = FALSE) +
    scale_y_continuous(limits = c(0, 54), breaks = seq(0, 54, 10)) +
    labs(x = "Timepoint", y = "N") +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
          legend.position = "right", legend.title = element_blank())
  pdf(file.path(out_dir, file), width = 12, height = 5); print(p); dev.off()
}

pam50_colors <- c("LumA" = "#6FC6FF", "LumB" = "#717FE9", "Her2" = "#A470FF",
                  "Basal" = "#46B3C7", "Normal" = "#FFD9C5", "NA" = "#999999")
tnbc_colors  <- c("BLIA" = "#E7969C", "BLIS" = "#5F9EA0", "MES" = "#B5CF6B",
                  "LAR" = "#BD9E39", "NA" = "#999999")

plot_alluvial(alluvial_data("PAM50_subtype",         names(pam50_colors)),
              pam50_colors, "SuppFig_S4_PAM50_alluvial.pdf")
plot_alluvial(alluvial_data("Burstein_TNBC_subtype", names(tnbc_colors)),
              tnbc_colors, "SuppFig_S4_Burstein_alluvial.pdf")
