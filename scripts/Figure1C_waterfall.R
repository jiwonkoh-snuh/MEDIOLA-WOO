# Figure 1C
# Waterfall plots of the early metabolic response (percentage change in PET SUV)
# at 2 weeks and 4 weeks, with each bar colored by the combination of the RAD51
# functional assay (fHRD/fHRP) and HRR gene mutation status, and pCR marked (*).
#
# Input  (place under data/):
#   sample_metadata.csv   (universal per-sample metadata; baseline rows used here)
# Output: Figure1C_waterfall_2wk.pdf, Figure1C_waterfall_4wk.pdf

library(ggplot2)

data_dir <- "data"
out_dir  <- "."

df <- read.csv(file.path(data_dir, "sample_metadata.csv"), check.names = FALSE)
df <- df[df$timepoint == "B", ]

# functional HRD (RAD51 assay) x HRR gene mutation status -> four groups
df$Waterfall_group <- factor(
  with(df, ifelse(RAD51 == "Y" & HRRmt == "N", "RAD51+ / HRRmt-",
           ifelse(RAD51 == "Y" & HRRmt == "Y", "RAD51+ / HRRmt+",
           ifelse(RAD51 == "N" & HRRmt == "N", "RAD51- / HRRmt-", "RAD51- / HRRmt+")))),
  levels = c("RAD51+ / HRRmt-", "RAD51+ / HRRmt+", "RAD51- / HRRmt-", "RAD51- / HRRmt+"))
df$pCR_annotation <- ifelse(df$pCR == "PCR", "*", "")

group_colors <- c("RAD51+ / HRRmt-" = "#689940", "RAD51+ / HRRmt+" = "#DDE8B2",
                  "RAD51- / HRRmt-" = "#886EB2", "RAD51- / HRRmt+" = "#C29DFF")
group_labels <- c("RAD51+ / HRRmt-" = "fHRD without HRR gene mutations",
                  "RAD51+ / HRRmt+" = "fHRD with HRR gene mutations",
                  "RAD51- / HRRmt-" = "fHRP without HRR gene mutations",
                  "RAD51- / HRRmt+" = "fHRP with HRR gene mutations")

waterfall <- function(suv_col, ylab) {
  d <- df
  d$SUV <- d[[suv_col]]
  d$ID  <- factor(d$ID, levels = d$ID[order(-d$SUV)])   # descending waterfall

  ggplot(d, aes(ID, SUV, fill = Waterfall_group, label = pCR_annotation)) +
    geom_bar(stat = "identity", color = "black") +
    geom_text(aes(y = SUV + ifelse(SUV > 0, 1, -1)),
              color = "black", size = 5, vjust = ifelse(d$SUV > 0, -0.5, 1.5)) +
    geom_hline(yintercept = -30, linetype = "dashed", color = "darkgrey", linewidth = 0.6) +
    scale_fill_manual(values = group_colors, labels = group_labels) +
    ylim(-100, 100) +
    annotate("text", x = 2, y = -90, label = "* pCR", hjust = 0, vjust = 0, size = 4) +
    labs(x = NULL, y = ylab, fill = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size = 10),
          panel.grid = element_blank(),
          legend.position = c(0.8, 0.8),
          legend.text = element_text(size = 9))
}

ggsave(file.path(out_dir, "Figure1C_waterfall_2wk.pdf"),
       waterfall("PET2wk_SUV", "SUV changes (%) at 2 weeks"), width = 9, height = 5)
ggsave(file.path(out_dir, "Figure1C_waterfall_4wk.pdf"),
       waterfall("PET4wk_SUV", "SUV changes (%) at 4 weeks"), width = 9, height = 5)
