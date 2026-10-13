# Figure 1B
# Baseline genomic and clinicopathologic landscape of the MEDIOLA WOO cohort:
# an oncoprint of somatic/germline driver mutations (14 genes x 54 patients) with
# clinical annotations on top and a per-gene mutation-frequency bar on the right.
#
# Input  (place under data/):
#   sample_metadata.csv   (universal per-sample metadata; baseline rows used here)
# Output: Figure1B_heatmap.pdf

library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(readr)

data_dir <- "data"
out_dir  <- "."

df <- read.csv(file.path(data_dir, "sample_metadata.csv"), check.names = FALSE,
               na.strings = c("", "NA"))   # empty gene cells -> NA (= no mutation)
df <- df[df$timepoint == "B", ]
rownames(df) <- df$ID

# --- mutation matrix ---
mutation_cols <- c("TP53", "gBRCA1", "sBRCA1", "gBRCA2", "gPALB2", "PIK3CA", "PTEN",
                   "PIK3R1", "NF1", "RB1", "ARID1A", "KMT2C", "SPEN", "CTCF")

mutation_matrix <- df[, mutation_cols] %>%
  mutate_all(~ recode(.,
                      "missense"   = "Missense",
                      "stopgain"   = "Nonsense",
                      "frameshift" = "Frameshift",
                      "splice"     = "Splicing",
                      "inframe"    = "Inframe",
                      "exondel"    = "Exonic deletion")) %>%
  as.matrix()
mutation_matrix[is.na(mutation_matrix)] <- "None"
rownames(mutation_matrix) <- rownames(df)

set.seed(1)   
age_levels <- c("20-29", "30-39", "40-49", "50-59", "60-69")
stopifnot(all(df$Age %in% age_levels))         
df$Age <- factor(df$Age, levels = age_levels)

# --- factor levels and display recoding for the annotations ---
df$pCR    <- factor(df$pCR,    levels = c("PCR", "nonPCR"))
df$PET2wk <- factor(df$PET2wk, levels = c("Y", "N"))
df$PET4wk <- factor(df$PET4wk, levels = c("Y", "N"))
df$CT4wk  <- factor(df$CT4wk,  levels = c("Y", "N"))
df$Recur  <- factor(df$Recur,  levels = c("Y", "N"))
df$RAD51  <- factor(df$RAD51,  levels = c("Y", "N"))
df$PDL1   <- factor(df$PDL1,   levels = c("P", "N"))

df$RAD51  <- recode(df$RAD51, "Y" = "fHRD", "N" = "fHRP")
df$ER     <- recode(df$ER,   "Low" = "Low", "N" = "Negative")
df$HER2   <- recode(df$HER2, "Low" = "Low", "N" = "Negative")
df$PDL1   <- recode(df$PDL1, "P" = "Positive", "N" = "Negative")
df$PET2wk <- recode(df$PET2wk, "Y" = "Yes", "N" = "No")
df$PET4wk <- recode(df$PET4wk, "Y" = "Yes", "N" = "No")
df$CT4wk  <- recode(df$CT4wk,  "Y" = "Yes", "N" = "No")
df$Recur  <- recode(df$Recur,  "Y" = "Yes", "N" = "No")
df$`Ki-67` <- df$`Ki-67` * 100

# order columns by response, then metabolic response, then RAD51 assay
sorted_samples  <- rownames(df)[order(df$pCR, df$PET2wk, df$RAD51, df$ID)]
mutation_matrix <- mutation_matrix[sorted_samples, ]

# --- top annotation ---
annotation_order <- c("pCR", "PET2wk", "PET4wk", "CT4wk", "Age", "Stage", "Recur",
                      "RAD51", "ER", "HER2", "PDL1", "Ki-67", "PAM50", "TNBC")
anno_combined <- df[sorted_samples, annotation_order]

anno_col_all <- list(
  pCR    = c("PCR" = "#8B8BFF", "nonPCR" = "#FF7373"),
  RAD51  = c("fHRD" = "#4DAF4A", "fHRP" = "#984EA3"),
  PET2wk = c("Yes" = "#9E9EFF", "No" = "#FF9090"),
  PET4wk = c("Yes" = "#BABAFF", "No" = "#FFACAC"),
  CT4wk  = c("Yes" = "#BABAFF", "No" = "#FFACAC"),
  Stage  = c("II" = "#6B94C1", "III" = "#2C3D50"),
  Recur  = c("Yes" = "#DD95E8", "No" = "#A9CD56"),
  PAM50  = c("LumA" = "#6FC6FF", "LumB" = "#717FE9", "Her2" = "#A470FF",
             "Basal" = "#46B3C7", "Normal" = "#FFD9C5"),
  TNBC   = c("BLIA" = "#E7969C", "BLIS" = "#5F9EA0", "MES" = "#B5CF6B", "LAR" = "#BD9E39"),
  ER     = c("Low" = "#FFB5B5", "Negative" = "#EFEFEF"),
  HER2   = c("Low" = "#E1C1FF", "Negative" = "#EFEFEF"),
  PDL1   = c("Positive" = "#B2DF8A", "Negative" = "#EFEFEF"),
  Age    = c("20-29" = "#435569", "30-39" = "#637488", "40-49" = "#8494A6",
             "50-59" = "#A4B4C5", "60-69" = "#C4D3E4"),
  `Ki-67` = colorRamp2(range(df$`Ki-67`, na.rm = TRUE), c("#E2C9C9", "darkred"))
)

top_anno <- HeatmapAnnotation(
  df = anno_combined, col = anno_col_all, na_col = "white",
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  annotation_legend_param = list(nrow = 2),
  annotation_label = c(pCR = "pCR", PET2wk = "Metabolic response OT1",
                       PET4wk = "Metabolic response OT2", CT4wk = "RECIST response OT2",
                       Age = "Age", Stage = "Stage", Recur = "Recurrence",
                       RAD51 = "RAD51 assay", ER = "ER", HER2 = "HER2",
                       PDL1 = "PD-L1 (22C3)", `Ki-67` = "Ki-67 (%)",
                       PAM50 = "PAM50", TNBC = "TNBC subtype"))

# --- per-gene mutation frequency (right bar) ---
mutation_pct <- apply(mutation_matrix, 2, function(col) 100 * sum(col != "None") / length(col))
right_bar <- rowAnnotation(
  "Mutation %" = anno_barplot(mutation_pct, gp = gpar(fill = "#9999FF", col = "white"),
                              border = FALSE, bar_width = 0.8, width = unit(2.5, "cm"),
                              ylim = c(0, 100),
                              axis_param = list(at = seq(0, 100, 25),
                                                labels = paste0(seq(0, 100, 25), "%"))),
  annotation_name_gp = gpar(fontsize = 9), annotation_label = NULL)

mut_colors <- c("Missense" = "#275D36", "Nonsense" = "#010001", "Frameshift" = "#3C519A",
                "Splicing" = "#E55398", "Inframe" = "#BF3438",
                "Exonic deletion" = "#E59947", "None" = "#EFEFEF")

ht <- Heatmap(t(mutation_matrix), name = "Mutation", col = mut_colors,
              top_annotation = top_anno, right_annotation = right_bar,
              cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
              row_names_gp = gpar(fontsize = 9, fontface = "bold.italic"),
              heatmap_legend_param = list(nrow = 2))

pdf(file.path(out_dir, "Figure1B_heatmap.pdf"), width = 14, height = 8)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "right")
decorate_heatmap_body("Mutation", {
  nr <- ncol(mutation_matrix)   # transposed: columns become rows
  for (i in 1:(nr - 1))
    grid.lines(y = unit(i / nr, "npc"), x = unit(c(0, 1), "npc"),
               gp = gpar(col = "white", lwd = 2))
})
dev.off()

writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "Figure1B_sessionInfo.txt"))
