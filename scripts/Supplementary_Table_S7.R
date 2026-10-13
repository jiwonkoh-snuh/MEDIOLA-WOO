# Supplementary Table S7
# Longitudinal changes in GSVA scores across treatment timepoints in the overall
# population (MEDIOLA WOO trial). For each retained Hallmark gene set, timepoint
# contrasts are tested two ways: a paired Wilcoxon signed-rank test (unadjusted
# sensitivity) and a tumor-area-adjusted linear mixed model (primary). Benjamini-
# Hochberg correction is applied across the complete set of 114 comparisons
# (38 SD-filtered gene sets x 3 timepoint contrasts).
#
# Inputs  (place under data/):
#   gsva_hallmark_scores.csv   gene set x sample GSVA scores
#   tissue_composition.csv           AI tissue composition (Tumor area = purity covariate)
# Output: gsva_overall_contrasts.csv  (Supplementary Table S7)

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(purrr)
library(lmerTest)
library(emmeans)

data_dir <- "data"
out_dir  <- "."
exclude    <- c("S25_B", "S38_B", "S40_B", "S44_B", "S60_B")
timepoints <- c("B", "OT1", "OT2")

gsva_result <- as.matrix(read.csv(file.path(data_dir, "gsva_hallmark_scores.csv"),
                                  row.names = 1, check.names = FALSE))
gsva_result <- gsva_result[, !colnames(gsva_result) %in% exclude]

ld <- read.csv(file.path(data_dir, "tissue_composition.csv"))
ld$sample <- paste(trimws(ld$Sample), trimws(ld$Class), sep = "_")
purity    <- setNames(ld$Tumor, ld$sample)

# retain gene sets with SD > 0.2 within at least one timepoint
tp_of <- sub(".*_", "", colnames(gsva_result))
keep  <- apply(gsva_result, 1, function(row)
  any(sapply(timepoints, function(tp) sd(row[tp_of == tp], na.rm = TRUE)) > 0.2))
fig3a_sets <- rownames(gsva_result)[keep]

gsva_long <- as.data.frame(gsva_result[fig3a_sets, , drop = FALSE]) %>%
  rownames_to_column("pathway") %>%
  pivot_longer(-pathway, names_to = "sample", values_to = "gsva") %>%
  mutate(timepoint  = str_extract(sample, "(B|OT1|OT2)$"),
         patient_id = str_replace(sample, "(_B|_OT1|_OT2)$", ""),
         Tumor      = purity[sample],
         timepoint  = factor(timepoint, levels = timepoints)) %>%
  filter(!is.na(Tumor))

comps <- list(c("B", "OT1"), c("B", "OT2"), c("OT1", "OT2"))

res <- map_dfr(fig3a_sets, function(ps) {
  d   <- gsva_long %>% filter(pathway == ps)
  fit <- tryCatch(suppressMessages(lmer(gsva ~ timepoint + Tumor + (1 | patient_id), data = d)),
                  error = function(e) NULL)
  pc  <- if (!is.null(fit)) as.data.frame(pairs(emmeans(fit, ~ timepoint), adjust = "none")) else NULL

  map_dfr(comps, function(cc) {
    w <- d %>% filter(timepoint %in% cc) %>%
      dplyr::select(patient_id, timepoint, gsva) %>%
      pivot_wider(names_from = timepoint, values_from = gsva) %>%
      filter(!is.na(.data[[cc[1]]]), !is.na(.data[[cc[2]]]))
    pv_w <- if (nrow(w) >= 2)
      wilcox.test(w[[cc[1]]], w[[cc[2]]], paired = TRUE, exact = FALSE)$p.value else NA_real_

    pv_adj <- NA_real_
    if (!is.null(pc)) {
      hit <- pc[pc$contrast %in% c(paste(cc[1], "-", cc[2]), paste(cc[2], "-", cc[1])), ]
      if (nrow(hit)) pv_adj <- hit$p.value[1]
    }
    tibble(pathway = ps, group1 = cc[1], group2 = cc[2], n_pairs = nrow(w),
           p_wilcox = pv_w, p_purityadj = pv_adj)
  })
}) %>%
  mutate(pathway        = sub("^HALLMARK_", "", pathway),
         p_wilcox_BH    = p.adjust(p_wilcox, "BH"),
         p_purityadj_BH = p.adjust(p_purityadj, "BH"))

write.csv(res, file.path(out_dir, "gsva_overall_contrasts.csv"), row.names = FALSE)
cat("gene sets:", dplyr::n_distinct(res$pathway),
    "| contrasts:", sum(!is.na(res$p_purityadj)), "(expect 38 / 114)\n")
