# Supplementary Table S8
# GSVA score changes over treatment, stratified by pCR status (MEDIOLA WOO trial).
# For every retained Hallmark gene set, within each response group (pCR, non-pCR),
# the three timepoint contrasts are tested two ways:
#   - primary   : tumor-area-adjusted linear mixed model (timepoint + Tumor area
#                 fixed effects, patient random intercept), timepoint contrasts by
#                 estimated marginal means
#   - sensitivity: paired Wilcoxon signed-rank test (unadjusted)
# Benjamini-Hochberg correction is applied across the complete set of 228
# comparisons (38 SD-filtered gene sets x 2 response groups x 3 contrasts) for
# each model separately.
#
# Inputs  (place under data/):
#   gsva_hallmark_scores.csv   gene set x sample GSVA scores
#   tissue_composition.csv           AI tissue composition (Tumor area, pCR)
# Output: GSVA_pCR_stratified_S8.csv  (Supplementary Table S8)

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
comps      <- list(c("B", "OT1"), c("B", "OT2"), c("OT1", "OT2"))

gsva <- as.matrix(read.csv(file.path(data_dir, "gsva_hallmark_scores.csv"),
                           row.names = 1, check.names = FALSE))
gsva <- gsva[, !colnames(gsva) %in% exclude]

ld <- read.csv(file.path(data_dir, "tissue_composition.csv"))
ld$sample <- paste(trimws(ld$Sample), trimws(ld$Class), sep = "_")

# retained gene sets: SD > 0.2 within at least one timepoint
tp_of <- sub(".*_", "", colnames(gsva))
keep  <- apply(gsva, 1, function(row)
  any(sapply(timepoints, function(tp) sd(row[tp_of == tp], na.rm = TRUE)) > 0.2))
sets  <- rownames(gsva)[keep]

long <- as.data.frame(gsva[sets, , drop = FALSE]) %>%
  rownames_to_column("pathway") %>%
  pivot_longer(-pathway, names_to = "sample", values_to = "score") %>%
  mutate(timepoint = factor(str_extract(sample, "(B|OT1|OT2)$"), levels = timepoints),
         patient   = str_replace(sample, "(_B|_OT1|_OT2)$", "")) %>%
  left_join(ld %>% dplyr::select(sample, Tumor, pCR), by = "sample") %>%
  mutate(pCR = trimws(pCR)) %>%
  filter(!is.na(Tumor), pCR %in% c("PCR", "nonPCR"))

# one gene set x response group: adjusted mixed model + unadjusted paired Wilcoxon
per_group <- function(d) {
  fit <- tryCatch(suppressMessages(lmer(score ~ timepoint + Tumor + (1 | patient), data = d)),
                  error = function(e) NULL)
  pc  <- if (!is.null(fit)) as.data.frame(pairs(emmeans(fit, ~ timepoint), adjust = "none")) else NULL
  map_dfr(comps, function(cc) {
    w <- d %>% filter(timepoint %in% cc) %>%
      dplyr::select(patient, timepoint, score) %>%
      pivot_wider(names_from = timepoint, values_from = score) %>%
      filter(!is.na(.data[[cc[1]]]), !is.na(.data[[cc[2]]]))
    p_unadj <- if (nrow(w) >= 2)
      wilcox.test(w[[cc[1]]], w[[cc[2]]], paired = TRUE, exact = FALSE)$p.value else NA_real_
    p_adj <- NA_real_
    if (!is.null(pc)) {
      hit <- pc[pc$contrast %in% c(paste(cc[1], "-", cc[2]), paste(cc[2], "-", cc[1])), ]
      if (nrow(hit)) p_adj <- hit$p.value[1]
    }
    tibble(timepoint1 = cc[1], timepoint2 = cc[2], p_purityadj = p_adj, p_unadj = p_unadj)
  })
}

out <- long %>%
  group_by(pathway, pCR) %>%
  group_modify(~ per_group(.x)) %>%
  ungroup() %>%
  mutate(pathway = sub("^HALLMARK_", "", pathway),
         response_group = recode(pCR, "PCR" = "pCR", "nonPCR" = "non-pCR"),
         comparison     = paste(timepoint1, "\u2192", timepoint2),
         p_purityadj_BH = p.adjust(p_purityadj, "BH"),
         p_unadj_BH     = p.adjust(p_unadj, "BH"),
         sig_purityadj  = !is.na(p_purityadj_BH) & p_purityadj_BH < 0.05,
         sig_unadj      = !is.na(p_unadj_BH)     & p_unadj_BH     < 0.05) %>%
  dplyr::select(pathway, response_group, comparison, timepoint1, timepoint2,
                p_purityadj, p_purityadj_BH, sig_purityadj,
                p_unadj, p_unadj_BH, sig_unadj)

write.csv(out, file.path(out_dir, "GSVA_pCR_stratified_S8.csv"), row.names = FALSE)
cat("rows:", nrow(out), "(expect 228 = 38 gene sets x 2 groups x 3 contrasts)\n")
