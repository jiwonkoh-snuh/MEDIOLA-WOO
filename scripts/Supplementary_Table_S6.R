# Supplementary Table S6
# Association between molecular subtype classification and pCR status at each
# timepoint (MEDIOLA WOO trial). For PAM50 and for the Burstein TNBC scheme, a
# Fisher's exact test of the subtype x pCR table is computed at baseline, OT1, and
# OT2, with Benjamini-Hochberg correction across the three timepoints within each
# scheme. Both subtype calls, the timepoint, and pCR are read from the deposited
# Supplementary Figure S4 Source Data, which already contains the analytic cohort
# (the five RNA-seq outlier baselines are absent; B = 49, OT1 = 54, OT2 = 52).
#
# Input  (place under data/):
#   Supplementary_Figure_S4.csv   sample, patient, timepoint, pCR,
#                                 PAM50_subtype, Burstein_TNBC_subtype
# Output: subtype_pCR_assoc_by_timepoint.csv  (Supplementary Table S6)

suppressMessages(library(dplyr))

data_dir <- "data"
out_dir  <- "."
tps      <- c("B", "OT1", "OT2")

s4 <- read.csv(file.path(data_dir, "Supplementary_Figure_S4.csv"), check.names = FALSE)

# exact test of the full subtype x pCR table (simulated fallback for large tables)
safe_fisher <- function(tab) {
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
  tryCatch(fisher.test(tab)$p.value,
           error = function(e) fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value)
}
test_one <- function(subtype_col, label, tp) {
  d   <- s4[s4$timepoint == tp & !is.na(s4[[subtype_col]]) & !is.na(s4$pCR), ]
  tab <- table(d[[subtype_col]], d$pCR)
  data.frame(scheme = label, timepoint = tp, n = sum(tab),
             p = signif(safe_fisher(tab), 3))
}

res <- bind_rows(
  lapply(tps, function(t) test_one("PAM50_subtype",         "PAM50",        t)),
  lapply(tps, function(t) test_one("Burstein_TNBC_subtype", "TNBC subtype", t))
) %>%
  group_by(scheme) %>%
  mutate(p_BH = signif(p.adjust(p, "BH"), 3)) %>%
  ungroup() %>% as.data.frame()

write.csv(res, file.path(out_dir, "subtype_pCR_assoc_by_timepoint.csv"), row.names = FALSE)
print(res)
