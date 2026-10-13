# Supplementary Table S4
# Linear mixed-effects model tests of response-dependent tissue-composition dynamics
# for the MEDIOLA WOO trial. For each digital-pathology compartment, a mixed model
# with timepoint, response status, and their interaction as fixed effects and a
# patient-level random intercept is fitted; the timepoint-by-response interaction
# tests whether the two response groups change differently over time.
#
# Run for three response endpoints (pCR, metabolic response at OT1, metabolic
# response at OT2). Benjamini-Hochberg correction is applied within each endpoint
# (across its four compartments). This is a supporting analysis for Figure 2 / S3,
# reported as Supplementary Table S4 (not a plotted panel).
#
# Input : tissue_composition.csv   (place under data/)
# Output: Supplementary_Table_S4.csv

library(dplyr)
library(lme4)
library(lmerTest)

data_dir <- "data"
out_dir  <- "."

df <- read.csv(file.path(data_dir, "tissue_composition.csv"))
df$Class    <- factor(df$Class, levels = c("B", "OT1", "OT2"))
df$pCR_f    <- factor(df$pCR,    levels = c("PCR", "nonPCR"), labels = c("pCR", "non-pCR"))
df$mRes_OT1 <- factor(df$PET2wk, levels = c("Y", "N"), labels = c("mRes", "non-mRes"))
df$mRes_OT2 <- factor(df$PET4wk, levels = c("Y", "N"), labels = c("mRes", "non-mRes"))

comp <- c("Tumor", "Stroma", "Necrosis", "Inflammatory")

interaction_test <- function(group_col, label) {
  res <- do.call(rbind, lapply(comp, function(v) {
    a    <- anova(lmer(as.formula(paste0(v, " ~ Class * ", group_col, " + (1 | Sample)")), data = df))
    term <- paste0("Class:", group_col)
    data.frame(Endpoint = label, Compartment = v,
               NumDF = a[term, "NumDF"], DenDF = a[term, "DenDF"],
               F = a[term, "F value"], p = a[term, "Pr(>F)"])
  }))
  res$p.adj <- p.adjust(res$p, method = "BH")   # within this endpoint (4 compartments)
  res
}

S4 <- rbind(
  interaction_test("pCR_f",    "Pathologic complete response"),
  interaction_test("mRes_OT1", "Metabolic response at OT1"),
  interaction_test("mRes_OT2", "Metabolic response at OT2")
)

write.csv(S4, file.path(out_dir, "Supplementary_Table_S4.csv"), row.names = FALSE)
print(S4, row.names = FALSE)
