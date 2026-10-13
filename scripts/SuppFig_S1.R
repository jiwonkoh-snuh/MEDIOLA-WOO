#!/usr/bin/env Rscript
# ============================================================
# Supplementary Figure S1 — genomic HRD & mutational signatures
# MEDIOLA WOO trial
#
# Panels:
#   A. Per-sample COSMIC v3.5 SBS contributions (named signatures + S-MSI + Others),
#      faceted by pCR, with a TMB top track and an fHRD / gBRCA / HRProfiler /
#      metabolic-response (OT1) annotation track. S58 (MSI-high) is marked.
#   C. HRProfiler HRD probability by fHRD status and by metabolic response (PET2wk)
#   D. SBS3 (HRD signature) proportion by fHRD status and by metabolic response
#   (Panel B, the MSI capillary-electrophoresis traces for S58, is an image and is
#    assembled separately.)
#
# Clinical annotation is read from the shared sample_metadata.csv (baseline rows);
# the HRProfiler probability/call comes from the HRProfiler output.
#
# Inputs  (place under data/):
#   sample_metadata.csv   ID, timepoint, RAD51, gBRCA, PET2wk, pCR (baseline rows used)
#   hrd_predictions.txt   HRProfiler output (samples, hrd.prob, prediction)
#   sbs_activities.txt    SigProfilerAssignment activities (Samples x SBS signatures)
#   tmb_by_sample.csv     per-sample TMB (ID, TMB_coding)
# Output: SuppFig_S1_signatures_HRD.pdf + genomicHRD_signature_group_pvalues.csv

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
  library(scales); library(readr); library(stringr)
})

data_dir <- "data"
out_dir  <- "."
OUT_PDF   <- file.path(out_dir, "SuppFig_S1_signatures_HRD.pdf")
OUT_STATS <- file.path(out_dir, "genomicHRD_signature_group_pvalues.csv")

## ---- palettes ---------------------------------------------
ann_pal <- c("fHRD" = "#4DAF4A", "fHRP" = "#7B4FA0",
             "Y" = "#333333", "N" = "#DDDDDD",
             "HRD" = "#3182BD", "HRP" = "#DEEBF7",
             "Responder" = "#000000", "Non-responder" = "#E0E0E0")

sig_named <- c("SBS1","SBS2","SBS3","SBS5","SBS8","SBS13",
               "SBS17a","SBS17b","SBS18","SBS30")
sig_smsi  <- c("SBS6","SBS14","SBS15","SBS20","SBS21","SBS26","SBS44")
get_cat <- function(s)
  ifelse(s %in% sig_named, s, ifelse(s %in% sig_smsi, "S-MSI", "Others"))

cat_levels <- c("SBS1","SBS2","SBS3","SBS5","SBS8","SBS13",
                "SBS17a","SBS17b","SBS18","SBS30","S-MSI","Others")
cat_pal <- c(SBS1 = "#A5C6E6", SBS2 = "#AE192C", SBS3 = "#F69D34",
             SBS5 = "#1E4C93", SBS8 = "#24B3E3", SBS13 = "#060606",
             SBS17a = "#89BE55", SBS17b = "#E87CA1", SBS18 = "#50317F",
             SBS30 = "#EA2D8D", "S-MSI" = "#188549", Others = "#CBCBCB")

## ---- canonical sample id (S01..S61) -----------------------
norm_id <- function(x) {
  d <- str_extract(x, "[0-9]+")
  ifelse(is.na(d), NA_character_, sprintf("S%02d", as.integer(d)))
}

## ---- 1. SigProfiler activities -> signature proportions ---
act <- read.delim(file.path(data_dir, "sbs_activities.txt"), check.names = FALSE)
act <- act[, !grepl("^$", colnames(act))]          # drop any trailing empty column
colnames(act)[1] <- "Sample"
act$ID   <- norm_id(act$Sample)
sig_cols <- setdiff(colnames(act), c("Sample", "ID"))
M <- as.matrix(sapply(act[sig_cols], as.numeric)); rownames(M) <- act$ID
prop <- M / rowSums(M); prop[is.na(prop)] <- 0

## ---- 2. clinical annotation (sample_metadata) + HRProfiler ----
sm <- read_csv(file.path(data_dir, "sample_metadata.csv"), show_col_types = FALSE)
names(sm) <- str_trim(names(sm))
sm <- sm %>% filter(timepoint == "B")

hp <- read.delim(file.path(data_dir, "hrd_predictions.txt"), check.names = FALSE)
hp <- hp[, colnames(hp) != "" & !is.na(colnames(hp))]   # drop empty row-index header
hp$ID <- norm_id(hp$samples)
hp <- hp %>% transmute(ID, hrd_prob = as.numeric(hrd.prob), hrd_call = prediction)

md <- sm %>%
  transmute(ID = norm_id(ID),
            fHRD  = ifelse(RAD51 == "Y", "fHRD", "fHRP"),
            gBRCA = gBRCA,
            mRes  = PET2wk,
            pCR   = pCR) %>%
  left_join(hp, by = "ID") %>%
  mutate(
    fHRD = factor(fHRD, levels = c("fHRD", "fHRP")),
    pCR  = factor(ifelse(grepl("non", pCR, ignore.case = TRUE), "non-pCR", "pCR"),
                  levels = c("pCR", "non-pCR")),
    mRes = factor(ifelse(mRes == "Y", "Responder", "Non-responder"),
                  levels = c("Responder", "Non-responder")),
    HRProfiler = factor(ifelse(hrd_call %in% c(1, "1"), "HRD", "HRP"),
                        levels = c("HRD", "HRP"))
  )
md$SBS3 <- if ("SBS3" %in% colnames(prop)) prop[md$ID, "SBS3"] else NA_real_

## ---- 3. signature-category long data ----------------------
long <- as.data.frame(prop) %>%
  mutate(ID = rownames(prop)) %>%
  pivot_longer(-ID, names_to = "Signature", values_to = "prop") %>%
  mutate(Category = get_cat(Signature)) %>%
  group_by(ID, Category) %>% summarise(prop = sum(prop), .groups = "drop") %>%
  left_join(md, by = "ID")

present <- intersect(cat_levels, unique(long$Category))
long$Category <- factor(long$Category, levels = present)

# order samples within each pCR x fHRD block by signature-composition clustering
etmat <- long %>% select(ID, Category, prop) %>%
  pivot_wider(names_from = Category, values_from = prop, values_fill = 0)
order_within <- function(ids) {
  if (length(ids) < 3) return(ids)
  m <- as.matrix(etmat[match(ids, etmat$ID), present, drop = FALSE])
  ids[hclust(dist(m), method = "ward.D2")$order]
}
ord <- md %>% group_by(pCR, fHRD) %>%
  summarise(ID = list(order_within(ID)), .groups = "drop") %>%
  arrange(pCR, fHRD) %>% tidyr::unnest(ID) %>% pull(ID)
long$ID <- factor(long$ID, levels = ord)
md$ID   <- factor(md$ID,   levels = ord)

## ---- TMB top track ----------------------------------------
tmb_tab <- read.csv(file.path(data_dir, "tmb_by_sample.csv"), stringsAsFactors = FALSE)
tmb_tab$ID <- norm_id(tmb_tab$ID)
burden <- md %>% transmute(ID = as.character(ID), pCR) %>%
  left_join(tmb_tab %>% select(ID, TMB_coding), by = "ID") %>%
  mutate(ID = factor(ID, levels = ord))
pTMB <- ggplot(burden, aes(ID, TMB_coding)) +
  geom_col(width = 1, fill = "grey25") +
  facet_grid(~pCR, scales = "free_x", space = "free_x") +
  scale_y_log10(expand = expansion(mult = c(0, 0.05)), labels = scales::label_number()) +
  labs(x = NULL, y = "TMB (mut/Mb)", tag = "A") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_blank(), panel.spacing = unit(4, "pt"))

## ---- Panel A: signature-stacked + S58 marker --------------
pA <- ggplot(long, aes(ID, prop, fill = Category)) +
  geom_col(width = 1) +
  facet_grid(~pCR, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = cat_pal[present], name = "Mutational\nsignatures", drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)), labels = percent) +
  labs(x = NULL, y = "Relative contribution") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.key.size = unit(9, "pt"), strip.background = element_blank(),
        strip.text = element_blank(), panel.spacing = unit(4, "pt"))
s58 <- md %>% filter(ID == "S58") %>% select(ID, pCR)
if (nrow(s58)) {
  pA <- pA +
    geom_text(data = s58, aes(x = ID, y = 1.05, label = "S58"),
              inherit.aes = FALSE, size = 2.6, fontface = "bold") +
    coord_cartesian(clip = "off")
}

## ---- Panel A annotation track -----------------------------
ann <- md %>%
  select(ID, pCR, fHRD, gBRCA, HRProfiler, mRes) %>%
  mutate(across(c(fHRD, gBRCA, HRProfiler, mRes), as.character)) %>%
  pivot_longer(c(fHRD, gBRCA, HRProfiler, mRes), names_to = "track", values_to = "val") %>%
  mutate(track = ifelse(track == "mRes", "mRes (OT1)", track),
         track = factor(track, levels = c("mRes (OT1)", "HRProfiler", "gBRCA", "fHRD")))
pAnn <- ggplot(ann, aes(ID, track, fill = val)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_grid(~pCR, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = ann_pal, name = NULL) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_blank(), panel.grid = element_blank(),
        strip.text = element_blank(), strip.background = element_blank(),
        panel.spacing = unit(4, "pt"), legend.key.size = unit(8, "pt"))

## ---- Panels C & D: group comparisons (Wilcoxon) -----------
make_group_panel <- function(d, yv, ylab, tag) {
  dd <- bind_rows(
    transmute(d, comp = "Functional HRD", grp = as.character(fHRD), y = .data[[yv]]),
    transmute(d, comp = "Metabolic response (OT1)", grp = as.character(mRes), y = .data[[yv]])
  ) %>% filter(!is.na(y), !is.na(grp))
  dd$comp <- factor(dd$comp, levels = c("Functional HRD", "Metabolic response (OT1)"))
  pv <- dd %>% group_by(comp) %>%
    summarise(p = wilcox.test(y ~ grp, exact = FALSE)$p.value,
              ymax = max(y, na.rm = TRUE), .groups = "drop") %>%
    mutate(lab = paste0("P = ", formatC(p, format = "f", digits = 3)))
  ggplot(dd, aes(grp, y)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, fill = "grey92") +
    geom_jitter(width = 0.12, size = 0.7, alpha = 0.6) +
    facet_wrap(~comp, scales = "free_x") +
    geom_text(data = pv, aes(x = 1.5, y = ymax * 1.06, label = lab),
              size = 2.8, inherit.aes = FALSE) +
    labs(x = NULL, y = ylab, tag = tag) +
    theme_classic(base_size = 9) +
    theme(strip.background = element_blank())
}
pC <- make_group_panel(md, "hrd_prob", "HRProfiler HRD probability", "C")
pD <- make_group_panel(md, "SBS3",     "SBS3 (HRD) proportion",      "D")

## ---- assemble & save --------------------------------------
fig <- (pTMB / pA / pAnn / pC / pD) + plot_layout(heights = c(1, 4, 1.2, 3, 3))
ggsave(OUT_PDF, fig, width = 8, height = 11)
message("Wrote ", OUT_PDF)

## ---- group-comparison stats (raw + BH) --------------------
stat_df <- data.frame(
  test = c("HRProfiler_prob ~ fHRD", "HRProfiler_prob ~ mRes",
           "SBS3 ~ fHRD", "SBS3 ~ mRes"),
  p_raw = c(wilcox.test(hrd_prob ~ fHRD, md, exact = FALSE)$p.value,
            wilcox.test(hrd_prob ~ mRes, md, exact = FALSE)$p.value,
            wilcox.test(SBS3 ~ fHRD, md, exact = FALSE)$p.value,
            wilcox.test(SBS3 ~ mRes, md, exact = FALSE)$p.value)
)
stat_df$p_BH <- p.adjust(stat_df$p_raw, "BH")
write.csv(stat_df, OUT_STATS, row.names = FALSE)

## ---- console summary --------------------------------------
n_hrd <- sum(md$HRProfiler == "HRD"); n_tot <- nrow(md)
cat("\n--- summary ---\n")
cat(sprintf("HRProfiler HRD calls: %d/%d (%.1f%%)\n", n_hrd, n_tot, 100 * n_hrd / n_tot))
cat(sprintf("Mean HRD prob  fHRD=%.3f  fHRP=%.3f\n",
            mean(md$hrd_prob[md$fHRD == "fHRD"]), mean(md$hrd_prob[md$fHRD == "fHRP"])))
cat(sprintf("Mean SBS3 prop fHRD=%.3f  fHRP=%.3f\n",
            mean(md$SBS3[md$fHRD == "fHRD"]), mean(md$SBS3[md$fHRD == "fHRP"])))
print(stat_df)
