# TNBC molecular subtyping (Burstein 80-gene signature)
# Samples are assigned to the four TNBC subtypes (BLIA, BLIS, MES, LAR) by
# nearest-centroid (Pearson) classification using the published 80-gene Burstein
# signature; subtype centroids are derived from the Burstein reference cohort
# (GEO accession GSE76124).
#
# Inputs  (place under data/):
#   rnaseq_tpm.txt         gene x sample TPM (RNA-seq; see GEO)
#   burstein_centroids.txt   centroid table (rows = genes, cols = BLIA/BLIS/MES/LAR)
# Output: tnbc_result_burstein80.csv
#   (all 160 libraries are classified; the 5 RNA-seq outliers are dropped later,
#    at the Supplementary Figure S4 / Table S6 stage.)

tpm_log2 <- log2(as.matrix(
  read.table(file.path("data", "rnaseq_tpm.txt"), header = TRUE,
             row.names = 1, check.names = FALSE)) + 1)

# 80-gene Burstein signature (Burstein et al., Clin Cancer Res 2015)
burstein80 <- c(
  "DHRS2","GABRP","AGR2","PIP","FOXA1","PROM1","TFF1","NAT1","BCL11A","ESR1",
  "FOXC1","CA12","TFF3","SCUBE2","SFRP1","ERBB4","SIDT1","PSAT1","CHI3L1","AR",
  "CD36","OGN","ABCA8","CFD","IGF1","HBB","CDH1","MEOX2","GPX3","SCARA5",
  "PDK4","ENPP2","AGTR1","LEP","LPL","DPT","TIMP4","FHL1","SRPX","EDNRB",
  "SERPINB5","SOX10","IRX1","MIA","DSC2","TTYH1","COL9A3","FGL2","RARRES3","PDE9A",
  "BST2","PTGER4","KCNK5","PSMB9","HLA-DMA","EPHB3","IGSF6","ST3GAL6","RHOH","SGPP1",
  "CXCL9","CXCL11","GBP5","GZMB","LAMP3","GBP1","ADAMDEC1","CCL5","SPON1","PBK",
  "STAT1","EZH2","PLAT","TAP2","SLAMF7","HERC5","SPOCK1","TAP1","CD2","AIM2")

# read the centroid table
lines  <- readLines(file.path("data", "burstein_centroids.txt"))
header <- strsplit(lines[1], "\t")[[1]]
rows   <- strsplit(lines[-1], "\t")
centroids <- as.data.frame(t(sapply(rows, function(x) as.numeric(x[2:5]))))
colnames(centroids) <- header
rownames(centroids) <- sapply(rows, `[`, 1)
centroids <- centroids[!grepl("///", rownames(centroids)), ]

# subset to the signature genes present in both the centroids and the expression
sig <- intersect(burstein80, rownames(centroids))
cat("signature genes in centroid file:", length(sig), "of 80\n")
common <- intersect(rownames(tpm_log2), rownames(centroids[sig, ]))
cat("matched to expression:", length(common), "genes\n")

cor_tnbc <- cor(tpm_log2[common, ], as.matrix(centroids[common, ]), method = "pearson")
tnbc_result <- data.frame(
  sample       = rownames(cor_tnbc),
  TNBC_subtype = apply(cor_tnbc, 1, function(x) colnames(cor_tnbc)[which.max(x)]),
  TNBC_cor     = round(apply(cor_tnbc, 1, max), 3),
  stringsAsFactors = FALSE)

write.csv(tnbc_result, file.path(".", "tnbc_result_burstein80.csv"), row.names = FALSE)
print(table(tnbc_result$TNBC_subtype))
