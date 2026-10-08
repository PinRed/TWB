# ==== Familial relatedness / highly similar genomes ====
library(GENESIS)
library(GWASTools)
library(SNPRelate)
library(Matrix)

## ======== TWB2 ========
path2 = "/SATA2TB/Taiwan_Biobank/TWB2_genotyped_GP_v20251211_n138615/King"
setwd(path2)



### ============ 轉換KingToMatrix格式 ============

kin <- read.table(
  "TWB2_relatedness_deg2plus.kin0",
  header = TRUE,
  comment.char = "",
  stringsAsFactors = FALSE
)
colnames(kin)[1:4] <- c("ID1", "ID2", "NSNP", "Kinship")

write.table(
  kin,
  "TWB2_relatedness_deg2plus.king",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


### ============ 1. Read KING output ============


# 38,920 samples
KINGmat <- kingToMatrix(
  "TWB2_relatedness_deg2plus.king",
  estimator = "Kinship",
  verbose = TRUE
)

cat("KING matrix dimensions:",
    paste(dim(KINGmat), collapse = " x "), "\n")

cat("Non-zero entries:",
    Matrix::nnzero(KINGmat), "\n")

### ============ 2. Read genotype ============

snpgdsBED2GDS(bed.fn = "TWB2_S04_ldpruned.bed", bim.fn = "TWB2_S04_ldpruned.bim", fam.fn = "TWB2_S04_ldpruned.fam", out.gdsfn = "TWB2_S04_ldpruned.gds")

geno <- GdsGenotypeReader(
  filename = "TWB2_S04_ldpruned.gds"
)

genoData <- GenotypeData(geno)
geno.ids <- getScanID(genoData)

rm(geno)


### ============ 3. Check sample IDs ============

cat("Genotype samples:", length(geno.ids), "\n")
cat("KING samples:", nrow(KINGmat), "\n")

missing.king <- setdiff(geno.ids, rownames(KINGmat))

cat("Genotype samples missing from KING:",
    length(missing.king), "\n")


### ============ 4. Reorder KING matrix ============

common.ids <- intersect(
  geno.ids,
  rownames(KINGmat)
)

KINGmat <- KINGmat[
  common.ids,
  common.ids,
  drop = FALSE
]


### ============ 5. PC-AiR ============

mypcair <- pcair(
  genoData,
  kinobj = KINGmat,
  divobj = NULL,
  kin.thresh = 0.0884,
  num.cores = 80,
  verbose = TRUE,
  maf = 0.05,
  missing.rate = 0.01,
  algorithm = "randomized",
  eigen.cnt = 32
)


### ============ 6. Output ============

colnames(mypcair$vectors) <-
  paste0("PCAiR", seq_len(ncol(mypcair$vectors)))

write.table(
  mypcair$vectors,
  file = "TWB2.PCAiR.tsv",
  col.names = NA,
  row.names = TRUE,
  sep = "\t",
  quote = FALSE
)

write.table(
  mypcair$rels,
  file = "TWB2.relatedsamples.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

write.table(
  mypcair$unrels,
  file = "TWB2.unrelatedsamples.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

write.table(
  mypcair$values,
  file = "TWB2.eigenvalue.txt",
  col.names = TRUE,
  row.names = TRUE,
  sep = "\t",
  quote = FALSE
)

write.table(
  mypcair$varprop,
  file = "TWB2.varprop.txt",
  col.names = TRUE,
  row.names = TRUE,
  sep = "\t",
  quote = FALSE
)