#!/bin/bash
# ===============================================================================
# Pipeline Step: Familial relatedness / highly similar genomes
# Target Dataset: TWB1 cohort
# Inputs:
#   - TWB1_S04_ldpruned.{bed,bim,fam}: LD-pruned genotype data for relatedness check
#   - TWB1_relatedness_deg2plus.kin0: KING relatedness check results for 2nd degree relatives
# Outputs:
#   - TWB1.PCAiR.tsv: PCA results from PC-AiR
#   - TWB1.relatedsamples.tsv: List of related samples identified by PC-AiR
#   - TWB1.unrelatedsamples.tsv: List of unrelated samples identified by PC-AiR
# Workflow:
#   1. Read KING output and convert to matrix format(.gds)
#   2. Check sample IDs and reorder KING matrix to match genotype data
#   3. Perform PC-AiR analysis to identify related and unrelated samples
#   4. Save PCA results and lists of related/unrelated samples
# ===============================================================================

# ==== Familial relatedness / highly similar genomes ====
library(GENESIS)
library(GWASTools)
library(SNPRelate)

## ==== TWB1 ====


kin <- read.table(
  "TWB1_relatedness_deg2plus.kin0",
  header = TRUE,
  comment.char = "",
  stringsAsFactors = FALSE
)

colnames(kin)[1:4] <- c("ID1", "ID2", "NSNP", "Kinship")

write.table(
  kin,
  "TWB1_relatedness_deg2plus.king",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


### ============ 1. Read KING output ============


KINGmat <- kingToMatrix(
  "TWB1_relatedness_deg2plus.king",
  estimator = "Kinship",
  verbose = TRUE
)

# 4042
cat("KING matrix dimensions:",
    paste(dim(KINGmat), collapse = " x "), "\n")

cat("Non-zero entries:",
    Matrix::nnzero(KINGmat), "\n")

### ============ 2. Read genotype ============

# snpgdsBED2GDS(bed.fn = "TWB1_S04_ldpruned.bed", bim.fn = "TWB1_S04_ldpruned.bim", fam.fn = "TWB1_S04_ldpruned.fam", out.gdsfn = "TWB2_S01_ldpruned.gds")

geno <- GdsGenotypeReader(
  filename = "TWB1_S04_ldpruned.gds"
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
  file = "TWB1_filter.PCAiR",
  col.names = NA,
  row.names = TRUE,
  sep = "\t",
  quote = FALSE
)

write.table(
  mypcair$rels,
  file = "TWB1_filter.relatedsamples",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

write.table(
  mypcair$unrels,
  file = "TWB1_filter.unrelatedsamples",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)


