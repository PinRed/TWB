
library(data.table)
library(Matrix)
library(GENESIS)
library(GWASTools)
library(gdsfmt)

setwd("/SATA2TB/Taiwan_Biobank/TWB2_genotyped_GP_v20251211_n138615/Population_structure_analysis")

# ============================================================
# 1. Read sample IDs
# ============================================================

ids <- fread(
  "TWB2_EAS.id",
  header = FALSE
)$V1

stopifnot(!anyDuplicated(ids))

# ============================================================
# 2. Read kinship file and construct K
# ============================================================

k <- fread("TWB2_EAS_pcair_kin.kin0")

setnames(k, sub("^#", "", names(k)))

id1_col <- intersect(c("IID1", "ID1"), names(k))[1]
id2_col <- intersect(c("IID2", "ID2"), names(k))[1]

stopifnot(!is.na(id1_col), !is.na(id2_col))
stopifnot("KINSHIP" %in% names(k))

i <- match(k[[id1_col]], ids)
j <- match(k[[id2_col]], ids)

stopifnot(!anyNA(i), !anyNA(j))

# 對稱 kinship matrix，對角線設為 0.5
K <- sparseMatrix(
  i = c(i, j, seq_along(ids)),
  j = c(j, i, seq_along(ids)),
  x = c(
    k$KINSHIP,
    k$KINSHIP,
    rep(0.5, length(ids))
  ),
  dims = c(length(ids), length(ids)),
  dimnames = list(ids, ids)
)

# ============================================================
# 3. Read genotype GDS
# ============================================================

gds_file <- "TWB2_EAS_pcair.gds"

gds_conn <- openfn.gds(gds_file)

gds_snp_ids <- read.gdsn(
  index.gdsn(gds_conn, "snp.id")
)

closefn.gds(gds_conn)


geno <- GdsGenotypeReader(filename = gds_file)
gdsobj <- GenotypeData(geno)


geno.ids <- getScanID(gdsobj)

stopifnot(!anyDuplicated(geno.ids))

cat("Number of genotype samples:", length(geno.ids), "\n")
cat("Number of GDS SNPs:", length(gds_snp_ids), "\n")


# ============================================================
# 4. Read LD-pruned SNP list
# ============================================================

pruned_snps <- fread(
  "TWB2_EAS_pcair_prune.prune.in",
  header = FALSE
)$V1

stopifnot(!anyDuplicated(pruned_snps))

cat("Number of SNPs in prune.in:", length(pruned_snps), "\n")

pcair_snps <- gds_snp_ids %in% pruned_snps

cat("Number of matched SNPs:", sum(pcair_snps), "\n")

missing_snps <- setdiff(pruned_snps, gds_snp_ids)

cat("Pruned SNPs missing from GDS:", length(missing_snps), "\n")

stopifnot(any(pcair_snps))
stopifnot(length(missing_snps) == 0)

# ============================================================
# 5. Check and align sample IDs
# ============================================================

stopifnot(setequal(ids, geno.ids))

# 將 K 重新排列成 genotype sample 順序
K <- K[geno.ids, geno.ids]

stopifnot(identical(rownames(K), geno.ids))
stopifnot(identical(colnames(K), geno.ids))

# ============================================================
# 6. Run PC-AiR
# ============================================================

pcair_fit <- pcair(
  gdsobj,
  kinobj = K,
  divobj = NULL,
  kin.thresh = 2^(-11/2),
  snp.include = gds_snp_ids
)

# ============================================================
# 7. Close GDS reader after analysis
# ============================================================

close(geno)

# 檢查結果
print(pcair_fit)

# 儲存結果
saveRDS(
  pcair_fit,
  file = "TWB2_EAS_pcair_fit.rds"
)
