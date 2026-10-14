# ===============================================================================
# Pipeline Step: Primary SNP QC 
# Target Dataset: TWB2 cohort
# Inputs:
#   - TWB2_S03_QC.bim: SNP information file (bim) from PLINK QC step
# Outputs:
#   - TWB2_PSNPQC01_structural_keep.snplist: List of SNPs that pass structural QC checks
# Workflow:
#   1. Extract chromosome, position, and allele information
#   2. Validate chromosome and base pair positions
#   3. Identify and exclude SNPs with duplicate IDs or genomic loci
#   4. Check for canonical biallelic SNPs and exclude non-canonical
#   5. Save the list of SNPs that pass all structural QC checks
# ===============================================================================

## ================== TWB2 ============================
setwd("/SATA2TB/Taiwan_Biobank/TWB2_genotyped_GP_v20251211_n138615/Primary_SNP_QC")
TWB2_bim <- "/SATA2TB/Taiwan_Biobank/TWB2_genotyped_GP_v20251211_n138615/Primary_SNP_QC/TWB2_S03_QC.bim"

# Read bim file
bim <- read.table(TWB2_bim, header=FALSE, stringsAsFactors=FALSE,
                  col.names=c("CHR","ID","CM","BP","A1","A2"),
                  colClasses="character")

chr <- toupper(trimws(bim$CHR))
# 判斷 chromosome 是否有效
valid_chr <- chr %in% c(as.character(1:22), "X","Y","XY","MT","M","23","24","25","26")
# BP = Base Pair position（鹼基對位置）
bp_num <- suppressWarnings(as.numeric(bim$BP))
invalid_coord <- !valid_chr | chr=="0" | is.na(bp_num) | bp_num<=0 | bp_num!=floor(bp_num)

# 檢查重複 SNP ID
dup_id <- duplicated(bim$ID) | duplicated(bim$ID, fromLast=TRUE)
# 建立 genomic locus key
loc_key <- ifelse(!invalid_coord, paste(chr, as.integer(bp_num), sep=":"), NA_character_)
# 重複 genomic locus
dup_locus <- !is.na(loc_key) & (duplicated(loc_key) | duplicated(loc_key, fromLast=TRUE))

# 檢查是否為 canonical biallelic SNP
a1 <- toupper(trimws(bim$A1)); a2 <- toupper(trimws(bim$A2))
canonical <- a1 %in% c("A","C","G","T") & a2 %in% c("A","C","G","T") & a1 != a2
structural_exclude <- invalid_coord | dup_id | dup_locus | !canonical

structural_keep <- bim$ID[!structural_exclude]

write.table(
  structural_keep,
  file = "TWB2_PSNPQC01_structural_keep.snplist",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)