# =================== SNP QC ==========================
## ================== TWB1 ============================
setwd("/SATA2TB/Taiwan_Biobank/TWB1_genotyped_GP_v20251211_n27652/genotyped_GP_v20251211_n27652/SNP_QC")
TWB1_bim <- "/SATA2TB/Taiwan_Biobank/TWB1_genotyped_GP_v20251211_n27652/genotyped_GP_v20251211_n27652/TWB1_S03_QC.bim"

# Read bim file
bim <- read.table(TWB1_bim, header=FALSE, stringsAsFactors=FALSE,
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
  file = "TWB1_PSNPQC01_structural_keep.snplist",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


