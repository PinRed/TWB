library(Rsamtools)

# Read 1000G Phase 3 GRCh38 SNP reference
phase3_1kgp <- read.table(
  "/SATA2TB/Taiwan_Biobank/1kGP_SGDP_prep/1kGP/reference/1000G_Phase3_GRCh38_SNP.tsv",
  header = FALSE,
  stringsAsFactors = FALSE,
  col.names = c("CHR", "BP", "ID", "REF", "ALT")
)

phase3_1kgp$CHR <- as.character(phase3_1kgp$CHR)
phase3_1kgp$BP  <- as.integer(phase3_1kgp$BP)

# GRCh38.fa
ref_fa <- "GRCh38_genome.fa"
fa <- FaFile(ref_fa)
open(fa)

fai <- scanFaIndex(fa)
fasta_chromosomes <- as.character(seqnames(fai))

setwd("/SATA2TB/Taiwan_Biobank/TWB2_genotyped_GP_v20251211_n138615/SNP_QC")
ARRAY <- "TWB2"
input_prefix <- paste0(ARRAY, "_F_SNPQC05_clean")
dir <- "SNP_QC_FILE/"
input_file <- paste0(dir,input_prefix)
ref_fa <- "GRCh38_genome.fa"
out_prefix <- paste0(ARRAY, "_F_SNPQC06")
out_dir = "Harmonization/"
out_file <- paste0(out_dir,out_prefix)
bim_file <- paste0(input_file, ".bim")

stopifnot(
  file.exists(paste0(input_file, ".bed")),
  file.exists(bim_file),
  file.exists(paste0(input_file, ".fam"))
)

bim <- read.table(
  bim_file,
  header = FALSE,
  stringsAsFactors = FALSE,
  col.names = c("CHR", "ID", "CM", "BP", "TWB_A1", "TWB_A2")
)

bim$CHR <- as.character(bim$CHR)
bim$BP  <- as.integer(bim$BP)

# 檢查重複 ID
bim$ID_DUPLICATE <- duplicated(bim$ID) |
  duplicated(bim$ID, fromLast = TRUE)

# 檢查重複座標
pos_key <- paste(bim$CHR, bim$BP, sep = ":")

bim$POSITION_DUPLICATE <- duplicated(pos_key) |
  duplicated(pos_key, fromLast = TRUE)

table(bim$ID_DUPLICATE)
table(bim$POSITION_DUPLICATE)


# Create chromosome-position keys
bim$POS_KEY <- paste(bim$CHR, bim$BP, sep = ":")
phase3_1kgp$POS_KEY <- paste(
  phase3_1kgp$CHR, phase3_1kgp$BP, sep = ":"
)

# Check whether each TWB SNP position exists in 1KGP
bim$IN_1KGP <- bim$POS_KEY %in% phase3_1kgp$POS_KEY
# Summary
table(bim$IN_1KGP)

# Number of unique positions in each dataset
length(unique(bim$POS_KEY))
length(unique(phase3_1kgp$POS_KEY))


# 確認 BIM 染色體都存在於 FASTA
missing_chr <- setdiff(unique(as.character(bim$CHR)), fasta_chromosomes)
if (length(missing_chr) > 0) {
  stop(
    "以下染色體在 FASTA 找不到：",
    paste(missing_chr, collapse = ", ")
  )
}

# 以 1-based 座標查詢
gr <- GRanges(
  seqnames = as.character(bim$CHR),
  ranges = IRanges(start = bim$BP, width = 1)
)

fasta_base <- as.character(scanFa(fa, param = gr))
fasta_base <- toupper(fasta_base)

close(fa)

bim$FASTA_BASE <- fasta_base


# ============================================================
# Classify TWB alleles against 1000G Phase 3 GRCh38 reference
# ============================================================

# Complementary alleles
complement_allele <- function(x) {
  chartr("ACGT", "TGCA", toupper(x))
}

# Keep reference columns and identify duplicate coordinates
ref <- phase3_1kgp

# IMPORTANT:
# Only merge positions that have exactly one reference record.
# Duplicate reference positions are excluded from automatic matching.
ref_unique <- ref[!(duplicated(ref$POS_KEY) |
                      duplicated(ref$POS_KEY, fromLast = TRUE)), ]

# Match by chromosome and position
ref_unique$CHR <- as.character(ref_unique$CHR)

result <- merge(
  bim,
  ref_unique,
  by = c("CHR", "BP"),
  all.x = TRUE,
  sort = FALSE
)

# Uppercase alleles
result$TWB_A1 <- toupper(result$TWB_A1)
result$TWB_A2 <- toupper(result$TWB_A2)
result$REF <- toupper(result$REF)
result$ALT <- toupper(result$ALT)

# Reference strand complements
result$REF_COMP <- complement_allele(result$REF)
result$ALT_COMP <- complement_allele(result$ALT)

# Classify each SNP
classify_alleles <- function(a1, a2, ref, alt) {
  
  # Missing reference or invalid allele
  if (any(is.na(c(a1, a2, ref, alt)))) {
    return("NOT_IN_REFERENCE")
  }
  
  # Only biallelic SNPs with single-base alleles
  valid <- c("A", "C", "G", "T")
  
  if (!all(c(a1, a2, ref, alt) %in% valid)) {
    return("NON_ACGT_OR_NON_SNP")
  }
  
  # Direct match
  if (a2== ref && a1 == alt) {
    return("MATCH")
  }
  
  # Same alleles, reversed order
  if (a2 == alt && a1 == ref) {
    return("ALLELE_ORDER_SWAP")
  }
  
  # Strand flip, preserving allele order
  if (a2 == complement_allele(ref) &&
      a1 == complement_allele(alt)) {
    return("STRAND_FLIP")
  }
  
  # Strand flip and reversed order
  if (a2 == complement_allele(alt) &&
      a1 == complement_allele(ref)) {
    return("FLIP_AND_SWAP")
  }
  
  return("UNRESOLVED_MISMATCH")
}

result$ALLELE_CLASS <- mapply(
  classify_alleles,
  result$TWB_A1,
  result$TWB_A2,
  result$REF,
  result$ALT,
  USE.NAMES = FALSE
)

# Summary counts and percentages
print(table(result$ALLELE_CLASS, useNA = "ifany"))

result$PALINDROMIC <- with(
  result,
  (TWB_A1 == "A" & TWB_A2 == "T") |
    (TWB_A1 == "T" & TWB_A2 == "A") |
    (TWB_A1 == "C" & TWB_A2 == "G") |
    (TWB_A1 == "G" & TWB_A2 == "C")
)

# Convert NA to FALSE
result$PALINDROMIC[is.na(result$PALINDROMIC)] <- FALSE
table(result$PALINDROMIC)

# Export each classification
write.table(
  result,
  file = paste0(out_file, "_allele_classification.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

