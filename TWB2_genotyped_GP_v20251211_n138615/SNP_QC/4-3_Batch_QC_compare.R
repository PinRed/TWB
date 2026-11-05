library(data.table)
##########################################
# Goal :
#   Compare Each two batch 
#   1. SNP missingness difference > 0.02
#   2. Allele-frequency difference > 0.10
##########################################

setwd("/SATA2TB/Taiwan_Biobank/TWB2_genotyped_GP_v20251211_n138615/SNP_QC")
ARRAY = "TWB2"
outdir = paste0(ARRAY, "_F_SNPQC05_batch_effect")
BATCH_LIST_FILE <- file.path(outdir, "processed_batches.txt")

MISSINGNESS_THRESHOLD <- 0.02
AF_THRESHOLD <- 0.10

# 2. Check and load processed batch list
if (!file.exists(BATCH_LIST_FILE)) {
  stop(paste("ERROR: Batch list file not found:", BATCH_LIST_FILE))
}

batches <- readLines(BATCH_LIST_FILE)
batches <- batches[batches != ""]
n_batches <- length(batches)

cat(sprintf("--> Found %d processed batches.\n", n_batches))

if (n_batches < 2) {
  stop("ERROR: At least two batches are required for batch-effect screening.")
}

# 3. Read and aggregate statistics across all batches
vmiss_list <- vector("list", n_batches)
afreq_list <- vector("list", n_batches)

for (i in seq_along(batches)) {
  b <- batches[i]
  prefix <- file.path(outdir, paste0("B", b), paste0(ARRAY, "_B", b, "_stats"))
  
  vmiss_file <- paste0(prefix, ".vmiss")
  afreq_file <- paste0(prefix, ".afreq")
  
  if (!file.exists(vmiss_file) || !file.exists(afreq_file)) {
    warning(sprintf("Warning: Missing stats files for batch B%s. Skipping.", b))
    next
  }
  
  # Read PLINK2 .vmiss (Columns: #CHROM, ID, MISS_CT, OBS_CT, F_MISS)
  dt_v <- fread(vmiss_file, select = c("#CHROM", "ID", "F_MISS"))
  setnames(dt_v, "#CHROM", "CHROM")
  dt_v[, BATCH := b]
  vmiss_list[[i]] <- dt_v
  
  # Read PLINK2 .afreq (Columns: #CHROM, ID, REF, ALT, ALT_FREQS, OBS_CT)
  dt_a <- fread(afreq_file, select = c("#CHROM", "ID", "ALT_FREQS"))
  setnames(dt_a, "#CHROM", "CHROM")
  dt_a[, BATCH := b]
  afreq_list[[i]] <- dt_a
}

# Clean null elements if any skipped
vmiss_list <- Filter(Negate(is.null), vmiss_list)
afreq_list <- Filter(Negate(is.null), afreq_list)

cat("--> Combining statistics across batches...\n")
all_vmiss <- rbindlist(vmiss_list)
all_afreq <- rbindlist(afreq_list)

# 4. Calculate max pairwise differences (Max - Min across batches)
cat("--> Screening for batch effects (Max difference calculation)...\n")

# Missingness summary
vmiss_summary <- all_vmiss[, .(
  MIN_MISS = min(F_MISS, na.rm = TRUE),
  MAX_MISS = max(F_MISS, na.rm = TRUE),
  DIFF_MISS = max(F_MISS, na.rm = TRUE) - min(F_MISS, na.rm = TRUE)
), by = .(CHROM, ID)]

vmiss_summary[, PROBLEMATIC_MISS := DIFF_MISS > MISSINGNESS_THRESHOLD]

# Allele Frequency summary
afreq_summary <- all_afreq[, .(
  MIN_AF = min(ALT_FREQS, na.rm = TRUE),
  MAX_AF = max(ALT_FREQS, na.rm = TRUE),
  DIFF_AF = max(ALT_FREQS, na.rm = TRUE) - min(ALT_FREQS, na.rm = TRUE)
), by = .(CHROM, ID)]

afreq_summary[, PROBLEMATIC_AF := DIFF_AF > AF_THRESHOLD]

# 5. Merge results and generate final flags
merged_summary <- merge(vmiss_summary, afreq_summary, by = c("CHROM", "ID"))
merged_summary[, PROBLEMATIC_ANY := PROBLEMATIC_MISS | PROBLEMATIC_AF]

# 6. Save outputs
out_summary_file <- file.path(outdir, "batch_effect_snp_qc_detailed_summary.tsv")
out_problem_file <- file.path(outdir, "problematic_variants_to_exclude.txt")

fwrite(merged_summary, out_summary_file, sep = "\t")

# Export plain list of SNP IDs flagged as problematic for subsequent --exclude in PLINK
problematic_ids <- merged_summary[PROBLEMATIC_ANY == TRUE, ID]
fwrite(data.table(ID = problematic_ids), out_problem_file, col.names = FALSE, sep = "\t")

# 7. Print summary statistics
cat("\n========================================\n")
cat(" F-SNPQC-05 Batch-Effect Screen Summary\n")
cat("========================================\n")
cat(sprintf("Total evaluated variants : %d\n", nrow(merged_summary)))
cat(sprintf("Flagged by Missingness diff (> %.2f): %d\n", MISSINGNESS_THRESHOLD, sum(merged_summary$PROBLEMATIC_MISS)))
cat(sprintf("Flagged by AF diff (> %.2f)         : %d\n", AF_THRESHOLD, sum(merged_summary$PROBLEMATIC_AF)))
cat(sprintf("Total unique problematic variants: %d\n", length(problematic_ids)))
cat(sprintf("Detailed report saved to : %s\n", out_summary_file))
cat(sprintf("Problematic SNP list saved: %s\n", out_problem_file))
cat("========================================\n")