# ==== TWB_Batch 檔案 (用於plink QC) ====

## ==== TWB2 ====
path = "/SATA2TB/Taiwan_Biobank/TWB2_genotyped_GP_v20251211_n138615"
setwd(path)

if(!file.exists("Sample_Batch")){
  dir.create("Sample_Batch")
}

TWB2 <- lab_info %>%
  filter(!is.na(TWB2_BATCH) & TWB2_BATCH != "") %>%
  select(c(TWB2_ID,TWB2_BATCH))

Batch <- unique(substr(TWB2$TWB2_BATCH,1,5))

TWB2_batch <- data.frame()
for (b in Batch){
  # Use filter() and select() or directly reference columns
  sample <- TWB2 %>% 
    filter(substr(TWB2_BATCH,1,5) == b)
  
  sample_ID <- data.frame(
    FID = sample$TWB2_ID,
    IID = sample$TWB2_ID
  )
  
  sample_batch <- data.frame(
    IID = sample$TWB2_ID,
    BATCH = substr(sample$TWB2_BATCH,1,5)
  )
  
  TWB2_batch <- rbind(TWB2_batch, sample_batch)
  
  file_name <- paste0("Sample_Batch/TWB2_", b, "_samples.keep")
  write.table(sample_ID, file = file_name, quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# SNPs Batch
path = "/SATA2TB/Taiwan_Biobank/TWB2_genotyped_GP_v20251211_n138615"
setwd(path)

TWB2_BR <- data.frame(read.csv("/SATA2TB/Taiwan_Biobank/Data_supply/TWB2_BatchRecommended_na36r2a2_20210913.csv"))
# Catch bim SNPs
bim_data <- read.table("TWB2.hg38.GP.20251211.bim", header = FALSE, stringsAsFactors = FALSE)
TWB2_snp<-bim_data$V2
rm(bim_data)

batchs <- colnames(TWB2_BR)[4:ncol(TWB2_BR)]

# 抓bim檔的SNPs
TWB2_snp_batch <- TWB2_BR %>%
  filter(SNP.ID %in% TWB2_snp) %>%
  rename(TWB2_snp = SNP.ID)

if(!file.exists("SNP_Batch")){
  dir.create("SNP_Batch")
}

for (b in batchs){
  snp <- TWB2_snp_batch %>%
    filter(.data[[b]] == 1) %>%
    pull(TWB2_snp)
  file_name <- paste0("SNP_Batch/TWB2_", b, "_recommended.snplist")
  write.table(snp, file = file_name, quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# ==== Batch GCP < 95的的 sample 彙整為 {ARRAY}_S01_GCR_exclude.txt ====
## ==== TWB2 ====
library(dplyr)
library(purrr)

path2 = "/SATA2TB/Taiwan_Biobank/TWB2_genotyped_GP_v20251211_n138615/Genotype_Call_Rate"
setwd(path2)

# List All .smiss files
smiss_files <- list.files(
  path = path,
  pattern = "\\.smiss$",
  full.names = TRUE
)
head(smiss_files)

filter_list <- list()

for (f in smiss_files) {
  TWB2_GCR <- read.table(f, header = FALSE)
  
  Filter_Sample <- TWB2_GCR %>%
    filter((1 - V6) < 0.95) %>%
    select(V1, V2)
  
  filter_list[[f]] <- Filter_Sample
}

exclude_samples <- bind_rows(filter_list) %>% distinct()

write.table(
  exclude_samples, 
  file = "TWB2_S01_GCR_exclude.txt", 
  row.names = FALSE, 
  col.names = FALSE, 
  quote = FALSE, 
  sep = "\t"
)

# ==== View Sex check ====
## ==== TWB2 ====
TWB2_sex <- read.table("/SATA2TB/Taiwan_Biobank/TWB2_genotyped_GP_v20251211_n138615/Primary_Sample_QC/Sex/TWB2_S02_sex_v2.sexcheck")
table(TWB2_sex$V4)

mismatch_or_na <- TWB2_sex[
  (TWB2_sex$V3 == 2 & TWB2_sex$V4 == 1) | 
    (TWB2_sex$V3 == 1 & TWB2_sex$V4 == 2) | 
    is.na(TWB2_sex$V3) | 
    is.na(TWB2_sex$V4), 
]

mismatch_or_na[is.na(mismatch_or_na$V4),]

# ==== Autosomal heterozygosity rate ====
run_het_qc <- function(path, het_path, batch_df, prefix) {
  #Set path
  setwd(path)
  
  # Read Data
  het_df <- read_delim(het_path, delim = '\t', show_col_types = FALSE)
  
  # Compute Heterozygosity Rate
  data <- het_df %>%
    inner_join(batch_df, by = "IID") %>%
    mutate(het_rate = (OBS_CT - `O(HOM)`) / OBS_CT) 
  
  # Plot and check batch effect
  p <- ggplot(data, aes(x = BATCH, y = het_rate, fill = BATCH)) +
    geom_boxplot(alpha = 0.7) +
    theme_minimal() +
    labs(title = paste("Heterozygosity Rate by Batch -", prefix), y = "Heterozygosity Rate")
  print(p)
  
  # Calculate the threshold based on the Mean ± 3 SD of each batch and label outliers.
  cleaned_data <- data %>%
    group_by(BATCH) %>%
    mutate(
      batch_mean = mean(het_rate, na.rm = TRUE),
      batch_sd = sd(het_rate, na.rm = TRUE),
      is_outlier = het_rate > (batch_mean + 3 * batch_sd) | het_rate < (batch_mean - 3 * batch_sd)
    ) %>%
    ungroup()
  
  # Output a list of samples to be removed
  outliers <- cleaned_data %>% 
    filter(is_outlier == TRUE) %>% 
    select('#FID', IID)
  
  dim(outliers)
  write_delim(outliers, paste0(prefix, "_batch_remove_samples_het.txt"), col_names = FALSE, delim = "\t")
  
  ## Without batch
  cleaned_data <- data %>%
    mutate(
      batch_mean = mean(het_rate, na.rm = TRUE),
      batch_sd = sd(het_rate, na.rm = TRUE),
      is_outlier = het_rate > (batch_mean + 3 * batch_sd) | het_rate < (batch_mean - 3 * batch_sd)
    ) %>%
    ungroup()
  
  # Output a list of samples to be removed
  outliers <- cleaned_data %>% 
    filter(is_outlier == TRUE) %>% 
    select('#FID', IID)
  
  dim(outliers)
  write_delim(outliers, paste0(prefix, "_remove_samples_het.txt"), col_names = FALSE, delim = "\t")
}
## ==== TWB2 ====
run_het_qc(
  path = "/SATA2TB/Taiwan_Biobank/TWB2_genotyped_GP_v20251211_n138615/Primary_Sample_QC/HET",
  het_path = "/SATA2TB/Taiwan_Biobank/TWB2_genotyped_GP_v20251211_n138615/Primary_Sample_QC/HET/TWB2_S03_het_pending.het",
  batch_df = TWB2_batch,
  prefix   = "TWB2"
)