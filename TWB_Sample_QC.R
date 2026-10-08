# ==== TWB_Batch 檔案 (用於plink QC) ====
## ==== TWB1 Sample ====
library(dplyr)

path = "/SATA2TB/Taiwan_Biobank/TWB1_genotyped_GP_v20251211_n27652"
setwd(path)

# 建立存放Sample Batch 的資料夾
if(!file.exists("Sample_Batch")){
  dir.create("Sample_Batch")
}

# 讀取sample batch 資訊
lab_info <- read.csv("../../lab_info/lab_info.csv")
TWB1 <- lab_info %>%
  filter(!is.na(TWB1_BATCH) & TWB1_BATCH != "") %>%
  select(c(TWB1_ID,TWB1_BATCH))

Batch <- unique(TWB1$TWB1_BATCH)

TWB1_batch <- data.frame()
for (b in Batch){
  # Use filter() and select() or directly reference columns
  sample <- TWB1 %>% 
    filter(TWB1_BATCH == b)
  
  sample_ID <- data.frame(
    FID = sample$TWB1_ID,
    IID = sample$TWB1_ID
  )
  
  sample_batch <- data.frame(
    IID = sample$TWB1_ID,
    BATCH = sample$TWB1_BATCH
  )
  
  TWB1_batch <- rbind(TWB1_batch, sample_batch)
  
  file_name <- paste0("Sample_Batch/TWB1_", b, "_samples.keep")
  write.table(sample_ID, file = file_name, quote = FALSE, row.names = FALSE, col.names = FALSE)
}


# SNPs Batch
library(readxl)
path = "/SATA2TB/Taiwan_Biobank/TWB1_genotyped_GP_v20251211_n27652"
setwd(path)

TWB1_BR <- data.frame(read_excel("/SATA2TB/Taiwan_Biobank/Data_supply/TWB1_BatchRecommended_na35v2_20180917.xlsx", sheet = 1))
# Catch bim SNPs
bim_data <- read.table("TWB1.hg19.GP.20251211.bim", header = FALSE, stringsAsFactors = FALSE)
TWB1_snp<-bim_data$V2
rm(bim_data)

batchs <- c("batch1", "batch2", "batch3", "batch4", "batch5", 
            "batch6", "batch7", "batch8", "batch9", "batch10", 
            "batch11", "batch12", "batch13", "batch14")

# 抓bim檔的SNPs
TWB1_snp_batch <- TWB1_BR %>%
  filter(Affy.SNP.ID %in% TWB1_snp | dbSNP.RS.ID %in% TWB1_snp) %>%
  mutate(TWB1_snp = if_else(
    Affy.SNP.ID %in% TWB1_snp,
    Affy.SNP.ID,
    dbSNP.RS.ID
  )) %>%
  select(TWB1_snp, all_of(batchs))

if(!file.exists("SNP_Batch")){
  dir.create("SNP_Batch")
}

for (b in batchs){
  snp <- TWB1_snp_batch %>%
    filter(.data[[b]] == 1) %>%
    pull(TWB1_snp)
  file_name <- paste0("SNP_Batch/TWB1_", b, "_recommended.snplist")
  write.table(snp, file = file_name, quote = FALSE, row.names = FALSE, col.names = FALSE)
}


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
## ==== TWB1 ====
path1 = "/SATA2TB/Taiwan_Biobank/TWB1_genotyped_GP_v20251211_n27652/Genotype_Call_Rate"
setwd(path1)

# List All .smiss files
library(dplyr)
library(purrr)

smiss_files <- list.files(
  path = path,
  pattern = "\\.smiss$",
  full.names = TRUE
)

filter_list <- list()

for (f in smiss_files) {
  TWB1_GCR <- read.table(f, header = FALSE)
  
  Filter_Sample <- TWB1_GCR %>%
    filter((1 - V6) < 0.95) %>%
    select(V1, V2)
  
  filter_list[[f]] <- Filter_Sample
}

exclude_samples <- bind_rows(filter_list) %>% distinct()

write.table(
  exclude_samples, 
  file = "TWB1_S01_GCR_exclude.txt", 
  row.names = FALSE, 
  col.names = FALSE, 
  quote = FALSE, 
  sep = "\t"
)

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
## ==== TWB1 ====
TWB1_sex <- read.table("/SATA2TB/Taiwan_Biobank/TWB1_genotyped_GP_v20251211_n27652/Primary_Sample_QC/Sex/TWB1_S02_sex.sexcheck")
table(TWB1_sex$V5,TWB1_sex$V4)

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
library(tidyverse)

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

## ===== TWB1 ===== 
run_het_qc(
  path = "/SATA2TB/Taiwan_Biobank/TWB1_genotyped_GP_v20251211_n27652/Primary_Sample_QC/HET",
  het_path = "/SATA2TB/Taiwan_Biobank/TWB1_genotyped_GP_v20251211_n27652/Primary_Sample_QC/HET/TWB1_S03_het_pending.het",
  batch_df = TWB1_batch,
  prefix   = "TWB1"
)

## ==== TWB2 ====
run_het_qc(
  path = "/SATA2TB/Taiwan_Biobank/TWB2_genotyped_GP_v20251211_n138615/Primary_Sample_QC/HET",
  het_path = "/SATA2TB/Taiwan_Biobank/TWB2_genotyped_GP_v20251211_n138615/Primary_Sample_QC/HET/TWB2_S03_het_pending.het",
  batch_df = TWB2_batch,
  prefix   = "TWB2"
)



