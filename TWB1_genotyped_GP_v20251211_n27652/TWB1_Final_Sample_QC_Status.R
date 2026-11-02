
# ============ TWB 1 ============
library(dplyr)

setwd("/SATA2TB/Taiwan_Biobank/TWB1_genotyped_GP_v20251211_n27652")
# Read all sample ID
sample <- read.table(
  "TWB1.hg19.GP.20251211.fam",
  header = FALSE,
  stringsAsFactors = FALSE
)
sample_id <- sample$V1
status <- c("GCR","SEX","HET","Kinship relation (> 0.0084)","Predict not EAS population")

# Create data frame
sample_status <- as.data.frame(matrix(0, nrow = length(sample_id), ncol = length(status)))
rownames(sample_status) <- sample_id
colnames(sample_status) <- status

read_exclude <- function(file) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  
  if (file.info(file)$size == 0) {
    return(data.frame(
      FID = character(),
      IID = character()
    ))
  }
  
  read.table(
    file,
    header = FALSE,
    col.names = c("FID", "IID"),
    colClasses = "character"
  )
}

GCR <- read_exclude("Primary_Sample_QC/GCR/TWB1_S01_GCR_exclude.txt")
SEX <- read_exclude("Primary_Sample_QC/Sex/TWB1_S02_sex.problems.txt")
HET <- read_exclude("Primary_Sample_QC/HET/TWB1_batch_remove_samples_het.txt")
King <-read_exclude("King/TWB1_filter.remove")
EAS <- read_exclude("Population_structure_analysis/TWB1_NOT_EAS.id")

# Create sample QC status
sample_status <- data.frame(
  GCR = as.integer(sample_id %in% GCR$IID),
  SEX = as.integer(sample_id %in% SEX$IID),
  HET = as.integer(sample_id %in% HET$IID),
  `Kinship relation (> 0.0884)` = as.integer(sample_id %in% King$IID),
  `Predict no EAS population` = as.integer(sample_id %in% EAS$IID),
  row.names = sample_id,
  check.names = FALSE
)

# Final QC status
sample_status <- sample_status %>%
  mutate(
    Final_Status = case_when(
      rowSums(across(everything())) > 0 ~ "Exclude",
      TRUE ~ "Pass"
    )
  )

# Summary
head(sample_status)

# Export
write.csv(
  sample_status,
  "TWB1_sample_QC_status.csv",
  quote = FALSE
)

