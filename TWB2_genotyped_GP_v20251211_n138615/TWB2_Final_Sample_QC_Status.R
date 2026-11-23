library(dplyr)
# ============ TWB 2 ============
setwd("/SATA2TB/Taiwan_Biobank/TWB2_genotyped_GP_v20251211_n138615")
# Read all sample ID
sample <- read.table(
  "TWB2.hg38.GP.20251211.fam",
  header = FALSE,
  stringsAsFactors = FALSE
)
sample_id <- sample$V1
status <- c("GCR","SEX","HET","Predict not EAS population","Kinship relation (> 0.0884)")

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
  
  x <- read.table(
    file,
    header = FALSE,
    colClasses = "character"
  )
  
  if (ncol(x) == 1) {
    x <- data.frame(
      FID = x[[1]],
      IID = x[[1]]
    )
  } else if (ncol(x) == 2) {
    colnames(x) <- c("FID", "IID")
  } else {
    stop("Unexpected number of columns: ", ncol(x))
  }
  
  return(x)
}

GCR <- read_exclude("Primary_Sample_QC/GCR/TWB2_S01_GCR_exclude.txt")
SEX <- read_exclude("Primary_Sample_QC/Sex/TWB2_S02_sex.problems.txt")
HET <- read_exclude("Primary_Sample_QC/HET/TWB2_batch_remove_samples_het.txt")
King <-read_exclude("Population_structure_analysis/TWB2.relatedsamples.tsv")
EAS <- read_exclude("Population_structure_analysis/TWB2_NOT_EAS.id")

# Create sample QC status
sample_status_2 <- data.frame(
  GCR = as.integer(sample_id %in% GCR$IID),
  SEX = as.integer(sample_id %in% SEX$IID),
  HET = as.integer(sample_id %in% HET$IID),
  `Predict no EAS population` = as.integer(sample_id %in% EAS$IID),
  `Kinship relation (> 0.0884)` = as.integer(sample_id %in% King$IID),
  row.names = sample_id,
  check.names = FALSE
)

# Final QC status
sample_status_2 <- sample_status_2 %>%
  mutate(
    Final_Status = case_when(
      rowSums(across(everything())) > 0 ~ "Exclude",
      TRUE ~ "Pass"
    )
  )

# Summary
head(sample_status_2)

# Export
write.csv(
  sample_status_2,
  "TWB2_sample_QC_status.csv",
  quote = FALSE
)
