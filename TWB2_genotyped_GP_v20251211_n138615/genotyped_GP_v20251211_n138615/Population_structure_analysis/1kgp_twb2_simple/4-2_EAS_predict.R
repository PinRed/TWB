setwd("/SATA2TB/Taiwan_Biobank/TWB2_genotyped_GP_v20251211_n138615/genotyped_GP_v20251211_n138615/Population_structure_analysis/1kgp_twb2_simple")

# ============================================================
# Random Forest classification
# Training : 1KGP
# Prediction : TWB2
#
# Classes:
#   AFR / AMR / EAS / EUR / SAS
#
# Rule:
#   P(EAS) > 0.9  -> EAS
#   otherwise     -> UNCLASSIFIED
# ============================================================


library(randomForest)


# ============================================================
# 1. 1KGP metadata
# ============================================================

pop_1kgp <- read.table(
  "/SATA2TB/Taiwan_Biobank/1kGP_SGDP_prep/metadata/igsr_samples.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)


# ============================================================
# 2. TWB2 FAM
# ============================================================

TWB2 <- read.delim(
  "/SATA2TB/Taiwan_Biobank/TWB2_genotyped_GP_v20251211_n138615/genotyped_GP_v20251211_n138615/King/TWB2_S05.fam",
  header = FALSE,
  sep = "",
  quote = ""
)

twb2fam <- data.frame(
  IID = TWB2$V2,
  stringsAsFactors = FALSE
)

rm(TWB2)


# ============================================================
# 3. PCA
# ============================================================

twb2pca <- read.table(
  "/SATA2TB/Taiwan_Biobank/TWB2_genotyped_GP_v20251211_n138615/genotyped_GP_v20251211_n138615/Population_structure_analysis/1kgp_twb2_simple/TWB2_pca.eigenvec",
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(twb2pca)[1:8] <- c(
  "FID",
  "IID",
  "PC1",
  "PC2",
  "PC3",
  "PC4",
  "PC5",
  "PC6"
)


# ============================================================
# 4. 建立 1KGP training data
# ============================================================

train <- merge(
  twb2pca[, c("IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6")],
  pop_1kgp[, c("Sample.name", "Superpopulation.code")],
  by.x = "IID",
  by.y = "Sample.name"
)


# 只保留五個 super-populations
train <- train[
  train$Superpopulation.code %in%
    c("AFR", "AMR", "EAS", "EUR", "SAS"),
]


# ============================================================
# 5. Training data 整理
# ============================================================

train$Population <- factor(
  train$Superpopulation.code,
  levels = c(
    "AFR",
    "AMR",
    "EAS",
    "EUR",
    "SAS"
  )
)

train_rf <- train[
  complete.cases(train[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "Population")]),
  c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "Population")
]


cat("========================================\n")
cat("Random Forest training data\n")
cat("========================================\n")
print(table(train_rf$Population))
cat("\n")
cat("Number of samples:", nrow(train_rf), "\n\n")


# ============================================================
# 6. Random Forest
# ============================================================

set.seed(26)

rf_model <- randomForest(
  Population ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6,
  data = train_rf,
  ntree = 500,
  importance = TRUE
)


# ============================================================
# 7. RF model summary
# ============================================================

cat("========================================\n")
cat("Random Forest model\n")
cat("========================================\n")

print(rf_model)

cat("\n========================================\n")
cat("Variable importance\n")
cat("========================================\n")

print(importance(rf_model))


# ============================================================
# 8. TWB2 prediction data
# ============================================================

twb2 <- twb2pca[
  twb2pca$IID %in% twb2fam$IID,
  c("IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6")
]

twb2 <- twb2[
  complete.cases(twb2[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")]),
]


# ============================================================
# 9. RF prediction probability
# ============================================================

rf_prob <- predict(
  rf_model,
  newdata = twb2[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")],
  type = "prob"
)


# ============================================================
# 10. RF predicted class
# ============================================================

rf_pred <- predict(
  rf_model,
  newdata = twb2[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")],
  type = "response"
)


classification <- ifelse(
  rf_prob[, "EAS"] > 0.9,
  "EAS",
  "UNCLASSIFIED"
)


# ============================================================
# 12. 建立 final result
# ============================================================

rf_result <- data.frame(
  IID = twb2$IID,
  
  P_AFR = rf_prob[, "AFR"],
  P_AMR = rf_prob[, "AMR"],
  P_EAS = rf_prob[, "EAS"],
  P_EUR = rf_prob[, "EUR"],
  P_SAS = rf_prob[, "SAS"],
  
  RF_prediction = as.character(rf_pred),
  
  Final_classification = classification,
  
  stringsAsFactors = FALSE
)

eas <- rf_result[rf_result$Final_classification == "EAS",]$IID

# ============================================================
# 13. Result summary
# ============================================================

cat("\n========================================\n")
cat("Final classification\n")
cat("========================================\n")

print(
  table(rf_result$Final_classification)
)


# ============================================================
# 14. Save result
# ============================================================

write.table(
  rf_result,
  file = "TWB2_RF_EAS_prediction.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  eas,
  file = "TWB2_EAS.id",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

