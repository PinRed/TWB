library(scales)

# =========== 1kGP metadata =========== 

pop_1kgp <- read.table(
  "/SATA2TB/Taiwan_Biobank/1kGP_SGDP_prep/metadata/igsr_samples.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

# EAS
iseas <- data.frame(
  IID = pop_1kgp$Sample.name[
    pop_1kgp$Superpopulation.code == "EAS"
  ]
)

# non-EAS
noteas <- data.frame(
  IID = pop_1kgp$Sample.name[
    pop_1kgp$Superpopulation.code != "EAS"
  ]
)


# =========== TWB1 FAM =========== 
setwd("/SATA2TB/Taiwan_Biobank/TWB1_genotyped_GP_v20251211_n27652/Population_structure_analysis")
TWB1 <- read.delim(
  "../King/TWB1_S05.fam",
  header = FALSE,
  sep = "",
  quote = ""
)

twb1fam <- data.frame(
  IID = TWB1$V2
)

rm(TWB1)

# =========== TWB1 PCA =========== 

twb1pca <- read.table(
  "meta_twb1_pca/TWB1_S05_pca.eigenvec",
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(twb1pca)[1:4] <- c(
  "FID",
  "IID",
  "PC1",
  "PC2"
)



# =========== colors ===========

colors <- c(
  AFR = "dodgerblue4",
  AMR = "gold",
  EAS = "forestgreen",
  EUR = "firebrick",
  SAS = "darkorchid"
)


# =========== function: assign 1kGP colors ===========

assign_1kg_color <- function(pca, pop_1kgp) {
  
  pca$color <- "black"
  
  pca$color[
    pca$IID %in%
      pop_1kgp$Sample.name[
        pop_1kgp$Superpopulation.code == "AFR"
      ]
  ] <- colors["AFR"]
  
  pca$color[
    pca$IID %in%
      pop_1kgp$Sample.name[
        pop_1kgp$Superpopulation.code == "AMR"
      ]
  ] <- colors["AMR"]
  
  pca$color[
    pca$IID %in%
      pop_1kgp$Sample.name[
        pop_1kgp$Superpopulation.code == "EAS"
      ]
  ] <- colors["EAS"]
  
  pca$color[
    pca$IID %in%
      pop_1kgp$Sample.name[
        pop_1kgp$Superpopulation.code == "EUR"
      ]
  ] <- colors["EUR"]
  
  pca$color[
    pca$IID %in%
      pop_1kgp$Sample.name[
        pop_1kgp$Superpopulation.code == "SAS"
      ]
  ] <- colors["SAS"]
  
  return(pca)
}


# =========== apply colors separately ===========
twb1pca <- assign_1kg_color(twb1pca, pop_1kgp)


# =========== TWB1 ===========

if (!dir.exists("plot")) {
  dir.create("plot")
}

idx_1kg <- twb1pca$IID %in% pop_1kgp$Sample.name
idx_twb <- twb1pca$IID %in% twb1fam$IID



## ---------- TWB1 color ----------
twb1pca$color <- "black"

twb1pca$color[
  twb1pca$IID %in%
    pop_1kgp$Sample.name[
      pop_1kgp$Superpopulation.code == "AFR"
    ]
] <- colors["AFR"]

twb1pca$color[
  twb1pca$IID %in%
    pop_1kgp$Sample.name[
      pop_1kgp$Superpopulation.code == "AMR"
    ]
] <- colors["AMR"]

twb1pca$color[
  twb1pca$IID %in%
    pop_1kgp$Sample.name[
      pop_1kgp$Superpopulation.code == "EAS"
    ]
] <- colors["EAS"]

twb1pca$color[
  twb1pca$IID %in%
    pop_1kgp$Sample.name[
      pop_1kgp$Superpopulation.code == "EUR"
    ]
] <- colors["EUR"]

twb1pca$color[
  twb1pca$IID %in%
    pop_1kgp$Sample.name[
      pop_1kgp$Superpopulation.code == "SAS"
    ]
] <- colors["SAS"]


## =========== TWB1 - plot 1KG 5 pop ===========

idx_1kg_twb1 <- twb1pca$IID %in% pop_1kgp$Sample.name

png(
  "plot/TWB1_1kg.png",
  width = 2000,
  height = 2000,
  units = "px",
  pointsize = 56
)

plot(
  twb1pca[idx_1kg_twb1, c("PC1", "PC2")],
  col = twb1pca$color[idx_1kg_twb1],
  pch = 1,
  lwd = 3
)

legend(
  "bottomleft",
  legend = c("AFR", "AMR", "EAS", "EUR", "SAS"),
  col = c(
    colors["AFR"],
    colors["AMR"],
    colors["EAS"],
    colors["EUR"],
    colors["SAS"]
  ),
  pch = 15
)

dev.off()

## =========== TWB1 - plot 1KG 5 pop + TWB1 ===========

idx_1kg_twb1 <- twb1pca$IID %in% pop_1kgp$Sample.name
idx_twb1 <- twb1pca$IID %in% twb1fam$IID

all_idx <- idx_1kg_twb1 | idx_twb1

png(
  "plot/TWB1_1kg_twb1.png",
  width = 2000,
  height = 2000,
  units = "px",
  pointsize = 56
)

plot(
  twb1pca[all_idx, c("PC1", "PC2")],
  type = "n",
  xlab = "PC1",
  ylab = "PC2"
)

# 先畫 1kGP
points(
  twb1pca[idx_1kg_twb1, c("PC1", "PC2")],
  col = twb1pca$color[idx_1kg_twb1],
  pch = 1,
  lwd = 3
)

# 再畫 TWB1
points(
  twb1pca[idx_twb1, c("PC1", "PC2")],
  col = "gray50",
  pch = 4,
  lwd = 3
)

legend(
  "bottomleft",
  legend = c(
    "TWB1",
    "AFR",
    "AMR",
    "EAS",
    "EUR",
    "SAS"
  ),
  col = c(
    "gray50",
    colors["AFR"],
    colors["AMR"],
    colors["EAS"],
    colors["EUR"],
    colors["SAS"]
  ),
  pch = c(4, 15, 15, 15, 15, 15)
)

dev.off()
