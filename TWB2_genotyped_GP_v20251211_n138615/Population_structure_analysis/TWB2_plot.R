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

# =========== TWB2 FAM ===========

TWB2 <- read.delim(
  "/SATA2TB/Taiwan_Biobank/TWB2_genotyped_GP_v20251211_n138615/King/TWB2_S05.fam",
  header = FALSE,
  sep = "",
  quote = ""
)

twb2fam <- data.frame(
  IID = TWB2$V2
)

rm(TWB2)

# =========== TWB2 PCA ===========

twb2pca <- read.table(
  "/SATA2TB/Taiwan_Biobank/TWB2_genotyped_GP_v20251211_n138615/Population_structure_analysis/TWB2_pca.eigenvec",
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(twb2pca)[1:4] <- c(
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
twb2pca <- assign_1kg_color(twb2pca, pop_1kgp)

# =========== TWB2 ===========

setwd("/SATA2TB/Taiwan_Biobank/TWB2_genotyped_GP_v20251211_n138615/Population_structure_analysis")

if (!dir.exists("plot")) {
  dir.create("plot")
}

idx_1kg <- twb2pca$IID %in% pop_1kgp$Sample.name
idx_twb <- twb2pca$IID %in% twb2fam$IID


## =========== TWB2 - color ===========

twb2pca$color <- "black"

twb2pca$color[
  twb2pca$IID %in%
    pop_1kgp$Sample.name[
      pop_1kgp$Superpopulation.code == "AFR"
    ]
] <- colors["AFR"]

twb2pca$color[
  twb2pca$IID %in%
    pop_1kgp$Sample.name[
      pop_1kgp$Superpopulation.code == "AMR"
    ]
] <- colors["AMR"]

twb2pca$color[
  twb2pca$IID %in%
    pop_1kgp$Sample.name[
      pop_1kgp$Superpopulation.code == "EAS"
    ]
] <- colors["EAS"]

twb2pca$color[
  twb2pca$IID %in%
    pop_1kgp$Sample.name[
      pop_1kgp$Superpopulation.code == "EUR"
    ]
] <- colors["EUR"]

twb2pca$color[
  twb2pca$IID %in%
    pop_1kgp$Sample.name[
      pop_1kgp$Superpopulation.code == "SAS"
    ]
] <- colors["SAS"]


## =========== TWB2 - plot 1KG 5 pop ===========

idx_1kg_twb2 <- twb2pca$IID %in% pop_1kgp$Sample.name

png(
  "plot/TWB2_1kg.png",
  width = 2000,
  height = 2000,
  units = "px",
  pointsize = 56
)

plot(
  twb2pca[idx_1kg_twb2, c("PC1", "PC2")],
  col = twb2pca$color[idx_1kg_twb2],
  pch = 1,
  lwd = 3
)

legend(
  "bottomleft",
  legend = c(
    "AFR",
    "AMR",
    "EAS",
    "EUR",
    "SAS"
  ),
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


## =========== TWB2 - plot 1KG 5 pop + TWB2 ===========

idx_1kg_twb2 <- twb2pca$IID %in% pop_1kgp$Sample.name
idx_twb2 <- twb2pca$IID %in% twb2fam$IID

all_idx <- idx_1kg_twb2 | idx_twb2

png(
  "plot/TWB2_1kg_twb2.png",
  width = 2000,
  height = 2000,
  units = "px",
  pointsize = 56
)

plot(
  twb2pca[all_idx, c("PC1", "PC2")],
  type = "n",
  xlab = "PC1",
  ylab = "PC2"
)

# 先畫 1kGP
points(
  twb2pca[idx_1kg_twb2, c("PC1", "PC2")],
  col = twb2pca$color[idx_1kg_twb2],
  pch = 1,
  lwd = 3
)

# 再畫 TWB2
points(
  twb2pca[idx_twb2, c("PC1", "PC2")],
  col = "gray50",
  pch = 4,
  lwd = 3
)

legend(
  "bottomleft",
  legend = c(
    "TWB2",
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
  pch = c(
    4,
    15,
    15,
    15,
    15,
    15
  )
)

dev.off()

