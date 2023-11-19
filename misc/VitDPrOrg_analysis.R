
# Vitamin D sufﬁciency enhances differentiation of patient-derived prostate epithelial organoids
# aka "VitDPrOrg"
# GSE159697 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159697
# https://www.sciencedirect.com/science/article/pii/S2589004220311718
setwd("/Users/liaoyj/Library/CloudStorage/OneDrive-UBC/CPSC545_Project/VitDPrOrg")

# Libraries ----
library(Seurat)

# Load raw data ----
# setwd("/Users/liaoyj/Library/CloudStorage/OneDrive-UBC/CPSC545_Project/VitDPrOrg/GSE159697_RAW")
# EtOH (vehicle)
debugonce(Read10X_Image)
img_EtOHVeh <- Read10X_Image(
  image.dir = "GSE159697_RAW",
  image.name = "GSM4837766_ETOH.png" # Doesn't work - wants specific PNG from spaceranger output
)

sobj_EtOHVeh <- Load10X_Spatial(
  data.dir = "GSE159697_RAW",
  filename = "GSM4837766_ETOH_filtered_feature_bc_matrix.h5",
  image = "GSM4837766_ETOH.png"
)

# D25