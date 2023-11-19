# Spatial resolution of cellular senescence dynamics in human colorectal liver metastasis
# aka "LiverMeta"
# GSE206552 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206552
# https://onlinelibrary.wiley.com/doi/10.1111/acel.13853
setwd("/Users/liaoyj/Library/CloudStorage/OneDrive-UBC/CPSC545_ProjectData/LiverMeta")

# Libraries ----
library(Seurat)
library(ggplot2)
library(patchwork)

# Load data ----
sobj_meta1 <- Load10X_Spatial(
  data.dir = "GSE206552_RAW/GSM6256810_meta1",
  filename = "filtered_feature_bc_matrix.h5"
)
sobj_meta2 <- Load10X_Spatial(
  data.dir = "GSE206552_RAW/GSM6256811_meta2",
  filename = "filtered_feature_bc_matrix.h5"
)
sobj_meta3 <- Load10X_Spatial(
  data.dir = "GSE206552_RAW/GSM6256812_meta3",
  filename = "filtered_feature_bc_matrix.h5"
)
sobj_meta4 <- Load10X_Spatial(
  data.dir = "GSE206552_RAW/GSM6256813_meta4",
  filename = "filtered_feature_bc_matrix.h5"
)

# Preliminary plots ----
m1_nct <- SpatialFeaturePlot(sobj_meta1, features = "nCount_Spatial")
m2_nct <- SpatialFeaturePlot(sobj_meta2, features = "nCount_Spatial")
m3_nct <- SpatialFeaturePlot(sobj_meta3, features = "nCount_Spatial")
m4_nct <- SpatialFeaturePlot(sobj_meta4, features = "nCount_Spatial")

pw_nct <- (m1_nct | m2_nct) / (m3_nct | m4_nct)
pw_nct
