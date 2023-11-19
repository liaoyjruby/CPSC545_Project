
# Spatial Gene Expression Analysis Reveals Characteristic Gene Expression Patterns
#   of De Novo Neuroendocrine Prostate Cancer Coexisting with Androgen Receptor Pathway Prostate Cancer
# aka "NEPCARPC"
# GSE230282 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230282
# https://www.mdpi.com/1422-0067/24/10/8955
setwd("/Users/liaoyj/Library/CloudStorage/OneDrive-UBC/CPSC545_ProjectData/NEPCARPC")

# Libraries ----
# BiocManager::install("DropletUtils")
library(DropletUtils)
library(Seurat)
library(dplyr)

# Load raw data ----
# setwd("/Users/liaoyj/Library/CloudStorage/OneDrive-UBC/CPSC545_ProjectData/NEPCARPC/GSE230282_RAW/")

# Create filtered_feature_bc_matrix
filt_mtx <- Seurat::Read10X("filtered_feature_bc_matrix/")
DropletUtils::write10xCounts(
  path = "data_SR/filtered_feature_bc_matrix.h5", x = filt_mtx, type = "HDF5",
  genome = "GRCh38", version = "3", overwrite = TRUE,
  gene.id = rownames(filt_mtx),
  gene.symbol = rownames(filt_mtx))

debugonce(Read10X_Image)
img_NEPCARPC <- Read10X_Image(
  image.dir = "data_SR/spatial"
)

coords_NEPCARPC <- read.csv("data_SR/spatial/tissue_positions.csv", row.names = 1) %>%
  mutate(
    imagerow = pxl_row_in_fullres,
    imagecol = pxl_col_in_fullres
  )

# https://github.com/satijalab/seurat/issues/3595
img_NEPCARPC <- new(
  Class = 'VisiumV1',
  image = png::readPNG("data_SR/spatial/tissue_hires_image.png"),
  scale.factors = scalefactors(spot = 1, fiducial = 1, hires = 1, lowres = 1),
  coordinates = coords_NEPCARPC
)

sobj_NEPCARPC <- Load10X_Spatial(
  data.dir = "data_SR",
  filename = "filtered_feature_bc_matrix.h5",
  image = img_NEPCARPC
)

Assays(sobj_NEPCARPC)

VlnPlot(sobj_NEPCARPC, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
SpatialFeaturePlot(sobj_NEPCARPC, features = "nCount_Spatial")
SpatialFeaturePlot(sobj_NEPCARPC, features = c("GREM1", "HMOX1"))


