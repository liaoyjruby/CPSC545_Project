
# Spatial resolution of cellular senescence dynamics in human colorectal liver metastasis
# aka "SCCOHT"; Sanders et al.
# GSE213699 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213699
# https://www.sciencedirect.com/science/article/pii/S2352578922001576
# https://github.com/kwells4/visium_ovarian_cancer
source("/Users/liaoyj/Documents/CPSC545_Project/SCCOHT_functions.R")
setwd("/Users/liaoyj/Library/CloudStorage/OneDrive-UBC/CPSC545_ProjectData/SCCOHT")

# Process sample D_GTFB1170 ----
samples
sample <- "D_GTFB1170_SmallCellOvarianCancer"
sobj <- process_sobj(sample) # Run SCT Seurat pipeline
sobj <- process_SPCA(sample) # Run SpatialPCA pipeline
sobj <- process_SingleR(sample) # Run SingleR cell + cluster type annotation
sobj <- get_sobj(sample, type = "SR") # Load sobj

# Cluster characterization
mrkrs_lst <- process_markers(sample)
names(mrkrs_lst)

# Clustifyr
# https://github.com/rnabioco/clustifyr
# SCTrr
DefaultAssay(sobj) <- "SCTrr"
mtx_SCTrr <- GetAssayData(sobj, layer = "data") %>%
  as.data.frame()
clust <- clustifyr::clustify(
  input = mtx_SCTrr,
  metadata = sobj@meta.data$SCTrr_cluster,
)

# make ref for clustifyr?
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147082

SpatialPlot(sobj, group.by = "cellanno_cl_SCTrr")
SpatialPlot(sobj, group.by = "cellanno_cl_spark") +
SpatialPlot(sobj, group.by = "cellanno_cl_custom")
