# Spatial resolution of cellular senescence dynamics in human colorectal liver metastasis
# aka "SCCOHT"; Sanders et al.
# GSE213699 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213699
# https://www.sciencedirect.com/science/article/pii/S2352578922001576
# https://github.com/kwells4/visium_ovarian_cancer
source("/Users/liaoyj/Documents/CPSC545_Project/SCCOHT_functions.R")
setwd("/Users/liaoyj/Library/CloudStorage/OneDrive-UBC/CPSC545_ProjectData/SCCOHT")

# Run once to generate "raw" data loaded straight from GEO
mclapply(
  samples,
  make_sobj,
  force = TRUE,
  mc.cores = cores
)

# Save all GEO sobjs as list ----
sobj_GEO_list <- lapply(
  samples,
  get_sobj,
  GEO = TRUE
)
names(sobj_GEO_list) <- samples
saveRDS(
  object = sobj_GEO_list,
  file = paste0(dir_data, "/sobj_GEO_list.rds")
)

# Preliminary plots ----
# Test plot
plot_SpatialPlot(
  sobj = sobj_GEO_list$A_955_OvarianTumor,
  feats = c("nCount_Spatial", "nFeature_Spatial"),
  assay = "SCT",
  savename = paste0(
    "plots/SpatialPlot_",
    Project(sobj_GEO_list$A_955_OvarianTumor), "_SCT_nCountnFeature.png")
)

# All sobj_GEOs SpatialPlots nCountnFeature ----
# (straight from GEO loading in)
# Spatial assay (raw counts)
mclapply(
  sobj_GEO_list,
  function(sobj, assay){
    plot_SpatialPlot(
      sobj = sobj,
      feats = c("nCount_Spatial", "nFeature_Spatial"),
      assay = assay,
      savename = paste0(
        "plots/SpatialPlot_",
        Project(sobj), "_", assay, "_nCountnFeature.png")
    )
  },
  assay = "Spatial",
  mc.cores = cores)
# SCT assay
mclapply(
  sobj_GEO_list,
  function(sobj, assay){
    plot_SpatialPlot(
      sobj = sobj,
      feats = c("nCount_Spatial", "nFeature_Spatial"),
      assay = assay,
      savename = paste0(
        "plots/SpatialPlot_",
        Project(sobj), "_", assay, "_nCountnFeature.png")
    )
  },
  assay = "SCT",
  mc.cores = cores)