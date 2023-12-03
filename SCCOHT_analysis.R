
# Spatial resolution of cellular senescence dynamics in human colorectal liver metastasis
# aka "SCCOHT"; Sanders et al.
# GSE213699 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213699
# https://www.sciencedirect.com/science/article/pii/S2352578922001576
# https://github.com/kwells4/visium_ovarian_cancer
setwd("/Users/liaoyj/Library/CloudStorage/OneDrive-UBC/CPSC545_ProjectData/SCCOHT")

# Libraries ----
library(data.table)
setDTthreads(4)
library(parallel)
cores <- detectCores() * 3/4 # overcomplicated way of saying "6"
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)

# Load samples ----
dir_data <- "/Users/liaoyj/Library/CloudStorage/OneDrive-UBC/CPSC545_ProjectData/SCCOHT/GSE213699_mod/"
samples <- list.files(dir_data)

# Run once to generate "raw" data loaded straight from GEO
make_sobj <- function(sample, force = FALSE){
  if(!("filtered_feature_bc_matrix.h5" %in% list.files(paste0(dir_data, "/", sample))) ||
     force == TRUE){
    # Load files
    cts_raw <- Rubrary::rread(
      paste0(dir_data, sample, "/raw_counts.csv"), row.names = 1, make.names = FALSE) %>%
      `colnames<-`(sub(".", "-", colnames(.), fixed = TRUE)) %>%
      as.matrix() %>%
      Matrix::Matrix(sparse = TRUE)
    # Create `filtered_feature_bc_matrix.h5` from this data
    # https://github.com/satijalab/seurat/issues/7157
    # https://rdrr.io/github/MarioniLab/DropletUtils/man/write10xCounts.html
    # https://samuel-marsh.github.io/scCustomize/reference/Create_10X_H5.html
    DropletUtils::write10xCounts(
      path = paste0(dir_data, sample, "/filtered_feature_bc_matrix/"),
      x = cts_raw,
      overwrite = TRUE
    )
    # Do h5 version
    DropletUtils::write10xCounts(
      path = paste0(dir_data, sample, "/filtered_feature_bc_matrix.h5"),
      x = cts_raw,
      overwrite = TRUE
    )
  }
  cts_SCT <- Rubrary::rread(
    paste0(dir_data, sample, "/sct_normalized_counts.csv"), row.names = 1, make.names = FALSE, to_df = TRUE) %>%
    `colnames<-`(sub(".", "-", colnames(.), fixed = TRUE))
  metadata <- Rubrary::rread(paste0(dir_data, sample, "/metadata.csv"), row.names = 1, make.names = FALSE)
  umap <- Rubrary::rread(paste0(dir_data, sample, "/umap_coords.csv"), row.names = 1, make.names = FALSE)
  # Read in & define image variables
  img <- Read10X_Image(image.dir = paste0(dir_data, sample, "/spatial/"))
  img@assay <- "Spatial"
  img@key <- "slice1_"
  # Create sobj from filtered_feature_bc_matrix
  sobj <- Load10X_Spatial(
    data.dir = paste0(dir_data, sample),
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "slice1"
  ) %>% # Add metadata + umap coords from Sanders et al.
    AddMetaData(metadata) %>%
    AddMetaData(umap)
  sobj@images$slice1 <- img # Patch in img
  # Add Sanders SCT transformed data
  sobj@assays[["SCT"]] <- CreateAssayObject(counts = cts_SCT) 
  Project(sobj) <- sample # Rename project
  Idents(sobj) <- Project(sobj) # Set ident
  saveRDS(
    sobj,
    file = paste0(dir_data, sample, "/sobj_GEO.rds")
  )
}

# mclapply(
#   samples,
#   make_sobj,
#   force = TRUE,
#   mc.cores = cores
# )

# Getter
get_sobj <- function(sample, GEO = FALSE){
  if(GEO){
    path <- paste0(dir_data, sample, "/sobj_GEO.rds")
  } else {
    path <- paste0(dir_data, sample, "/sobj.rds")
  }
  return(readRDS(path))
}

# Retrieve all GEO sobjs as list
sobj_GEO_list <- lapply(
  samples,
  get_sobj,
  GEO = TRUE
)
names(sobj_GEO_list) <- samples

# function to spit out preliminary nCount_Spatial + nFeature_Spatial spatial tissue plots
plot_SpatialPlot <- function(
    sobj, feats, assay = "SCT", savename = NULL, width = 12, height = 6){
  DefaultAssay(sobj) <- assay
  plt <- lapply(feats, function(feats) {
    SpatialPlot(
      sobj, features = feats, slot = "counts") +
      labs(title = feats) +
      theme(legend.position = "right")}) %>%
    wrap_plots() + plot_annotation(
      title = Project(sobj),
      subtitle = paste0("Assay: ", DefaultAssay(sobj)))
  if(!is.null(savename)){
    ggsave(
      plot = plt,
      filename = savename,
      width = width, height = height
    )
  }
  return(plt)
}

# Test plot
plot_SpatialPlot(
  sobj = sobj_GEO_list$A_955_OvarianTumor,
  feats = c("nCount_Spatial", "nFeature_Spatial"),
  assay = "SCT",
  savename = paste0(
    "plots/SpatialPlot_",
    Project(sobj_GEO_list$A_955_OvarianTumor), "_SCT_nCountnFeature.png")
)

# Generate all nCountnFeature SpatialPlots for all sobj_GEOs
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

# "D_GTFB1170_SmallCellOvarianCancer" is the one in the publication
sobj <- sobj_GEO_list$D_GTFB1170_SmallCellOvarianCancer

SpatialPlot(sobj, features = "SCT_cluster")

