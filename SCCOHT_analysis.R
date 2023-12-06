
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
samples <- list.files(dir_data)[grepl("^[A-D]_", list.files(dir_data))]

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
  img <- Read10X_Image(
    image.dir = paste0(dir_data, sample, "/spatial/"),
    filter.matrix = TRUE)
  img@assay <- "Spatial"
  img@key <- "slice1_"
  # Create sobj from filtered_feature_bc_matrix
  sobj <- Load10X_Spatial(
    data.dir = paste0(dir_data, sample),
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    filter.matrix = TRUE,
    slice = "slice1"
  ) %>% # Add metadata + umap coords from Sanders et al.
    AddMetaData(metadata) %>%
    AddMetaData(umap)
  sobj@images$slice1 <- img # Patch in img
  # Add Sanders SCT transformed data
  sobj@assays[["SCT"]] <- CreateAssayObject(counts = cts_SCT, key = "sct_") 
  Project(sobj) <- sample # Rename project
  Idents(sobj) <- Project(sobj) # Set ident
  saveRDS(
    sobj,
    file = paste0(dir_data, sample, "/sobj_GEO.rds")
  )
}

mclapply(
  samples,
  make_sobj,
  force = TRUE,
  mc.cores = cores
)

# Getter
get_sobj <- function(sample, GEO = FALSE){
  if(GEO){
    path <- paste0(dir_data, sample, "/sobj_GEO.rds")
  } else {
    path <- paste0(dir_data, sample, "/sobj.rds")
  }
  return(readRDS(path))
}

# Retrieve all GEO sobjs as list ----
# sobj_GEO_list <- lapply(
#   samples,
#   get_sobj,
#   GEO = TRUE
# )
# names(sobj_GEO_list) <- samples
# saveRDS(
#   object = sobj_GEO_list,
#   file = paste0(dir_data, "/sobj_GEO_list.rds")
# )

sobj_GEO_list <- readRDS(paste0(dir_data, "/sobj_GEO_list.rds"))

# Preliminary plots ----
# function to spit out nCount_Spatial + nFeature_Spatial spatial tissue plots
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

# sobj processing functions ----
make_sobj <- function(sobj_GEO){
  SpatialPlot(sobj, group.by = "SCT_cluster") # UMAP clusters not given?
  ggsave(
    plot = sobj
  )
}

# Sample D_GTFB1170 ----
# "D_GTFB1170_SmallCellOvarianCancer" is the one in the publication
sobj_GEO <- sobj_GEO_list$D_GTFB1170_SmallCellOvarianCancer
Keys(sobj_GEO)
# Rerub generic Seurat pipeline for this sample
# https://github.com/kwells4/visium_ovarian_cancer/tree/main/src/scripts/D_GTFB1170_SmallCellOvarianCancer
sobj_SCT <- sobj_GEO %>% SCTransform(assay = "Spatial", new.assay.name = "SCTrr")
sobj <- sobj_SCT
# sobj <- readRDS("GSE213699_mod/D_GTFB1170_SmallCellOvarianCancer/sobj.rds")
# sobj <- sobj
## Initial processing ----
# https://github.com/kwells4/visium_ovarian_cancer/blob/main/src/scripts/D_GTFB1170_SmallCellOvarianCancer/01_Initial_processing.R
assay = "Spatial"
DefaultAssay(sobj) <- assay
RNA_nPCs = 33
resolution = 0.5
set.seed(0) # from GH code
sobj <- sobj %>%
  NormalizeData(normalization.method = "LogNormalize") %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = NULL)%>%
  RunPCA(assay = "Spatial", reduction.name = "pca",
         reduction.key = paste0(assay, "PC_")) %>%
  FindNeighbors(reduction = "pca", dims = 1:RNA_nPCs) %>%
  FindClusters(resolution = resolution) %>%
  RunUMAP(
    metric = "correlation",
    dims = 1:RNA_nPCs,
    reduction = "pca",
    reduction.name = "umap",
    reduction.key = paste0(assay, "UMAP_"),
    assay = assay
  )
sobj[["Spatial_cluster"]] <- Idents(sobj)
# https://github.com/kwells4/visium_ovarian_cancer/blob/main/src/scripts/D_GTFB1170_SmallCellOvarianCancer/02_PCA_UMAP.R
## SCT assay ----
assay = "SCTrr"
DefaultAssay(sobj) <- assay
sobj <- sobj %>%
  RunPCA(
    assay = assay,
    reduction.name = "sctpca",
    reduction.key = paste0(assay, "PC_")) %>%
  FindNeighbors(reduction = "sctpca", dims = 1:RNA_nPCs) %>%
  FindClusters(resolution = resolution) %>%
  RunUMAP(
    metric = "correlation",
    dims = 1:RNA_nPCs,
    reduction = "sctpca",
    reduction.name = "sctumap",
    reduction.key = paste0(assay, "UMAP_"),
    assay = assay
  )
sobj[["SCTrr_cluster"]] <- Idents(sobj)

## SpatialPCA ----
# Linux paths
sobj_GEO <- readRDS("~/data/sobj_GEO.rds")
n_cores <- ceiling(parallel::detectCores()/2)
# BiocManager::install("xzhoulab/SPARK")
# BiocManager::install("shangll123/SpatialPCA")
library(Seurat)
library(SpatialPCA)
library(dplyr)
sobj_GEO <- sobj_GEO %>%
  NormalizeData(normalization.method = "LogNormalize") %>%
  FindVariableFeatures() # nfeatures = 2000

cts <- sobj_GEO@assays$Spatial$counts
# View(sobj_GEO@meta.data)
loc <- sobj_GEO@meta.data %>%
  select(row, col) %>%
  as.matrix()

### Create initial objects ----
SPCA_obj <- SpatialPCA::CreateSpatialPCAObject(
  counts = cts,
  location = loc,
  project = Project(sobj_GEO),
  gene.type = "spatial",
  numCores_spark = n_cores
)
saveRDS(
  SPCA_obj,
  file = "~/data/SPCA_obj.rds"
)

SPCA_obj_genelist <- SpatialPCA::CreateSpatialPCAObject(
  counts = cts,
  location = loc,
  project = Project(sobj_GEO),
  numCores_spark = n_cores,
  gene.type = "custom",
  customGenelist = VariableFeatures(sobj_GEO)
)
saveRDS(
  SPCA_obj_genelist,
  file = "~/data/SPCA_obj_genelist.rds"
)

SPCA_obj <- readRDS("~/data/SPCA_obj.rds")
SPCA_obj_genelist <- readRDS("~/data/SPCA_obj_genelist.rds")

### Estimate spatial PCs ----
SPCA_obj_SPCs <- SPCA_obj %>%
  SpatialPCA_buildKernel( # Defaults, for smaller datasets
    kerneltype = "gaussian",
    bandwidthtype = "SJ",
    bandwidth.set.by.user = NULL,
    sparseKernel = TRUE,
    sparseKernel_tol = 1e-20,
    sparseKernel_ncore = n_cores) %>%
  SpatialPCA_EstimateLoading(
    maxiter = 300,
    initial_tau = 1,
    fast = TRUE,
    SpatialPCnum = 20) %>%
  SpatialPCA_SpatialPCs(fast = TRUE)
saveRDS(
  SPCA_obj_SPCs,
  file = "~/data/SPCA_obj_SPCs.rds"
)

SPCA_obj_genelist_SPCs <- SPCA_obj_genelist %>%
  SpatialPCA_buildKernel( # Defaults, for smaller datasets
    kerneltype = "gaussian",
    bandwidthtype = "SJ",
    bandwidth.set.by.user = NULL,
    sparseKernel = TRUE,
    sparseKernel_tol = 1e-20,
    sparseKernel_ncore = n_cores) %>%
  SpatialPCA_EstimateLoading(
    maxiter = 300,
    initial_tau = 1,
    fast = TRUE,
    SpatialPCnum = 20) %>%
  SpatialPCA_SpatialPCs(fast = TRUE)
saveRDS(
  SPCA_obj_genelist_SPCs,
  file = "~/data/SPCA_obj_genelist_SPCs.rds"
)

### Write SPCs to text file ----
SPCA_obj_SPCs_df <- SPCA_obj_SPCs@SpatialPCs %>%
  t() %>%
  as.data.frame() %>%
  `rownames<-`(SPCA_obj_SPCs@counts@Dimnames[[2]])
Rubrary::rwrite(
  SPCA_obj_SPCs_df,
  "~/data/SPCA_obj_SPCs_df.txt"
)

SPCA_obj_genelist_SPCs_df <- SPCA_obj_genelist_SPCs@SpatialPCs %>%
  t() %>%
  as.data.frame() %>%
  `rownames<-`(SPCA_obj_genelist_SPCs@counts@Dimnames[[2]])
Rubrary::rwrite(
  SPCA_obj_genelist_SPCs_df,
  "~/data/SPCA_obj_genelist_SPCs_df.txt"
)

# Overwrite RDS ----
saveRDS(
  object = sobj,
  file = "GSE213699_mod/D_GTFB1170_SmallCellOvarianCancer/sobj.rds"
)
sobj <- get_sobj("D_GTFB1170_SmallCellOvarianCancer")

# Visualization ----
# View(sobj@meta.data)
SpatialDimPlot(
  sobj, group.by = "SCTrr_cluster",
  label = TRUE, label.size = 3)

## Spatial plot ----
# SCT reclustering largely similar
plt_SCTcl <- SpatialDimPlot(
  sobj, group.by = "SCT_cluster",
  label = TRUE, label.size = 3) +
  labs(
    title = "SCT Cluster - GEO"
  )
plt_SCTrrcl <- SpatialDimPlot(
  sobj, group.by = "SCTrr_cluster",
  label = TRUE, label.size = 3) +
  labs(
    title = "SCT Cluster - Rerun"
  )

plt_SCTcl + plt_SCTrrcl

## UMAP ----
Reductions(sobj)
DimPlot(
  sobj, group.by = "SCT_cluster"
)

# Plot UMAP
sobj@meta.data$SCT_cluster <- factor(
  sobj@meta.data$SCT_cluster,
  levels = sort(unique(sobj@meta.data$SCT_cluster))
)
ggplot(
  sobj@meta.data,
  aes(x = UMAP_1, y = UMAP_2, color = SCT_cluster)
) +
  geom_point() +
  theme_classic()

DimPlot(sobj)

# SingleR ----
# https://bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html
# https://bioconductor.org/books/devel/SingleRBook/
library(celldex) # https://bioconductor.org/packages/3.18/data/experiment/html/celldex.html
library(SingleR)
library(BiocParallel)
# https://bioconductor.org/packages/3.18/data/experiment/vignettes/celldex/inst/doc/userguide.html
# Human Primary Cell Atlas
hpca.se <- celldex::HumanPrimaryCellAtlasData()
# Get numeric log-expression normalized df
library(scRNAseq)
smp <- samples[8]
sobj_SR <- get_sobj(smp)
sobj_SR[["barcode"]] <- rownames(sobj_SR@meta.data)
View(sobj_SR@meta.data)
rns <- rownames(sobj_SR)
sce <- Seurat::as.SingleCellExperiment(sobj_SR, assay = "Spatial")
scater::plotUMAP(sce, color_by = "SCTrr_cluster")
# Per cell annotation
cellanno <- SingleR::SingleR(
  test = sce,
  ref = hpca.se,
  labels = hpca.se$label.main,
  BPPARAM = MulticoreParam()
)
saveRDS(
  cellanno,
  file = "GSE213699_mod/D_GTFB1170_SmallCellOvarianCancer/cellanno_SR_percell.rds"
)
sobj_SR[["cellanno"]] <- cellanno$labels
sobj_SR[["cellanno_pruned"]] <- cellanno$pruned.labels
View(sobj_SR@meta.data)

cellanno <- SingleR::SingleR(
  test = sce,
  ref = hpca.se,
  labels = hpca.se$label.main,
  BPPARAM = MulticoreParam()
)
# Cluster level annotation
# https://bioconductor.org/books/devel/SingleRBook/advanced-options.html#cluster-level-annotation
colLabels(sce) <- colData(sce)$SCTrr_cluster
cellanno_cl <- SingleR::SingleR(
  test = sce,
  ref = hpca.se,
  clusters = colLabels(sce),
  labels = hpca.se$label.main,
  # BPPARAM = MulticoreParam()
) 

saveRDS(
  cellanno_cl,
  file = "GSE213699_mod/D_GTFB1170_SmallCellOvarianCancer/cellanno_SR_perSCTrrcluster.rds"
)

cellanno_cl_df <- cellanno_cl %>%
  as.data.frame() %>%
  select(labels, delta.next, pruned.labels) %>%
  `colnames<-`(c("cellanno_cl", "cellanno_cl_deltanext", "cellanno_cl_pruned")) %>%
  tibble::rownames_to_column("SCTrr_cluster")
# Add to metadata
sobj_SR@meta.data <- sobj_SR@meta.data %>%
  left_join(cellanno_cl_df, by = "SCTrr_cluster")
rownames(sobj_SR@meta.data) <- sobj_SR$barcode
View(sobj_SR@meta.data)

plt_spatial_ca <- SpatialPlot(
  sobj_SR, group.by = "cellanno_cl",
)
plt_spatial_ca
Reductions(sobj_SR)
plt_umap_SCTrrcl <- DimPlot(
  sobj_SR, group.by = "SCTrr_cluster",
  reduction = "umap"
)
plt_umap_SCTrrcl
plt_umap_ca <- DimPlot(
  sobj_SR, group.by = "cellanno_cl",
  reduction = "umap"
)
plt_umap_ca

plt_spatial_ca + plt_umap_ca
