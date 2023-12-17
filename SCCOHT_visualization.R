# Spatial resolution of cellular senescence dynamics in human colorectal liver metastasis
# aka "SCCOHT"; Sanders et al.
# GSE213699 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213699
# https://www.sciencedirect.com/science/article/pii/S2352578922001576
# https://github.com/kwells4/visium_ovarian_cancer
source("/Users/liaoyj/Documents/CPSC545_Project/SCCOHT_functions.R")
setwd("/Users/liaoyj/Library/CloudStorage/OneDrive-UBC/CPSC545_ProjectData/SCCOHT")

# Load sobj ----
samples
sample <- "D_GTFB1170_SmallCellOvarianCancer"
sobj <- get_sobj(sample, type = "SR") # Load sobj

# PCA dimplot ----
Reductions(sobj)
plt_SpatialPCA <- DimPlot(
  sobj, reduction = "pca", group.by = "Spatial_cluster") +
  labs(title = "Spatial PCA") +
  theme(legend.position = "bottom")
plt_SCTrrPCA <- DimPlot(
  sobj, reduction = "sctpca", group.by = "SCTrr_cluster") +
  labs(title = "SCT PCA (Rerun)") +
  theme(legend.position = "bottom")
plt_SPCAsparkPCA <- DimPlot(
  sobj, reduction = "spcasparkpca", group.by = "spcaspark_clusters") +
  labs(title = "Spatial PCA (SPARK Genes)") +
  theme(legend.position = "bottom")
plt_SPCAcustomPCA <- DimPlot(
  sobj, reduction = "spcacustompca", group.by = "spcacustom_clusters") +
  labs(title = "Spatial PCA (Seurat Variable Genes)") +
  theme(legend.position = "bottom")

# Patchworks
pw_PCA <- (plt_SpatialPCA | plt_SCTrrPCA | plt_SPCAsparkPCA | plt_SPCAcustomPCA) +
  plot_annotation(title = paste0(sample, " PCA"))
ggsave(
  plot = pw_PCA,
  filename = paste0("plots/DimPlotPCA_", sample, "_pw_spatial_sct_spcaspark_spcacustom.png"),
  width = 20, height = 7
)
pw_PCA <- (plt_SpatialPCA | plt_SCTrrPCA | plt_SPCAsparkPCA | plt_SPCAcustomPCA) +
  plot_annotation(title = paste0(sample, " PCA"))
# UMAP dimplot ----
Reductions(sobj)
plt_SpatialUMAP <- DimPlot(
  sobj, reduction = "umap", group.by = "Spatial_cluster") +
  labs(title = "Spatial UMAP") +
  theme(legend.position = "bottom")
plt_SCTrrUMAP <- DimPlot(
  sobj, reduction = "sctumap", group.by = "SCTrr_cluster") +
  labs(title = "SCT UMAP (Rerun)") +
  theme(legend.position = "bottom")
plt_SPCAsparkUMAP <- DimPlot(
  sobj, reduction = "spcasparkumap", group.by = "spcaspark_clusters") +
  labs(title = "Spatial UMAP (SPARK Genes)") +
  theme(legend.position = "bottom")
plt_SPCAcustomUMAP <- DimPlot(
  sobj, reduction = "spcacustomumap", group.by = "spcacustom_clusters") +
  labs(title = "Spatial UMAP (Seurat Variable Genes)") +
  theme(legend.position = "bottom")
# Patchworks
pw_UMAP <- (plt_SpatialUMAP | plt_SCTrrUMAP | plt_SPCAsparkUMAP | plt_SPCAcustomUMAP) +
  plot_annotation(title = paste0(sample, " UMAP"))
ggsave(
  plot = pw_UMAP,
  filename = paste0("plots/DimPlotUMAP_", sample,"_pw_spatial_sct_spcaspark_spcacustom.png"),
  width = 20, height = 7
)

# PCA + UMAP PW ----
pw_PCAUMAP <- pw_PCA / pw_UMAP
ggsave(
  plot = pw_PCAUMAP,
  filename = paste0("plots/DimPlotPCAUMAP_", sample,"_pw_spatial_sct_spcaspark_spcacustom.png"),
  width = 20, height = 12
)

## Clustering spatial plot ----
names(sobj@meta.data)
plt_SpatialCl <- SpatialDimPlot(
  sobj, group.by = "Spatial_cluster",
  label = TRUE, label.size = 3) +
  labs(
    title = "Spatial Cluster"
  )
plt_SCTrrCl <- SpatialDimPlot(
  sobj, group.by = "SCTrr_cluster",
  label = TRUE, label.size = 3) +
  labs(title = "SCT Cluster (Rerun)") +
  theme(legend.position = "bottom")
plt_SPCAspark <- SpatialDimPlot(
  sobj, group.by = "spcaspark_clusters",
  label = FALSE, label.size = 3) +
  labs(title = "Spatial PCA Cluster (SPARK Genes)") +
  theme(legend.position = "bottom")
plt_SPCAcustom <- SpatialDimPlot(
  sobj, group.by = "spcacustom_clusters",
  label = FALSE, label.size = 3) +
  labs(title = "Spatial PCA Cluster (Seurat Variable Genes)") +
  theme(legend.position = "bottom")

### Save ----
# Individual
ggsave(
  plot = plt_SCTrrCl,
  filename = paste0("plots/SpatialDimPlot_", sample,"_clusters_sctrr.png"),
  width = 6, height = 6
)

# Patchworks
pw_Spatial <- (plt_SpatialCl + plt_SCTrrCl) / (plt_SPCAspark + plt_SPCAcustom) +
  plot_annotation(title = paste0(sample, " Spatial Clustering"))
ggsave(
  plot = pw_Spatial,
  filename = paste0("plots/SpatialDimPlot_", sample,"_pw_clusters_spatial_sct_spatialpca.png"),
  width = 15, height = 15
)

pw2_Spatial <- (plt_SCTrrCl + plt_SPCAspark + plt_SPCAcustom) +
  plot_annotation(title = paste0(sample, " Spatial Clustering"))
ggsave(
  plot = pw2_Spatial,
  filename = paste0("plots/SpatialDimPlot_", sample,"_pw_clusters_sctrr_spcaspark_spcacustom.png"),
  width = 15, height = 8
)

pw_spca <- plt_SPCAspark + plt_SPCAcustom +
  plot_annotation(title = paste0(sample, " Spatial PCA Clustering"))
ggsave(
  plot = pw_spca,
  filename = paste0("plots/SpatialDimPlot_", sample,"_pw_clusters_spcaspark_spcacustom.png"),
  width = 12, height = 6
)

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

sample <- samples[8]
sample
ggsave(
  plot = plt_SCTrrcl,
  filename = paste0("plots/SpatialDimPlot_", sample, "_SCTrrcl.png"),
  width = 8, height = 6
)

## SingleR plots ----
names(sobj@meta.data)[grepl("cluster", names(sobj@meta.data))]
names(sobj@meta.data)[grepl("cellanno", names(sobj@meta.data))]

plt_SR_cell <- SpatialPlot(sobj, group.by = "cellanno") +
    labs(title = paste0("Cell SingleR"))
plt_SR_cell
# SCTrr cluster + anno
plt_clSCTrr <- SpatialPlot(sobj, group.by = "SCTrr_cluster") +
    labs(title = paste0("SCT Rerun Clusters"))
plt_SR_clSCTrr <- SpatialPlot(sobj, group.by = "cellanno_cl_SCTrr") +
    labs(title = paste0("SCT Rerun Cluster SingleR"))
plt_clSCTrr + plt_SR_clSCTrr
# SPCA SPARK cluster + anno
plt_clSPARK <- SpatialPlot(sobj, group.by = "spcaspark_clusters") +
    labs(title = paste0("SpatialPCA SPARK Clusters"))
plt_SR_clSPARK <- SpatialPlot(sobj, group.by = "cellanno_cl_spark") +
    labs(title = paste0("SpatialPCA SPARK Clusters SingleR"))
plt_clSPARK + plt_SR_clSPARK
# SPCA custom cluster + anno
plt_clcustom <- SpatialPlot(sobj, group.by = "spcacustom_clusters") +
    labs(title = paste0("SpatialPCA SVG Clusters"))
plt_SR_clcustom <- SpatialPlot(sobj, group.by = "cellanno_cl_custom") +
    labs(title = paste0("SpatialPCA SVG Clusters SingleR"))
plt_clcustom + plt_SR_clcustom