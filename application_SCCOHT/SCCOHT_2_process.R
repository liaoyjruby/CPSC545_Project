
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

# Visualization ----
SpatialPlot(sobj, group.by = "cellanno_cl_SCTrr")
SpatialPlot(sobj, group.by = "cellanno_cl_spark") +
SpatialPlot(sobj, group.by = "cellanno_cl_custom")

## SCT Clusters GEO vs SCTrr ----
sobj@meta.data$SCTrr_cluster <- factor(
    sobj@meta.data$SCTrr_cluster,
    levels = as.character(0:13))
plt_SCTcl <- SpatialDimPlot(
  sobj, group.by = "SCT_cluster",
  label = TRUE, label.size = 3) +
  labs(title = "SCT Cluster - GEO") +
  theme(legend.position = "bottom")
plt_SCTrrcl <- SpatialDimPlot(
  sobj, group.by = "SCTrr_cluster",
  label = TRUE, label.size = 3) +
  labs(title = "SCT Cluster - Rerun") +
  theme(legend.position = "bottom")

plt_SCTcl_GEO_rr <- plt_SCTcl + plt_SCTrrcl

ggsave(
  plot = plt_SCTcl_GEO_rr,
  filename = paste0("plots/SpatialDimPlot_", sample, "_SCTrrcl.png"),
  width = 10, height = 6
)

## PCA dimplot ----
Reductions(sobj)
plt_SpatialPCA <- DimPlot(
  sobj, reduction = "pca", group.by = "Spatial_cluster") +
  labs(title = "No SCT PCA") +
  theme(legend.position = "bottom")
plt_SCTrrPCA <- DimPlot(
  sobj, reduction = "sctpca", group.by = "SCTrr_cluster") +
  labs(title = "SCTrr PCA") +
  theme(legend.position = "bottom")
plt_SPCAsparkPCA <- DimPlot(
  sobj, reduction = "spcasparkpca", group.by = "spcaspark_clusters") +
  labs(title = "Spatial PCA SPARK") +
  theme(legend.position = "bottom")
plt_SPCAcustomPCA <- DimPlot(
  sobj, reduction = "spcacustompca", group.by = "spcacustom_clusters") +
  labs(title = "Spatial PCA SVF") +
  theme(legend.position = "bottom")

# Patchworks
pw_PCA <- (plt_SpatialPCA | plt_SCTrrPCA | plt_SPCAsparkPCA | plt_SPCAcustomPCA) +
  plot_annotation(title = paste0(sample, " PCA"))
ggsave(
  plot = pw_PCA,
  filename = paste0("plots/DimPlotPCA_", sample, "_pw_spatial_sct_spcaspark_spcacustom.png"),
  width = 20, height = 7
)

pw_PCA2 <- (plt_SCTrrPCA | plt_SPCAsparkPCA | plt_SPCAcustomPCA)
ggsave(
  plot = pw_PCA2,
  filename = paste0("plots/DimPlotPCA_", sample, "_pw_strip.png"),
  width = 15, height = 6
)

pw_PCA_sq <- (plt_SpatialPCA + plt_SCTrrPCA + plt_SPCAsparkPCA + plt_SPCAcustomPCA) +
    plot_layout(ncol = 2)
pw_PCA_sq
ggsave(
  plot = pw_PCA_sq,
  filename = paste0("plots/DimPlotPCA_", sample, "_pw_sq.png"),
  width = 9, height = 10
)

## UMAP dimplot ----
Reductions(sobj)
plt_SpatialUMAP <- DimPlot(
  sobj, reduction = "umap", group.by = "Spatial_cluster") +
  labs(title = "No SCT UMAP") +
  theme(legend.position = "bottom")
plt_SCTrrUMAP <- DimPlot(
  sobj, reduction = "sctumap", group.by = "SCTrr_cluster") +
  labs(title = "SCTrr UMAP") +
  theme(legend.position = "bottom")
plt_SPCAsparkUMAP <- DimPlot(
  sobj, reduction = "spcasparkumap", group.by = "spcaspark_clusters") +
  labs(title = "Spatial UMAP SPARK") +
  theme(legend.position = "bottom")
plt_SPCAcustomUMAP <- DimPlot(
  sobj, reduction = "spcacustomumap", group.by = "spcacustom_clusters") +
  labs(title = "Spatial UMAP SVF") +
  theme(legend.position = "bottom")
# Patchworks
pw_UMAP <- (plt_SpatialUMAP | plt_SCTrrUMAP | plt_SPCAsparkUMAP | plt_SPCAcustomUMAP) +
  plot_annotation(title = paste0(sample, " UMAP"))
ggsave(
  plot = pw_UMAP,
  filename = paste0("plots/DimPlotUMAP_", sample, "_pw_spatial_sct_spcaspark_spcacustom.png"),
  width = 20, height = 7
)

pw_UMAP2 <- (plt_SCTrrUMAP | plt_SPCAsparkUMAP | plt_SPCAcustomUMAP)
ggsave(
  plot = pw_UMAP2,
  filename = paste0("plots/DimPlotUMAP_", sample, "_pw_strip.png"),
  width = 15, height = 6
)

pw_UMAP_sq <- (plt_SpatialUMAP + plt_SCTrrUMAP + plt_SPCAsparkUMAP + plt_SPCAcustomUMAP) + plot_layout(ncol = 2)
pw_UMAP_sq
ggsave(
  plot = pw_UMAP_sq,
  filename = paste0("plots/DimPlotUMAP_", sample, "_pw_sq.png"),
  width = 9, height = 10
)

## PCA + UMAP PW ----
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
  labs(title = "No SCT Cluster")
plt_SCTrrCl <- SpatialDimPlot(
  sobj, group.by = "SCTrr_cluster",
  label = TRUE, label.size = 3) +
  labs(title = "SCTrr Cluster") +
  theme(legend.position = "bottom")
plt_SPCAspark <- SpatialDimPlot(
  sobj, group.by = "spcaspark_clusters",
  label = FALSE, label.size = 3) +
  labs(title = "SpatialPCA Cluster SPARK") +
  theme(legend.position = "bottom")
plt_SPCAsvf <- SpatialDimPlot(
  sobj, group.by = "spcacustom_clusters",
  label = FALSE, label.size = 3) +
  labs(title = "SpatialPCA Cluster SVF") +
  theme(legend.position = "bottom")

### Save ----
# Individual
ggsave(
  plot = plt_SCTrrCl,
  filename = paste0("plots/SpatialDimPlot_", sample,"_clusters_sctrr.png"),
  width = 6, height = 6
)

# Patchworks
pw_Spatial <- (plt_SpatialCl + plt_SCTrrCl) / (plt_SPCAspark + plt_SPCAsvf) +
  plot_annotation(title = paste0(sample, " Spatial Clustering"))
pw_Spatial
ggsave(
  plot = pw_Spatial,
  filename = paste0("plots/SpatialDimPlot_", sample,"_pw_clusters_spatial_sct_spatialpca.png"),
  width = 15, height = 15
)

pw2_Spatial <- (plt_SPCAspark + plt_SPCAsvf) # +
  # plot_annotation(title = paste0(sample, " Spatial Clustering"))
pw2_Spatial
ggsave(
  plot = pw2_Spatial,
  filename = paste0("plots/SpatialDimPlot_", sample, "_pw_clusters_spcasparksvf.png"),
  width = 10, height = 7
)

pw_spca <- plt_SPCAspark + plt_SPCAcustom +
  plot_annotation(title = paste0(sample, " Spatial PCA Clustering"))
ggsave(
  plot = pw_spca,
  filename = paste0("plots/SpatialDimPlot_", sample,"_pw_clusters_spcaspark_spcacustom.png"),
  width = 12, height = 6
)

## SingleR plots ----
names(sobj@meta.data)[grepl("cluster", names(sobj@meta.data))]
names(sobj@meta.data)[grepl("cellanno", names(sobj@meta.data))]

plt_SR_cell <- SpatialPlot(sobj, group.by = "cellanno") +
    labs(title = paste0("Capture Spot SingleR")) +
    guides(fill = "none")
# plt_SR_cell
# SCTrr cluster + anno
plt_clSCTrr <- SpatialPlot(sobj, group.by = " ") +
    labs(title = paste0("SCT Rerun Clusters"))
plt_SR_clSCTrr <- SpatialPlot(sobj, group.by = "cellanno_cl_SCTrr") +
    labs(title = paste0("SCTrr Cluster SingleR"))
# plt_clSCTrr + plt_SR_clSCTrr
# SPCA SPARK cluster + anno
plt_clSPARK <- SpatialPlot(sobj, group.by = "spcaspark_clusters") +
    labs(title = paste0("SpatialPCA SPARK Clusters"))
plt_SR_clSPARK <- SpatialPlot(sobj, group.by = "cellanno_cl_spark") +
    labs(title = paste0("SpatialPCA SPARK Clusters SingleR"))
# plt_clSPARK + plt_SR_clSPARK
# SPCA custom cluster + anno
plt_clcustom <- SpatialPlot(sobj, group.by = "spcacustom_clusters") +
    labs(title = paste0("SpatialPCA SVF Clusters"))
plt_SR_clcustom <- SpatialPlot(sobj, group.by = "cellanno_cl_custom") +
    labs(title = paste0("SpatialPCA SVF Clusters SingleR"))
# plt_clcustom + plt_SR_clcustom

# Cluster composition of capture spot cell type anno
comp_SCTrr <- Rubrary::plot_comp_barplot(
  sobj = sobj,
  breaks = "SCTrr_cluster",
  xlabel = "SCTrr Cluster",
  group = "cellanno",
  stack = TRUE,
  counts = TRUE
)
comp_spark <- Rubrary::plot_comp_barplot(
  sobj = sobj,
  breaks = "spcaspark_clusters",
  xlabel = "SpatialPCA SPARK Cluster",
  group = "cellanno",
  stack = TRUE,
  counts = TRUE
)

comp_SVF <- Rubrary::plot_comp_barplot(
  sobj = sobj,
  breaks = "spcacustom_clusters",
  xlabel = "SpatialPCA SVF Cluster",
  group = "cellanno",
  stack = TRUE,
  counts = TRUE
)
# Per capture spot PW
plt_SR_cell
pw_spot <- (plt_SR_cell | comp_SCTrr) / (comp_spark | comp_SVF) +
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom')
pw_spot
ggsave(
  plot = pw_spot,
  file = paste0("plots/SRSpot_", sample, "_pw_sq.png"),
  width = 8, height = 9
)

# SingleR cluster PW
pw_SR <- (plt_SR_clSCTrr | plt_SR_clSPARK | plt_SR_clcustom) &
  theme(
    legend.position = "right",
    legend.title = element_blank()
  )
pw_SR
ggsave(
  plot = pw_SR,
  file = paste0("plots/SRCluster_", sample, "_pw_sq.png"),
  width = 18, height = 6
)
