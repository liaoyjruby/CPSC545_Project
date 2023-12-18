
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

# Genes of interest expression ---
Reductions(sobj)
names(sobj@meta.data)
gois <- c("SMARCA4", "SMARCB1", "BRCA2", "EZH2")
# SCTrr clusters expression
goi_exp_SCTrr <- lapply(
    gois,
    plot_exp,
    sobj = sobj,
    reduction = "sctumap",
    ident = "SCTrr_cluster"
)
pw_goi_exp_SCTrr <- wrap_plots(
    goi_exp_SCTrr,
    ncol = 2
)
pw_goi_exp_SCTrr
ggsave(
    plot = pw_goi_exp_SCTrr,
    filename = paste0("plots/GOI_exp_", sample, "_clSCTrr", ".png"),
    width = 17, height = 8
)
# SPCAspark clusters expression
goi_exp_spark <- lapply(
    gois,
    plot_exp,
    sobj = sobj,
    reduction = "spcasparkumap",
    ident = "spcaspark_clusters"
)
pw_goi_exp_spark <- wrap_plots(
    goi_exp_spark,
    ncol = 2
)
pw_goi_exp_spark
ggsave(
    plot = pw_goi_exp_spark,
    filename = paste0("plots/GOI_exp_", sample, "_clSPARK", ".png"),
    width = 17, height = 8
)

# SPCA SVF clusters expression
goi_exp_svf <- lapply(
    gois,
    plot_exp,
    sobj = sobj,
    reduction = "spcacustomumap",
    ident = "spcacustom_clusters"
)
pw_goi_exp_svf <- wrap_plots(
    goi_exp_svf,
    ncol = 2
)
pw_goi_exp_svf
ggsave(
    plot = pw_goi_exp_svf,
    filename = paste0("plots/GOI_exp_", sample, "_clSVF", ".png"),
    width = 17, height = 8
)

# Spatial plot expression
SpatialPlot(
    sobj,
    features = gois
) + viridis::scale_fill_viridis(option = "magma")

plot_exp_spatial(gois, sobj)

# Find per-cluster markers ----
mrkrs_lst <- process_markers(sample)
mrkrs_lst_pos <- process_markers(sample, pos = TRUE)
names(mrkrs_lst)
