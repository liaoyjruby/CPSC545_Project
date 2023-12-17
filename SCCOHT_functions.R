# Libraries ----
library(parallel)
library(BiocParallel) # Parallelization
library(data.table) # Efficient tables
library(dplyr) # Wrangling
library(Seurat) # Single cell analysis
library(ggplot2) # Plotting
library(patchwork) # Plotting as grid
# BiocManager::install("xzhoulab/SPARK")
# BiocManager::install("shangll123/SpatialPCA")
library(SpatialPCA) # Spatially aware dimensional reduction
library(celldex) # https://bioconductor.org/packages/3.18/data/experiment/html/celldex.html
library(SingleR) # Single cell type annotation from reference
library(scRNAseq)
n_cores <- detectCores() - 2 # overcomplicated way of saying "6"
setDTthreads(n_cores)

# Variables ----
dir_data <- "/Users/liaoyj/Library/CloudStorage/OneDrive-UBC/CPSC545_ProjectData/SCCOHT/GSE213699_mod/"
samples <- list.files(dir_data)[grepl("^[A-D]_", list.files(dir_data))]

# Get Seurat object from path ----
# sample: string; name of sample
# type: string; c("GEO", "SPCA", "SR")
get_sobj <- function(sample = "D_GTFB1170_SmallCellOvarianCancer", type = NULL) {
  if (is.null(type)) {
    path <- paste0(dir_data, sample, "/sobj.rds")
  } else if (type == "GEO") {
    path <- paste0(dir_data, sample, "/sobj_GEO.rds")
  } else if (type == "SPCA") {
    path <- paste0(dir_data, sample, "/sobj_SPCA.rds")
  } else if (type == "SR") {
    path <- paste0(dir_data, sample, "/sobj_SPCA_SR.rds")
  }
  return(readRDS(path))
}

# Create Seurat object ----
# Create Seurat object for "raw" data loaded straight from GEO
# force = TRUE will recreate `filtered_feature_bc_matrix.h5`
make_sobj <- function(sample, force = FALSE) {
  if (!("filtered_feature_bc_matrix.h5" %in% list.files(paste0(dir_data, "/", sample))) ||
     force == TRUE) {
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

# Rerun Seurat pipeline based on Sanders et al. parameters
# https://github.com/kwells4/visium_ovarian_cancer/tree/main/src/scripts/D_GTFB1170_SmallCellOvarianCancer
# RNA_nPCs <- 33
# Number of PCs for UMAP from 02_PCA_UMAP.R files from Sanders GH repo
RNA_nPCs_smp <- c(
    "A_955_OvarianTumor" = 21,
    "A_GTFB1154_OvarianCancerTumor" = 23,
    "B_1180_Omentum" = 12,
    "B_GTFB1191_OvarianCancerTumor" = 26,
    "C_GTFB1170_SmallCellOvarianCancer" = 33,
    "C_GTFB_1191_OmentumTumor" = 19,
    "D_GTFB1170_SmallCellOvarianCancer" = 33,
    "D_GTFB_1230_OvarianTumor" = 30)

# Replicate SCTransform ----
process_sobj <- function(sample, cl_resolution = 0.6, save = FALSE) {
    message(paste0("Sample: ", sample))
    RNA_nPCs <- RNA_nPCs_smp[[sample]]
    sobj_GEO <- get_sobj(sample, type = "GEO")

    # Rerun generic Seurat pipeline for multiple assays ----
    # Resolution + seed from Sanders GH repo
    ## "Spatial" assay: raw counts ----
    assay <- "Spatial"
    DefaultAssay(sobj_GEO) <- assay
    RNA_nPCs <- RNA_nPCs_smp[[sample]]
    set.seed(0)
    sobj <- sobj_GEO %>%
        NormalizeData(normalization.method = "LogNormalize") %>%
        FindVariableFeatures() %>%
        ScaleData(vars.to.regress = NULL) %>%
        RunPCA(assay = "Spatial", reduction.name = "pca",
                reduction.key = paste0(assay, "PC_")) %>%
        FindNeighbors(reduction = "pca", dims = 1:RNA_nPCs) %>%
        FindClusters(resolution = cl_resolution) %>%
        RunUMAP(
            metric = "correlation",
            dims = 1:RNA_nPCs,
            reduction = "pca",
            reduction.name = "umap",
            reduction.key = paste0(assay, "UMAP_"),
            assay = assay
        )
    sobj[["Spatial_cluster"]] <- sobj@meta.data$seurat_clusters
    ## "SCTrr" assay: rerunning SCtransform ----
    sobj <- sobj %>% SCTransform(assay = "Spatial", new.assay.name = "SCTrr")
    assay <- "SCTrr"
    DefaultAssay(sobj) <- assay
    sobj <- sobj %>%
        RunPCA(
            assay = assay,
            reduction.name = "sctpca",
            reduction.key = paste0(assay, "PC_")) %>%
        FindNeighbors(reduction = "sctpca", dims = 1:RNA_nPCs) %>%
        FindClusters(resolution = cl_resolution) %>%
        RunUMAP(
            metric = "correlation",
            dims = 1:RNA_nPCs,
            reduction = "sctpca",
            reduction.name = "sctumap",
            reduction.key = paste0(assay, "UMAP_"),
            assay = assay
        )
    sobj[["SCTrr_cluster"]] <- sobj@meta.data$seurat_clusters

    # Write sobj to RDS ----
    if (save) {
        saveRDS(
            object = sobj,
            file = paste0(dir_data, sample, "/sobj.rds")
        )
        message(paste0("* ", sample, "/sobj.rds", " saved!"))
    }
    return(sobj)
}

# Run SpatialPCA ----
# sample: string, name of sample to process
# n_cores: integer, number of cores for parallelization where applicable
# spark: logical, process with variable genes derived from `spark`? (SLOW)
# save: logical, save objects as RDS?
process_SPCA <- function(
    sample, n_cores = parallel::detectCores() - 2,
    spark = TRUE, save = TRUE) {
    message(paste0("Sample: ", sample))
    dir.create(
        paste0(dir_data, sample, "/SpatialPCA/"),
        showWarnings = FALSE
    )
    # Load data ----
    # Use raw counts from GEO sobj
    sobj <- get_sobj(sample)
    DefaultAssay(sobj) <- "Spatial"

    cts <- sobj@assays$Spatial$counts
    loc <- sobj@meta.data %>%
        select(row, col) %>%
        as.matrix()

    ## Spark gene selection ----
    if (spark) {
    # Create SpatialPCA object
    # small sample size data = higher detection power of spatial genes
        SPCA_obj <- SpatialPCA::CreateSpatialPCAObject(
            counts = cts,
            location = loc,
            project = Project(sobj),
            gene.type = "spatial",
            sparkversion = "spark",
            numCores_spark = n_cores
        )

        # Estimate spatial PCs ----
        SPCA_obj <- SPCA_obj %>%
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

        # Extract data from SPCAobj
        # PCs
        SPCA_SPCs_spark <- SPCA_obj@SpatialPCs %>%
            t() %>%
            as.data.frame() %>%
            `rownames<-`(SPCA_obj@counts@Dimnames[[2]]) %>%
            `colnames<-`(paste0("SpPC", 1:20)) %>%
            as.matrix()
        # Raw counts
        SPCA_cts_spark <- SPCA_obj@counts
        # Norm counts (SCT)
        SPCA_ctsnorm_spark <- SPCA_obj@normalized_expr

        if (save) {
            saveRDS(
                SPCA_obj,
                file = paste0(dir_data, sample, "/SpatialPCA/SPCA_obj.rds"))
            message(paste0("* ", sample, "/SpatialPCA/SPCA_obj.rds", " saved!"))
            Rubrary::rwrite(
                tibble::rownames_to_column(as.data.frame(SPCA_SPCs_spark), "barcode"),
                paste0(dir_data, sample, "/SpatialPCA/SPCA_obj_SPCs.txt"))
            Rubrary::rwrite(
                tibble::rownames_to_column(as.data.frame(SPCA_cts_spark), "gene"),
                paste0(dir_data, sample, "/SpatialPCA/SPCA_obj_cts.txt"))
            Rubrary::rwrite(
                tibble::rownames_to_column(as.data.frame(SPCA_ctsnorm_spark), "gene"),
                paste0(dir_data, sample, "/SpatialPCA/SPCA_obj_cts_normexpr.txt"))
        }
    }

    # Using Seurat VariableFeatures as gene list
    SPCA_obj_genelist <- SpatialPCA::CreateSpatialPCAObject(
        counts = cts,
        location = loc,
        project = Project(sobj),
        gene.type = "custom",
        customGenelist = VariableFeatures(sobj),
        numCores_spark = n_cores
    )
    # Get spatial PCs
    SPCA_obj_genelist <- SPCA_obj_genelist %>%
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

    # Extract data from SPCAobj
    # PCs
    SPCA_SPCs_custom <- SPCA_obj_genelist@SpatialPCs %>%
        t() %>%
        as.data.frame() %>%
        `rownames<-`(SPCA_obj_genelist@counts@Dimnames[[2]]) %>%
        `colnames<-`(paste0("SpPC", 1:20)) %>%
        as.matrix()
    # Counts (should be same as sobj)
    SPCA_cts_custom <- SPCA_obj_genelist@counts
    # Norm counts (SCT)
    SPCA_ctsnorm_custom <- SPCA_obj_genelist@normalized_expr

    if (save) {
        saveRDS(
            SPCA_obj_genelist,
            file = paste0(dir_data, sample, "/SpatialPCA/SPCA_obj_genelist.rds"))
        message(paste0("* ", sample, "/SpatialPCA/SPCA_obj_genelist.rds", " saved!"))
        Rubrary::rwrite(
            tibble::rownames_to_column(as.data.frame(SPCA_SPCs_custom), "barcode"),
            paste0(dir_data, sample, "/SpatialPCA/SPCA_obj_genelist_SPCs.txt"))
        Rubrary::rwrite(
            tibble::rownames_to_column(as.data.frame(SPCA_cts_custom), "gene"),
            paste0(dir_data, sample, "/SpatialPCA/SPCA_obj_genelist_cts.txt"))
        Rubrary::rwrite(
            tibble::rownames_to_column(as.data.frame(SPCA_ctsnorm_custom), "gene"),
            paste0(dir_data, sample, "/SpatialPCA/SPCA_obj_genelist_cts_normexpr.txt"))
    }

    # Incorporate SPCA into sobj ----
    if (spark) {
        sobj@assays[["spcaspark"]] <- CreateAssayObject(
            # counts = SPCA_cts_spark,
            data = as.matrix(SPCA_ctsnorm_spark),
            key = "spcaspark_"
        )
        sobj[["spcasparkpca"]] <- CreateDimReducObject(
            embeddings = SPCA_SPCs_spark,
            key = "spcasparkpca_",
            assay = "spcaspark"
        )
    }
    sobj@assays[["spcacustom"]] <- CreateAssayObject(
        # counts = SPCA_cts_custom,
        data = as.matrix(SPCA_ctsnorm_custom),
        key = "spcacustom_"
    )
    sobj[["spcacustompca"]] <- CreateDimReducObject(
        embeddings = SPCA_SPCs_custom,
        key = "spcacustompca_",
        assay = "spcacustom"
    )

    # Louvain (re)clustering ----
    if (spark) {
        # spcaspark
        DefaultAssay(sobj) <- "spcaspark"
        sobj <- sobj %>%
            FindNeighbors(
                reduction = "spcasparkpca",
                assay = "spcaspark",
                dims = 1:20) %>%
            FindClusters(resolution = 0.8) %>%
            RunUMAP(
                metric = "correlation",
                dims = 1:20,
                reduction = "spcasparkpca",
                reduction.name = "spcasparkumap",
                reduction.key = "spcasparkUMAP_",
                assay = "spcaspark"
            )
        sobj[["spcaspark_clusters"]] <- sobj@meta.data$seurat_clusters
    }

    # spcacustom
    DefaultAssay(sobj) <- "spcacustom"
    sobj <- sobj %>%
        FindNeighbors(
            reduction = "spcacustompca",
            assay = "spcacustom",
            dims = 1:20) %>%
        FindClusters(resolution = 0.8) %>%
        RunUMAP(
            metric = "correlation",
            dims = 1:20,
            reduction = "spcacustompca",
            reduction.name = "spcacustomumap",
            reduction.key = "spcacustomUMAP_",
            assay = "spcacustom"
        )
    sobj[["spcacustom_clusters"]] <- sobj@meta.data$seurat_clusters

    # Write sobj to RDS ----
    if (save) {
        saveRDS(
            object = sobj,
            file = paste0(dir_data, sample, "/sobj_SPCA.rds")
        )
        message(paste0("* ", sample, "/sobj_SPCA.rds", " saved!"))
    }
    return(sobj)
}

# Run SingleR ----
# https://bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html
# https://bioconductor.org/books/devel/SingleRBook/
# should be ran after SPCA analysis is performed
process_SingleR <- function(sample) {
    dir.create(
        paste0(dir_data, sample, "/SingleR/"),
        showWarnings = FALSE
    )
    # Human Primary Cell Atlas
    hpca.se <- celldex::HumanPrimaryCellAtlasData()
    # Get numeric log-expression normalized df
    sobj <- get_sobj(sample, type = "SPCA")
    sobj[["barcode"]] <- rownames(sobj@meta.data)
    rns <- rownames(sobj)
    sce <- Seurat::as.SingleCellExperiment(sobj, assay = "SCTrr")
    # Per-cell annotation
    cellanno <- SingleR::SingleR(
        test = sce,
        ref = hpca.se,
        labels = hpca.se$label.main,
        BPPARAM = MulticoreParam(workers = n_cores)
    )
    saveRDS(
        cellanno,
        file = paste0(dir_data, sample, "/SingleR/cellanno_SR_percell.rds")
    )
    sobj[["cellanno"]] <- cellanno$labels
    sobj[["cellanno_pruned"]] <- cellanno$pruned.labels
    # Cluster level annotation ----
    # https://bioconductor.org/books/devel/SingleRBook/advanced-options.html#cluster-level-annotation
    # SCTrr cluster
    cellanno_cl_SCTrr <- SingleR::SingleR(
        test = sce,
        ref = hpca.se,
        clusters = SingleCellExperiment::colData(sce)$SCTrr_cluster,
        labels = hpca.se$label.main,
    )
    saveRDS(
        cellanno_cl_SCTrr,
        file = paste0(dir_data, sample, "/SingleR/cellanno_SR_perSCTrrcluster.rds")
    )
    cellanno_cl_SCTrr_df <- cellanno_cl_SCTrr %>%
        as.data.frame() %>%
        select(labels, delta.next, pruned.labels) %>%
        `colnames<-`(c("cellanno_cl_SCTrr", "cellanno_cl_SCTrr_deltanext", "cellanno_cl_SCTrr_pruned")) %>%
        tibble::rownames_to_column("SCTrr_cluster")
    # Add cluster annotation to metadata
    sobj@meta.data <- sobj@meta.data %>%
        left_join(cellanno_cl_SCTrr_df, by = "SCTrr_cluster")

    # spcaspark cluster
    message("* Annotating SPCA SPARK clusters...")
    if ("spcaspark" %in% Seurat::Assays(sobj)) {
        sce_spcaspark <- Seurat::as.SingleCellExperiment(sobj, assay = "spcaspark")
        cellanno_cl_spcaspark <- SingleR::SingleR(
            test = sce_spcaspark,
            ref = hpca.se,
            clusters = SingleCellExperiment::colData(sce_spcaspark)$spcaspark_clusters,
            labels = hpca.se$label.main,
        )
        saveRDS(
            cellanno_cl_spcaspark,
            file = paste0(dir_data, sample, "/SingleR/cellanno_SR_perSPCAsparkcluster.rds")
        )
        cellanno_cl_spcaspark_df <- cellanno_cl_spcaspark %>%
            as.data.frame() %>%
            select(labels, delta.next, pruned.labels) %>%
            `colnames<-`(c("cellanno_cl_spark", "cellanno_cl_spark_deltanext", "cellanno_cl_spark_pruned")) %>%
            tibble::rownames_to_column("spcaspark_clusters")
        # Add cluster annotation to metadata
        sobj@meta.data <- sobj@meta.data %>%
            left_join(cellanno_cl_spcaspark_df, by = "spcaspark_clusters")
    }

    # spcacustom clusters
    message("* Annotating SPCA custom clusters...")
    sce_spcacustom <- Seurat::as.SingleCellExperiment(sobj, assay = "spcacustom")
    cellanno_cl_spcacustom <- SingleR::SingleR(
        test = sce_spcacustom,
        ref = hpca.se,
        clusters = SingleCellExperiment::colData(sce_spcacustom)$spcacustom_clusters,
        labels = hpca.se$label.main,
    )
    saveRDS(
        cellanno_cl_spcacustom,
        file = paste0(dir_data, sample, "/SingleR/cellanno_SR_perSPCAcustomcluster.rds")
    )
    cellanno_cl_spcacustom_df <- cellanno_cl_spcacustom %>%
        as.data.frame() %>%
        select(labels, delta.next, pruned.labels) %>%
        `colnames<-`(c("cellanno_cl_custom", "cellanno_cl_custom_deltanext", "cellanno_cl_custom_pruned")) %>%
        tibble::rownames_to_column("spcacustom_clusters")
    # Add cluster annotation to metadata
    sobj@meta.data <- sobj@meta.data %>%
        left_join(cellanno_cl_spcacustom_df, by = "spcacustom_clusters")

    rownames(sobj@meta.data) <- sobj$barcode
    saveRDS(
        object = sobj,
        file = paste0(dir_data, sample, "/sobj_SPCA_SR.rds")
    )
    message(paste0("* ", sample, "/sobj_SPCA_SR.rds", " saved!"))
    return(sobj)
}

# Find markers ----
process_markers <- function(sample) {
    # BiocManager::install('immunogenomics/presto')
    dir.create(
        paste0(dir_data, sample, "/Markers/"), showWarnings = FALSE)
    sobj <- get_sobj(sample, type = "SR")
    # Always use RNA/Spatial assay for DE
    DefaultAssay(sobj) <- "Spatial"
    names(sobj@meta.data)
    # SCTrr cluster markers
    message("* SCTrr cluster markers")
    Idents(sobj) <- "SCTrr_cluster"
    mrkrs_SCTrr <- FindAllMarkers(sobj, only.pos = TRUE)
    Rubrary::rwrite(
        mrkrs_SCTrr,
        paste0(dir_data, sample, "/Markers/clSCTrr_mrkrs.tsv")
    )
    # SPARK cluster markers
    message("* SPARK cluster markers")
    Idents(sobj) <- "spcaspark_clusters"
    mrkrs_SPARK <- FindAllMarkers(sobj, only.pos = TRUE)
    Rubrary::rwrite(
        mrkrs_SPARK,
        paste0(dir_data, sample, "/Markers/clSPARK_mrkrs.tsv")
    )
    # SVF cluster markers
    message("* SVF cluster markers")
    Idents(sobj) <- "spcacustom_clusters"
    mrkrs_SVF <- FindAllMarkers(sobj, only.pos = TRUE)
    Rubrary::rwrite(
        mrkrs_SVF,
        paste0(dir_data, sample, "/Markers/clSVF_mrkrs.tsv")
    )
    return(list(
        SCTrr = mrkrs_SCTrr,
        SPARK = mrkrs_SPARK,
        SVF = mrkrs_SVF
    ))
}

# Visualization ----
# function to spit out nCount_Spatial + nFeature_Spatial spatial tissue plots
plot_SpatialPlot <- function(
    sobj, feats, assay = "SCT", savename = NULL, width = 12, height = 6) {
  DefaultAssay(sobj) <- assay
  plt <- lapply(feats, function(feats) {
    SpatialPlot(
      sobj, features = feats, slot = "counts") +
      labs(title = feats) +
      theme(legend.position = "right")}) %>%
    wrap_plots() + plot_annotation(
      title = Project(sobj),
      subtitle = paste0("Assay: ", DefaultAssay(sobj)))
  if (!is.null(savename)) {
    ggsave(
      plot = plt,
      filename = savename,
      width = width, height = height
    )
  }
  return(plt)
}