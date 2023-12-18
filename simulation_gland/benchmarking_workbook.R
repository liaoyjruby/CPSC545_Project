# Benchmarking scratch script
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(magrittr)
  library(ggplot2)
  library(data.table)
  library(viridis)
  library(Seurat)
  library(SpatialPCA)
  library(patchwork)
  library(SRTsim)
  library(lisi)
  library(viridisLite)
  library(cluster)
  library(patchwork)
})



#############################################
#############################################
#############################################
#############################################
################ SRTsim Part ################

# To create the simulation dataset we used the UI as it provided more ability to
# control the spatial layout of cells (at least for a new user of that package)
# After many many tries, a dataset was exported and saved with the second two commands
shinySRT <- SRTsim_shiny()
simSRT <- Shiny2SRT(shinySRT)
glandular1 <- saveRDS(simSRT, file = '~/Documents/9th_year_2023/CPSC545/term_project/data/simulated_data_1_simsrtobj.rds')




######## SRTsim  Testing section Section
######## Learning to use this package
# str(exampleLIBD)
# ## Set a seed for reproducible simulation
# set.seed(1)

# counts_hcc <- my_object@assays[["Spatial"]]@layers[["counts"]]
# rownames(counts_hcc) <- Features(my_object)
# info_hcc <- my_object@images[["Slice1"]]@coordinates
# info_hcc <- info_hcc[,c("imagecol","imagerow","tissue")]
# colnames(counts_hcc) <- rownames(info_hcc)
# colnames(info_hcc) <- c("x","y","label")
# simSRT  <- createSRT(count_in = counts_hcc, loc_in = info_hcc)

# ## Estimate model parameters for data generation
# simSRT1 <- srtsim_fit(simSRT,sim_schem="tissue")
# ## Generate synthetic data with estimated parameters
# simSRT2 <- srtsim_count(simSRT1)
# simSRT3 <- srtsim_fit(simSRT,sim_schem="domain")
# simSRT4 <- srtsim_count(simSRT3)
# simSRT2   <- compareSRT(simSRT2)
# simSRT4   <- compareSRT(simSRT4)
# visualize_metrics(simSRT2)
# visualize_metrics(simSRT4)

# # SRTsim example
# str(exampleLIBD)
# example_count   <- exampleLIBD$count
# example_loc     <- exampleLIBD$info[,c("imagecol","imagerow","layer")]
# colnames(example_loc) <- c("x","y","label")
# e_simSRT  <- createSRT(count_in=example_count,loc_in =example_loc)
# e_simSRT1 <- srtsim_fit(e_simSRT,sim_schem="tissue")
# e_simSRT2 <- srtsim_count(e_simSRT1)
# e_simSRT2 <- compareSRT(e_simSRT2)
# visualize_metrics(e_simSRT2)




#############################################
#############################################
#############################################
#############################################
############ SpatialPCA Part ################

# Created custom glandular dataset but unable to set celltype proportions using UI
# Also unable to set shape other than a square
# Built glandular 2-D layout inside of square with all extraneous cells
# belonging to two groups: A and D
glandular1 <- readRDS('~/Documents/9th_year_2023/CPSC545/term_project/data/simulated_data_1_simsrtobj.rds')
celltypes_df <- data.table::fread(input = '~/Documents/9th_year_2023/CPSC545/term_project/data/celltype_proportions.csv',
                                  sep = ',', header = TRUE)
rownames(celltypes_df) <- celltypes_df$V1
celltypes_df$V1 <- NULL

# Subset glandular1 to glandular1_red (remove A and D groups)
mask <- glandular1@simcolData@listData[["group"]] %in% c('B', 'C', 'E', 'F')
xy_coords <- cbind(glandular1@simcolData@listData[["x"]],
                   glandular1@simcolData@listData[["y"]])
xy_coords_red <- xy_coords[mask,]
glandular1_counts_red <- glandular1@simCounts[,mask]  # The count matrix
rownames(xy_coords_red) <- colnames(glandular1_counts_red)
groups_red <- glandular1@simcolData@listData[["group"]][mask]
locations <- cbind(xy_coords_red, groups_red)

# Attempts to modify cell-proportions post-hoc
# Didn't work super well
# Create SRT object out of reduced data
# simSRT  <- createSRT(count_in=glandular1_counts_red,loc_in = locations)
# glandular1_red_domain <- srtsim_fit(simSRT, sim_schem="domain")
# glandular1_red_tissue <- srtsim_fit(simSRT, sim_schem="tissue")
# glandular1_red_shiftcells <- srtsim_cci_ref(
#           EstParam = glandular1_red_tissue@EstParam,
#           numGene = 5000,
#           location_in  = locations,
#           region_cell_map = celltypes_df,
#           sim_seed = 123,
#           fc = 3)



# Reduced dataset
glandular_sp_red <- CreateSpatialPCAObject(counts=glandular1_counts_red,
                                       location=xy_coords_red,
                                       project = "glandular1",
                                       gene.type="spatial",
                                       sparkversion="spark",
                                       numCores_spark=5,
                                       gene.number=3000,
                                       customGenelist=NULL,
                                       min.loctions = 20,
                                       min.features=20)
glandular_sp_red <- SpatialPCA_buildKernel(glandular_sp_red,
                                           kerneltype="gaussian",
                                           bandwidthtype="SJ",
                                           bandwidth.set.by.user=NULL)
glandular_sp_red <- SpatialPCA_EstimateLoading(glandular_sp_red,
                                               fast=FALSE,
                                               SpatialPCnum=20)
glandular_sp_red <- SpatialPCA_SpatialPCs(glandular_sp_red, fast=FALSE)
clusterlabel3 <- walktrap_clustering(clusternum = 3,
                                     latent_dat = glandular_sp_red@SpatialPCs,
                                     knearest = 70 )
clusterlabel4 <- walktrap_clustering(clusternum = 4,
                                    latent_dat = glandular_sp_red@SpatialPCs,
                                    knearest = 70 )
clusterlabel5 <- walktrap_clustering(clusternum = 5,
                                     latent_dat = glandular_sp_red@SpatialPCs,
                                     knearest = 70 )
cbp <- c("#9C9EDE" ,"#5CB85C" ,"#E377C2", "#4DBBD5" ,"#FED439" ,"#FF9896", "#FFDC91")



### Plot the SpatialPCA Domains
# Re-map labels
# P1
char_int_mapping <- c("C" = "A", "B" = "B", "E" = "C", "F" = "D")
new_labels <- char_int_mapping[groups_red]
p1 <- plot_cluster(location = xy_coords_red,
             clusterlabel = new_labels,    # groups_red clusterlabel
             pointsize = 2,
             text_size = 12,
             title_in = paste0("Ground truth simulated data"),
             color_in = cbp,
             legend = 'bottom') + theme(legend.position = "none")
# P2
char_int_mapping <- c('3' = "A", '1' = "B", '2' = "C")
new_labels <- char_int_mapping[clusterlabel3]
p2 <- plot_cluster(location = xy_coords_red,
             clusterlabel = new_labels,    # groups_red clusterlabel
             pointsize = 2,
             text_size = 12,
             title_in = paste0("Walktrap 3 clusters"),
             color_in = cbp,
             legend = 'bottom') + theme(legend.position = "none")
# P3
char_int_mapping <- c('3' = "A", '4' = "B", '1' = "C", "2" = 'D')
new_labels <- char_int_mapping[clusterlabel4]
p3 <- plot_cluster(location = xy_coords_red,
             clusterlabel = new_labels,    # groups_red clusterlabel
             pointsize = 2,
             text_size = 12,
             title_in = paste0("Walktrap 4 clusters"),
             color_in = cbp,
             legend = 'bottom') + theme(legend.position = "none")
# P4
char_int_mapping <- c('2' = "A", '3' = "B", '5' = "C", "1" = 'D', "4" = 'E')
new_labels <- char_int_mapping[clusterlabel5]
p4 <- plot_cluster(location = xy_coords_red,
             clusterlabel = new_labels,    # groups_red clusterlabel
             pointsize = 2,
             text_size = 12,
             title_in = paste0("Walktrap 5 clusters"),
             color_in = cbp,
             legend = 'right')

legend <- cowplot::get_legend(p4)
pp <- (p1 + p2) / (p3 + (p4 + theme(legend.position = 'none')))

# Arrange plots and legend
final_plot <- cowplot::plot_grid(
  pp,                       # Combined plots without the legend
  legend,                   # Legend as a separate element
  ncol = 2,                 # Number of columns in the grid
  rel_widths = c(6, 1)      # Adjust the width ratio between plots and legend
)

final_plot
ggsave(final_plot, filename = '~/Documents/9th_year_2023/CPSC545/term_project/figures/clusters_plot.png',
       dpi = 300, width = 22, height = 18, units = 'cm')

# Save output:
saveRDS(glandular_sp_red, file = '~/Documents/9th_year_2023/CPSC545/term_project/data/simulated_data_1_reduced_spatialpca.rds')


# Full dataset
sample_names <- c("glandular_sim")
count_sub <- glandular1@simCounts  # The count matrix
xy_coords <- cbind(glandular1@simcolData@listData[["x"]],
                      glandular1@simcolData@listData[["y"]])
rownames(xy_coords) <- glandular1@simcolData@rownames
glandular_sp <- CreateSpatialPCAObject(counts=count_sub,
                              location=xy_coords,
                              project = "glandular1",
                              gene.type="spatial",
                              sparkversion="spark",
                              numCores_spark=5,
                              gene.number=3000,
                              customGenelist=NULL,
                              min.loctions = 20,
                              min.features=20)
glandular_sp <- SpatialPCA_buildKernel(glandular_sp,
                              kerneltype="gaussian",
                              bandwidthtype="SJ",
                              bandwidth.set.by.user=NULL)
glandular_sp <- SpatialPCA_EstimateLoading(glandular_sp,
                                          fast=FALSE,
                                          SpatialPCnum=20)
glandular_sp <- SpatialPCA_SpatialPCs(glandular_sp, fast=FALSE)
clusterlabel <- walktrap_clustering(clusternum = 5,
                                    latent_dat = glandular_sp@SpatialPCs,
                                    knearest = 70 )
cbp <- c("#9C9EDE" ,"#5CB85C" ,"#E377C2", "#4DBBD5" ,"#FED439" ,"#FF9896", "#FFDC91")
plot_cluster(location = xy_coords,
             clusterlabel = glandular1@simcolData@listData[["group"]],
             pointsize = 1.5,
             text_size = 20,
             title_in = paste0("SpatialPCA"),
             color_in = cbp,
             legend = 'bottom')

saveRDS(glandular_sp, file = '~/Documents/9th_year_2023/CPSC545/term_project/data/simulated_data_1_full_spatialpca.rds')


#############################################
#############################################
#############################################
#############################################
############ Metrics Portion ################

# We will calculate cLISI, Silhouette, and CHAOS scores
glandular_sp_red <- readRDS(file = '~/Documents/9th_year_2023/CPSC545/term_project/data/simulated_data_1_reduced_spatialpca.rds')

# Calculate Silhouette Score

# We will use the Minkowski distance to calculate the distance matrix
mink_dist <- stats::dist(t(glandular_sp_red@SpatialPCs),
                         method = 'minkowski',
                         p = 3)
# Calculate silhouette scores considering spatial information
# Convert characters to integers using factors

# Reminder of the mappings we are using for cluster/domain IDs
# char_int_mapping <- c("C" = "A", "B" = "B", "E" = "C", "F" = "D")
# char_int_mapping <- c('3' = "A", '1' = "B", '2' = "C")
# char_int_mapping <- c('3' = "A", '4' = "B", '1' = "C", "2" = 'D')
# char_int_mapping <- c('2' = "A", '3' = "B", '5' = "C", "1" = 'D', "4" = 'E')

integer_vector <- as.integer(as.factor(groups_red))
sil_vals_truth <- silhouette(integer_vector, mink_dist)
sil_vals_3 <- silhouette(as.integer(clusterlabel3), mink_dist)
sil_vals_4 <- silhouette(as.integer(clusterlabel4), mink_dist)
sil_vals_5 <- silhouette(as.integer(clusterlabel5), mink_dist)
summary(sil_vals_truth)$clus.avg.widths
summary(sil_vals_3)$clus.avg.widths
summary(sil_vals_4)$clus.avg.widths
summary(sil_vals_5)$clus.avg.widths

# Calculate CHAOS Score
# Low CHAOS score indicates that clusters/domains are continuous and smooth
fx_CHAOS(groups_red, xy_coords_red)
fx_CHAOS(clusterlabel3, xy_coords_red)
fx_CHAOS(clusterlabel4, xy_coords_red)
fx_CHAOS(clusterlabel5, xy_coords_red)




# Calculate LISI
clusters <- cbind(groups_red, clusterlabel3, clusterlabel4, clusterlabel5)
colnames(clusters) <- c('ground_truth', 'walktrap_3_clusters',
                        'walktrap_4_clusters', 'walktrap_5_clusters')

as.data.frame(xy_coords_red) %>%
  cbind(clusters) %>%
  sample_frac(1L, FALSE) %>%
  gather(key, val, ground_truth, walktrap_3_clusters,
         walktrap_4_clusters, walktrap_5_clusters) %>%
  ggplot(aes(V1, V2, color = val)) +
  geom_point(shape = 19) +
  facet_wrap(~key)
lisi_score <- compute_lisi(xy_coords_red, clusters, c('ground_truth', 'walktrap_3_clusters',
                                                      'walktrap_4_clusters', 'walktrap_5_clusters'))
as.data.frame(xy_coords_red) %>%
  cbind(lisi_score) %>%
  sample_frac(1L, FALSE) %>%
  gather(key, val, ground_truth, walktrap_3_clusters,
         walktrap_4_clusters, walktrap_5_clusters) %>%
  ggplot(aes(V1, V2, color = val)) +
  geom_point(shape = 19, size = 3) +
  facet_wrap(~key) + labs(color = "LISI Score") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        strip.text = element_text(size = 12)) +
  scale_color_viridis(option = "magma")

ggsave(filename = '~/Documents/9th_year_2023/CPSC545/term_project/figures/LISI_plot.png',
       dpi = 300, width = 22, height = 18, units = 'cm')



#############################################
#############################################
#############################################
#############################################
################# Extra code ################

# fx_lisi = function(clusterlabel_in, location){
#   dat = data.frame("clusterlabel"=clusterlabel_in)
#   lisi=compute_lisi(location, dat, c("clusterlabel"))
#   return("lisi"=lisi$clusterlabel)
# }


# Example
# head(lisi::meta_data)
# lisi::X %>%
#   cbind(lisi::meta_data) %>%
#   sample_frac(1L, FALSE) %>%
#   gather(key, val, label1, label2) %>%
#   ggplot(aes(X1, X2, color = val)) +
#   geom_point(shape = 21) +
#   facet_wrap(~key)
# lisi_res <- compute_lisi(lisi::X, lisi::meta_data, c('label1', 'label2'))
# head(lisi_res)
# lisi::X %>%
#   cbind(lisi_res) %>%
#   sample_frac(1L, FALSE) %>%
#   gather(key, lisi_value, label1, label2) %>%
#   ggplot(aes(X1, X2, color = lisi_value)) +
#   geom_point(shape = 21) +
#   facet_wrap(~key)


#############################################
#############################################
#############################################
#############################################
####### Abandoned analysis directions #######

# HCC_pathC1 <- '~/repos/CPSC545_Project/data/Control1'
# HCC_pathC2 <- '~/repos/CPSC545_Project/data/Control2/'
# HCC_pathT1 <- '~/repos/CPSC545_Project/data/Tumor1/'
# HCC_pathT2 <- '~/repos/CPSC545_Project/data/Tumor2/'
#  gene.column = 2 by default and that was causing an error, update to 1 for last version with gziped
# my_object <- CreateSeuratObject(
#   counts = Read10X( data.dir = paste0(HCC_pathT1, '/filtered_feature_bc_matrix')),
#   assay = 'Spatial'
# )
# the_image <- Read10X_Image(HCC_pathT2)
# image <- the_image[Cells(x = my_object)]
# DefaultAssay(object = image) <- 'Spatial'
# my_object[['Slice1']] <- image
# plot1 <- VlnPlot(my_object, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
# plot2 <- SpatialFeaturePlot(my_object, features = "nCount_Spatial") + theme(legend.position = "right")
# wrap_plots(plot1, plot2)

# Old tonsil stuff
# tonsil_obj <- readRDS("~/repos/CPSC545_Project/data/spatial_transcriptomics/20220527_tonsil_atlas_spatial_seurat_obj.rds")
# plot1 <- VlnPlot(tonsil_obj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
# plot2 <- SpatialFeaturePlot(tonsil_obj, features = "nCount_Spatial") + theme(legend.position = "right")
# wrap_plots(plot1, plot2)

