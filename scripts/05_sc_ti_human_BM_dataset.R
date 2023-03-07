#------------------Testing human bone marrow data set script-------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Description: R script to test Totem on CD34+ human bone marrow cells (replicate 1) 
#used in the Palantir MS 
# Date: 27/02/2023
# Last update: 27/02/2023
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
# Load packages
library("dplyr")
library("Totem")
library("ggplot2")
library("dynwrap")
library("dyndimred")
if (!"zellkonverter" %in% installed.packages()) BiocManager::install("zellkonverter")
if (!"scater" %in% installed.packages()) BiocManager::install("scater")
library("scater")
#library("zellkonverter")
library("SingleCellExperiment")
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Download CD34+ human bone marrow single-cell (replicate 1) data set used in Palantir's MS: 
# https://data.humancellatlas.org/explore/projects/091cf39b-01bc-42e5-9437-f419a66c8a45
#links below expire after some minutes
data.set <- "human_bm_palantir"
data.dir <- file.path("data", data.set)
file2down <- file.path(data.dir, "cd34_human_bm_palantir_rep1.h5ad")
if (!file.exists(file2down)) {
  url2down <- "https://storage.googleapis.com/datarepo-4ef3f5a3-bucket/0e5e329e-2709-4ceb-bfe3-97d23a652ac0/3428f967-3376-4051-b6f7-8dd84580ca5b/human_cd34_bm_rep1.h5ad?X-Goog-Algorithm=GOOG4-RSA-SHA256&X-Goog-Credential=datarepo-jade-api%40terra-datarepo-production.iam.gserviceaccount.com%2F20230227%2Fauto%2Fstorage%2Fgoog4_request&X-Goog-Date=20230227T122200Z&X-Goog-Expires=900&X-Goog-SignedHeaders=host&requestedBy=azul-ucsc-0-public-prod%40platform-hca-prod.iam.gserviceaccount.com&userProject=datarepo-1dbff5cd&X-Goog-Signature=07bd5186f33602b05a1cd198b6ae4b8a16a6a106f738b1a1b15d8df588078090651e9bd0b0fd8e84913e4cc550caeed18430b3e259f176eb1a35065843834ffbfe569b072edf217aefc612cab42ff32d36033d3c9b36f9c07e64018fcf76ac5a92ee95661e2e2c985e2a16d794a82f78813a9b19217f655276f7237b0542c02888e8100490d1f3eca870d7c3c76fa7c2f67e62433157b701eac806819077fac105708ec10c6cd7db4c0f0fff7728e5b1e2ca448eb92299b67b5ed4a01b19c0ba27c3ee42c723756ae469e89fc3be6cd1042f8e9081eca422899c908b3976bd0d9664e88722a96c2818aebc1de0fba768a6e9e889974e2e61acbbe02e4bcc7377"
  if (!dir.exists(data.dir)) dir.create(data.dir, recursive = TRUE)
  download.file(url = url2down, destfile = file2down)
}
sce <- zellkonverter::readH5AD(file = file2down)

## Add cell annotations
meta2down <- file.path(data.dir, "HaematopoieticProfiling-10x_cell_type_2020-03-12.csv")
if (!file.exists(meta2down)) {
  meta.url <- "https://storage.googleapis.com/datarepo-4ef3f5a3-bucket/0e5e329e-2709-4ceb-bfe3-97d23a652ac0/b4fcd6f8-0ec5-4608-bfe6-fa76406ace9f/HaematopoieticProfiling-10x_cell_type_2020-03-12.csv?X-Goog-Algorithm=GOOG4-RSA-SHA256&X-Goog-Credential=datarepo-jade-api%40terra-datarepo-production.iam.gserviceaccount.com%2F20230227%2Fauto%2Fstorage%2Fgoog4_request&X-Goog-Date=20230227T132451Z&X-Goog-Expires=900&X-Goog-SignedHeaders=host&requestedBy=azul-ucsc-0-public-prod%40platform-hca-prod.iam.gserviceaccount.com&userProject=datarepo-1dbff5cd&X-Goog-Signature=6ef99a32aeff49ecef62b24a57762691b5f6a27d2df9ddfa8e62c5a605b8bc338b4cfb4102ef38f9a38250ab7ed22d014c34fd03aae225e0c4544183ed23e6ed115d9a8eb675e019bdc1855d26a2f92ae6fc480800c971a32eaecfe3931aac8c21bbcb7d71c8be272344b3732278df24969683a9af47088271cf08da5fe70a2ee6c350e0a5030440438c965ec913bbfed49e9e4af051f3692f61c16b22c96f5f2274fd0809b54acf0b903620100b2f3d482ca8f67ccd97e62b817af92cb51f62cb75e94bd2d2f6d864017130fa60907b14d8d2fd2d1529948c011f29b2b3822c0b006ac46b83da95fe85c4a415328b8747937870b07a11c789319ba5a2ba539b"
  download.file(meta.url, destfile = meta2down)
}
cell.annot <- read.table(file = meta2down, header = TRUE, sep = ",")
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Prepare SCE object for Totem
# Remove palantir related reduced dimensional results - let only 'tsne'
reducedDim(sce, type = "MAGIC_imputed_data") <- NULL
reducedDim(sce, type = "palantir_branch_probs") <- NULL
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Totem's workflow
output <- file.path("results", data.set, paste0(data.set, ".pdf"))
if (!dir.exists(dirname(output))) dir.create(dirname(output), recursive = TRUE)

## Clustering scRNA-seq data w/ CLARA (k-medoids), MSTs & smoothing 
#w/ principal curves algorithm
set.seed(123)
sce <- RunClustering(object = sce) %>% 
  SelectClusterings(object = .) %>% 
  RunSmoothing(object = .)

## Define root
root.cluster <- 15
select.cluster <- names(metadata(sce)$totem$slingshot_trajectory)
sce <- ChangeTrajRoot(object = sce, traj.name = select.cluster, 
                      root.cluster = root.cluster)

## Parse data
# Format names of the cluster colors of original clusters
names(metadata(sce)$cluster_colors) <- levels(sce$clusters) 
assay(sce, "logcounts") <- assay(sce, "X") # create a fake 'logcounts' assay just for plotting

## Visualize pseudotime
## Export the results as PDF
pdf(output, width = 18, height=5)
# Project original data
print(
  VizCellConnectivity(object = sce,
                      dim.red.type = "tsne",
                      viz.dim.red = reducedDim(sce))
)
print(
  cowplot::plot_grid(
    VizSmoothedTraj(object = sce,
                    traj.names = select.cluster,
                    viz.dim.red = reducedDim(sce),
                    plot.pseudotime = TRUE), 
    plotReducedDim(sce, dimred = "tsne", colour_by = "palantir_pseudotime") + 
      ggtitle("Palantir's pseudotime"),
    plotReducedDim(sce, dimred = "tsne", colour_by = "clusters") + 
      scale_color_manual(values=unlist(metadata(sce)$cluster_colors)) + 
      ggtitle("Original clusters"),
    ncol = 3)
)
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## R packages and versions used in these analyses
sessionInfo()
#
#------------------------------------------------------------------------------#
