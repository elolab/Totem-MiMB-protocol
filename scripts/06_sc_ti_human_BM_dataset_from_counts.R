#------------------Testing human bone marrow data set script-------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Description: R script to test Totem on CD34+ human bone marrow cells (replicate 1) 
#used in the Palantir MS, but starting from counts
# Date: 01/03/2023
# Last update: 01/03/2023
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
### Import data 
## Download CD34+ human bone marrow single-cell counts used in Palantir's MS available on github: 
# https://github.com/dpeerlab/Palantir/blob/master/notebooks/Palantir_sample_notebook.ipynb
data.set <- "human_bm_palantir"
data.dir <- file.path("data", data.set)
file2down <- file.path(data.dir, "cd34_human_bm_palantir_counts_rep1.h5ad")
if (!file.exists(file2down)) {
  url2down <- "https://dp-lab-data-public.s3.amazonaws.com/palantir/marrow_sample_scseq_counts.h5ad"
  if (!dir.exists(data.dir)) dir.create(data.dir, recursive = TRUE)
  download.file(url = url2down, destfile = file2down)
}
sce <- zellkonverter::readH5AD(file = file2down)

## Download CD34+ human bone marrow single-cell (replicate 1) data set used in Palantir's MS: 
# https://data.humancellatlas.org/explore/projects/091cf39b-01bc-42e5-9437-f419a66c8a45
file2down <- file.path(data.dir, "cd34_human_bm_palantir_rep1.h5ad")
if (!file.exists(file2down)) { # link from HCA expires after some minutes 
  url2down <- "https://storage.googleapis.com/datarepo-4ef3f5a3-bucket/0e5e329e-2709-4ceb-bfe3-97d23a652ac0/3428f967-3376-4051-b6f7-8dd84580ca5b/human_cd34_bm_rep1.h5ad?X-Goog-Algorithm=GOOG4-RSA-SHA256&X-Goog-Credential=datarepo-jade-api%40terra-datarepo-production.iam.gserviceaccount.com%2F20230227%2Fauto%2Fstorage%2Fgoog4_request&X-Goog-Date=20230227T122200Z&X-Goog-Expires=900&X-Goog-SignedHeaders=host&requestedBy=azul-ucsc-0-public-prod%40platform-hca-prod.iam.gserviceaccount.com&userProject=datarepo-1dbff5cd&X-Goog-Signature=07bd5186f33602b05a1cd198b6ae4b8a16a6a106f738b1a1b15d8df588078090651e9bd0b0fd8e84913e4cc550caeed18430b3e259f176eb1a35065843834ffbfe569b072edf217aefc612cab42ff32d36033d3c9b36f9c07e64018fcf76ac5a92ee95661e2e2c985e2a16d794a82f78813a9b19217f655276f7237b0542c02888e8100490d1f3eca870d7c3c76fa7c2f67e62433157b701eac806819077fac105708ec10c6cd7db4c0f0fff7728e5b1e2ca448eb92299b67b5ed4a01b19c0ba27c3ee42c723756ae469e89fc3be6cd1042f8e9081eca422899c908b3976bd0d9664e88722a96c2818aebc1de0fba768a6e9e889974e2e61acbbe02e4bcc7377"
  if (!dir.exists(data.dir)) dir.create(data.dir, recursive = TRUE)
  download.file(url = url2down, destfile = file2down)
}
sce.hca <- zellkonverter::readH5AD(file = file2down)

## Add cell annotations
# WARNING: can't match cell barcodes used in CSV file annotation with cell barcodes used in the anndata objects
# meta2down <- file.path(data.dir, "HaematopoieticProfiling-10x_cell_type_2020-03-12.csv")
# if (!file.exists(meta2down)) {
#   meta.url <- "https://storage.googleapis.com/datarepo-4ef3f5a3-bucket/0e5e329e-2709-4ceb-bfe3-97d23a652ac0/b4fcd6f8-0ec5-4608-bfe6-fa76406ace9f/HaematopoieticProfiling-10x_cell_type_2020-03-12.csv?X-Goog-Algorithm=GOOG4-RSA-SHA256&X-Goog-Credential=datarepo-jade-api%40terra-datarepo-production.iam.gserviceaccount.com%2F20230227%2Fauto%2Fstorage%2Fgoog4_request&X-Goog-Date=20230227T132451Z&X-Goog-Expires=900&X-Goog-SignedHeaders=host&requestedBy=azul-ucsc-0-public-prod%40platform-hca-prod.iam.gserviceaccount.com&userProject=datarepo-1dbff5cd&X-Goog-Signature=6ef99a32aeff49ecef62b24a57762691b5f6a27d2df9ddfa8e62c5a605b8bc338b4cfb4102ef38f9a38250ab7ed22d014c34fd03aae225e0c4544183ed23e6ed115d9a8eb675e019bdc1855d26a2f92ae6fc480800c971a32eaecfe3931aac8c21bbcb7d71c8be272344b3732278df24969683a9af47088271cf08da5fe70a2ee6c350e0a5030440438c965ec913bbfed49e9e4af051f3692f61c16b22c96f5f2274fd0809b54acf0b903620100b2f3d482ca8f67ccd97e62b817af92cb51f62cb75e94bd2d2f6d864017130fa60907b14d8d2fd2d1529948c011f29b2b3822c0b006ac46b83da95fe85c4a415328b8747937870b07a11c789319ba5a2ba539b"
#   download.file(meta.url, destfile = meta2down)
# }
# cell.annot <- read.table(file = meta2down, header = TRUE, sep = ",")
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
### Parse data
## Add cell metadata from .h5ad object made available through HCA ("sce.hca") to the 
#counts .h5ad object ("sce" - both made available by the dpeerlab) 
stopifnot(all(colnames(sce) %in% colnames(sce.hca)))
colData(sce) <- colData(sce.hca)[colnames(sce),]

## Add tSNE projection
reducedDim(sce, type="tsne") <- reducedDim(sce.hca, type="tsne")[colnames(sce),]

## Fix 'cluster_colors' names
metadata(sce)$cluster_colors <- metadata(sce.hca)$cluster_colors
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Totem's workflow
output <- file.path("results", data.set, paste0(data.set, "_counts.pdf"))
if (!dir.exists(dirname(output))) dir.create(dirname(output), recursive = TRUE)

## Set seed
set.seed(1024)

## Log-normalize
log.norm <- apply(assay(sce, "X"), 2, function(x) log1p(x/sum(x)*10000))
assay(sce, "logcounts") <- as(log.norm, "sparseMatrix")
rm(log.norm); gc(); # rm env to free memory

## Preparing data to Totem
sce <- PrepareTotem(sce)

## Selection of HVG w/ scran R package
var.genes <- scran::modelGeneVar(sce)
hvg <- scran::getTopHVGs(var.genes, n = 2000)

## Dimensional reduction: 
# PCA: clustering
## Export the results as PDF
pdf(output, width=18, height=9)
sce <- RunDimRed(object = sce,
                 dim.red.method = "pca",
                 dim.red.features = hvg,
                 dim.reduction.par.list = list(ndim=50))

## Elbow, PCA & UMAP
# Format names of the cluster colors of original clusters
names(metadata(sce)$cluster_colors) <- levels(sce$clusters) 
elbow.plt <- reducedDim(sce, "pca") %>% 
  apply(X = ., MARGIN = 2, FUN = function(x) sd(x)) %>%
  as.data.frame(.) %>% 
  `colnames<-`("Standard Deviation") %>% 
  mutate("PCs" = factor(1:nrow(.), levels = 1:nrow(.))) %>% 
  ggplot(data = ., mapping = aes(x = PCs, y = `Standard Deviation`)) + 
  geom_point() + 
  theme_bw()
pca.plt <- reducedDim(sce, "pca") %>% 
  as.data.frame(.) %>%
  mutate("Cell_ID" = row.names(.)) %>% 
  cbind(., colData(sce)) %>% 
  ggplot(data = ., mapping = aes(x=comp_1, y=comp_2, color=clusters)) + 
  geom_point() + 
  labs(x = "PC1", y = "PC2") + 
  scale_color_manual(values=unlist(metadata(sce)$cluster_colors)) + 
  theme_bw()
pca.scores.plt <- reducedDim(sce, "pca") %>% 
  as.data.frame(.) %>% 
  mutate("Cell_ID"=row.names(.)) %>% 
  cbind(., colData(sce)) %>% 
  tidyr::pivot_longer(., cols=comp_1:comp_20, names_to="PCs", values_to ="Scores") %>% 
  mutate("PCs" = factor(PCs, levels = paste0("comp_", 1:20))) %>% 
  ggplot(data= ., aes(x=PCs, y=Scores)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=clusters), size=0.25) + 
  scale_color_manual(values=unlist(metadata(sce)$cluster_colors)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))
print(
  cowplot::plot_grid(cowplot::plot_grid(elbow.plt, pca.plt, ncol=2), pca.scores.plt, ncol=1)
)

## Set the no. of PCs
pcs <- 3

## Select only the first three PCs 
#(which hold biological variance from the viz above)
reducedDim(sce, "pca") <- reducedDim(sce, "pca")[,1:pcs]

# Remove original tSNE & save it for later
orig.dimred <- reducedDim(sce, "tsne") # original tSNE - used for comparison
reducedDim(sce, "tsne") <- NULL

## Compute UMAP for visualization
umap.dimred <- dimred_umap(t(logcounts(sce[hvg,])), ndim=2, pca_components = pcs)

## Clustering scRNA-seq data w/ CLARA (k-medoids), MSTs & smoothing 
#w/ principal curves algorithm
set.seed(123)
sce <- RunClustering(object = sce) %>% 
  SelectClusterings(object = .) %>% 
  RunSmoothing(object = .)

## Get name of the best cluster result
select.cluster <- names(metadata(sce)$totem$slingshot_trajectory)

## Visualize cluster & MST
print(
  VizMST(sce, clustering.names = select.cluster, viz.dim.red = umap.dimred)
)

## Define root
root.cluster <- 2
sce <- ChangeTrajRoot(object = sce, traj.name = select.cluster, 
                      root.cluster = root.cluster)

## Visualize pseudotime
# Project original data
print(
  VizCellConnectivity(object = sce,
                      viz.dim.red = umap.dimred)
)
reducedDim(sce, "tsne") <- orig.dimred
reducedDim(sce, "umap") <- umap.dimred
print(
  cowplot::plot_grid(
    VizSmoothedTraj(object = sce,
                    traj.names = select.cluster,
                    viz.dim.red = umap.dimred,
                    plot.pseudotime = TRUE),
    plotReducedDim(sce, dimred = "umap", colour_by = "palantir_pseudotime") + 
      ggtitle("Palantir's pseudotime"),
    plotReducedDim(sce, dimred = "umap", colour_by = "clusters") + 
      scale_color_manual(values=unlist(metadata(sce)$cluster_colors)) + 
      ggtitle("Original clusters"),
    VizSmoothedTraj(object = sce,
                    traj.names = select.cluster,
                    viz.dim.red = orig.dimred,
                    plot.pseudotime = TRUE),
    plotReducedDim(sce, dimred = "tsne", colour_by = "palantir_pseudotime") + 
      ggtitle("Palantir's pseudotime"),
    plotReducedDim(sce, dimred = "tsne", colour_by = "clusters") + 
      scale_color_manual(values=unlist(metadata(sce)$cluster_colors)) + 
      ggtitle("Original clusters"),
    ncol = 3)
)
dev.off()

## Pseudotime correlation: Palantir vs Totem
cor(sce$palantir_pseudotime, dynwrap::calculate_pseudotime(metadata(sce)$totem$dynwrap_trajectory[[select.cluster]])) # 0.88
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## R packages and versions used in these analyses
sessionInfo()
#
#------------------------------------------------------------------------------#
