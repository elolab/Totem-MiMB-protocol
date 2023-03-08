# Load packages
library("dplyr")
library("Totem")
library("scater")
library("ggplot2")
library("SingleCellExperiment")

# Set seed
set.seed(1204)

## Import the scRNA-seq SCE object
sce <- readRDS(file = "../data/human_cd34_bm_rep1.rds")

## ## Log normalization
## counts <- assay(altExp(sce, "raw"), "X")
## log.counts <- apply(counts, 2, function(x) log1p(x/sum(x)*10000)) # log1p normalization w/ 10K scaling factor

## Remove tSNE from SCE object
reducedDim(sce, "tsne") <- NULL

## Prepare data for Totem: remove non-expressed genes
sce <- PrepareTotem(object = sce)

## Selection of HVG w/ scran R package
var.genes <- scran::modelGeneVar(sce)
hvg <- scran::getTopHVGs(var.genes, n = 2000)

## Dimensionality reduction
sce <- RunDimRed(object = sce,
                 dim.red.method = "pca",
                 dim.red.features = hvg,
                 dim.reduction.par.list = list(ndim=50))

## ## Add low-dimensional representation to SCE object
## own_dim_red <- reducedDim(sce) # substitute this line by importing your own dimensional result (define class as matrix (rows x cols: cells x latent dimensions))
## reducedDim(sce, type = "pca") <- own_dim_red # type can be 'pca', 'umap' whatever you want - this is the name given to the dimensional result

## Inspect PCA variance

# Elbow plot
elbow.plt <- reducedDim(sce, "pca") %>% 
  apply(X = ., MARGIN = 2, FUN = function(x) sd(x)) %>%
  as.data.frame(.) %>% 
  `colnames<-`("Standard Deviation") %>% 
  mutate("PCs" = factor(1:nrow(.), levels = 1:nrow(.))) %>% 
  ggplot(data = ., mapping = aes(x = PCs, y = `Standard Deviation`)) + 
  geom_point() + 
  theme_bw()

# PCA plot
pca.plt <- reducedDim(sce, "pca") %>% 
  as.data.frame(.) %>%
  mutate("Cell_ID" = row.names(.)) %>% 
  cbind(., colData(sce)) %>% 
  ggplot(data = ., mapping = aes(x=comp_1, y=comp_2, color=cell_types_short)) + 
  geom_point() + 
  labs(x = "PC1", y = "PC2") + 
  scale_color_manual(values=as.character(unlist(metadata(sce)$cluster_colors)[!duplicated(levels(sce$cell_types_short))])) + 
  theme_bw()

# Box plot: top 20 PCs w/ cell types highlight by PC
pca.scores.plt <- reducedDim(sce, "pca") %>% 
  as.data.frame(.) %>% 
  mutate("Cell_ID"=row.names(.)) %>% 
  cbind(., colData(sce)) %>% 
  tidyr::pivot_longer(., cols=comp_1:comp_20, names_to="PCs", values_to ="Scores") %>% 
  mutate("PCs" = factor(PCs, levels = paste0("comp_", 1:20))) %>% 
  ggplot(data= ., aes(x=PCs, y=Scores)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=cell_types_short), size=0.25) + 
  scale_color_manual(values=as.character(unlist(metadata(sce)$cluster_colors)[!duplicated(levels(sce$cell_types_short))])) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))

# Plot altogether 
cowplot::plot_grid(cowplot::plot_grid(elbow.plt, pca.plt, ncol=2), pca.scores.plt, ncol=1)

## Pick PCs
pick.pcs <- 1:6
reducedDim(sce, "pca") <- reducedDim(sce, "pca")[, pick.pcs ] 

## UMAP dimensional reduction for visualization
set.seed(123)
sce <- RunDimRed(object = sce, 
                 dim.red.method = "umap", 
                 dim.red.features = hvg, 
                 dim.reduction.par.list = list(ndim=2, pca_components = 6))
dim_red <- reducedDim(sce, "umap")

## Visualize 'Group' in UMAP projection
plotReducedDim(object = sce, dimred = "umap", colour_by = "cell_types_short", point_size=0.5) + 
  scale_color_manual(name= "Cell Types", 
                     values=as.character(unlist(metadata(sce)$cluster_colors)[!duplicated(levels(sce$cell_types_short))])) + 
  theme_void()
reducedDim(sce, "uma") <- NULL # remove UMAP from 'sce' object

## Clustering scRNA-seq data w/ CLARA (k-medoids)
set.seed(123)
sce <- RunClustering(object = sce, k.range = 3:20,
                     min.cluster.size = 5, N.clusterings = 10000)

## Visualization of cell connectivity
VizCellConnectivity(object = sce, viz.dim.red = dim_red)

## Select best clusters for MST calculation
sce <- SelectClusterings(sce, selection.method = 3,
                         selection.N.models = 6,
                         selection.stratified = FALSE,
                         prior.clustering = NULL)

## Visualize selected clusters
select.clusters <- ReturnTrajNames(sce)
VizMST(object = sce, clustering.names = select.clusters, viz.dim.red = dim_red)

## Smoothing MSTs selected w/ Slingshot
sce <- RunSmoothing(sce)

## Visualize smoothed MSTs
smooth.msts.names <- ReturnTrajNames(sce)
VizSmoothedTraj(object = sce,
                traj.names = smooth.msts.names,
                viz.dim.red = dim_red,plot.pseudotime = FALSE)

## Define the root of the cell trajectory
select.traj <- "8.135"
root.cluster <- 2
sce <- ChangeTrajRoot(object = sce, traj.name = select.traj, root.cluster = root.cluster)

cowplot::plot_grid(
  VizSmoothedTraj(object = sce,
                  traj.names = select.traj,
                  viz.dim.red = dim_red, plot.pseudotime = FALSE),
  VizSmoothedTraj(object = sce,
                  traj.names = select.traj,
                  viz.dim.red = dim_red, plot.pseudotime = TRUE), 
  ncol=2
)

## Compare Totem pseudotime against ground-truth
reducedDim(sce, "umap") <- dim_red
cowplot::plot_grid(
  (VizSmoothedTraj(object = sce,
                traj.names = select.traj,
                viz.dim.red = dim_red, plot.pseudotime = TRUE) + 
     ggtitle("Totem")),
  (plotReducedDim(sce, dimred = "umap", colour_by = "palantir_pseudotime") + 
      ggtitle("Palantir's pseudotime")) + 
    theme_void() + 
    theme(legend.position = "bottom"), 
   ncol=2
   ) 

## R packages and versions used in these analyses
sessionInfo()
