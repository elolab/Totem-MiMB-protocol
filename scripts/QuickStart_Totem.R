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

## Prepare data for Totem: remove non-expressed genes
sce <- PrepareTotem(object = sce)

## Totem's TI workflow: 
# (1) clustering dimensional reduction w/ CLARA (k-medoids), MSTs 
# (2) smoothing best clustering/MSTs w/ principal curves algorithm
set.seed(123) # keep reproducibility
sce <- RunClustering(object = sce) %>% 
  SelectClusterings(object = .) %>% 
  RunSmoothing(object = .)

## ## Totem's TI workflow:
## # (1) clustering dimensional reduction w/ CLARA (k-medoids), MSTs
## # (2) smoothing best clustering/MSTs w/ principal curves algorithm
## set.seed(123) # keep reproducibility
## sce <- RunClustering(object = sce)
## sce <- SelectClusterings(object = sce)
## sce <- RunSmoothing(object = sce)

## Visualization of cell connectivity
dim_red <- reducedDim(sce, "tsne") # retrieve tSNE
VizCellConnectivity(object = sce, viz.dim.red = dim_red) # plot

## Visualization of clustering/MST
select.cluster <- names(metadata(sce)$totem$slingshot_trajectory) # retrieve the name of the best cluster
cowplot::plot_grid(
  plotReducedDim(sce, dimred = "tsne", color_by = "clusters") + 
    scale_color_manual(name= "Clusters", values=unlist(metadata(sce)$cluster_colors)) + 
    ggtitle("Palantir's clusters") + 
    theme_void() + 
    theme(legend.position = "bottom"),
  plotReducedDim(sce, dimred = "tsne", color_by = "cell_types_short") +
    scale_color_manual(name= "Cell Types", 
                       values=as.character(unlist(metadata(sce)$cluster_colors)[!duplicated(levels(sce$cell_types_short))])) +
    ggtitle("Cell types") + 
    theme_void() + 
    theme(legend.position = "bottom"),
  VizMST(object = sce, clustering.names = select.cluster, viz.dim.red = dim_red), 
  ncol=3, align = "v"
)

## Define the root of the cell trajectory
root.cluster <- 15
sce <- ChangeTrajRoot(object = sce, traj.name = select.cluster, root.cluster = root.cluster)

cowplot::plot_grid(
  VizSmoothedTraj(object = sce,
                  traj.names = select.cluster,
                  viz.dim.red = dim_red, plot.pseudotime = FALSE),
  VizSmoothedTraj(object = sce,
                  traj.names = select.cluster,
                  viz.dim.red = dim_red, plot.pseudotime = TRUE), 
  ncol=2
)

## Compare Totem pseudotime against ground-truth
cowplot::plot_grid(
  (VizSmoothedTraj(object = sce,
                traj.names = select.cluster,
                viz.dim.red = dim_red, plot.pseudotime = TRUE) + 
     ggtitle("Totem")),
  (plotReducedDim(sce, dimred = "tsne", colour_by = "palantir_pseudotime") + 
      ggtitle("Palantir's pseudotime")) + 
    theme_void() + 
    theme(legend.position = "bottom"), 
   ncol=2
   ) 

## R packages and versions used in these analyses
sessionInfo()
