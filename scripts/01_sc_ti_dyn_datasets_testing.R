#-----------------------Testing dynverse data sets script----------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Description: R script to testing dynverse real data sets  
# How to execute the R script in the command-line ('~/data' & '~/results' folder 
#directories need to exist - example): 
# Rscript ~/scripts/01_sc_ti_dyn_datasets_testing.R \
#   --input embronic-mesenchyme-neuron-differentiation_mca \
#   --output ~/results/embronic-mesenchyme-neuron-differentiation_mca.pdf \
#   --data ~/data --type silver
# Date: 20/02/2023
# Last update: 20/02/2023
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#
# Load packages
library("dplyr")
library("Totem")
library("ggplot2")
library("dynwrap")
library("optparse")
library("dyndimred")
library("SingleCellExperiment")

## Parse input parameters
opts = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="dynverse real data set to download and analyse", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file path PDF name to save results", metavar="character"), 
  make_option(c("-d", "--data"), type="character", default=NULL, 
              help="data folder to save the rds dynverse data set", metavar="character"), 
  make_option(c("-t", "--type"), type="character", default=NULL, 
              help="type of data set - one of 'gold' or 'silver'", metavar="character")
)
opts_parser <- OptionParser(option_list=opts)
opts_parsed <- parse_args(opts_parser)


# Download data sets
url2down <- paste0("https://zenodo.org/record/1443566/files/real/",opts_parsed$type, "/", opts_parsed$input, ".rds?download=1")
file2proc <- file.path(opts_parsed$data, paste0(opts_parsed$input, ".rds"))
download.file(url = url2down, destfile = file2proc)

# Set seed
set.seed(1204)

# Import data sets
dyn <- readRDS(file = file2proc)

## Export the results as PDF
pdf(opts_parsed$output)

## Create a SCE object for Totem

# Retrieve counts and log normalized and groupings data from 'dyn' object 
#or use your own
counts <- get_expression(dataset = dyn, expression_source = "counts") # get 'counts'
log.counts <- get_expression(dataset = dyn, expression_source = "expression") # get 'counts'
groups <- get_grouping(dataset = dyn)
cell_ids <- names(groups)

# Transpose (rows x cols): cells x genes --> genes x cells 
# dynwrap assumes: cells x genes
# SCE object assumes: genes x cells 
counts <- t(counts)
log.counts <- t(log.counts)

# Create SCE object
sce <- SingleCellExperiment(assays = list("counts" = counts, "logcounts" = log.counts), 
                            colData = data.frame("Group" = groups, row.names = cell_ids))


## Prepare data for Totem: remove non-expressed genes
sce <- PrepareTotem(object = sce)


## Selection of HVG w/ scran R package
var.genes <- scran::modelGeneVar(sce)
hvg <- scran::getTopHVGs(var.genes, n = 2000)


## Dimensionality reduction
sce <- RunDimRed(object = sce,
                 dim.red.method = "lmds",
                 dim.red.features = hvg,
                 dim.reduction.par.list = list(ndim=5))


## Select 2D dimensional reduction for visualization
dim_red <- dimred_mds(t(log.counts), ndim=2)


## Visualize 'Group' in UMAP projection
print(
  cbind(dim_red, "Cell_type"=sce$Group) %>% 
    as.data.frame(.) %>% 
    mutate_at(c("comp_1", "comp_2"), as.numeric) %>% 
    ggplot(data = ., mapping = aes(x = comp_1, y = comp_2, color = Cell_type)) + 
    geom_point() + 
    theme_bw()
)

## Clustering scRNA-seq data w/ CLARA (k-medoids)
set.seed(123)
sce <- RunClustering(object = sce, k.range = 3:20,
                     min.cluster.size = 5, N.clusterings = 10000)


## Visualization of cell connectivity
print(
  VizCellConnectivity(object = sce, viz.dim.red = dim_red)
)


## Select best clusters for MST calculation
sce <- SelectClusterings(sce, selection.method = 3,
                         selection.N.models = 5,
                         selection.stratified = FALSE,
                         prior.clustering = NULL)


## Visualize selected clusters
select.clusters <- ReturnTrajNames(sce)
print(
  VizMST(object = sce, clustering.names = select.clusters, viz.dim.red = dim_red)
)


## Smoothing MSTs selected w/ Slingshot
sce <- RunSmoothing(sce)


## Visualize smoothed MSTs
smooth.msts.names <- ReturnTrajNames(sce)
print(
  VizSmoothedTraj(object = sce,
                  traj.names = smooth.msts.names,
                  viz.dim.red = dim_red,plot.pseudotime = FALSE)
)


## Define the root of the cell trajectory
# root.cluster <- 3
# sce <- ChangeTrajRoot(object = sce, traj.name = smooth.msts.names[1], root.cluster = root.cluster)


## Visualize pseudotime
print(
  VizSmoothedTraj(object = sce,
                  traj.names = smooth.msts.names[1],
                  viz.dim.red = dim_red, plot.pseudotime = FALSE)
)
print(
  VizSmoothedTraj(object = sce,
                  traj.names = smooth.msts.names[1],
                  viz.dim.red = dim_red, plot.pseudotime = TRUE) 
)


## Compare Totem pseudotime against ground-truth
print(
  cowplot::plot_grid(
    (VizSmoothedTraj(object = sce,
                     traj.names = smooth.msts.names[1],
                     viz.dim.red = dim_red, plot.pseudotime = TRUE) + 
       ggtitle("Totem")),
    (dynplot::plot_dimred(dyn, color_cells = "pseudotime", dimred = dim_red) + 
       ggtitle("Ground-truth")), 
    ncol=2
  ) 
)
dev.off()

## R packages and versions used in these analyses
sessionInfo()
#
#------------------------------------------------------------------------------#
