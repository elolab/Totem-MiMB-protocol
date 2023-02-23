#---------------------------Optimizing params script---------------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Description: R script to optimize Totem HVG & dim parameters for dyn real data 
#set 'mesoderm-development_loh'
# Date: 22/02/2023
# Last update: 22/02/2023
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#
# Load packages
library("dplyr")
library("Totem")
library("ggplot2")
library("dynwrap")
library("dyndimred")
library("SingleCellExperiment")
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
# Params to optimize
input <- "data/mesoderm-development_loh.rds"
n.hvg <- c(100, 1000, 2000, 4000)
n.dims <- c(3, 5, 9)
params.comp <- expand.grid(n.hvg, n.dims)
colnames(params.comp) <- c("hvg", "dims")
output.dir <- "results/02-mesoderm-development_loh-optimization"
if (!dir.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Upstream Workflow
# Import data sets
dyn <- readRDS(file = input)
# Create a SCE object for Totem
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
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Downstream workflow
RunTotemWkf <- function(input, n.hvg, n.dims, output) {
  # Comp
  comp.name <- paste(n.hvg, n.dims, sep="_")
  outs[[comp.name]] <- outs <- list()
  
  # Set seed
  set.seed(1204)
  
  ## Export the results as PDF
  pdf(output)

  ## Selection of HVG w/ scran R package
  var.genes <- scran::modelGeneVar(sce)
  hvg <- scran::getTopHVGs(var.genes, n = n.hvg)
  
  
  ## Dimensionality reduction
  sce <- RunDimRed(object = sce,
                   dim.red.method = "lmds",
                   dim.red.features = hvg,
                   dim.reduction.par.list = list(ndim = n.dims))
  
  
  ## Select 2D dimensional reduction for visualization
  dim_red <- dimred_mds(t(log.counts[hvg,]), ndim=2)
  
  
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
  outs[[comp.name]][["plot"]] <- cowplot::plot_grid(
    (VizSmoothedTraj(object = sce,
                     traj.names = smooth.msts.names[1],
                     viz.dim.red = dim_red, plot.pseudotime = TRUE) + 
       ggtitle("Totem")),
    (dynplot::plot_dimred(dyn, color_cells = "pseudotime", dimred = dim_red) + 
       ggtitle("Ground-truth")), 
    ncol=2
  ) 
  print(
    outs[[comp.name]][["plot"]] 
  )
  totem.pseudo <- dynwrap::calculate_pseudotime(metadata(sce)$totem$dynwrap_trajectory[[smooth.msts.names[1]]])
  true.pseudo <- dynwrap::calculate_pseudotime(dyn)
  outs[[comp.name]][["pseudo.corr"]] <- cor(x=true.pseudo, y=totem.pseudo, method="spearman")
  plot.new()
  text(0.5, 0.5, paste0("Pseudotime's Spearman correlation: ", outs[[comp.name]][["pseudo.corr"]]))
  dev.off()
  rm(sce); gc(); 
  return(outs)
}
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Test params
comp <- list()
for (i in 1:nrow(params.comp)) {
  i.hvg <- params.comp[i, "hvg"]
  i.dim <- params.comp[i, "dims"]
  output <- file.path(output.dir, paste0("mesoderm-development_loh-", i.hvg, "_hvg-",  i.dim, "_dims.pdf"))
  comp[[i]] <- RunTotemWkf(input = input, n.hvg = i.hvg, n.dims = i.dim, output = output)
}
spear.corr <- lapply(comp, function(x) x[[1]]$pseudo.corr) %>% unlist() %>% round(., 2)
label.plts <- lapply(comp, names) %>% unlist %>% 
  gsub("_", " hvg & ", .) %>% 
  paste0(., " dims - Pseudotime's Spearman: ", spear.corr)
pdf(file.path(output.dir, "mesoderm-development_loh_pseudotime_comp_all.pdf"), width=22, height=22)
cowplot::plot_grid(comp[[1]][[1]]$plot, comp[[5]][[1]]$plot, comp[[9]][[1]]$plot,
                   comp[[2]][[1]]$plot, comp[[6]][[1]]$plot, comp[[10]][[1]]$plot,  
                   comp[[3]][[1]]$plot, comp[[7]][[1]]$plot, comp[[11]][[1]]$plot,
                   comp[[4]][[1]]$plot, comp[[8]][[1]]$plot, comp[[12]][[1]]$plot,
                   ncol=3,
                   labels = label.plts[c(1, 5, 9, 2, 6, 10, 3, 7, 11, 4, 8, 12)])
dev.off()
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## R packages and versions used in these analyses
sessionInfo()
#
#------------------------------------------------------------------------------#
