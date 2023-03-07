#--------------Testing packer_embryo Monocle3's data set script----------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Description: R script to test Monocle3 packer_embryo's data with Totem 
# Date: 24/02/2023
# Last update: 24/02/2023
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
## Download 'packer_embryo' data used in Monocle3 tutorial: 
# https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/
# urls2down <- list("gene_exp" = "https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_expression.rds", 
#                   "cell_meta" = "https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_colData.rds", 
#                   "gene_meta" = "https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_rowData.rds")
data.set <- "packer_embryo_monocle3"
# scrna <- list()
# for (u in names(urls2down)) {
#   file2import <- file.path("data", data.set, basename(urls2down[[u]]))
#   if (!dir.exists(dirname(file2import))) dir.create(dirname(file2import), recursive = TRUE)
#   #download.file(url = urls2down[[u]], destfile = file2import)
#   scrna[[u]] <- readRDS(file2import) 
# }
out.dir <- file.path("results", data.set) 
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
file2import <- "results/repro_monocle3_TI_tutorial/packer_embryo_cds.rds"
cds <- readRDS(file = file2import)
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Upstream Workflow
# Create a SCE object for Totem
# stopifnot(all(colnames(scrna$gene_exp)==row.names(scrna$cell_meta)))
# stopifnot(all(row.names(scrna$gene_exp)==row.names(scrna$gene_meta)))
counts <- assay(cds, "counts")
row.names(counts) <- rowData(cds)[,"gene_short_name"]
log.counts <- apply(counts, 2, function(x) log1p(x/sum(x)*10000)) # log1p normalization w/ 10K scaling factor

# Create SCE object
sce <- SingleCellExperiment(assays = list("counts" = counts, "logcounts" = log.counts), 
                            colData = colData(cds), 
                            reducedDims = list("umap"=reducedDim(cds, type = "UMAP")))

## Prepare data for Totem: remove non-expressed genes
sce <- PrepareTotem(object = sce)
#rm(cds); gc(); 
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Downstream workflow
output <- file.path("results", data.set, paste0(data.set, ".pdf"))
if (!dir.exists(dirname(output))) dir.create(dirname(output), recursive = TRUE)

# Set seed
set.seed(1204)

## Clustering scRNA-seq data w/ CLARA (k-medoids), MSTs & smoothing 
#w/ principal curves algorithm
set.seed(123)
sce <- RunClustering(object = sce) %>% 
  SelectClusterings(object = .) %>% 
  RunSmoothing(object = .)

## Define root
root.cluster <- 2
select.cluster <- names(metadata(sce)$totem$slingshot_trajectory)
sce <- ChangeTrajRoot(object = sce, traj.name = select.cluster, 
                      root.cluster = root.cluster)

## Visualize pseudotime
## Export the results as PDF
pdf(output, width = 16)
print(
  VizCellConnectivity(object = sce, 
                      viz.dim.red = reducedDim(sce))
)
print(
  VizSmoothedTraj(object = sce,
                  traj.names = select.cluster,
                  viz.dim.red = reducedDim(sce), 
                  plot.pseudotime = FALSE)
)
print(
  cowplot::plot_grid(
    VizSmoothedTraj(object = sce,
                  traj.names = select.cluster,
                  viz.dim.red = reducedDim(sce),
                  plot.pseudotime = TRUE), 
    monocle3::plot_cells(cds,
               color_cells_by = "pseudotime",
               label_cell_groups=FALSE,
               label_leaves=FALSE,
               label_branch_points=FALSE,
               graph_label_size=1.5),
    ncol = 2)
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
