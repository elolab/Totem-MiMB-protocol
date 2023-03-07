# Load packages
library("dplyr")
library("Totem")
library("scater")
library("ggplot2")
library("zellkonverter")
library("SingleCellExperiment")

# Set seed
set.seed(1204)

## Download the anndata object 'human_cd34_bm_rep1.h5ad' from HCA Portal: 
# https://data.humancellatlas.org/explore/projects/091cf39b-01bc-42e5-9437-f419a66c8a45
data.dir <- "../data/quickstart" # create folder to save the data 
dir.create(path = data.dir, recursive = TRUE)
file2down <- file.path(data.dir, "human_cd34_bm_rep1.h5ad") # path to save the file

# Check if the data exists locally before attempting to download it
if (!file.exists(file2down)) { # the link to download expires, thus it needs to be updated
  url2down <- "https://storage.googleapis.com/datarepo-4ef3f5a3-bucket/0e5e329e-2709-4ceb-bfe3-97d23a652ac0/3428f967-3376-4051-b6f7-8dd84580ca5b/human_cd34_bm_rep1.h5ad?X-Goog-Algorithm=GOOG4-RSA-SHA256&X-Goog-Credential=datarepo-jade-api%40terra-datarepo-production.iam.gserviceaccount.com%2F20230302%2Fauto%2Fstorage%2Fgoog4_request&X-Goog-Date=20230302T141550Z&X-Goog-Expires=900&X-Goog-SignedHeaders=host&requestedBy=azul-ucsc-0-public-prod%40platform-hca-prod.iam.gserviceaccount.com&userProject=datarepo-1dbff5cd&X-Goog-Signature=2c5c19abfe88e1d2d2fd7e50c3135e878f757183ddf383fa43e12285338998f1419dfb73dd4185a89ea9135a696feb65042fee38c4d833f7cf255af4b2e4be53a6c0ff8ab2f653c7ca1843fd31e850ad7d8532d598ca44b074c47d1dc80eb508a4ae8f843f8ba12869cc88f51037c50cd802c02061f9bd8ff1de3bd082b61135f635468c643bdf73d04b8f2ede122d5871848b519b10cc6d56bb3b03cdcd7bfc8996e23f4fe3b8685b51f5770d5c62466388cd1f77eaf2b29242faae3612f21e949600fc03cce94b1c97c050f9fc53aedd500966de0221a429a4541c62b3a81f8024182dfe177d78d1ea008757ca839a726dda6f975c5d34c0ac9ff3034d9b18"
  download.file(url = url2down, destfile = file2down)
}

## Import the anndata '.h5ad' object as SCE object
sce <- readH5AD(file = file2down, X_name = "scaled", raw = TRUE)

## Remove unnecessary reductions
reducedDims(sce)[c("MAGIC_imputed_data", "palantir_branch_probs")] <- NULL

## logcounts
stopifnot(all(colnames(altExp(sce, "raw"))==colnames(sce)))
stopifnot(all(row.names(altExp(sce, "raw"))==row.names(sce)))
counts <- assay(altExp(sce, "raw"), "X")
log.counts <- apply(counts, 2, function(x) log1p(x/sum(x)*10000)) # log1p normalization w/ 10K scaling factor
assay(sce, "logcounts") <- as(log.counts, "sparseMatrix")
rm(list = c("counts", "log.counts")) # remove unnecessary objects

## Cell cluster annotations
# Format names of the cluster colors of original clusters
names(metadata(sce)$cluster_colors) <- levels(sce$clusters) 

# Add cell annotations based on github issue: 
# https://github.com/dpeerlab/Palantir/issues/40
cluster2annot.long <- c("0" = "Hematopoietic stem cells", "1" = "Hematopoietic multipotent progenitors", 
                        "2" = "Erythroid progenitors", "3" = "Monocyte progenitors", 
                        "4" = "Myeloid progenitors", "5" = "Common lymphoid progenitors", 
                        "6" = "Monocyte progenitors", "7" = "Dendritic cell progenitors", 
                        "8" = "Erythroid progenitors", "9" = "Megakaryocyte progenitors")
cluster2annot.short <- c("0" = "HSC", "1" = "HMP", "2" = "EP", "3" = "MoP", 
                         "4" = "MyP", "5" = "CLP", "6" = "MP", "7" = "DCP", 
                         "8" = "EP", "9" = "MKP")
colData(sce)$cell_types_long <- colData(sce)$clusters
levels(colData(sce)$cell_types_long) <- cluster2annot.long
colData(sce)$cell_types_short <- colData(sce)$clusters
levels(colData(sce)$cell_types_short) <- cluster2annot.short

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
                       values=as.character(unlist(metadata(sce)$cluster_colors)[!duplicated(cluster2annot.short)])) +
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
