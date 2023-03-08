#----------------download, parse & convert anndata to SCE rds------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Description: R script to download 'human_cd34_bm_rep1.h5ad' anndata h5ad from 
# HCA Portal (reference project: 
# https://data.humancellatlas.org/explore/projects/091cf39b-01bc-42e5-9437-f419a66c8a45),
# parse it (log-normalization, remove unnecessary dimensional reduction, fix cluster 
# color names, add cell annotations - https://github.com/dpeerlab/Palantir/issues/40), 
# and save it as SCE RDS file format. 
# Date: 08/03/2023
# Last update: 08/03/2023
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Import packages
library("zellkonverter")
library("SingleCellExperiment")
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Download and import anndata h5ad file as SCE object 

# Download the anndata object 'human_cd34_bm_rep1.h5ad' from HCA Portal: 
# https://data.humancellatlas.org/explore/projects/091cf39b-01bc-42e5-9437-f419a66c8a45
# WARNING: the following data file url was used on the 08/03/2023 - 13:20 UTC time
# the link expires after some time and, thus, the 'file2down' character needs to be 
# updated by the user before running it 
file2down <- "https://storage.googleapis.com/datarepo-4ef3f5a3-bucket/0e5e329e-2709-4ceb-bfe3-97d23a652ac0/3428f967-3376-4051-b6f7-8dd84580ca5b/human_cd34_bm_rep1.h5ad?X-Goog-Algorithm=GOOG4-RSA-SHA256&X-Goog-Credential=datarepo-jade-api%40terra-datarepo-production.iam.gserviceaccount.com%2F20230308%2Fauto%2Fstorage%2Fgoog4_request&X-Goog-Date=20230308T132036Z&X-Goog-Expires=900&X-Goog-SignedHeaders=host&requestedBy=azul-ucsc-0-public-prod%40platform-hca-prod.iam.gserviceaccount.com&userProject=datarepo-1dbff5cd&X-Goog-Signature=9f146db3f8977359f51915a01c7eebfd4877e3ee981e015c286cf2545d3d305c5156538d8384b2b17eff2c91f243585b5f790433c67b0d393dca395eced58b4ac2148b8c592292c88712395e981fbc0dec470a8e9eae2de26e6cbee3c1f6be6daac2b6b897d96d65686e61a575ecea68e6c12179eb4dcf598c05de6f43d29c539a4a49425ca10f91c68e8a3b7ee1a47bbc20cbf3c108bed894de82c5377a306088d3a26e8ebaa514305584e69bab0f2728d4679862a9b273476007654f1c245ac5a17dd9973bde5e329fd834961437ed3134321ddba083a7153d5d6d05aaeba93dd33613a28eb8ec7712f56b5acb216c4188cac4825713b749427e6d7c3837ca" # url to download file
file2save <- "data/human_cd34_bm_rep1.h5ad" # path to save the file
download.file(url = file2down, destfile = file2save)

# Import the anndata '.h5ad' object as SCE object
sce <- readH5AD(file = file2save, X_name = "scaled", raw = TRUE)
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Parse SCE object

# Add 'counts' to the main expression object & Perform log1p normalization w/ 10K scaling factor
stopifnot(all(colnames(altExp(sce, "raw"))==colnames(sce))) # check that cell names match
stopifnot(all(row.names(altExp(sce, "raw"))==row.names(sce))) # check that gene names match
assay(sce, "counts") <- assay(altExp(sce, "raw"), "X") # add counts to main exp obj
altExp(sce, "raw") <- NULL # delete alternative exp obj 
log.counts <- apply(assay(sce, "counts"), 2, function(x) log1p(x/sum(x)*10000)) # log1p normalization w/ 10K scaling factor
assay(sce, "logcounts") <- as(log.counts, "sparseMatrix")
rm(log.counts); gc();  # remove unnecessary objects
assays(sce) <- assays(sce)[c("counts", "logcounts", "scaled")] # reorder assays

# Change unnecessary reductions to 'metadata(sce)'
metadata(sce)[["MAGIC_imputed_data"]] <- reducedDim(sce, "MAGIC_imputed_data")
metadata(sce)[["palantir_branch_probs"]] <- reducedDim(sce, "palantir_branch_probs")
reducedDims(sce)[c("MAGIC_imputed_data", "palantir_branch_probs")] <- NULL

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
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Export SCE RDS file
saveRDS(object = sce, file = "data/human_cd34_bm_rep1.rds")
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## R packages and versions used in these analyses
sessionInfo()
#
#------------------------------------------------------------------------------#