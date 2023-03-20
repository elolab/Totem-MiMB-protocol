### Script to generate MS figures

# Load packages
library("dplyr")
library("Totem")
library("scater")
library("ggplot2")
library("SingleCellExperiment")

# General params: 
out.folder <- "results/ms_figs"
if (!dir.exists(out.folder)) dir.create(out.folder, recursive = TRUE)

## QuickStart protocol figs: 

# Params: 
quickstart <- "results/sce_quickstart.rds"
sce <- readRDS(quickstart)

# Figure 1: 
dim_red <- reducedDim(sce, "tsne")
pdf(file = file.path(out.folder, "figure_1.pdf"), width=4, height=4)
VizCellConnectivity(object = sce, viz.dim.red = dim_red)
dev.off()

# Figure 2: 
select.cluster <- names(metadata(sce)$totem$slingshot_trajectory) 
pdf(file = file.path(out.folder, "figure_2.pdf"), width=12, height=4)
print(
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
)
dev.off()

# Figure 3: 
pdf(file = file.path(out.folder, "figure_3.pdf"), width=8, height=4)
print(
  cowplot::plot_grid(
    VizSmoothedTraj(object = sce,
                    traj.names = select.cluster,
                    viz.dim.red = dim_red, plot.pseudotime = FALSE),
    VizSmoothedTraj(object = sce,
                    traj.names = select.cluster,
                    viz.dim.red = dim_red, plot.pseudotime = TRUE), 
    ncol=2
  )
)
dev.off()

# Figure 4: 
pdf(file = file.path(out.folder, "figure_4.pdf"), width=8, height=4)
print(
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
)
dev.off()

# Remove object
rm(sce); gc(); 

## GuidedStart protocol figs: 

# Params: 
guidedstart <- "results/sce_guidedstart.rds"
sce <- readRDS(guidedstart)

# Figure 5: 
# Selection of HVG w/ scran R package
var.genes <- scran::modelGeneVar(sce)
hvg <- scran::getTopHVGs(var.genes, n = 2000)
# PCA
sce <- RunDimRed(object = sce,
                 dim.red.method = "pca",
                 dim.red.features = hvg,
                 dim.reduction.par.list = list(ndim=50))
# Elbow plot
elbow.plt <- reducedDim(sce, "pca") %>% 
  apply(X = ., MARGIN = 2, FUN = function(x) sd(x)) %>%
  as.data.frame(.) %>% 
  `colnames<-`("Standard Deviation") %>% 
  mutate("PCs" = factor(1:nrow(.), levels = 1:nrow(.))) %>% 
  ggplot(data = ., mapping = aes(x = PCs, y = `Standard Deviation`)) + 
  geom_point() + 
  theme_bw() + 
  theme(axis.text = element_text(size=6))
# PCA plot
pca.plt <- reducedDim(sce, "pca") %>% 
  as.data.frame(.) %>%
  mutate("Cell_ID" = row.names(.)) %>% 
  cbind(., colData(sce)) %>% 
  ggplot(data = ., mapping = aes(x=comp_1, y=comp_2, color=cell_types_short)) + 
  geom_point(size=0.25) + 
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
pdf(file = file.path(out.folder, "figure_5.pdf"), width=9.5, height=5)
print(
  cowplot::plot_grid(cowplot::plot_grid(elbow.plt, pca.plt, ncol=2, rel_widths = c(0.55, 0.45)), pca.scores.plt, ncol=1)
)
dev.off()

# Figure 6: 
dim_red <- reducedDim(sce, "umap")
pdf(file = file.path(out.folder, "figure_6.pdf"), width=4, height=4)
plotReducedDim(object = sce, dimred = "umap", colour_by = "cell_types_short", point_size=0.5) + 
  scale_color_manual(name= "Cell Types", 
                     values=as.character(unlist(metadata(sce)$cluster_colors)[!duplicated(levels(sce$cell_types_short))])) + 
  theme_void()
dev.off()

# Figure 7:
pdf(file = file.path(out.folder, "figure_7.pdf"), width=4, height=4)
VizCellConnectivity(object = sce, viz.dim.red = dim_red)
dev.off()

# Figure 8:
select.clusters <- ReturnTrajNames(sce)
pdf(file = file.path(out.folder, "figure_8.pdf"), width=12, height=12)
VizMST(object = sce, clustering.names = select.clusters, viz.dim.red = dim_red)
dev.off()

# Figure 9: 
smooth.msts.names <- ReturnTrajNames(sce)
pdf(file = file.path(out.folder, "figure_9.pdf"), width=12, height=12)
VizSmoothedTraj(object = sce,
                traj.names = smooth.msts.names,
                viz.dim.red = dim_red,plot.pseudotime = FALSE)
dev.off()

# Figure 10: 
select.traj <- "8.135"
pdf(file = file.path(out.folder, "figure_10.pdf"), width=8, height=4)
print(
  cowplot::plot_grid(
    VizSmoothedTraj(object = sce,
                    traj.names = select.traj,
                    viz.dim.red = dim_red, plot.pseudotime = FALSE),
    VizSmoothedTraj(object = sce,
                    traj.names = select.traj,
                    viz.dim.red = dim_red, plot.pseudotime = TRUE), 
    ncol=2
  )
)
dev.off()

# Figure 11: 
pdf(file = file.path(out.folder, "figure_11.pdf"), width=8, height=4)
print(
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
)
dev.off()

# Remove object
rm(sce); gc();
