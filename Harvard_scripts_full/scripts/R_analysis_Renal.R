library(cowplot)
theme_set(theme_cowplot())
library(RColorBrewer)
library(viridis)

umap_theme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_blank()
)

pdf(paste0(fig_dir, "qc_violin_plot.pdf"), width=10, height=10)
# png(paste0(fig_dir, "pngs/qc_violin_plot.png"), width=10, height=10, res=250, units='in')
VlnPlot(seurat_obj, group.by="SampleID", features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size=0, )
dev.off()

#seurat_obj <- subset(seurat_obj, nFeature_RNA > 200 & nFeature_RNA < 10000 )
#seurat_obj <- subset(seurat_obj, SampleID != 'Sample-101')

# compute PCA:
seurat_obj <- RunPCA(seurat_obj)

# plot pca heatmap
png(paste0(fig_dir, "pngs/pca_heatmap.png"), width=10, height=10, res=300, units='in')
DimHeatmap(seurat_obj, dims = 1:12, cells = 500, balanced = TRUE)
dev.off()

# plot variance explained by PC:
png(paste0(fig_dir, "pngs/pca_elbow_plot.png"), width=5, height=3, res=300, units='in')
ElbowPlot(seurat_obj)
dev.off()

# UMAP and clustering with top PCs
seurat_obj <- RunUMAP(seurat_obj, reduction='pca', dims = 1:30)
seurat_obj <- FindNeighbors(seurat_obj, reduction='pca')
seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)