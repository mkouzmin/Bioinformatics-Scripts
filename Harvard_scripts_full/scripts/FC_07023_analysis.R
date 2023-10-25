
library(Seurat)
library(tidyverse)
library(cowplot)
library(harmony)

in_data_dir <- "/n/groups/flyrnai/mikhail/FC_07023_Liz/run_count_all/"

#for SCENIC
dbFiles <- "./Databases/"

# create seurat object:

samples <- dir(in_data_dir)


seurat_list <- lapply(samples, function(sample){
  cur_data <- Read10X(paste0(in_data_dir,sample,'/outs/filtered_feature_bc_matrix/'))
  cur_seurat <- CreateSeuratObject(
    counts = cur_data,
    min.cells=3,
    min.features=200,
    project='RENAL2'
  )
  cur_seurat$SampleID <- sample
  return(cur_seurat)
})
#rename cells to  Loupe naming convention - make sure this is in same order as aggr csv
seurat_list[[2]] <- RenameCells(object = seurat_list[[2]], new.names = paste0(substring(Cells(x = seurat_list[[2]]),1,16),"-2"))
seurat_list[[3]] <- RenameCells(object = seurat_list[[3]], new.names = paste0(substring(Cells(x = seurat_list[[3]]),1,16),"-3"))
seurat_list[[4]] <- RenameCells(object = seurat_list[[4]], new.names = paste0(substring(Cells(x = seurat_list[[4]]),1,16),"-4"))

# merge seurat object
seurat_obj <- merge(x=seurat_list[[1]], y=seurat_list[2:length(seurat_list)])

# clean up
rm(seurat_list)
gc()

obj<-Seurat::NormalizeData(seurat_obj,verbose = FALSE)%>%
FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>% 
ScaleData(verbose = FALSE)%>%
RunPCA(npcs=30)

View(obj@assays$RNA@scale.data)

expr_mat <- as.matrix(GetAssayData(object = seurat_obj, slot = "counts"))
write.csv(expr_mat,"Liz_count_matrix_raw.csv")

obj<-Seurat::NormalizeData(seurat_obj,verbose = FALSE)%>%
  ScaleData(verbose = FALSE,do.center = FALSE)

write.csv(obj@assays$RNA@scale.data,"Liz_count_matrix_scaled_normalized.csv")


top10 <- head(VariableFeatures(obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#check elbow plot
ElbowPlot(obj, ndims = 30)

options(repr.plot.height = 10, repr.plot.width = 12)
p1 <- DimPlot(object = obj, reduction = "pca", pt.size = .1, group.by = "SampleID")
p2 <- VlnPlot(object = obj, features = "PC_1", group.by = "SampleID", pt.size = .1)
plot_grid(p1,p2)

VlnPlot(obj, group.by="SampleID", features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size=0, )
#obj <- subset(obj, subset = nFeature_RNA < 2500 & nCount_RNA < 10000)
obj <- subset(obj, subset = nFeature_RNA < 4000)
obj <- subset(obj, subset = nCount_RNA < 40000)
#Run Harmony
options(repr.plot.height = 2.5, repr.plot.width = 6)
obj <- obj %>% 
  RunHarmony("SampleID", plot_convergence = TRUE)

#check embeddings
harmony_embeddings <- Embeddings(obj, 'harmony')
harmony_embeddings[1:5, 1:5]
write.csv(harmony_embeddings,"Liz_harmony_pca_embeddings.csv")

#plot of Samples with Harmony
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = obj, reduction = "harmony", pt.size = .1, group.by = "SampleID")
p2 <- VlnPlot(object = obj, features = "harmony_1", group.by = "SampleID", pt.size = .1)
plot_grid(p1,p2)
#umap_no_harmony
obj <- obj %>% 
  RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()

options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(obj, reduction = "umap", group.by = "SampleID", pt.size = .1, split.by = 'SampleID')
DimPlot(obj, reduction = "umap", group.by = "SampleID", pt.size = .1)
DimPlot(obj, reduction = "umap", pt.size = .1, split.by = 'SampleID')
DimPlot(obj, reduction = "umap", pt.size = .1)
FeaturePlot(obj,reduction= "umap", split.by = 'SampleID', features = "Octalpha2R")

obj <- obj %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()

options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(obj, reduction = "umap", group.by = "SampleID", pt.size = .1, split.by = 'SampleID')
DimPlot(obj, reduction = "umap", group.by = "SampleID", pt.size = .1)
obj <- obj %>% 
  RunTSNE(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()

options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(obj, reduction = "tsne", group.by = "SampleID", pt.size = .1, split.by = 'SampleID')


options(repr.plot.height = 4, repr.plot.width = 6)
DimPlot(obj, reduction = "umap", label = TRUE, pt.size = .1)

#remove clusters
sub_obj <- subset(obj, idents = c(2))
#re-do clustering
DimPlot(sub_obj, reduction = "umap", label = TRUE, pt.size = .1)
sub_obj <- sub_obj %>% 
  RunUMAP(reduction = "harmony",dims = 1:10) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:10) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()


DimPlot(sub_obj, reduction = "umap", pt.size = .5, split.by = 'SampleID')
DimPlot(sub_obj, reduction = "umap", group.by = "SampleID", pt.size = .5)
FeaturePlot(sub_obj,reduction= "umap", features = "Octalpha2R")

sub_obj <- RenameIdents(object = sub_obj, `0` = "main segment principal cells", `1` = "upper ureter principal cells", `2` = "renal stem cells", `3` = "lower tubule principal cells", `4` = "initial and transitional principal cells", `5` = "lower segment principal cells", `6` = "stellate cells", `7` = "lower ureter principal cells", `8` = "main segment principal cells")
options(repr.plot.height = 4, repr.plot.width = 6)
DimPlot(sub_obj, reduction = "umap", label = TRUE, pt.size = .1)
#sort by positive fold change

#find Markers
obj.markers <- FindAllMarkers(obj, min.pct = 0.25)
marker_data<-obj.markers %>%
  group_by(cluster)
write.csv(marker_data,"marker_data_Liz.csv",row.names = FALSE)

obj.markers_Liz_7_v_3 <- FindMarkers(obj, min.pct = 0.25,ident.1 = '7', ident.2 = "3")
write.csv(obj.markers_Liz_7_v_3,"marker_data_Liz_no_batch_correct_7_v_3.csv")

marker_data_top<-obj.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
#marker_data_top
write.csv(marker_data_top,"marker_data_top_FC.csv",row.names = FALSE)


#ElbowPlot(obj, ndims = 30)

DimHeatmap(sub_obj, dims = 1:15, cells = 500, balanced = TRUE, ncol = 4)

umap = cbind("Barcode" = rownames(Embeddings(object = obj, reduction = "umap")), Embeddings(object = obj, reduction = "umap"))
write.table(umap, file="./umap_FC_07023.csv", sep = ",", quote = F, row.names = F, col.names = T)


umap_names = as.data.frame(as.matrix(obj$seurat_clusters))
write.table(umap_names, file="./FC_07023_names.csv", sep = ",", quote = F, row.names = T, col.names = F)


umap = cbind("Barcode" = rownames(Embeddings(object = sub_obj, reduction = "umap")), Embeddings(object = sub_obj, reduction = "umap"))
write.table(umap, file="./umap_FC_07023_sub_obj.csv", sep = ",", quote = F, row.names = F, col.names = T)


umap_names = as.data.frame(as.matrix(sub_obj$seurat_clusters))
write.table(umap_names, file="./FC_07023_names_sub_obj.csv", sep = ",", quote = F, row.names = T, col.names = F)


#get SampleIDs
ids= as.data.frame(as.matrix(obj$SampleID))
write.table(ids, file="./umap_renal_sampleIDs.csv", sep = ",", quote = F, row.names = T, col.names = F)

#umap for subset object
umap = cbind("Barcode" = rownames(Embeddings(object = sub_obj, reduction = "umap")), Embeddings(object = sub_obj, reduction = "umap"))
write.table(umap, file="./sub_umap_renal.csv", sep = ",", quote = F, row.names = F, col.names = T)


umap_names = as.data.frame(as.matrix(Idents(object = sub_obj)))
write.table(umap_names, file="./sub_umap_renal_names_3.csv", sep = ",", quote = F, row.names = T, col.names = F)

counts = as.data.frame(as.matrix(as.matrix(sub_obj@assays$RNA@counts)))
write.table(counts, file="./subset_counts_matrix_2", sep = ",", quote = F, row.names = T, col.names = T)

#group by Sample without replicate number
group_metadata=paste0(substring(obj$SampleID,1,2),substring(obj$SampleID,nchar(obj$SampleID),nchar(obj$SampleID)))
group_metadata=obj$SampleID
obj<-AddMetaData(obj, group_metadata, 'Sample')

obj.markers_2_v_1_by_clust_list <- lapply(0:(length(unique(Idents(obj)))-1), function(id){
  obj.markers_2_v_1_cur <- FindMarkers(obj,logfc.threshold = 0, ident.1 = 'L2', ident.2 = "L1", group.by='SampleID',subset.ident = id)
  obj.markers_2_v_1_cur$SampleID = id
  obj.markers_2_v_1_cur <- cbind(rownames(obj.markers_2_v_1_cur), data.frame(obj.markers_2_v_1_cur, row.names=NULL))
  return(obj.markers_2_v_1_cur)
})
obj.markers_2_v_1_by_clust <- bind_rows(obj.markers_2_v_1_by_clust_list)

write.csv(obj.markers_2_v_1_by_clust,"marker_data_Liz_2_v_1_full.csv",row.names = FALSE)

obj.markers_3_v_1_by_clust_list <- lapply(0:(length(unique(Idents(obj)))-1), function(id){
  obj.markers_3_v_1_cur <- FindMarkers(obj,logfc.threshold = 0, ident.1 = 'L3', ident.2 = "L1", group.by='SampleID',subset.ident = id)
  obj.markers_3_v_1_cur$SampleID = id
  obj.markers_3_v_1_cur <- cbind(rownames(obj.markers_3_v_1_cur), data.frame(obj.markers_3_v_1_cur, row.names=NULL))
  return(obj.markers_3_v_1_cur)
})
obj.markers_3_v_1_by_clust <- bind_rows(obj.markers_3_v_1_by_clust_list)

write.csv(obj.markers_3_v_1_by_clust,"marker_data_Liz_3_v_1_full.csv",row.names = FALSE)

obj.markers_4_v_2_by_clust_list <- lapply(0:(length(unique(Idents(obj)))-1), function(id){
  obj.markers_4_v_2_cur <- FindMarkers(obj,logfc.threshold = 0, ident.1 = 'L4', ident.2 = "L2", group.by='SampleID',subset.ident = id)
  obj.markers_4_v_2_cur$SampleID = id
  obj.markers_4_v_2_cur <- cbind(rownames(obj.markers_4_v_2_cur), data.frame(obj.markers_4_v_2_cur, row.names=NULL))
  return(obj.markers_4_v_2_cur)
})
obj.markers_4_v_2_by_clust <- bind_rows(obj.markers_4_v_2_by_clust_list)

write.csv(obj.markers_4_v_2_by_clust,"marker_data_Liz_4_v_2_full.csv",row.names = FALSE)

obj.markers_4_v_3_by_clust_list <- lapply(0:(length(unique(Idents(obj)))-1), function(id){
  obj.markers_4_v_3_cur <- FindMarkers(obj,logfc.threshold = 0, ident.1 = 'L4', ident.2 = "L3", group.by='SampleID',subset.ident = id)
  obj.markers_4_v_3_cur$SampleID = id
  obj.markers_4_v_3_cur <- cbind(rownames(obj.markers_4_v_3_cur), data.frame(obj.markers_4_v_3_cur, row.names=NULL))
  return(obj.markers_4_v_3_cur)
})
obj.markers_4_v_3_by_clust <- bind_rows(obj.markers_4_v_3_by_clust_list)

write.csv(obj.markers_4_v_3_by_clust,"marker_data_Liz_4_v_3_full.csv",row.names = FALSE)








#Use SCENIC to do regulatory network analysis
cellInfo <- data.frame(seuratCluster=Idents(sub_obj))
expr_mat <- as.matrix(GetAssayData(object = sub_obj, slot = "counts"))
write.csv(expr_mat,"Expr_Mat_for_SCENIC")
library(SCENIC)
scenicOptions <- initializeScenic(org="dmel", dbDir=dbFiles, nCores=10)
names(cellInfo)[1]<-"CellType"
cellInfo$nGene <- sub_obj@meta.data$nFeature_RNA
cellInfo$nUMI <- sub_obj@meta.data$nCount_RNA
saveRDS(cellInfo, file="int/cellInfo.Rds")
#scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
genekept<-geneFiltering(expr_mat,scenicOptions = scenicOptions)
exprMat_filtered <- expr_mat[genekept,]
exprMat_filtered_log <- log2(exprMat_filtered+1) 
#running Genie takes a long time - do it as part of a batch script
#write.csv(exprMat_filtered,"Expr_Mat_for_GENIE3")
#exportsForArboreto(exprMat_filtered, scenicOptions)
#runGenie3(exprMat_filtered_log, scenicOptions)
runCorrelation(exprMat_filtered, scenicOptions)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
#resource intensive:
#scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
#scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, expr_mat, skipHeatmap=TRUE, skipTsne=TRUE)
scenicOptions <- readRDS("int/scenicOptions.Rds")
