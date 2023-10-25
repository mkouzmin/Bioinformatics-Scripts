
library(Seurat)
library(tidyverse)
library(cowplot)
library(harmony)

in_data_dir <- "/n/groups/flyrnai/mikhail/RAS_tumor/run_count_all/"


#for SCENIC
dbFiles <- "./Databases/"



samples <- dir(in_data_dir)

# create seurat object:
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

seurat_list[[2]] <- RenameCells(object = seurat_list[[2]], new.names = paste0(substring(Cells(x = seurat_list[[2]]),1,16),"-2"))
seurat_list[[3]] <- RenameCells(object = seurat_list[[3]], new.names = paste0(substring(Cells(x = seurat_list[[3]]),1,16),"-3"))
seurat_list[[4]] <- RenameCells(object = seurat_list[[4]], new.names = paste0(substring(Cells(x = seurat_list[[4]]),1,16),"-4"))
seurat_list[[5]] <- RenameCells(object = seurat_list[[5]], new.names = paste0(substring(Cells(x = seurat_list[[5]]),1,16),"-5"))
seurat_list[[6]] <- RenameCells(object = seurat_list[[6]], new.names = paste0(substring(Cells(x = seurat_list[[6]]),1,16),"-6"))

# merge seurat object
seurat_obj <- merge(x=seurat_list[[1]], y=seurat_list[2:length(seurat_list)])

# clean up
rm(seurat_list)
gc()

obj<-Seurat::NormalizeData(seurat_obj,verbose = FALSE)%>%
FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
ScaleData(verbose = FALSE)%>%
RunPCA(npcs=35)

top10 <- head(VariableFeatures(obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 

#check elbow plot
ElbowPlot(obj, ndims = 35)

options(repr.plot.height = 10, repr.plot.width = 12)
p1 <- DimPlot(object = obj, reduction = "pca", pt.size = .1, group.by = "SampleID")
p2 <- VlnPlot(object = obj, features = "PC_1", group.by = "SampleID", pt.size = .1)
plot_grid(p1,p2)

VlnPlot(seurat_obj, group.by="SampleID", features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size=0, )
obj <- subset(obj, subset = nFeature_RNA < 3000)


#Run Harmony
options(repr.plot.height = 2.5, repr.plot.width = 6)
obj <- obj %>% 
  RunHarmony("SampleID", plot_convergence = TRUE)

#check embeddings
harmony_embeddings <- Embeddings(obj, 'harmony')
harmony_embeddings[1:5, 1:5]

#plot of Samples with Harmony
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = obj, reduction = "harmony", pt.size = .1, group.by = "SampleID")
p2 <- VlnPlot(object = obj, features = "harmony_1", group.by = "SampleID", pt.size = .1)
plot_grid(p1,p2)

#no harmony
obj <- obj %>% 
  RunUMAP(dims = 1:20) %>% 
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()

options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(obj, reduction = "umap", group.by = "SampleID", pt.size = .1, split.by = 'SampleID')
DimPlot(obj, reduction = "umap", group.by = "SampleID", pt.size = .1)

#plot with harmony
obj <- obj %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()


options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(obj, reduction = "umap", group.by = "SampleID", pt.size = .1, split.by = 'SampleID')
DimPlot(obj, reduction = "umap", split.by = "SampleID", pt.size = .1)



#find Markers
obj.markers <- FindAllMarkers(obj, min.pct = 0.25)
marker_data<-obj.markers %>%
  group_by(cluster)
write.csv(marker_data,"marker_data_RAS.csv",row.names = FALSE)

sub_obj.markers <- FindAllMarkers(sub_obj, min.pct = 0.25)
marker_data<-sub_obj.markers %>%
  group_by(cluster)
write.csv(marker_data,"marker_data_final2.csv",row.names = FALSE)

marker_data_top<-obj.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
#marker_data_top
write.csv(marker_data_top,"marker_data_top_FC.csv",row.names = FALSE)


#ElbowPlot(obj, ndims = 30)

DimHeatmap(sub_obj, dims = 1:15, cells = 500, balanced = TRUE, ncol = 4)

umap = cbind("Barcode" = rownames(Embeddings(object = obj, reduction = "umap")), Embeddings(object = obj, reduction = "umap"))
write.table(umap, file="./umap_renal.csv", sep = ",", quote = F, row.names = F, col.names = T)

tsne = cbind("Barcode" = rownames(Embeddings(object = obj, reduction = "tsne")), Embeddings(object = obj, reduction = "tsne"))
write.table(tsne, file="./RAS_Tumor_T-SNE_v2.csv", sep = ",", quote = F, row.names = F, col.names = c("Barcode","t-SNE-1","t-SNE-2"))


umap_names = as.data.frame(as.matrix(obj$RNA_snn_res.0.4))
write.table(umap_names, file="./umap_renal_names.csv", sep = ",", quote = F, row.names = T, col.names = F)

t_SNE_names = as.data.frame(as.matrix(Idents(object = obj)))
t_SNE_names$Barcode<-rownames(t_SNE_names)
write.table(t_SNE_names[c(2,1)], file="./RAS_Tumor_T-SNE_v2_names.csv", sep = ",", quote = F, row.names = F, col.names = c("barcode","id"))

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

clusters <-read.table(file = "./RAS_and_8_names.csv", sep = ",", row.names = 1, header = T)

Idents(obj) <- clusters

#group by Sample without replicate number
group_metadata=paste0(substring(obj$SampleID,1,2),substring(obj$SampleID,nchar(obj$SampleID),nchar(obj$SampleID)))
group_metadata<-str_replace_all(group_metadata,c("8W1" = "W1","8W2"="W1","8Y3"="8","8Y4"="8","8Y5"="8","8Y6"="8","RA1"="W2","RA2"="W2","RA3"="R","RA4"="R","RA5"="R","RA6"="R"))
obj<-AddMetaData(obj, group_metadata, 'R_vs_8')
#find markers grouped by SampleID for cluster 0
obj.markers_8 <- FindMarkers(obj, min.pct = 0.25,ident.1 = '8', ident.2 = "W", group.by='R_vs_8')
write.csv(obj.markers_8,"marker_data_8_vs_W.csv")

obj.markers_R <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'R', ident.2 = "W", group.by='R_vs_8')
write.csv(obj.markers_R,"marker_data_R_vs_W.csv")

obj.markers_R_v_8 <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'R', ident.2 = "8", group.by='R_vs_8')
write.csv(obj.markers_R_v_8,"marker_data_R_vs_8.csv")

unique(Idents(obj))

obj.markers_8_c <- FindMarkers(obj, min.pct = 0.25,ident.1 = '8', ident.2 = "W", group.by='R_vs_8',subset.ident = 28)
obj.markers_8_c$SampleID = 1 
max(as.numeric(unique(Idents(obj))))
unique(Idents(obj))
length(unique(Idents(obj)))
obj.markers_8_by_clust_list <- lapply(0:(length(unique(Idents(obj)))-1), function(id){
  obj.markers_8_cur <- FindMarkers(obj,logfc.threshold = 0, ident.1 = '8', ident.2 = "W1", group.by='R_vs_8',subset.ident = id)
  obj.markers_8_cur$SampleID = id
  return(obj.markers_8_cur)
})
obj.markers_8_by_clust <- bind_rows(obj.markers_8_by_clust_list)

write.csv(obj.markers_8_by_clust,"marker_data_8_vs_W_by_clust_full.csv")


obj.markers_8_c <- FindMarkers(obj, min.pct = 0.25,ident.1 = '8', ident.2 = "W1", group.by='R_vs_8',subset.ident = 0)
obj.markers_8_c$SampleID = 1 
unique(Idents(obj))
length(unique(Idents(obj)))
obj.markers_8_by_clust_list <- lapply(0:(length(unique(Idents(obj)))-1), function(id){
  obj.markers_8_cur <- FindMarkers(obj, logfc.threshold = 0,ident.1 = '8', ident.2 = "W1", group.by='R_vs_8',subset.ident = id)
  obj.markers_8_cur$SampleID = id
  obj.markers_8_cur <- cbind(rownames(obj.markers_8_cur), data.frame(obj.markers_8_cur, row.names=NULL))
  return(obj.markers_8_cur)
})
obj.markers_8_by_clust <- bind_rows(obj.markers_8_by_clust_list)

write.csv(obj.markers_8_by_clust,"marker_data_8_vs_W_by_clust_full.csv",row.names = FALSE)

obj.markers_R_by_clust_list <- lapply(0:(length(unique(Idents(obj)))-1), function(id){
  obj.markers_R_cur <- FindMarkers(obj, logfc.threshold = 0, ident.1 = 'R', ident.2 = "W2", group.by='R_vs_8',subset.ident = id)
  obj.markers_R_cur$SampleID = id
  obj.markers_R_cur <- cbind(rownames(obj.markers_R_cur), data.frame(obj.markers_R_cur, row.names=NULL))
  return(obj.markers_R_cur)
})
obj.markers_R_by_clust <- bind_rows(obj.markers_R_by_clust_list)

write.csv(obj.markers_R_by_clust,"marker_data_R_vs_W_by_clust_full.csv", row.names = FALSE)

obj.markers_R_v_8_by_clust_list <- lapply(0:(length(unique(Idents(obj)))-1), function(id){
  obj.markers_R_v_8_cur <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'R', ident.2 = "8", group.by='R_vs_8',subset.ident = id)
  obj.markers_R_v_8_cur$SampleID = id
  obj.markers_R_v_8_cur <- cbind(rownames(obj.markers_R_v_8_cur), data.frame(obj.markers_R_v_8_cur, row.names=NULL))
  return(obj.markers_R_v_8_cur)
})
obj.markers_R_v_8_by_clust <- bind_rows(obj.markers_R_v_8_by_clust_list)

write.csv(obj.markers_R_v_8_by_clust,"marker_data_R_vs_8_by_clust.csv", row.names = FALSE)


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
