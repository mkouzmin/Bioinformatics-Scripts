
library(Seurat)
library(tidyverse)
library(cowplot)
library(harmony)

in_data_dir <- "/n/groups/flyrnai/mikhail/FC_07486_JA/run_count_all/"
fig_dir <- "/n/groups/flyrnai/mikhail/RAS_tumor_and_8_day/batch_correction_figs/"
out_data_dir <- "/n/groups/flyrnai/mikhail/RAS_tumor_and_8_day/batch_correction_data/"

#for SCENIC
dbFiles <- "./Databases/"



samples <- dir(in_data_dir)

# create seurat object w/o SoupX:
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

#run SoupX
seurat_list <- lapply(samples, function(sample){
  cur_filtered_data <- Read10X(paste0(in_data_dir,sample,'/outs/filtered_feature_bc_matrix/'))
  cur_raw_data <- Read10X(paste0(in_data_dir,sample,'/outs/raw_feature_bc_matrix/'))
  cur_sc=SoupChannel(cur_raw_data, cur_filtered_data)
  cur_meta = read.table(paste0(in_data_dir,sample,"/outs/analysis/clustering/graphclust/clusters.csv"),header = T, row.names = 1, sep=",")
  cur_DR = read.table(paste0(in_data_dir,sample,"/outs/analysis/umap/2_components/projection.csv"),header = T, row.names = 1, sep=",")
  cur_meta<-merge(cur_meta,cur_DR, by="row.names")
  row.names(cur_meta)<-cur_meta$Row.names
  cur_meta<-cur_meta[2:length(cur_meta)]
  cur_sc = setClusters(cur_sc, setNames(cur_meta$Cluster, row.names(cur_meta)))
  cur_sc = setDR(cur_sc, cur_meta[colnames(cur_sc$toc), c("UMAP.1", "UMAP.2")])
  cur_sc = autoEstCont(cur_sc)
  cur_out = adjustCounts(cur_sc)
  cur_seurat <- CreateSeuratObject(
    counts = cur_out,
    min.cells=3,
    min.features=200,
    project='RENAL2'
  )
  cur_seurat$SampleID <- sample
  return(cur_seurat)
})


#run SoupX - estimate with Ras85D / Gal4
seurat_list <- lapply(samples[1:4], function(sample){
  cur_filtered_data <- Read10X(paste0(in_data_dir,sample,'/outs/filtered_feature_bc_matrix/'))
  cur_raw_data <- Read10X(paste0(in_data_dir,sample,'/outs/raw_feature_bc_matrix/'))
  cur_sc=SoupChannel(cur_raw_data, cur_filtered_data)
  cur_meta = read.table(paste0(in_data_dir,sample,"/outs/analysis/clustering/graphclust/clusters.csv"),header = T, row.names = 1, sep=",")
  cur_DR = read.table(paste0(in_data_dir,sample,"/outs/analysis/umap/2_components/projection.csv"),header = T, row.names = 1, sep=",")
  cur_meta<-merge(cur_meta,cur_DR, by="row.names")
  row.names(cur_meta)<-cur_meta$Row.names
  cur_meta<-cur_meta[2:length(cur_meta)]
  cur_sc = setClusters(cur_sc, setNames(cur_meta$Cluster, row.names(cur_meta)))
  cur_sc = setDR(cur_sc, cur_meta[colnames(cur_sc$toc), c("UMAP.1", "UMAP.2")])
  cur_sc = autoEstCont(cur_sc)
  igGenes = c("Ras85D","gal4")
  head(cur_sc$soupProfile[order(cur_sc$soupProfile$est, decreasing = TRUE), ], n = 20)
  plotMarkerDistribution(cur_sc)
  useToEst = estimateNonExpressingCells(cur_sc, nonExpressedGeneList = list(Ras = c("Ras85D","gal4")), 
                                        clusters = FALSE)
  plotMarkerMap(cur_sc, geneSet = igGenes, useToEst = useToEst)
  cur_sc = calculateContaminationFraction(cur_sc, list(IG = igGenes), useToEst = useToEst)
  cur_out = adjustCounts(cur_sc)
  cur_seurat <- CreateSeuratObject(
    counts = cur_out,
    min.cells=3,
    min.features=200,
    project='RENAL2'
  )
  cur_seurat$SampleID <- sample
  return(cur_seurat)
})

gg = plotMarkerMap(cur_sc, "Ras85D")
plot(gg)

#use this to run only on "Wild Type" Samples

#seurat_list <- lapply(samples[5:8], function(sample){
  #run SoupX
  #cur_filtered_data <- Read10X(paste0(in_data_dir,sample,'/outs/filtered_feature_bc_matrix/'))
  #cur_raw_data <- Read10X(paste0(in_data_dir,sample,'/outs/raw_feature_bc_matrix/'))
  #cur_sc=SoupChannel(cur_raw_data, cur_filtered_data)
  #cur_meta = read.table(paste0(in_data_dir,sample,"/outs/analysis/clustering/graphclust/clusters.csv"),header = T, row.names = 1, sep=",")
  #cur_DR = read.table(paste0(in_data_dir,sample,"/outs/analysis/umap/2_components/projection.csv"),header = T, row.names = 1, sep=",")
  #cur_meta<-merge(cur_meta,cur_DR, by="row.names")
  #row.names(cur_meta)<-cur_meta$Row.names
  #cur_meta<-cur_meta[2:length(cur_meta)]
  #cur_sc = setClusters(cur_sc, setNames(cur_meta$Cluster, row.names(cur_meta)))
  #cur_sc = setDR(cur_sc, cur_meta[colnames(cur_sc$toc), c("UMAP.1", "UMAP.2")])
  #cur_sc = autoEstCont(cur_sc)
  #cur_out = adjustCounts(cur_sc)
  #cur_seurat <- CreateSeuratObject(
  #  counts = cur_out,
  #  min.cells=3,
  #  min.features=200,
  #  project='RENAL2'
  #)
  #cur_seurat$SampleID <- sample
  #return(cur_seurat)
#})

seurat_list <- lapply(samples[5:8], function(sample){
  cur_filtered_data <- Read10X(paste0(in_data_dir,sample,'/outs/filtered_feature_bc_matrix/'))
  cur_raw_data <- Read10X(paste0(in_data_dir,sample,'/outs/raw_feature_bc_matrix/'))
  cur_sc=SoupChannel(cur_raw_data, cur_filtered_data)
  cur_meta = read.table(paste0(in_data_dir,sample,"/outs/analysis/clustering/graphclust/clusters.csv"),header = T, row.names = 1, sep=",")
  cur_DR = read.table(paste0(in_data_dir,sample,"/outs/analysis/umap/2_components/projection.csv"),header = T, row.names = 1, sep=",")
  cur_meta<-merge(cur_meta,cur_DR, by="row.names")
  row.names(cur_meta)<-cur_meta$Row.names
  cur_meta<-cur_meta[2:length(cur_meta)]
  cur_sc = setClusters(cur_sc, setNames(cur_meta$Cluster, row.names(cur_meta)))
  cur_sc = setDR(cur_sc, cur_meta[colnames(cur_sc$toc), c("UMAP.1", "UMAP.2")])
  cur_sc = autoEstCont(cur_sc)
  igGenes = c("nub","pdm2","ptc")
  head(cur_sc$soupProfile[order(cur_sc$soupProfile$est, decreasing = TRUE), ], n = 20)
  plotMarkerDistribution(cur_sc)
  useToEst = estimateNonExpressingCells(cur_sc, nonExpressedGeneList = list(wt = c("nub")), 
                                        clusters = FALSE)
  plotMarkerMap(cur_sc, geneSet = igGenes, useToEst = useToEst)
  cur_sc = calculateContaminationFraction(cur_sc, list(IG = igGenes), useToEst = useToEst)
  cur_out = adjustCounts(cur_sc)
  cur_seurat <- CreateSeuratObject(
    counts = cur_out,
    min.cells=3,
    min.features=200,
    project='RENAL2'
  )
  cur_seurat$SampleID <- sample
  return(cur_seurat)
})

#rename cell ids for ease of use with Loupe Browser - attach number to end of each sample

seurat_list[[2]] <- RenameCells(object = seurat_list[[2]], new.names = paste0(substring(Cells(x = seurat_list[[2]]),1,16),"-2"))
seurat_list[[3]] <- RenameCells(object = seurat_list[[3]], new.names = paste0(substring(Cells(x = seurat_list[[3]]),1,16),"-3"))
seurat_list[[4]] <- RenameCells(object = seurat_list[[4]], new.names = paste0(substring(Cells(x = seurat_list[[4]]),1,16),"-4"))
seurat_list[[5]] <- RenameCells(object = seurat_list[[5]], new.names = paste0(substring(Cells(x = seurat_list[[5]]),1,16),"-5"))
seurat_list[[6]] <- RenameCells(object = seurat_list[[6]], new.names = paste0(substring(Cells(x = seurat_list[[6]]),1,16),"-6"))
seurat_list[[7]] <- RenameCells(object = seurat_list[[7]], new.names = paste0(substring(Cells(x = seurat_list[[7]]),1,16),"-7"))
seurat_list[[8]] <- RenameCells(object = seurat_list[[8]], new.names = paste0(substring(Cells(x = seurat_list[[8]]),1,16),"-8"))

# merge seurat object
seurat_obj <- merge(x=seurat_list[[1]], y=seurat_list[2:length(seurat_list)])

# clean up
rm(seurat_list)
gc()

obj<-Seurat::NormalizeData(seurat_obj,verbose = FALSE)%>%
FindVariableFeatures(selection.method = "vst", nfeatures = 5000) %>% 
ScaleData(verbose = FALSE)%>%
RunPCA(npcs=35)


expr_mat <- as.matrix(GetAssayData(object = seurat_obj, slot = "counts"))
write.csv(expr_mat,"Justin_Ankita_count_matrix_raw.csv")

obj<-Seurat::NormalizeData(seurat_obj,verbose = FALSE)%>%
  ScaleData(verbose = FALSE,do.center = FALSE)

write.csv(obj@assays$RNA@scale.data,"Justin_Ankita_count_matrix_scaled_normalized.csv")



top10 <- head(VariableFeatures(obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#check elbow plot
ElbowPlot(obj, ndims = 35)

options(repr.plot.height = 10, repr.plot.width = 12)
p1 <- DimPlot(object = obj, reduction = "pca", pt.size = .1, group.by = "SampleID")
p2 <- VlnPlot(object = obj, features = "PC_1", group.by = "SampleID", pt.size = .1)
plot_grid(p1,p2)

VlnPlot(obj, group.by="SampleID", features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size=0, )
#obj <- subset(obj, subset = nFeature_RNA < 2500 & nCount_RNA < 10000)
obj <- subset(obj, subset = nFeature_RNA < 3000)
obj <- subset(obj, subset = nCount_RNA < 10000)
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

obj <- obj %>% 
  RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()

options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(obj, reduction = "umap", group.by = "SampleID", pt.size = .1, split.by = 'SampleID')
DimPlot(obj, reduction = "umap", group.by = "SampleID", pt.size = .1)
FeaturePlot(obj,reduction= "umap", split.by = 'SampleID', features = "Ras85D")

obj <- obj %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()

options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(obj, reduction = "umap", group.by = "SampleID", pt.size = .1, split.by = 'SampleID')
DimPlot(obj, reduction = "umap", group.by = "SampleID", pt.size = .1)
FeaturePlot(obj,reduction= "umap", split.by = 'SampleID', features = "Ras85D")



options(repr.plot.height = 4, repr.plot.width = 6)
DimPlot(obj, reduction = "umap", label = TRUE, pt.size = .1)

#remove clusters
sub_obj <- subset(obj, idents = c(5))
#re-do clustering
DimPlot(sub_obj, reduction = "umap", label = TRUE, pt.size = .1)
sub_obj <- sub_obj %>% 
  RunUMAP(reduction = "harmony", dims = 1:10) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()


options(repr.plot.height = 4, repr.plot.width = 6)
DimPlot(sub_obj, reduction = "umap", label = TRUE, pt.size = .1)

obj <- RenameIdents(object = obj, `0` = "Leg Disc", `1` = "Fat", `2` = "Glia/Epidermis", `3` = "Neurons", `4` = "hemocytes", `5` = "wing disc", `6` = "gut", `7` = "trachea", `8` = "Fat like",`9` = "Eye Disk",`10` = "oenocytes",`11` = "Salivary Glands",`12` = "Gonads",`13` = "Neurons",`14` = "Muscle",`15` = "Malphigian Tubules")
options(repr.plot.height = 4, repr.plot.width = 6)
DimPlot(obj, reduction = "umap", label = TRUE, pt.size = .1)
#sort by positive fold change

#find Markers
obj.markers <- FindAllMarkers(obj, min.pct = 0.25,group.by())
marker_data<-obj.markers %>%
  group_by(cluster)
write.csv(marker_data,"marker_data_JA.csv",row.names = FALSE)


marker_data_top<-obj.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
#marker_data_top
write.csv(marker_data_top,"marker_data_top_FC.csv",row.names = FALSE)

#group by Sample without replicate number
group_metadata=substring(obj$SampleID,1,1)
obj<-AddMetaData(obj, group_metadata, 'Sample_short')
#find markers grouped by SampleID for cluster 0
obj.markers_T0 <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'T', group.by='Sample_short', subset.ident= '0')
obj.markers_T0$zeroes<-integer(nrow(obj.markers_T0))
obj.markers_T1 <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'T', group.by='Sample_short', subset.ident= '1')
obj.markers_T1$zeroes<-integer(nrow(obj.markers_T1))
obj.markers_T2 <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'T', group.by='Sample_short', subset.ident= '2')
obj.markers_T2$zeroes<-integer(nrow(obj.markers_T2))
obj.markers_T3 <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'T', group.by='Sample_short', subset.ident= '3')
obj.markers_T3$zeroes<-integer(nrow(obj.markers_T3))
obj.markers_T1 <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'T', group.by='Sample_short', subset.ident= '1')
obj.markers_T1$zeroes<-integer(nrow(obj.markers_T1))
obj.markers_T1 <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'T', group.by='Sample_short', subset.ident= '1')
obj.markers_T1$zeroes<-integer(nrow(obj.markers_T1))
obj.markers_T1 <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'T', group.by='Sample_short', subset.ident= '1')
obj.markers_T1$zeroes<-integer(nrow(obj.markers_T1))
obj.markers_T1 <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'T', group.by='Sample_short', subset.ident= '1')
obj.markers_T1$zeroes<-integer(nrow(obj.markers_T1))
obj.markers_T1 <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'T', group.by='Sample_short', subset.ident= '1')
obj.markers_T1$zeroes<-integer(nrow(obj.markers_T1))
obj.markers_T1 <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'T', group.by='Sample_short', subset.ident= '1')
obj.markers_T1$zeroes<-integer(nrow(obj.markers_T1))
obj.markers_T1 <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'T', group.by='Sample_short', subset.ident= '1')
obj.markers_T1$zeroes<-integer(nrow(obj.markers_T1))
obj.markers_T1 <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'T', group.by='Sample_short', subset.ident= '1')
obj.markers_T1$zeroes<-integer(nrow(obj.markers_T1))
obj.markers_T1 <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'T', group.by='Sample_short', subset.ident= '1')
obj.markers_T1$zeroes<-integer(nrow(obj.markers_T1))
obj.markers_T1 <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'T', group.by='Sample_short', subset.ident= '1')
obj.markers_T1$zeroes<-integer(nrow(obj.markers_T1))
obj.markers_T1 <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'T', group.by='Sample_short', subset.ident= '1')
obj.markers_T1$zeroes<-integer(nrow(obj.markers_T1))
obj.markers_T1 <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'T', group.by='Sample_short', subset.ident= '1')
obj.markers_T1$zeroes<-integer(nrow(obj.markers_T1))
obj.markers_W0 <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'W', group.by='Sample_short', subset.ident= '0')

obj.markers_T <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'T', group.by='Sample_short')


write.csv(marker_data,"marker_data_JA_g0.csv",row.names = FALSE)



#ElbowPlot(obj, ndims = 30)

DimHeatmap(sub_obj, dims = 1:15, cells = 500, balanced = TRUE, ncol = 4)

umap = cbind("Barcode" = rownames(Embeddings(object = obj, reduction = "umap")), Embeddings(object = obj, reduction = "umap"))
write.table(umap, file="./umap_JA.csv", sep = ",", quote = F, row.names = F, col.names = T)



umap_names = as.data.frame(as.matrix(obj$seurat_clusters))
write.table(umap_names, file="./JA_names.csv", sep = ",", quote = F, row.names = T, col.names = F)

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




#Use SCENIC to do regulatory network analysis
cellInfo <- data.frame(seuratCluster=Idents(obj))
expr_mat <- as.matrix(GetAssayData(object = obj, slot = "counts"))
write.csv(expr_mat,"Expr_Mat_for_JA_SCENIC")
library(SCENIC)
data(defaultDbNames)
dbs <- defaultDbNames[["dmel"]]
scenicOptions <- initializeScenic(org="dmel", dbDir=dbFiles, datasetTitle="Justin_Ankita_SCENIC_analysis", nCores=10)
names(cellInfo)[1]<-"CellType"
cellInfo$nGene <- obj@meta.data$nFeature_RNA
cellInfo$nUMI <- obj@meta.data$nCount_RNA
saveRDS(cellInfo, file="int/cellInfo.Rds")
#scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
genekept<-geneFiltering(expr_mat,scenicOptions = scenicOptions)
exprMat_filtered <- expr_mat[genekept,]
exprMat_filtered_log <- log2(exprMat_filtered+1) 
#running Genie takes a long time - do it as part of a batch script
write.csv(exprMat_filtered,"Expr_Mat_JA_for_GENIE3")
#exportsForArboreto(exprMat_filtered, scenicOptions)
#runGenie3(exprMat_filtered_log, scenicOptions)
runCorrelation(exprMat_filtered, scenicOptions)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
#resource intensive:
#scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
#scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, expr_mat, skipHeatmap=TRUE, skipTsne=TRUE)
scenicOptions <- readRDS("int/scenicOptions.Rds")
