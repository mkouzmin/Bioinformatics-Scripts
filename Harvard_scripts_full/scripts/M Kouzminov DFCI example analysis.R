# Note: this was created to run in R Studio - visualizations are not automatically saved by running code as further visualization was often done in the Loupe Browser - csv files are saved instead for visualizations to be transfered into loupe
library(Seurat)
library(tidyverse)
library(cowplot)
library(harmony)
library(SoupX)
library(remotes)


in_data_dir <- "/n/groups/flyrnai/mikhail/FC_07486_pedro/run_count_all/"


#for SCENIC
dbFiles <- "./Databases/"


set.seed(123)
samples <- dir(in_data_dir)

# create seurat object - no SoupX:
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

#run SoupX - as cell free mRNA fraction is suspected to be high in this sample.
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

#run from here for both versions
#rename cell ids for ease of use with Loupe Browser - attach number to end of each sample
seurat_list[[2]] <- RenameCells(object = seurat_list[[2]], new.names = paste0(substring(Cells(x = seurat_list[[2]]),1,16),"-2"))
seurat_list[[3]] <- RenameCells(object = seurat_list[[3]], new.names = paste0(substring(Cells(x = seurat_list[[3]]),1,16),"-3"))
seurat_list[[4]] <- RenameCells(object = seurat_list[[4]], new.names = paste0(substring(Cells(x = seurat_list[[4]]),1,16),"-4"))
seurat_list[[5]] <- RenameCells(object = seurat_list[[5]], new.names = paste0(substring(Cells(x = seurat_list[[5]]),1,16),"-5"))
seurat_list[[6]] <- RenameCells(object = seurat_list[[6]], new.names = paste0(substring(Cells(x = seurat_list[[6]]),1,16),"-6"))


# merge seurat object
seurat_obj <- merge(x=seurat_list[[1]], y=seurat_list[2:length(seurat_list)])

#following code runs soupX manually for one sample for exploratative analysis and comparison of contamination fraction - already done above automatically
#cur_filtered_data <- Read10X(paste0(in_data_dir,samples[1],'/outs/filtered_feature_bc_matrix/'))
#cur_raw_data <- Read10X(paste0(in_data_dir,samples[1],'/outs/raw_feature_bc_matrix/'))
#sc1=SoupChannel(cur_raw_data, cur_filtered_data)
#meta = read.table(paste0(in_data_dir,samples[1],"/outs/analysis/clustering/graphclust/clusters.csv"),header = T, row.names = 1, sep=",")
#DR = read.table(paste0(in_data_dir,samples[1],"/outs/analysis/umap/2_components/projection.csv"),header = T, row.names = 1, sep=",")
#meta<-merge(meta,DR, by="row.names")
#row.names(meta)<-meta$Row.names
#meta<-meta[2:length(meta)]

#sc1 = setClusters(sc1, setNames(meta$Cluster, row.names(meta)))
#sc1 = setDR(sc1, meta[colnames(sc1$toc), c("UMAP.1", "UMAP.2")])
#library(ggplot2)
#dd = meta[colnames(sc1$toc), ]
#mids = aggregate(cbind(UMAP.1, UMAP.2) ~ Cluster, data = dd, FUN = mean)
#gg = ggplot(dd, aes(UMAP.1, UMAP.2)) + geom_point(aes(colour = Cluster), siz= 0.2) + 
#  geom_label(data = mids, aes(label = Cluster)) + ggtitle("default_clusters") + 
#  guides(colour = guide_legend(override.aes = list(size = 1)))
#plot(gg)
#sc1 = setContaminationFraction(sc1, 0.2)
## can auto set contmination fraction based on subset of genes
##nonExpressedGeneList = list(HB = c("HBB", "HBA2"), IG = c("IGKC"))
##sc1 = autoEstCont(sc1,priorRhoStdDev = 0.3)
#sc1 = autoEstCont(sc1)
##manual version
##check most highly expressed genes in background
##head(sc1$soupProfile[order(sc1$soupProfile$est, decreasing = TRUE), ], n = 20)
#plotMarkerDistribution(sc1)
#out1 = adjustCounts(sc1)
#cntSoggy = rowSums(sc1$toc > 0)
#cntStrained = rowSums(out1 > 0)
#mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
#mostZeroed
#tail(sort(rowSums(sc1$toc > out1)/rowSums(sc1$toc > 0)), n = 20)

#plotChangeMap(sc1, out1, "Obp19d")


#possibly run for other samples later - second seurat_list function does the above automatically

# clean up
rm(seurat_list)
gc()

obj<-Seurat::NormalizeData(seurat_obj,verbose = FALSE)%>%
FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
ScaleData(verbose = FALSE)%>%
RunPCA(npcs=30)

expr_mat <- as.matrix(GetAssayData(object = seurat_obj, slot = "counts"))
write.csv(expr_mat,"Pedro_count_matrix_raw_no_SoupX.csv")

obj<-Seurat::NormalizeData(seurat_obj,verbose = FALSE)%>%
  ScaleData(verbose = FALSE,do.center = FALSE)

write.csv(obj@assays$RNA@scale.data,"Pedro_count_matrix_scaled_normalized.csv")


top10 <- head(VariableFeatures(obj), 10)


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 

#check elbow plot
ElbowPlot(obj, ndims = 30)

#make pca plots
options(repr.plot.height = 10, repr.plot.width = 12)
p1 <- DimPlot(object = obj, reduction = "pca", pt.size = .1, group.by = "SampleID")
p2 <- VlnPlot(object = obj, features = "PC_1", group.by = "SampleID", pt.size = .1)
plot_grid(p1,p2)

#subset by counts & features
VlnPlot(seurat_obj, group.by="SampleID", features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size=0, )
obj <- subset(obj, subset = nFeature_RNA < 800)
obj <- subset(obj, subset = nCount_RNA < 1500)
VlnPlot(obj, group.by="SampleID", features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size=0, )

#View(obj@assays$RNA@counts)
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

#create umap without harmony
obj <- obj %>% 
  RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()

#make umap visualizations
options(repr.plot.height = 4, repr.plot.width = 11)
DimPlot(obj, reduction = "umap", group.by = "SampleID", pt.size = .1, split.by = 'SampleID')
DimPlot(obj, reduction = "umap", group.by = "SampleID", pt.size = .1)
DimPlot(obj, reduction = "umap", pt.size = .1)
FeaturePlot(obj,reduction= "umap", split.by = 'SampleID', features = "Hsp26")


#run Seurat with Harmony
obj <- obj %>% 
  RunUMAP(reduction = "harmony", dims = 1:25) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:25) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()

options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(obj, reduction = "umap", group.by = "SampleID", pt.size = .1, split.by = 'SampleID')
DimPlot(obj, reduction = "umap", group.by = "SampleID", pt.size = .1)
DimPlot(obj, reduction = "umap", pt.size = .1)


options(repr.plot.height = 4, repr.plot.width = 6)
DimPlot(obj, reduction = "umap", label = TRUE, pt.size = .1)



#find Markers
obj.markers <- FindAllMarkers(obj, min.pct = 0.25)
marker_data<-obj.markers %>%
  group_by(cluster)
write.csv(marker_data,"marker_data_SoupX_pedro.csv",row.names = FALSE)


marker_data_top<-obj.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
#marker_data_top
write.csv(marker_data_top,"marker_data_top_FC.csv",row.names = FALSE)


#ElbowPlot(obj, ndims = 30)

DimHeatmap(obj, dims = 1:15, cells = 500, balanced = TRUE, ncol = 4)

obj.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
marker_data_top<-obj.markers %>%
  group_by()
DoHeatmap(obj, features = marker_data_top$gene) + NoLegend()

umap = cbind("Barcode" = rownames(Embeddings(object = obj, reduction = "umap")), Embeddings(object = obj, reduction = "umap"))
write.table(umap, file="./umap_pedro.csv", sep = ",", quote = F, row.names = F, col.names = T)

umap_names = as.data.frame(as.matrix(obj$RNA_snn_res.0.4))
write.table(umap_names, file="./umap_pedro_names.csv", sep = ",", quote = F, row.names = T, col.names = F)

#combine samples by type for finding DE Genes
group_metadata=obj$SampleID
group_metadata<-str_replace_all(group_metadata,c("pedro1" = "wt","pedro2"="wt","pedro3"="foxo","pedro4"="foxo","pedro5"="reptor","pedro6"="reptor"))
obj<-AddMetaData(obj, group_metadata, 'sample_type')


umap_sample_type = as.data.frame(as.matrix(obj$sample_type))
write.table(umap_sample_type, file="./umap_pedro_sample_type.csv", sep = ",", quote = F, row.names = T, col.names = F)

group_metadata=obj$SampleID
group_metadata<-str_replace_all(group_metadata,c("pedro1" = "wt","pedro2"="wt","pedro3"="foxo","pedro4"="foxo","pedro5"="reptor","pedro6"="reptor"))
obj<-AddMetaData(obj, group_metadata, 'sample_type')

obj.markers_foxo_by_clust_list <- lapply(0:(length(unique(Idents(obj)))-1), function(id){
  obj.markers_foxo_cur <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'foxo', ident.2 = "wt", group.by='sample_type',subset.ident = id)
  obj.markers_foxo_cur$clusterID = id
  obj.markers_foxo_cur <- cbind(rownames(obj.markers_foxo_cur), data.frame(obj.markers_foxo_cur, row.names=NULL))
  return(obj.markers_foxo_cur)
})
obj.markers_foxo_by_clust <- bind_rows(obj.markers_foxo_by_clust_list)
write.csv(obj.markers_foxo_by_clust,"marker_data_foxo_vs_W_by_clust.csv",row.names = FALSE)


obj.markers_reptor_by_clust_list <- lapply(0:(length(unique(Idents(obj)))-1), function(id){
  obj.markers_reptor_cur <- FindMarkers(obj, min.pct = 0.25,ident.1 = 'reptor', ident.2 = "wt", group.by='sample_type',subset.ident = id)
  obj.markers_reptor_cur$clusterID = id
  obj.markers_reptor_cur <- cbind(rownames(obj.markers_reptor_cur), data.frame(obj.markers_reptor_cur, row.names=NULL))
  return(obj.markers_reptor_cur)
})
obj.markers_reptor_by_clust <- bind_rows(obj.markers_reptor_by_clust_list)
write.csv(obj.markers_reptor_by_clust,"marker_data_reptor_vs_W_by_clust.csv",row.names = FALSE)

obj.markers_foxo_by_clust %>%
  group_by(clusterID) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10_foxo

obj.markers_reptor_by_clust %>%
  group_by(clusterID) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10_reptor

obj.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10_markers

#load selected genes
pedro_genes_0<-read.csv("genes_for_hetamap_pedro_data.csv")
pedro_genes_2<-read.csv("genes_for_hetamap_pedro_data_2.csv")

#scale all data - not just sub-cluster
obj <- ScaleData(object = obj, features = rownames(obj))

#Heatmap for selected genes
DoHeatmap(obj,cells = WhichCells(obj, idents = 0), features = pedro_genes_0$Gene, group.by = "sample_type") + NoLegend()
DoHeatmap(obj,cells = WhichCells(obj, idents = 2), features = pedro_genes_2$Gene, group.by = "sample_type") + NoLegend()

#If not likely to run each heatmap individually, implement for loop - this format is useful for running and viewing each heatmap separately before importing into loupe browser - may not be best practice for readability
#cluster 0
DoHeatmap(obj,cells = WhichCells(obj, idents = 0), features = c(top10_foxo$`rownames(obj.markers_foxo_cur)`[1:10], top10_reptor$`rownames(obj.markers_reptor_cur)`[1:10]), group.by = "sample_type") + NoLegend()
#cluster 1
DoHeatmap(obj,cells = WhichCells(obj, idents = 1), features = c(top10_foxo$`rownames(obj.markers_foxo_cur)`[11:20], top10_reptor$`rownames(obj.markers_reptor_cur)`[11:20]), group.by = "sample_type") + NoLegend()
#cluster 2
DoHeatmap(obj,cells = WhichCells(obj, idents = 2), features = c(top10_foxo$`rownames(obj.markers_foxo_cur)`[21:30], top10_reptor$`rownames(obj.markers_reptor_cur)`[21:30]), group.by = "sample_type") + NoLegend()
#cluster 3
DoHeatmap(obj,cells = WhichCells(obj, idents = 3), features = c(top10_foxo$`rownames(obj.markers_foxo_cur)`[31:40], top10_reptor$`rownames(obj.markers_reptor_cur)`[31:40]), group.by = "sample_type") + NoLegend()
#cluster 4
DoHeatmap(obj,cells = WhichCells(obj, idents = 4), features = c(top10_foxo$`rownames(obj.markers_foxo_cur)`[41:50], top10_reptor$`rownames(obj.markers_reptor_cur)`[41:50]), group.by = "sample_type") + NoLegend()
#cluster 5
DoHeatmap(obj,cells = WhichCells(obj, idents = 5), features = c(top10_foxo$`rownames(obj.markers_foxo_cur)`[51:60], top10_reptor$`rownames(obj.markers_reptor_cur)`[51:60]), group.by = "sample_type") + NoLegend()
#cluster 6
DoHeatmap(obj,cells = WhichCells(obj, idents = 6), features = c(top10_foxo$`rownames(obj.markers_foxo_cur)`[61:70], top10_reptor$`rownames(obj.markers_reptor_cur)`[61:70]), group.by = "sample_type") + NoLegend()
#cluster 7
DoHeatmap(obj,cells = WhichCells(obj, idents = 7), features = c(top10_foxo$`rownames(obj.markers_foxo_cur)`[71:80], top10_reptor$`rownames(obj.markers_reptor_cur)`[71:80]), group.by = "sample_type") + NoLegend()
#cluster 8
DoHeatmap(obj,cells = WhichCells(obj, idents = 8), features = c(top10_foxo$`rownames(obj.markers_foxo_cur)`[81:90], top10_reptor$`rownames(obj.markers_reptor_cur)`[81:90]), group.by = "sample_type") + NoLegend()
#cluster 9
DoHeatmap(obj,cells = WhichCells(obj, idents = 9), features = c(top10_foxo$`rownames(obj.markers_foxo_cur)`[91:100], top10_reptor$`rownames(obj.markers_reptor_cur)`[91:100]), group.by = "sample_type") + NoLegend()
#cluster 10
DoHeatmap(obj,cells = WhichCells(obj, idents = 10), features = c(top10_foxo$`rownames(obj.markers_foxo_cur)`[101:110], top10_reptor$`rownames(obj.markers_reptor_cur)`[101:110]), group.by = "sample_type") + NoLegend()
#cluster 11
DoHeatmap(obj,cells = WhichCells(obj, idents = 11), features = c(top10_foxo$`rownames(obj.markers_foxo_cur)`[111:120], top10_reptor$`rownames(obj.markers_reptor_cur)`[111:120]), group.by = "sample_type") + NoLegend()
#cluster 12
DoHeatmap(obj,cells = WhichCells(obj, idents = 12), features = c(top10_foxo$`rownames(obj.markers_foxo_cur)`[121:130], top10_reptor$`rownames(obj.markers_reptor_cur)`[121:130]), group.by = "sample_type") + NoLegend()
#cluster 13
DoHeatmap(obj,cells = WhichCells(obj, idents = 13), features = c(top10_foxo$`rownames(obj.markers_foxo_cur)`[131:140], top10_reptor$`rownames(obj.markers_reptor_cur)`[131:140]), group.by = "sample_type") + NoLegend()

#heatmap for marker genes
#cluster 0
DoHeatmap(obj,cells = WhichCells(obj, idents = 0), features = top10_markers$gene[1:10], group.by = "sample_type") + NoLegend()
#cluster 1
DoHeatmap(obj,cells = WhichCells(obj, idents = 1), features = top10_markers$gene[11:20], group.by = "sample_type") + NoLegend()
#cluster 2
DoHeatmap(obj,cells = WhichCells(obj, idents = 2), features = top10_markers$gene[21:30], group.by = "sample_type") + NoLegend()
#cluster 3
DoHeatmap(obj,cells = WhichCells(obj, idents = 3), features = top10_markers$gene[31:40], group.by = "sample_type") + NoLegend()
#cluster 4
DoHeatmap(obj,cells = WhichCells(obj, idents = 4), features = top10_markers$gene[41:50], group.by = "sample_type") + NoLegend()
#cluster 5
DoHeatmap(obj,cells = WhichCells(obj, idents = 5), features = top10_markers$gene[51:60], group.by = "sample_type") + NoLegend()
#cluster 6
DoHeatmap(obj,cells = WhichCells(obj, idents = 6), features = top10_markers$gene[61:70], group.by = "sample_type") + NoLegend()
#cluster 7
DoHeatmap(obj,cells = WhichCells(obj, idents = 7), features = top10_markers$gene[71:80], group.by = "sample_type") + NoLegend()
#cluster 8
DoHeatmap(obj,cells = WhichCells(obj, idents = 8), features = top10_markers$gene[81:90], group.by = "sample_type") + NoLegend()
#cluster 9
DoHeatmap(obj,cells = WhichCells(obj, idents = 9), features = top10_markers$gene[91:100], group.by = "sample_type") + NoLegend()
#cluster 10
DoHeatmap(obj,cells = WhichCells(obj, idents = 10), features = top10_markers$gene[101:110], group.by = "sample_type") + NoLegend()
#cluster 11
DoHeatmap(obj,cells = WhichCells(obj, idents = 11), features = top10_markers$gene[111:120], group.by = "sample_type") + NoLegend()
#cluster 12
DoHeatmap(obj,cells = WhichCells(obj, idents = 12), features = top10_markers$gene[121:130], group.by = "sample_type") + NoLegend()
#cluster 13
DoHeatmap(obj,cells = WhichCells(obj, idents = 13), features = top10_markers$gene[131:140], group.by = "sample_type") + NoLegend()


# run monocle3
library(SeuratWrappers)
library(monocle3)
obj.cds <- as.cell_data_set(obj)
obj.cds <- cluster_cells(cds = obj.cds, reduction_method = "UMAP")
obj.cds <- learn_graph(obj.cds, use_partition = TRUE)
obj.cds <- order_cells(obj.cds, reduction_method = "UMAP")
plot_cells(
  cds = obj.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

#rename cluster identities for SCENIC
#clusters identified separately - with some input from biologists - read paper report for specific marker genes.
rename_obj <- RenameIdents(object = obj, `0` = "Indirect flight muscle", `1` = "Adipose tissue", `2` = "Leg and somatic muscle", `3` = "ventral nervous system", `4` = "unknown1", `5` = "Adult epidermis", `6` = "unknown2", `7` = "Hemocyte", `8` = "Peripheral nervous system", `9` = "Trachea", `10` = "unknown3", `11` = "Glial cell", `12` = "unknown4", `13` = "Adult oenocyte" )
saveRDS(rename_obj, file="Scenic_obj_pedro.Rds")


#GENIE 3 & SCENIC for gene network visualization
#BiocManager::install(c("AUCell", "RcisTarget"))
#BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost
#BiocManager::install(c("zoo", "mixtools", "rbokeh"))
# For various visualizations and perform t-SNEs:
#BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
# To support paralell execution in tSNE (not available in Windows):
#BiocManager::install(c("doMC", "doRNG"))
# To export/visualize in http://scope.aertslab.org
#if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
#install.packages("hdf5r")
#devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
#devtools::install_github("aertslab/SCENIC") 
packageVersion("SCENIC")

#create cell Info matrix - saves to file
cellInfo <- data.frame(seuratCluster=Idents(rename_obj))
names(cellInfo)[1]<-"CellType"
rownames(cellInfo)<-gsub("-",".",rownames(cellInfo))
#differentiate cells based on sample
cellInfo$CellType <- paste(cellInfo$CellType, " ", rename_obj$SampleID)
cellInfo$nGene <- rename_obj@meta.data$nFeature_RNA
cellInfo$nUMI <- rename_obj@meta.data$nCount_RNA
saveRDS(cellInfo, file="int/cellInfo.Rds")

#Use SCENIC to do regulatory network analysis
expr_mat <- as.matrix(GetAssayData(object = rename_obj, slot = "counts"))
write.csv(expr_mat,"Expr_Mat_for_SCENIC")
library(SCENIC)
library(AUCell)
scenicOptions <- initializeScenic(org="dmel", dbDir=dbFiles, nCores=10)
#if dbFiles not found, load databases

#run from here to read prepared RDS
cellInfo <- readRDS("int/cellInfo.Rds")
#scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
genekept<-geneFiltering(expr_mat,scenicOptions = scenicOptions)
exprMat_filtered <- expr_mat[genekept,]
exprMat_log <- log2(expr_mat+1)
#running Genie takes a long time - do it as part of a batch script -
#see GENIE3_pedro.R and run_GENIE_pedro.sh


write.csv(exprMat_filtered,"Expr_Mat_pedro_for_GENIE3")
#NOTICE: run run_GENIE3_pedro.sh here
exportsForArboreto(exprMat_filtered, scenicOptions, dir = "int")
#runGenie3(exprMat_filtered, scenicOptions)
runCorrelation(exprMat_filtered, scenicOptions)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
#resource intensive:
#NOTICE: Run SCENIC_find_regulons_pedro.sh and SCENIC_score_cells_pedro.sh here

#scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
#scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, expr_mat, skipHeatmap=TRUE, skipTsne=TRUE)
scenicOptions <- readRDS("int/scenicOptions.Rds")

nPcs <- c(15)
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)
par(mfrow=c(1,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)



#set pca dims
scenicOptions@settings$defaultTsne$dims <- 15

motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="achi"]
viewMotifs(motifEnrichment_selfMotifs_wGenes) 

regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="achi" & highConfAnnot==TRUE]
viewMotifs(tableSubset) 
viewMotifs(regulonTargetsInfo[highConfAnnot==TRUE])

#write regulons to file
regulons <- loadInt(scenicOptions, "regulons")
capture.output(regulons, file = "SCENIC_regulons.txt")

#make heatmap
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, column_order = order(as.numeric(gsub("column", "", colnames(regulonActivity_byCellType_Scaled)))), name="Regulon activity")

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"], )
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)
