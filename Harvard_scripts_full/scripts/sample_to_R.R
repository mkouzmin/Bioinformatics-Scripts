
library(Seurat)
library(tidyverse)
library(cowplot)
library(harmony)

in_data_dir <- "/n/groups/flyrnai/mikhail/Renal_yki_Jun/cellranger_count_results_fixed/"


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
    project='RENAL'
  )
  cur_seurat$SampleID <- sample
  return(cur_seurat)
})

seurat_list[[2]] <- RenameCells(object = seurat_list[[2]], new.names = paste0(substring(Cells(x = seurat_list[[2]]),1,16),"-2"))

# merge seurat object
seurat_obj <- merge(x=seurat_list[[1]], y=seurat_list[2:length(seurat_list)])

# clean up
rm(seurat_list)
gc()

obj<-Seurat::NormalizeData(seurat_obj,verbose = FALSE)%>%
FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
ScaleData(verbose = FALSE)%>%
RunPCA(npcs=15)

expr_mat <- as.matrix(GetAssayData(object = seurat_obj, slot = "counts"))
write.csv(expr_mat,"Jun_count_matrix_raw.csv")

obj<-Seurat::NormalizeData(seurat_obj,verbose = FALSE)%>%
  ScaleData(verbose = FALSE,do.center = FALSE)

write.csv(obj@assays$RNA@scale.data,"Jun_count_matrix_scaled_normalized.csv")


ElbowPlot(obj, ndims = 15)

top10 <- head(VariableFeatures(obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 
#check elbow plot
#ElbowPlot(obj)

options(repr.plot.height = 10, repr.plot.width = 12)
p1 <- DimPlot(object = obj, reduction = "pca", pt.size = .1, group.by = "SampleID")
#p2 <- VlnPlot(object = obj, features = "PC_1", group.by = "SampleID", pt.size = .1)
plot_grid(p1)

VlnPlot(seurat_obj, group.by="SampleID", features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size=0, )


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
plot_grid(p1)

obj <- obj %>% 
  RunUMAP(dims = 1:15) %>% 
  FindNeighbors(dims = 1:15) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()


obj <- obj %>% 
  RunUMAP(reduction = "harmony", dims = 1:15) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:15) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()

options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(obj, reduction = "umap", group.by = "SampleID", pt.size = .1, split.by = 'SampleID')
DimPlot(obj, reduction = "umap", group.by = "SampleID", pt.size = .1)


options(repr.plot.height = 4, repr.plot.width = 6)
DimPlot(obj, reduction = "umap", label = TRUE, pt.size = .1)

#remove clusters
sub_obj <- subset(obj, idents = c(0,1,6,10,11,12,14), invert = TRUE)
#re-do clustering
DimPlot(sub_obj, reduction = "umap", label = TRUE, pt.size = .1)
sub_obj <- sub_obj %>% 
  RunUMAP(reduction = "harmony", dims = 1:10) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()


sub_obj <- RenameIdents(object = sub_obj, `0` = "main segment principal cells", `1` = "upper ureter principal cells", `2` = "renal stem cells", `3` = "lower tubule principal cells", `4` = "initial and transitional principal cells", `5` = "lower segment principal cells", `6` = "stellate cells", `7` = "lower ureter principal cells", `8` = "main segment principal cells")
options(repr.plot.height = 4, repr.plot.width = 8)
DimPlot(sub_obj, reduction = "umap", label = TRUE, pt.size = .1)
#sort by positive fold change

#find Markers
obj.markers <- FindAllMarkers(obj, min.pct = 0.25)
marker_data<-obj.markers %>%
  group_by(cluster)
write.csv(marker_data,"marker_data_fixed_2.csv",row.names = FALSE)

sub_obj.markers <- FindAllMarkers(sub_obj, min.pct = 0.25)
marker_data<-sub_obj.markers %>%
  group_by(cluster)
write.csv(marker_data,"marker_data_final2.csv",row.names = FALSE)

marker_data_top<-obj.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
#marker_data_top
write.csv(marker_data_top,"marker_data_top_FC.csv",row.names = FALSE)


ElbowPlot(obj, ndims = 15)

DimHeatmap(sub_obj, dims = 1:15, cells = 500, balanced = TRUE, ncol = 4)

umap = cbind("Barcode" = rownames(Embeddings(object = obj, reduction = "umap")), Embeddings(object = obj, reduction = "umap"))
write.table(umap, file="./umap_renal.csv", sep = ",", quote = F, row.names = F, col.names = T)


umap_names = as.data.frame(as.matrix(obj$RNA_snn_res.0.4))
write.table(umap_names, file="./umap_renal_names.csv", sep = ",", quote = F, row.names = T, col.names = F)

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
expr_mat <- as.matrix(GetAssayData(object = sub_obj, slot = "counts"))
write.csv(expr_mat,"Expr_Mat_for_SCENIC")
library(SCENIC)
library(AUCell)
scenicOptions <- initializeScenic(org="dmel", dbDir=dbFiles, nCores=10)
#create cell Info matrix
cellInfo <- data.frame(seuratCluster=Idents(sub_obj))
names(cellInfo)[1]<-"CellType"
rownames(cellInfo)<-gsub("-",".",rownames(cellInfo))
#differentiate cells based on sample
cellInfo$CellType <- paste(cellInfo$CellType, " ", sub_obj$SampleID)
cellInfo$nGene <- sub_obj@meta.data$nFeature_RNA
cellInfo$nUMI <- sub_obj@meta.data$nCount_RNA
saveRDS(cellInfo, file="int/cellInfo.Rds")
#scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
genekept<-geneFiltering(expr_mat,scenicOptions = scenicOptions)
exprMat_filtered <- expr_mat[genekept,]
exprMat_log <- log2(expr_mat+1)
#running Genie takes a long time - do it as part of a batch script
#write.csv(exprMat_filtered,"Expr_Mat_for_GENIE3")
#exportsForArboreto(exprMat_filtered, scenicOptions)
#runGenie3(exprMat_filtered, scenicOptions)
runCorrelation(exprMat_filtered, scenicOptions)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
#resource intensive:
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
