
library(Seurat)
library(tidyverse)
library(cowplot)
library(harmony)
library(SoupX)
library(remotes)

dbFiles <- "./Databases/"

#create cell Info matrix - saves to file - this was run in Pedro_Tumor_full_body_analysis.R
cellInfo <- data.frame(seuratCluster=Idents(rename_obj))
names(cellInfo)[1]<-"CellType"
rownames(cellInfo)<-gsub("-",".",rownames(cellInfo))
#differentiate cells based on sample
#cellInfo$CellType <- paste(cellInfo$CellType, " ", rename_obj$SampleID)
cellInfo$nGene <- rename_obj@meta.data$nFeature_RNA
cellInfo$nUMI <- rename_obj@meta.data$nCount_RNA
saveRDS(cellInfo, file="int/cellInfo.Rds")


#SCENIC expects cell Info to be in int folder

rename_obj <- readRDS(file="Scenic_obj_pedro.Rds")



#Use SCENIC to do regulatory network analysis
expr_mat <- as.matrix(GetAssayData(object = rename_obj, slot = "counts"))
write.csv(expr_mat,"Expr_Mat_for_SCENIC.csv")
library(SCENIC)
library(AUCell)
scenicOptions <- initializeScenic(org="dmel", dbDir=dbFiles, nCores=8)
#if dbFiles not found, load databases

#run from here to read prepared RDS
cellInfo <- readRDS("int/cellInfo.Rds")
#scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
genekept<-geneFiltering(expr_mat,scenicOptions = scenicOptions, minCountsPerGene=3*.01*ncol(expr_mat),
                        minSamples=ncol(expr_mat)*.01)
exprMat_filtered <- expr_mat[genekept,]
exprMat_log <- log2(expr_mat+1)
#running Genie takes a long time - do it as part of a batch script -
#see GENIE3_pedro.R and run_GENIE_pedro.sh


write.csv(exprMat_filtered,"Expr_Mat_pedro_for_GENIE3")
#NOTICE: run run_GENIE3_pedro.sh here or uncomment lines
exportsForArboreto(exprMat_filtered, scenicOptions, dir = "int")
#runGenie3(exprMat_filtered, scenicOptions)
runCorrelation(exprMat_filtered, scenicOptions)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
#resource intensive:
#NOTICE: Run SCENIC_find_regulons_pedro.sh and SCENIC_score_cells_pedro.sh here

#scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
#scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, expr_mat, skipHeatmap=TRUE, skipTsne=TRUE)
# Might be error here - may need to re-run with expr_mat_log
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

scenicOptions <- readRDS("int/scenicOptions.Rds")

nPcs <- c(15)
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
#plot Tsne to compare normal and only high confidence gene Tsnes
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

#create table of top regulator genes
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)

#create rss plot
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"], )
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)
