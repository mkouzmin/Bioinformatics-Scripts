
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager",Dependencies=TRUE)
#BiocManager::install("locfit")
#BiocManager::install(version = "3.12")
#BiocManager::install("DESeq2")
library("BiocManager")
library("DESeq2")
library("dplyr")
library("ggplot2")
library(data.table)

dt = data.frame(fread("./gene_rpkm_matrix_fb_2022_02.tsv.gz"), row.names=1)
dt = data.frame(fread("./related_dataset_Mary-Lee_late_point_FPKM.txt"), row.names=1)


#Fly Atlas Counts Matrix Prep
dt = data.frame(fread("./FlyAtlas_Features.txt"), row.names=1)
flyAtlas_names = data.frame(fread("./flyatlas2_bam2tissue_folder1.txt", header = FALSE), row.names=1)
colnames(dt) <- c(colnames(dt)[1:5],flyAtlas_names[colnames(dt)[6:47],])

#subset samples from flybase
#128-169 FlyAtlas
#101-127 modEncode Cell Line
#17-75 - tissue developmental
dt<-dt[,101:127]
flysample<-names(dt)
zeroes<-integer(length(flysample))
flysample<-cbind(flysample,zeroes,flysample,zeroes)
rownames(flysample)<-flysample[,1]
flysample<-flysample[,2:4]
colnames(flysample)<-c("SampleID","CellLine","Triplicate")
flysample<-flysample[6:47,]
#read column Data of project data
colData <-as.matrix(read.csv("/n/groups/flyrnai/mikhail/bulk_rna_seq_1/BAM_run1/coldata.txt",sep="\t",row.names=1, header = TRUE))
#append flybase columns
colData <- rbind(colData,flysample)
#read count files
cts1 <- as.matrix(read.csv("/n/groups/flyrnai/mikhail/bulk_rna_seq_1/BAM_run1/featurecounts_1.txt",sep="\t",row.names="Geneid",skip = 1, header = TRUE))
cts2 <- as.matrix(read.csv("/n/groups/flyrnai/mikhail/bulk_rna_seq_1/BAM_run2/featurecounts_2.txt",sep="\t",row.names="Geneid",skip = 1, header = TRUE))
cts3 <- as.matrix(read.csv("/n/groups/flyrnai/mikhail/bulk_rna_seq_1/BAM_run3/featurecounts_3.txt",sep="\t",row.names="Geneid",skip = 1, header = TRUE))
tf <- data.frame(fread("./Fly_TF.txt"))

cts1 <- as.data.frame(cts1)
cts2 <- as.data.frame(cts2)
cts3 <- as.data.frame(cts3)
cts_only1 <- cts1[,6:41]
cts_only2 <- cts2[,6:41]
cts_only3 <- cts3[,6:41]
#convert to numeric
cts_only1 <-as.data.frame(data.matrix(cts_only1))
cts_only2 <-as.data.frame(data.matrix(cts_only2))
cts_only3 <-as.data.frame(data.matrix(cts_only3))

# create a new variable from the rownames
cts_only1$rn <- rownames(cts_only1)
cts_only2$rn <- rownames(cts_only2)
cts_only3$rn <- rownames(cts_only3)


# bind the two dataframes together by row and aggregate
cts_all <- aggregate(. ~ rn, rbind(cts_only1,cts_only2), sum)
cts_all <- aggregate(. ~ rn, rbind(cts_all,cts_only3), sum)

# assign the rownames again
rownames(cts_all) <- cts_all$rn

# get rid of the 'rn' column
cts_all <- cts_all[, -1]


#separate sample count
#colnames(cts_only1) <- paste(colnames(cts_only1),"1",sep="_")
#colnames(cts_only2) <- paste(colnames(cts_only2),"2",sep="_")
#colnames(cts_only3) <- paste(colnames(cts_only3),"3",sep="_")
#cts_all <- merge(cts_only1,cts_only2,by = 'row.names', all = TRUE)
#rownames(cts_all) <- cts_all[,1]
#cts_all[,1] <- NULL
#cts_all <- merge(cts_all,cts_only3,by = 'row.names', all = TRUE)
#rownames(cts_all) <- cts_all[,1]
#cts_all[,1] <- NULL


#join with flybase data
dt1<-dt[,6:47]
cts_all <- merge(cts_all,dt1,by = 'row.names')
rownames(cts_all) <- cts_all[,1]
cts_all[,1] <- NULL
cts_all[is.na(cts_all)] = 0

#use only flyTF
cts_all<-cts_all[tf[,1],]
dt<-dt[tf[,1],]
#write counts matrix
#write.table(cts_all, file="./counts_matrix__aggregate.csv", sep = ",", quote = F, row.names = T, col.names = T)


#get rows containg FBgn in gene name
cts_subset <- cts_all[grep("FBgn", rownames(cts_all)), ]
cts_length <-dt[rownames(cts_subset), ]$Length
cts_matrix <-data.matrix(cts_subset)
DESeqobj <- DESeqDataSetFromMatrix(cts_matrix,colData,design = ~CellLine)
#median of ratios normalization
DESeqobj <- estimateSizeFactors(DESeqobj)
sizeFactors(DESeqobj)
normalized_counts <- counts(DESeqobj, normalized=TRUE)
#write.csv(normalized_counts,"FlyAtlas_normalized_counts")
mcols(DESeqobj)$basepairs<-as.numeric(cts_length)
fpkm_table <- fpkm(DESeqobj)
write.csv(fpkm_table,"Fly_Atlas_and_New_Cell_Lines_FPKM.csv")

keep <- rowSums(counts(DESeqobj) >= 10) >= 3
DESeqobj <- DESeqobj[keep,]
fpkm_table <- fpkm_table[keep,]

vsd <- vst(DESeqobj, blind = FALSE)
 head(assay(vsd), 3)
rld <- rlog(DESeqobj, blind = FALSE)

df <- bind_rows(
  as_data_frame(log2(counts(DESeqobj, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  
lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation) 


library("pheatmap")
library("RColorBrewer")


#join with flybase data for heatmap
colnames(fpkm_table)<-DESeqobj$CellLine
fpkm_table <- merge(fpkm_table,dt,by = 'row.names')#ignore column name warning
rownames(fpkm_table) <- fpkm_table[,1]
fpkm_table[,1] <- NULL
fpkm_table[is.na(fpkm_table)] = 0
sampleDists<- dist(t(fpkm_table))
sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

write.csv(sampleDistMatrix,"Fly_Atlas_TF_only.csv")

#distance heatmap for normalized counts
assay(vsd)
sampleDists <- dist(t(assay(vsd)))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$CellLine, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

plotPCA(vsd, intgroup = c("CellLine")) + geom_point(size = 0.1)
colData(vsd)
write.table(normalized_counts, file="./normalized_counts.txt", sep="\t", quote=F)
write.table(fpkm_table, file="./fpkm_values.txt", sep="\t", quote=F)






#get Differentially expressed genes in Sample

colData1 <-as.matrix(read.csv("./sample_backgrounds.txt",sep="\t",row.names=1, header = TRUE))
DESeqobj <- DESeqDataSetFromMatrix(cts_matrix,colData,design = ~CellLine)
#create factors where all samples not of a specific Background group are 0 - to compare sample against all other background groups
DESeqobj$group1 <- factor(replace(DESeqobj$Sample.Group,DESeqobj$Background.Group != 1,))
DESeqobj$group1 <- factor(replace(replace(DESeqobj$Triplicate,DESeqobj$Triplicate != 0,"CL"),DESeqobj$Triplicate == 0,"ME"))

#change object design and re-run DESeq
design(DESeqobj) <- ~group1
DESeqobj<-DESeq(DESeqobj)
#output names of resulting group comparisons
resultsNames(DESeqobj)
comparison_323A1 <- results(DESeqobj, contrast=c("group1", "323A1", "0"))
#order by p value
order_323A1<-comparison_323A1[order(comparison_323A1$pvalue),]
plotMA(order_323A1)
write.csv(order_323A1,"323A1_vs_Background")
comparison_324A1 <- results(DESeqobj, contrast=c("group1", "324A1", "0"))
order_324A1<-comparison_324A1[order(comparison_324A1$pvalue),]
plotMA(comparison_324A1)
write.csv(order_324A1,"324A1_vs_Background")
comparison_325A1 <- results(DESeqobj, contrast=c("group1", "325A1", "0"))
order_325A1<-comparison_325A1[order(comparison_323A1$pvalue),]
write.csv(order_325A1,"325A1_vs_Background")

#DESeqobj$group1 <- factor(replace(replace(DESeqobj$Triplicate,DESeqobj$Triplicate != 0,"CL"),DESeqobj$Triplicate == 0,"FA"))
comparison_v_Fly_Atlas <- results(DESeqobj, contrast=c("group1", "CL", "ME"))
#order by p value
order_v_Fly_Atlas<-comparison_v_Fly_Atlas[order(comparison_v_Fly_Atlas$pvalue),]
plotMA(order_v_Fly_Atlas)
write.csv(order_v_Fly_Atlas,"CellLine_vs_modEncode.csv")





DESeqobj$group2 <- factor(replace(DESeqobj$Sample.Group,DESeqobj$Background.Group != 2,0))
design(DESeqobj) <- ~group2
DESeqobj<-DESeq(DESeqobj)
resultsNames(DESeqobj)
comparison_329A1 <- results(DESeqobj, contrast=c("group2", "329A1", "0"))
order_329A1<-comparison_329A1[order(comparison_329A1$pvalue),]
write.csv(order_329A1,"329A1_vs_Background")
comparison_330A1 <- results(DESeqobj, contrast=c("group2", "330A1", "0"))
order_330A1<-comparison_330A1[order(comparison_330A1$pvalue),]
write.csv(order_330A1,"330A1_vs_Background")
comparison_331A1 <- results(DESeqobj, contrast=c("group2", "331A1", "0"))
order_331A1<-comparison_331A1[order(comparison_331A1$pvalue),]
write.csv(order_331A1,"331A1_vs_Background")


DESeqobj$group3 <- factor(replace(DESeqobj$Sample.Group,DESeqobj$Background.Group != 3,0))
design(DESeqobj) <- ~group3
DESeqobj<-DESeq(DESeqobj)
resultsNames(DESeqobj)
comparison_285A1 <- results(DESeqobj, contrast=c("group3", "285A1", "0"))
order_285A1<-comparison_285A1[order(comparison_285A1$pvalue),]
write.csv(order_285A1,"285A1_vs_Background")
comparison_286A1 <- results(DESeqobj, contrast=c("group3", "286A1", "0"))
order_286A1<-comparison_286A1[order(comparison_286A1$pvalue),]
write.csv(order_286A1,"286A1_vs_Background")
comparison_332A1 <- results(DESeqobj, contrast=c("group3", "332A1", "0"))
order_332A1<-comparison_332A1[order(comparison_332A1$pvalue),]
write.csv(order_332A1,"332A1_vs_Background")


DESeqobj$group4 <- factor(replace(DESeqobj$Sample.Group,DESeqobj$Background.Group != 4,0))
design(DESeqobj) <- ~group4
DESeqobj<-DESeq(DESeqobj)
resultsNames(DESeqobj)
comparison_326A1 <- results(DESeqobj, contrast=c("group4", "326A1", "0"))
order_326A1<-comparison_326A1[order(comparison_326A1$pvalue),]
write.csv(order_326A1,"326A1_vs_Background")
comparison_327A1 <- results(DESeqobj, contrast=c("group4", "327A1", "0"))
order_327A1<-comparison_327A1[order(comparison_327A1$pvalue),]
write.csv(order_327A1,"327A1_vs_Background")
comparison_328A1 <- results(DESeqobj, contrast=c("group4", "328A1", "0"))
order_328A1<-comparison_328A1[order(comparison_328A1$pvalue),]
write.csv(order_328A1,"328A1_vs_Background")

#comparison of all samples
deseq2VST <- vst(DESeqobj)
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)
#DE_results <- lfcShrink(DESeqobj, coef=2)
Deseq2VST_long <- melt(deseq2VST, id.vars=c("Gene"))
heatmap <- ggplot(Deseq2VST_long, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap
write.csv(deseq2VST,"FlyAtlasVST.csv")

# Compute a distance calculation on both dimensions of the matrix
distanceGene <- dist(deseq2VST)
distanceSample <- dist(t(deseq2VST))
clusterGene <- hclust(distanceGene, method="average")
clusterSample <- hclust(distanceSample, method="average")
Deseq2VST_long$variable <- factor(Deseq2VST_long$variable, levels=clusterSample$labels[clusterSample$order])
Deseq2VST_long$Gene <- factor(Deseq2VST_long$Gene, levels=clusterGene$labels[clusterGene$order])

deseq2VSTMatrix <- dcast(Deseq2VST_long, Gene ~ variable)

heatmap <- ggplot(Deseq2VST_long, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap

deseq2VSTMatrix <- dcast(Deseq2VST_long, Gene ~ variable)
rownames(deseq2VSTMatrix)<-deseq2VSTMatrix$Gene
deseq2VSTMatrix<-deseq2VSTMatrix[,2:43]
LFC_matrix <- deseq2VSTMatrix - rowMeans(deseq2VSTMatrix)
LFC_mat<-as.matrix(LFC_matrix)
LFC_max<-rowMax(abs(LFC_mat))
LFC_matrix$abs_max<-LFC_max
LFC_sig <- LFC_mat[which(rowMax(abs(LFC_mat))>2),]
write.csv(LFC_matrix, "FlyAtlas_LFC_Matrix")
#DE_results <- lfcShrink(DESeqobj, coef=2)


write.csv(deseq2VSTMatrix, "clustered_DESEQ.csv")
