library(tidyverse)
library(SCENIC)

dbFiles <- "./Databases/"
expr <- read.csv("Expr_Mat_for_GENIE3")
rownames(expr) <- expr$X
expr<-expr[-1]
expr_mat <- as.matrix(expr)
scenicOptions <- readRDS("int/scenicOptions.Rds")
genekept<-geneFiltering(expr_mat,scenicOptions = scenicOptions,minCountsPerGene=3*.01*ncol(expr_mat),
                        minSamples=ncol(expr_mat)*.01)
exprMat_filtered <- expr_mat[genekept,]
#exprMat_filtered <- log2(exprMat_filtered+1) 

runCorrelation(exprMat_filtered, scenicOptions)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)

saveRDS(scenicOptions, file="int/scenicOptions.Rds")
