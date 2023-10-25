library(tidyverse)
library(SCENIC)

dbFiles <- "./Databases/"

scenicOptions <- initializeScenic(org="dmel", dbDir=dbFiles, datasetTitle="Pedro_SCENIC_analysis", nCores=12)
expr <- read.csv("Expr_Mat_pedro_for_GENIE3")#already filtred, but check
rownames(expr) <- expr$X
expr<-expr[-1]
expr_mat <- as.matrix(expr)
genekept<-geneFiltering(expr_mat,scenicOptions = scenicOptions, minCountsPerGene=3*.01*ncol(expr_mat),
                        minSamples=ncol(expr_mat)*.01)
exprMat_filtered <- expr_mat[genekept,]
#exprMat_filtered <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered, scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
