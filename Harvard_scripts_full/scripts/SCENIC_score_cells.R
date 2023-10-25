library(SCENIC)

dbFiles <- "./Databases/"

scenicOptions <- readRDS("int/scenicOptions.Rds")
expr <- read.csv("Expr_Mat_for_SCENIC")
rownames(expr) <- expr$X
expr<-expr[-1]
expr_mat <- as.matrix(expr)
#exprMat_log <- log2(expr_mat+1)
#dim(expr_Mat_log)
#scenicOptions@settings$seed <- 123 
#scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, expr_mat)
#scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, expr_mat, skipHeatmap=TRUE, skipTsne=TRUE)

saveRDS(scenicOptions, file="int/scenicOptions.Rds")
