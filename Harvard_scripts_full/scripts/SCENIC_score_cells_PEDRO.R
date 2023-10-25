library(SCENIC)

dbFiles <- "./Databases/"

scenicOptions <- readRDS("int/scenicOptions.Rds")
expr <- read.csv("Expr_Mat_pedro_for_GENIE3")
rownames(expr) <- expr$X
expr<-expr[-1]
#colnames(expr)<-paste0(substr(colnames(expr),1,16),"-", substr(colnames(expr),18,18))
expr_mat <- as.matrix(expr)
#used expr_mat to run scenic - this is not the log2 one from the tutorial, but the lognormalized and scaled one from seurat.
#The code to log normalize the expression matrix is thus commented out

#exprMat_log <- log2(expr_mat+1)
#dim(expr_Mat_log)
#scenicOptions@settings$seed <- 123 
#scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, expr_mat)
#scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, expr_mat, skipHeatmap=TRUE, skipTsne=TRUE)

saveRDS(scenicOptions, file="int/scenicOptions.Rds")

int <- readRDS("int/cellInfo.Rds")
