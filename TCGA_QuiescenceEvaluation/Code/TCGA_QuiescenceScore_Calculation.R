###############################################
###TCGA Quiescence Score calculation
###############################################

#Load required pacakges:
library(GSVA)
library(biomaRt)


############################################
#Common Quiescence score

#Load RNA-seq data:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
load("combined_scaled_expr_data.RData")
#This a limited example dataframe with only 100 entries 
#To download the full dataset follow instructions in the "TCGA_DataDownload" folder

#Load quiescence biomarker genes:
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("downregulated_common_ENSG.RData")
load("upregulated_common_ENSG.RData")

#Select only genes involved in quiescence to reduce the size of the dataframe:
CancerType <- combined.scaled$CancerType 
expr.data <- combined.scaled
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated_common_ENSG, upregulated_common_ENSG))]

#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
gene_lists <- list(upregulated_common_ENSG, downregulated_common_ENSG)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
z_score$X1 <- NULL
z_score$X2 <- NULL
z_score$CancerType <- CancerType

#Save the quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
save(z_score, file = "TCGA_common_QS_purity_scaled.RData")




############################################
#Common Quiescence score with reduced signature

#Load RNA-seq data:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
load("combined_scaled_expr_data.RData")
#This a limited example dataframe with only 100 entries 
#To download the full dataset follow instructions in the "TCGA_DataDownload" folder

#Load quiescence biomarker genes:
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("downregulated_common_ENSG.RData")
load("upregulated_common_ENSG.RData")
load("RefinedGenePanel_35.RData")
#Select genes from the reduced signature:
human = useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "useast.ensembl.org")
results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=GenePanel35, mart=human, filters = "hgnc_symbol")
GenePanel35 <- results$ensembl_gene_id
downregulated_common_ENSG <- downregulated_common_ENSG[downregulated_common_ENSG %in% GenePanel35]
upregulated_common_ENSG <- upregulated_common_ENSG[upregulated_common_ENSG %in% GenePanel35]
#Select only genes involved in quiescence to reduce the size of the dataframe:
CancerType <- combined.scaled$CancerType 
expr.data <- combined.scaled
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated_common_ENSG, upregulated_common_ENSG))]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
gene_lists <- list(upregulated_common_ENSG, downregulated_common_ENSG)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
z_score$X1 <- NULL
z_score$X2 <- NULL
z_score$CancerType <- CancerType
#Save the quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
save(z_score, file = "TCGA_common_QS_purity_scaled_35GenePanel.RData")





############################################
#CDK4/6 inhibition Quiescence score

#Load RNA-seq data:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
load("combined_scaled_expr_data.RData")
#This a limited example dataframe with only 100 entries 
#To download the full dataset follow instructions in the "TCGA_DataDownload" folder

#Load quiescence biomarker genes:
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("cdk_downregulated_ENSG.RData")
load("cdk_upregulated_ENSG.RData")

#Select only genes involved in quiescence to reduce the size of the dataframe:
CancerType <- combined.scaled$CancerType 
expr.data <- combined.scaled
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated_ENSG, upregulated_ENSG))]

#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
gene_lists <- list(upregulated_ENSG, downregulated_ENSG)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
z_score$X1 <- NULL
z_score$X2 <- NULL
z_score$CancerType <- CancerType

#Save the quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
save(z_score, file = "TCGA_CDK_inhibition_QS_purity_scaled.RData")



############################################
#Contact inhibition Quiescence score

#Load RNA-seq data:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
load("combined_scaled_expr_data.RData")
#This a limited example dataframe with only 100 entries 
#To download the full dataset follow instructions in the "TCGA_DataDownload" folder

#Load quiescence biomarker genes:
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("contact_downregulated_ENSG.RData")
load("contact_upregulated_ENSG.RData")

#Select only genes involved in quiescence to reduce the size of the dataframe:
CancerType <- combined.scaled$CancerType 
expr.data <- combined.scaled
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated_ENSG, upregulated_ENSG))]

#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
gene_lists <- list(upregulated_ENSG, downregulated_ENSG)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
z_score$X1 <- NULL
z_score$X2 <- NULL
z_score$CancerType <- CancerType

#Save the quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
save(z_score, file = "TCGA_contact_inhibition_QS_purity_scaled.RData")





############################################
#MEK inhibition Quiescence score

#Load RNA-seq data:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
load("combined_scaled_expr_data.RData")
#This a limited example dataframe with only 100 entries 
#To download the full dataset follow instructions in the "TCGA_DataDownload" folder

#Load quiescence biomarker genes:
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("mek_downregulated_ENSG.RData")
load("mek_upregulated_ENSG.RData")

#Select only genes involved in quiescence to reduce the size of the dataframe:
CancerType <- combined.scaled$CancerType 
expr.data <- combined.scaled
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated_ENSG, upregulated_ENSG))]

#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
gene_lists <- list(upregulated_ENSG, downregulated_ENSG)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
z_score$X1 <- NULL
z_score$X2 <- NULL
z_score$CancerType <- CancerType

#Save the quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
save(z_score, file = "TCGA_MEK_inhibition_QS_purity_scaled.RData")


############################################
#Serum Starvation Quiescence score

#Load RNA-seq data:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
load("combined_scaled_expr_data.RData")
#This a limited example dataframe with only 100 entries 
#To download the full dataset follow instructions in the "TCGA_DataDownload" folder

#Load quiescence biomarker genes:
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("serum_downregulated_ENSG.RData")
load("serum_upregulated_ENSG.RData")

#Select only genes involved in quiescence to reduce the size of the dataframe:
CancerType <- combined.scaled$CancerType 
expr.data <- combined.scaled
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated_ENSG, upregulated_ENSG))]

#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
gene_lists <- list(upregulated_ENSG, downregulated_ENSG)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
z_score$X1 <- NULL
z_score$X2 <- NULL
z_score$CancerType <- CancerType

#Save the quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
save(z_score, file = "TCGA_serum_starvation_QS_purity_scaled.RData")



############################################
#Spontaneous Quiescence score

#Load RNA-seq data:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
load("combined_scaled_expr_data.RData")
#This a limited example dataframe with only 100 entries 
#To download the full dataset follow instructions in the "TCGA_DataDownload" folder

#Load quiescence biomarker genes:
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("spontaneous_downregulated_ENSG.RData")
load("spontaneous_upregulated_ENSG.RData")


#Select only genes involved in quiescence to reduce the size of the dataframe:
CancerType <- combined.scaled$CancerType 
expr.data <- combined.scaled
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated_ENSG, upregulated_ENSG))]

#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
gene_lists <- list(upregulated_ENSG, downregulated_ENSG)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
z_score$X1 <- NULL
z_score$X2 <- NULL
z_score$CancerType <- CancerType

#Save the quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
save(z_score, file = "TCGA_spotaneous_quiescence_QS_purity_scaled.RData")


