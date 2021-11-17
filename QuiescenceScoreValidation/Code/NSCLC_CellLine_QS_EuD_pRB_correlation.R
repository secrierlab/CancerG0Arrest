####################################################################################
#Quiescence Scores vs Eud and PhosphoRB estimates in lung adenocarcinoma cell lines
####################################################################################

#Load required packages:
library(GSVA)
library(biomaRt)
library(ggplot2)
library(ggpubr)
library(reshape)


###################################
#Load required data
###################################

#Load pRB and EuD estimates for NSCLC cell lines
setwd("~/Documents/GitHub/CancerDormancy/Data/NSCLC_CellLine_pRB_EuD_estimates")
load("NSCLC_CellLine_pRB_EuD_estimates.RData")

#CCLE expression data for NSCLC cell lines
setwd("~/Documents/GitHub/CancerDormancy/Data/NSCLC_CellLine_CCLE_expression/")
load("NSCLC_CellLine_CCLE_expr.Rdata")


###################################
#Calculating the quiescence scores
###################################

refined_expr <- as.matrix(refined_expr)
refined_expr <- data.frame(t(refined_expr))
#Load list of genes associated with quiescence
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("downregulated_common_ENSG.RData")
load("upregulated_common_ENSG.RData")
#Look at which for which genes out of the list of genes differentially expressed in quiescence we have expression data for:
all_genes <- rownames(refined_expr)
upregulated <- upregulated_common_ENSG[upregulated_common_ENSG %in% all_genes]
downregulated <- downregulated_common_ENSG[downregulated_common_ENSG %in% all_genes]
#calculate quiescence scores
gene_lists <- list(upregulated, downregulated)
refined_expr <- as.matrix(refined_expr)
z_score <- gsva(refined_expr, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
Validation_dataset$Quiescence <- z_score$z_score

###################################
#Calculating G1 score
###################################

#Load list of genes associated with G1
setwd("~/Documents/GitHub/CancerDormancy/Data/OtherGeneLists/")
markers <- read.table("REACTOME_G1.txt", header = FALSE,sep = "\t") #Obtained from msigdb (REACTOME_G1_PHASE)
markers <- as.character(markers$V1)
human = useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "useast.ensembl.org")
results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=markers, mart=human, filters = "hgnc_symbol")
markers <- results$ensembl_gene_id
gene_lists <- list(markers)
refined_expr <- as.matrix(refined_expr)
z_score <- gsva(refined_expr, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
Validation_dataset$G1 <- z_score$z_score




#########################
#Correlation with EuD and PhosphoRB measurements:
#Plot the correlation:
melt.data <- melt(Validation_dataset, id.vars = c("Cell_line","p53_status","Percentage_G0_p21","Percentage_G0_EuD","Percentage_G0_phosphoRB"))
colnames(melt.data) <- c("Cell_line","p53_status","Percentage_G0_p21","Percentage_G0_EuD","Percentage_G0_phosphoRB","Z_score","value")
pdf("NSCLC_CellLine_QS_GQ_EuD_correlation.pdf", width = 5, height = 6)

p <- ggscatter(melt.data, x = "Percentage_G0_EuD", y = "value",
               add = "reg.line",                         # Add regression line
               conf.int = TRUE,                          # Add confidence interval
               color = "Z_score", palette = "jco",           # Color by groups "cyl"
               shape = "Z_score"                             # Change point shape by groups "cyl"
)+
  stat_cor(aes(color = Z_score), label.x = 3)           # Add correlation coefficient
p + labs(x = "AVG % Qui EuD", y = "Z-score value")          # Add correlation coefficient
dev.off()



pdf("NSCLC_CellLine_QS_GQ_pRB_correlation.pdf", width = 5, height = 6)

p <- ggscatter(melt.data, x = "Percentage_G0_phosphoRB", y = "value",
               add = "reg.line",                         # Add regression line
               conf.int = TRUE,                          # Add confidence interval
               color = "Z_score", palette = "jco",           # Color by groups "cyl"
               shape = "Z_score"                             # Change point shape by groups "cyl"
)+
  stat_cor(aes(color = Z_score), label.x = 3)           # Add correlation coefficient
p + labs(x = "AVG % Qui PRb", y = "Z-score value")          # Add correlation coefficient
dev.off()


