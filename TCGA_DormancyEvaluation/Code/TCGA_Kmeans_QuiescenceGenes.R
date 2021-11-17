##############################################################################################
###K-means clustering on the expression of quiescence genes (combat corrected for cancer type)
###############################################################################################

#Load required packages:
library(data.table)
library(ggplot2)
library(cluster)
library(factoextra)

#Select genes of interest
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("upregulated_common.RData")
load("downregulated_common.RData")

#Load expression data
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
combined_data <- fread("TCGA_combat_tumor_type_correction.txt", sep = "\t") 
#This is a toy example data with 100 entries only 

#Load quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
load("TCGA_common_QS_purity_scaled.RData")


#Select samples of interest (which have quiescence score annotation)
combined_data <- data.frame(combined_data)
combined_data <- combined_data[combined_data$genes %in% c(upregulated_common, downregulated_common),]
rownames(combined_data) <- combined_data$genes
combined_data$genes <- NULL
colnames(combined_data) <- gsub('\\.', '-', colnames(combined_data))
colnames(combined_data)  <- sapply(colnames(combined_data) , function(x)
  paste(strsplit(x,"-")[[1]][2:8],collapse="-"))
combined_data <- data.frame(t(combined_data))
combined_data <- combined_data[rownames(combined_data) %in% rownames(z_score),]


#Clustering based on the expression data:
set.seed(123)
clustering <- kmeans(combined_data, 2, iter.max = 100, nstart = 10)
combined_data$cluster <- clustering$cluster
#Add this information to z_score data.frame
cluster1 <- combined_data[combined_data$cluster %in% 1,]
cluster1 <- as.character(rownames(cluster1))
cluster2 <- combined_data[combined_data$cluster %in% 2,]
cluster2 <- as.character(rownames(cluster2))
z_score$cluster <- sapply(rownames(z_score), function(x)
  ifelse(x %in% cluster1,"cluster1",
         ifelse(x %in% cluster2,"cluster2","other")))
table(z_score$cluster)
#Save quiescence scores with cluster annotations:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
save(z_score, file = "QS_with_clusters.RData")

