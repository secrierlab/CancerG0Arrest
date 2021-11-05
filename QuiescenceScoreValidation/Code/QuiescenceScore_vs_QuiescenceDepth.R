##############################################
#QS scores at different quiescence depths:
##############################################

#Load required packages:
library("biomaRt") 
library(dplyr)
library(GSVA)
library(ggplot2)


################################################################
#Load expression data from cells at different quiescence depths:
expr.data <- read.table("GSE124109_genes.processed.fpkm_table.txt", header = TRUE, sep = "\t") #Data can be downloaded from GEO using accession code GSE124109


#Convert rat to human gene names:
Gene <- unique(as.character(expr.data$Gene))
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
rat = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
genesV2 = getLDS(attributes = c("rgd_symbol"), 
                 filters = "rgd_symbol", 
                 values = Gene , 
                 mart = rat, 
                 attributesL = c("hgnc_symbol"), 
                 martL = human, 
                 uniqueRows=T)
merged.data <- merge(expr.data, genesV2,
                     by.x = "Gene", by.y = "RGD.symbol")
merged.data$Gene <- NULL
merged.data <- merged.data %>% group_by(HGNC.symbol) %>% mutate_each(funs(mean)) %>% distinct
merged.data <- data.frame(merged.data)
rownames(merged.data) <- merged.data$HGNC.symbol
merged.data$HGNC.symbol <- NULL
#Log transform the data:
merged.data <- as.matrix(merged.data)
merged.data <- log2(merged.data + 1)
merged.data <- data.frame(merged.data)


#######################
#Quiescence scoring:
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("upregulated_common.RData")
load("downregulated_common.RData")
merged.data <- as.matrix(merged.data)
gene_lists <- list(upregulated_common, downregulated_common)
es.dif <- gsva(merged.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
es.dif <- t(es.dif)
es.dif <- data.frame(es.dif)
es.dif$GSVA_score <- es.dif$X1 - es.dif$X2
z_score <- es.dif$GSVA_score
z_score <- data.frame(z_score)
rownames(z_score) <- colnames(merged.data)


##########################
#Plot the information:
z_score$CellType <- 1:30
z_score$Day <- c(rep(0,3),rep(2,3),rep(3,3),rep(4,3),rep(6,3),rep(8,3),rep(10,3),rep(12,3),rep(14,3),rep(16,3))
z_score$Day <- factor(z_score$Day, levels = c(0,2,3,4,6,8,10,12,14,16))
# Barplot
p <- ggplot(z_score, aes(x=CellType, y=z_score, fill = Day)) + 
  geom_bar(stat = "identity") + theme_classic()
setwd("~/Documents/GitHub/CancerDormancy/QuiescenceScoreValidation/Figures/")
pdf("QS_vs_QuiescenceDepth.pdf",height = 10,width = 5)
p
dev.off()

