##################################################################
###GSE149224 - Colorectal cancer cell line 5-FU treatment response
##################################################################


###Load required packages:
library(Seurat)
library(dplyr)
library(cowplot)
library(GSVA)
library(data.table)
library(ggplot2)


##Load expression data and annotation (can be downloaded from GEO using GSE149224 accession code)
anno <-read.csv("GSE149224_meta.information.csv", header = TRUE)
expr.data <- read.table("GSE149224_RSH.all.txt", header = TRUE,sep = " ") 
#(data is already log normalised)


#Load quiescence biomarker genes
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("upregulated_common.RData")
load("downregulated_common.RData")


################################
##Quiescence Score Calculation:
################################
gene_lists <- list(upregulated_common, downregulated_common)
expr.data <- as.matrix(expr.data)
es.dif <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
es.dif <- t(es.dif)
es.dif <- data.frame(es.dif)
es.dif$Common_score <- es.dif$X1 - es.dif$X2
z_scores <- es.dif
z_scores$X1 <- NULL
z_scores$X2 <- NULL


#Merge with annotation data:
z_scores$SampleID <- rownames(z_scores)
z_scores$SampleID <- gsub('\\.', '-', z_scores$SampleID)
z_scores$SampleID <- gsub('\\X', '', z_scores$SampleID)
all(z_scores$SampleID == anno$X)
anno$QS <- z_scores$Common_score
rownames(anno) <- rownames(z_scores)
anno$SampleID <- rownames(anno)
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Results/")
save(anno, file = "GSE149224_QuiescenceScores.RData")





###########################################
####UMAP ANALYSIS - RKO cell line

#Select RKO cells
anno_RKO <- anno[anno$df.gid %in% "RKO",]
anno_RKO <- as.character(rownames(anno_RKO))
expr.data_RKO <- expr.data[,colnames(expr.data) %in% anno_RKO]

#Dose 0
anno_RKO <- anno[anno$df.gid %in% "RKO",]
anno_RKO <- anno_RKO[anno_RKO$dose %in% 0,]
anno_RKO <- as.character(rownames(anno_RKO))
umap.expr <- expr.data_RKO[,colnames(expr.data_RKO) %in% anno_RKO]
umap.expr <- umap.expr[rownames(umap.expr) %in% c(downregulated_common, upregulated_common),]
umap.expr <- as.matrix(t(umap.expr))
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_D0 <- merge(UMAP_coordinates, anno,
                             by.x = "Sample", by.y = "SampleID")

#Dose 10
anno_RKO <- anno[anno$df.gid %in% "RKO",]
anno_RKO <- anno_RKO[anno_RKO$dose %in% 10,]
anno_RKO <- as.character(rownames(anno_RKO))
umap.expr <- expr.data_RKO[,colnames(expr.data_RKO) %in% anno_RKO]
umap.expr <- umap.expr[rownames(umap.expr) %in% c(downregulated_common, upregulated_common),]
umap.expr <- as.matrix(t(umap.expr))
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_D10 <- merge(UMAP_coordinates, anno,
                              by.x = "Sample", by.y = "SampleID")

#Dose 50
anno_RKO <- anno[anno$df.gid %in% "RKO",]
anno_RKO <- anno_RKO[anno_RKO$dose %in% 50,]
anno_RKO <- as.character(rownames(anno_RKO))
umap.expr <- expr.data_RKO[,colnames(expr.data_RKO) %in% anno_RKO]
umap.expr <- umap.expr[rownames(umap.expr) %in% c(downregulated_common, upregulated_common),]
umap.expr <- as.matrix(t(umap.expr))
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_D50 <- merge(UMAP_coordinates, anno,
                              by.x = "Sample", by.y = "SampleID")


#Dose 200
anno_RKO <- anno[anno$df.gid %in% "RKO",]
anno_RKO <- anno_RKO[anno_RKO$dose %in% 200,]
anno_RKO <- as.character(rownames(anno_RKO))
umap.expr <- expr.data_RKO[,colnames(expr.data_RKO) %in% anno_RKO]
umap.expr <- umap.expr[rownames(umap.expr) %in% c(downregulated_common, upregulated_common),]
umap.expr <- as.matrix(t(umap.expr))
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
##This is the matrix with cooridnates
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_D200 <- merge(UMAP_coordinates, anno,
                               by.x = "Sample", by.y = "SampleID")


###Plot combined UMAP plots:
UMAP_coordinates <- rbind(UMAP_coordinates_D0, UMAP_coordinates_D10, UMAP_coordinates_D50, UMAP_coordinates_D200)
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Figures")
pdf("GSE149224_UMAP_RKO.pdf", height = 3, width = 9)
ggplot(UMAP_coordinates, aes(x=UMAP1, y=UMAP2, colour = QS)) +
  geom_point() +
  scale_color_gradient2(low = "#59ac53", midpoint = 0,mid = "grey95", high = "#8b5aa8") + theme_classic() + facet_wrap(~dose,nrow = 1)
dev.off()



#Percentage of cells in quiescence at each dose:
#Summarise the data:
anno$CellStatus <- sapply(anno$QS, function(x)
  ifelse(x < 0, "Proliferating","Quiescent"))
Dose <- c(0,0,10,10,50,50,200,200)
CellStatus <- c("Proliferating","Quiescent","Proliferating","Quiescent","Proliferating","Quiescent","Proliferating","Quiescent")
N <- NULL
for (i in c(0,10,50,200)) {
  
  print(i)
  test <- anno[anno$df.gid %in% "RKO",]
  test <- test[test$dose %in% i,]
  test <- table(test$CellStatus)
  n <- test[1]
  N <- c(N,n)
  n <- test[2]
  N <- c(N,n)
  
}
Summary <- data.frame(Dose, CellStatus, N)
Summary$Dose <- factor(Summary$Dose, levels = c(0,10,50,200))
pdf("GSE149224_barplot_cell_composition_RKO.pdf", height = 5, width = 4)
p <- ggplot(Summary, aes(fill=CellStatus, y=N, x=Dose, width = 0.75)) + 
  geom_bar(stat = "identity", position = "fill") + theme_classic()
p + rotate_x_text(45) + scale_fill_manual(values = c("#666666", "#D95F02"))
dev.off()










###########################################
####UMAP ANALYSIS - SW480 cell line

#Select SW480 cells
anno_SW480 <- anno[anno$df.gid %in% "SW480",]
anno_SW480 <- as.character(rownames(anno_SW480))
expr.data_SW480 <- expr.data[,colnames(expr.data) %in% anno_SW480]

#Dose 0
anno_SW480 <- anno[anno$df.gid %in% "SW480",]
anno_SW480 <- anno_SW480[anno_SW480$dose %in% 0,]
anno_SW480 <- as.character(rownames(anno_SW480))
umap.expr <- expr.data_SW480[,colnames(expr.data_SW480) %in% anno_SW480]
umap.expr <- umap.expr[rownames(umap.expr) %in% c(downregulated_common, upregulated_common),]
umap.expr <- as.matrix(t(umap.expr))
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_D0 <- merge(UMAP_coordinates, anno,
                             by.x = "Sample", by.y = "SampleID")

#Dose 10
anno_SW480 <- anno[anno$df.gid %in% "SW480",]
anno_SW480 <- anno_SW480[anno_SW480$dose %in% 10,]
anno_SW480 <- as.character(rownames(anno_SW480))
umap.expr <- expr.data_SW480[,colnames(expr.data_SW480) %in% anno_SW480]
umap.expr <- umap.expr[rownames(umap.expr) %in% c(downregulated_common, upregulated_common),]
umap.expr <- as.matrix(t(umap.expr))
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_D10 <- merge(UMAP_coordinates, anno,
                              by.x = "Sample", by.y = "SampleID")

#Dose 50
anno_SW480 <- anno[anno$df.gid %in% "SW480",]
anno_SW480 <- anno_SW480[anno_SW480$dose %in% 50,]
anno_SW480 <- as.character(rownames(anno_SW480))
umap.expr <- expr.data_SW480[,colnames(expr.data_SW480) %in% anno_SW480]
umap.expr <- umap.expr[rownames(umap.expr) %in% c(downregulated_common, upregulated_common),]
umap.expr <- as.matrix(t(umap.expr))
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_D50 <- merge(UMAP_coordinates, anno,
                              by.x = "Sample", by.y = "SampleID")

#Dose 200
anno_SW480 <- anno[anno$df.gid %in% "SW480",]
anno_SW480 <- anno_SW480[anno_SW480$dose %in% 200,]
anno_SW480 <- as.character(rownames(anno_SW480))
umap.expr <- expr.data_SW480[,colnames(expr.data_SW480) %in% anno_SW480]
umap.expr <- umap.expr[rownames(umap.expr) %in% c(downregulated_common, upregulated_common),]
umap.expr <- as.matrix(t(umap.expr))
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_D200 <- merge(UMAP_coordinates, anno,
                               by.x = "Sample", by.y = "SampleID")




###Plot combined UMAP plots:
UMAP_coordinates <- rbind(UMAP_coordinates_D0, UMAP_coordinates_D10, UMAP_coordinates_D50, UMAP_coordinates_D200)
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Figures")
pdf("GSE149224_UMAP_SW480.pdf", height = 3, width = 9)
ggplot(UMAP_coordinates, aes(x=UMAP1, y=UMAP2, colour = QS)) +
  geom_point() +
  scale_color_gradient2(low = "#59ac53", midpoint = 0,mid = "grey95", high = "#8b5aa8") + theme_classic() + facet_wrap(~dose,nrow = 1)
dev.off()



#Percentage of cells in quiescence at each dose:
#Summarise the data:
anno$CellStatus <- sapply(anno$QS, function(x)
  ifelse(x < 0, "Proliferating","Quiescent"))
Dose <- c(0,0,10,10,50,50,200,200)
CellStatus <- c("Proliferating","Quiescent","Proliferating","Quiescent","Proliferating","Quiescent","Proliferating","Quiescent")
N <- NULL
for (i in c(0,10,50,200)) {
  
  print(i)
  test <- anno[anno$df.gid %in% "SW480",]
  test <- test[test$dose %in% i,]
  test <- table(test$CellStatus)
  n <- test[1]
  N <- c(N,n)
  n <- test[2]
  N <- c(N,n)
  
}
Summary <- data.frame(Dose, CellStatus, N)
Summary$Dose <- factor(Summary$Dose, levels = c(0,10,50,200))
pdf("GSE149224_barplot_cell_composition_SW480.pdf", height = 5, width = 4)
p <- ggplot(Summary, aes(fill=CellStatus, y=N, x=Dose, width = 0.75)) + 
  geom_bar(stat = "identity", position = "fill") + theme_classic()
p + rotate_x_text(45) + scale_fill_manual(values = c("#666666", "#D95F02"))
dev.off()


