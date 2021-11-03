###########################################
###GSE134839 treatment response analysis
###########################################


###Load required packages:
library(Seurat)
library(dplyr)
library(cowplot)
library(GSVA)
library(umap)
library(ggplot2)
library(ggpubr)


################
#Load expression data (can be downloaded from GEO using the GSE149383 accession code):
D0 <- read.table("GSM3972657_D0.dge.txt", sep = "\t", header = TRUE)
D1 <- read.table("GSM3972658_D1.dge.txt", sep = "\t", header = TRUE)
D2 <- read.table("GSM3972659_D2.dge.txt", sep = "\t", header = TRUE)
D4 <- read.table("GSM3972660_D4.dge.txt", sep = "\t", header = TRUE)
D9 <- read.table("GSM3972661_D9.dge.txt", sep = "\t", header = TRUE)
D11 <- read.table("GSM3972662_D11.dge.txt", sep = "\t", header = TRUE)

##########
##Change colnames to reflect the treatment condition:
#D0:
names <- as.character(1:2500)
names <- paste("D0_", names,sep = "")
names <- c("GENE",names)
colnames(D0) <- names
#D1
names <- as.character(1:650)
names <- paste("D1_", names,sep = "")
names <- c("GENE",names)
colnames(D1) <- names
#D2
names <- as.character(1:790)
names <- paste("D2_", names,sep = "")
names <- c("GENE",names)
colnames(D2) <- names
#D4
names <- as.character(1:1343)
names <- paste("D4_", names,sep = "")
names <- c("GENE",names)
colnames(D4) <- names
#D9
names <- as.character(1:494)
names <- paste("D9_", names,sep = "")
names <- c("GENE",names)
colnames(D9) <- names
#D11
names <- as.character(1:731)
names <- paste("D11_", names,sep = "")
names <- c("GENE",names)
colnames(D11) <- names


##################
#Merge the datasets:
count.data <- merge(D0,D1,
                    by.x = "GENE",by.y = "GENE")
count.data <- merge(count.data,D2,
                    by.x = "GENE",by.y = "GENE")
count.data <- merge(count.data,D4,
                    by.x = "GENE",by.y = "GENE")
count.data <- merge(count.data,D9,
                    by.x = "GENE",by.y = "GENE")
count.data <- merge(count.data,D11,
                    by.x = "GENE",by.y = "GENE")
rownames(count.data) <- count.data$GENE
count.data$GENE <- NULL


#######################
#SEURAT ANALYSIS
#######################

#Set up a seurat object:
lc.data <- CreateSeuratObject(counts = count.data, min.cells = 3, min.features  = 200, project = "PC9_longitudinal", assay = "RNA")
#QS and selecting cells for further analysis 
mito.genes <- grep(pattern = "^MT-", x = rownames(lc.data@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(lc.data@assays[["RNA"]][mito.genes, ])/Matrix::colSums(lc.data@assays[["RNA"]])
lc.data <- AddMetaData(object = lc.data, metadata = percent.mito, col.name = "percent.mito") 
VlnPlot(object = lc.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
par(mfrow = c(1, 2))
FeatureScatter(object = lc.data, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object = lc.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
lc.data <- subset(x = lc.data, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mito >  -Inf & percent.mito < 0.15 )
#Data normalization:
lc.data <- NormalizeData(object = lc.data, normalization.method = "LogNormalize", scale.factor = 10000)
#Extract normalized data
expr.data <- as.matrix(GetAssayData(lc.data, slot = "data"))



##########################################################################
###Quiescence scoring
##########################################################################

#Load quiescence biomarker genes:
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList")
load("upregulated_common.RData")
load("downregulated_common.RData")
quiescence_genes <- c(downregulated_common, upregulated_common)
expr.data <- expr.data[rownames(expr.data) %in% quiescence_genes,]
gene_lists <- list(upregulated_common, downregulated_common)
expr.data <- as.matrix(expr.data)
es.dif <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
es.dif <- t(es.dif)
es.dif <- data.frame(es.dif)
es.dif$Common_score <- es.dif$X1 - es.dif$X2
QS <- es.dif
QS$X1 <- NULL
QS$X2 <- NULL
QS$SampleID <- rownames(QS)
QS$Day <- sapply(QS$SampleID, function(x)
  strsplit(x,"_")[[1]][1])
table(QS$Day)
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Results/")
save(QS, file = "GSE134839_QuiescenceScores.RData")


###########################################
####UMAP analysis 
###########################################

expr.data <- data.frame(expr.data)

######
#DAY0
######
D0_samples <- QS[QS$Day %in% "D0",]
D0_samples <- as.character(D0_samples$SampleID)
umap.expr <- expr.data[,colnames(expr.data) %in% D0_samples]
umap.expr <- as.matrix(t(umap.expr))
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_D0 <- merge(UMAP_coordinates, QS,
                          by.x = "Sample", by.y = "SampleID")

######
#DAY1
#####
D1_samples <- QS[QS$Day %in% "D1",]
D1_samples <- as.character(D1_samples$SampleID)
umap.expr <- expr.data[,colnames(expr.data) %in% D1_samples]
umap.expr <- as.matrix(t(umap.expr))
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_D1 <- merge(UMAP_coordinates, QS,
                          by.x = "Sample", by.y = "SampleID")

######
#DAY2
######
D2_samples <- QS[QS$Day %in% "D2",]
D2_samples <- as.character(D2_samples$SampleID)
umap.expr <- expr.data[,colnames(expr.data) %in% D2_samples]
umap.expr <- as.matrix(t(umap.expr))
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_D2 <- merge(UMAP_coordinates, QS,
                          by.x = "Sample", by.y = "SampleID")

######
#DAY4
######
D4_samples <- QS[QS$Day %in% "D4",]
D4_samples <- as.character(D4_samples$SampleID)
umap.expr <- expr.data[,colnames(expr.data) %in% D4_samples]
umap.expr <- as.matrix(t(umap.expr))
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_D4 <- merge(UMAP_coordinates, QS,
                          by.x = "Sample", by.y = "SampleID")

######
#DAY9
######
D9_samples <- QS[QS$Day %in% "D9",]
D9_samples <- as.character(D9_samples$SampleID)
umap.expr <- expr.data[,colnames(expr.data) %in% D9_samples]
umap.expr <- as.matrix(t(umap.expr))
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_D9 <- merge(UMAP_coordinates, QS,
                          by.x = "Sample", by.y = "SampleID")

######
#DAY11
######
D11_samples <- QS[QS$Day %in% "D11",]
D11_samples <- as.character(D11_samples$SampleID)
umap.expr <- expr.data[,colnames(expr.data) %in% D11_samples]
umap.expr <- as.matrix(t(umap.expr))
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_D11 <- merge(UMAP_coordinates, QS,
                          by.x = "Sample", by.y = "SampleID")

######Plot combined UMAP coordinates:
UMAP_coordinates <- rbind(UMAP_coordinates_D0, UMAP_coordinates_D1, UMAP_coordinates_D2,UMAP_coordinates_D4, UMAP_coordinates_D9, UMAP_coordinates_D11)
UMAP_coordinates$Day <- factor(UMAP_coordinates$Day, levels = c("D0","D1","D2","D4","D9","D11"))
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Figures")
pdf("GSE134839_UMAP.pdf", height = 3, width = 9)
ggplot(UMAP_coordinates, aes(x=UMAP1, y=UMAP2, colour = Common_score)) +
  geom_point() +
  scale_color_gradient2(low = "#59ac53", midpoint = 0,mid = "grey95", high = "#8b5aa8") + theme_classic() + facet_wrap(~Day,nrow = 1)
dev.off()





#################################################
#Percentage of quiescent cells at each timepoint
#################################################
QS$CellStatus <- sapply(QS$Common_score, function(x)
  ifelse(x < 0, "Proliferating","Quiescent"))
table(QS$Day)
Days <- c(0,0,1,1,2,2,4,4,9,9,11,11)
CellStatus <- c("Proliferating","Quiescent","Proliferating","Quiescent","Proliferating","Quiescent","Proliferating","Quiescent","Proliferating","Quiescent","Proliferating","Quiescent")
N <- NULL
for (i in c("D0","D1","D2","D4","D9","D11")) {
  
  print(i)
  test <- QS[QS$Day %in% i,]
  test <- table(test$CellStatus)
  n <- test[1]
  N <- c(N,n)
  n <- test[2]
  N <- c(N,n)
  
}
Summary <- data.frame(Days, CellStatus, N)
Summary$Days <- factor(Summary$Days, levels = c(0,1,2,4,9,11))
pdf("GSE134839_barplot_cell_composition.pdf", height = 5, width = 4)
p <- ggplot(Summary, aes(fill=CellStatus, y=N, x=Days, width = 0.75)) + 
  geom_bar(stat = "identity", position = "fill") + theme_classic()
p + rotate_x_text(45) + scale_fill_manual(values = c("#666666", "#D95F02"))
dev.off()

