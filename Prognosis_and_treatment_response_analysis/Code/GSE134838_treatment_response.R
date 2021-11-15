#########################################################################
###GSE134838 treatment response analysis (melanoma and braf inhibitors)
#########################################################################

###Load required packages:
library(Seurat)
library(dplyr)
library(cowplot)
library(GSVA)
library(umap)
library(ggplot2)
library(ggpubr)

################
#Load data:
D0 <- read.table("GSM3972655_M14Day0.dge.txt", sep = "\t", header = TRUE)
D3 <- read.table("GSM3972656_M14Day3_Vem.dge.txt", sep = "\t", header = TRUE)
#Data can be obtained from GEO using accession code GSE134838

##########
##Change colnames to reflect condition:
#D0:
names <- as.character(1:5260)
names <- paste("D0_", names,sep = "")
names <- c("GENE",names)
colnames(D0) <- names
#D3
names <- as.character(1:7031)
names <- paste("D3_", names,sep = "")
names <- c("GENE",names)
colnames(D3) <- names
##################
#Merge the datasets:
count.data <- merge(D0,D3,
                    by.x = "GENE",by.y = "GENE")
#####
#Clean up:
rownames(count.data) <- count.data$GENE
count.data$GENE <- NULL
################
#Set up a seurat object:
lc.data <- CreateSeuratObject(counts = count.data, min.cells = 3, min.features  = 200, project = "PC9_longitudinal", assay = "RNA")
#####################
#####QS and selecting cells for further analysis 
mito.genes <- grep(pattern = "^MT-", x = rownames(lc.data@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(lc.data@assays[["RNA"]][mito.genes, ])/Matrix::colSums(lc.data@assays[["RNA"]])
lc.data <- AddMetaData(object = lc.data, metadata = percent.mito, col.name = "percent.mito") 
VlnPlot(object = lc.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
par(mfrow = c(1, 2))
FeatureScatter(object = lc.data, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object = lc.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# Filter out cells that have unique gene counts over 6,000 or mitochondrial
# content over 15%
lc.data <- subset(x = lc.data, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mito >  -Inf & percent.mito < 0.15 )
#Normalisign the data:
lc.data <- NormalizeData(object = lc.data, normalization.method = "LogNormalize", scale.factor = 10000)
expr.data <- as.matrix(GetAssayData(lc.data, slot = "data"))




##########################################################################
###Quiescence scoring

##Load quiescence marker genes
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("upregulated_common.RData")
load("downregulated_common.RData")
quiescence_genes <- c(downregulated_common, upregulated_common)
expr.data <- data.frame(expr.data)
expr.data <- expr.data[rownames(expr.data) %in% quiescence_genes,]
##Calculting quiescence scores:
gene_lists <- list(upregulated_common, downregulated_common)
expr.data <- as.matrix(expr.data)
es.dif <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
es.dif <- t(es.dif)
es.dif <- data.frame(es.dif)
es.dif$Common_score <- es.dif$X1 - es.dif$X2
z_scores <- es.dif
z_scores$X1 <- NULL
z_scores$X2 <- NULL
QS <- z_scores
QS$X1 <- NULL
QS$X2 <- NULL
#Save quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Results/")
save(z_scores, file = "GSE134838_QS.RData")


###########################################
####UMAP:
QS$SampleID <- rownames(QS)
QS$Day <- sapply(QS$SampleID, function(x)
  strsplit(x,"_")[[1]][1])

######
#DAY0
D0_samples <- QS[QS$Day %in% "D0",]
D0_samples <- as.character(D0_samples$SampleID)
umap.expr <- expr.data[,colnames(expr.data) %in% D0_samples]
umap.expr <- umap.expr[rownames(umap.expr) %in% c(downregulated_common, upregulated_common),]
umap.expr <- as.matrix(t(umap.expr))
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
##This is the matrix with cooridnates
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_D0 <- UMAP_coordinates
UMAP_coordinates_D0$Day <- "Day 0"

######
#DAY3
D3_samples <- QS[QS$Day %in% "D3",]
D3_samples <- as.character(D3_samples$SampleID)
umap.expr <- expr.data[,colnames(expr.data) %in% D3_samples]
umap.expr <- umap.expr[rownames(umap.expr) %in% c(downregulated_common, upregulated_common),]
umap.expr <- as.matrix(t(umap.expr))
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
##This is the matrix with cooridnates
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_D3 <- UMAP_coordinates
UMAP_coordinates_D3$Day <- "Day 3"


######
#Plot UMAP coordiantes together:
UMAP_coordinates <- rbind(UMAP_coordinates_D0, UMAP_coordinates_D3)
UMAP_coordinates <- merge(UMAP_coordinates, QS,
                          by.x = "Sample", by.y = "SampleID")
UMAP_coordinates$Day <- UMAP_coordinates$Day.x
UMAP_coordinates$QuiescenceScore <- UMAP_coordinates$Common_score
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Figures")
pdf("GSE134838_UMAP.pdf", height = 4, width = 8)
ggplot(UMAP_coordinates, aes(x=UMAP1, y=UMAP2, colour = QuiescenceScore)) +
  geom_point() +
  scale_color_gradient2(low = "#59ac53", midpoint = 0,mid = "grey95", high = "#8b5aa8") + theme_classic() + facet_wrap(~Day)
dev.off()





#############################
#Percentage of quiescent cells in each condition:
QS$CellStatus <- sapply(QS$Common_score, function(x)
  ifelse(x < 0, "Proliferating","Quiescent"))
table(QS$Day)
Days <- c(0,0,3,3)
CellStatus <- c("Proliferating","Quiescent","Proliferating","Quiescent")
N <- NULL
for (i in c('D0',"D3")) {
  
  print(i)
  test <- QS[QS$Day %in% i,]
  test <- table(test$CellStatus)
  n <- test[1]
  N <- c(N,n)
  n <- test[2]
  N <- c(N,n)
  
}
Summary <- data.frame(Days, CellStatus, N)
Summary$Days <- as.factor(Summary$Days)
pdf("GSE134838_barplot_cell_composition.pdf", height = 5, width = 4)
p <- ggplot(Summary, aes(fill=CellStatus, y=N, x=Days, width = 0.75)) + 
  geom_bar(stat = "identity", position = "fill") + theme_classic()
p + rotate_x_text(45) + scale_fill_manual(values = c("#666666", "#D95F02"))
dev.off()





