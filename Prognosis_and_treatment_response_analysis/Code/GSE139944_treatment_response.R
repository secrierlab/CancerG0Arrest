##################################################################
#GSE139944 - Quiescence dynamics upon various treatment modalities
##################################################################

#Load required R packages:
library(Seurat)
library(dplyr)
library(cowplot)
library(GSVA)
library(ggpubr)


######################
#A549 cells
######################


#Load rna-seq data
combined.data <- readRDS("GSM4150378_sciPlex3_A549_24hrs.rds") #Data can be obtained from GEO using accession code GSE139944
###This is annotation data
anno <- combined.data@colData
anno <- data.frame(anno)
#This is expression data:
expr.data <- combined.data@assays$data$counts
#Set up Seurat object
lc.data <- CreateSeuratObject(counts = expr.data, min.cells = 3, min.features  = 200, project = "ScreenStudy", assay = "RNA")
mito.genes <- grep(pattern = "^MT-", x = rownames(lc.data@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(lc.data@assays[["RNA"]][mito.genes, ])/Matrix::colSums(lc.data@assays[["RNA"]])
lc.data <- AddMetaData(object = lc.data, metadata = percent.mito, col.name = "percent.mito") 
VlnPlot(object = lc.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
par(mfrow = c(1, 2))
FeatureScatter(object = lc.data, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object = lc.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# Filter out cells that have unique gene counts over 6,000 or mitochondrial
# content over 15%
lc.data <- subset(x = lc.data, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mito >  -Inf & percent.mito < 0.15 )
#Normalising expression data:
lc.data <- NormalizeData(object = lc.data, normalization.method = "LogNormalize", scale.factor = 10000)
#Select expression data for quiescence gens only
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("downregulated_common.RData")
load("upregulated_common.RData")
genes <- c(upregulated_common, downregulated_common)
gene.anno <- read.table("GSM4150378_sciPlex3_A549_MCF7_K562_screen_gene.annotations.txt", header = TRUE, sep = " ")#Data can be obtained from GEO using accession code GSE139944
gene.anno <- gene.anno[gene.anno$gene_short_name %in% genes,]
genes <- as.character(gene.anno$id)
subset.seurat <- subset(lc.data, features = genes)
expr.data <- as.matrix(GetAssayData(subset.seurat, slot = "data"))

#############################
#Quiescence Scoring
upregulated_common <- gene.anno[gene.anno$gene_short_name %in% upregulated_common,]
upregulated_common <- upregulated_common$id
downregulated_common <- gene.anno[gene.anno$gene_short_name %in% downregulated_common,]
downregulated_common <- downregulated_common$id
gene_lists <- list(upregulated_common, downregulated_common)
expr.data <- as.matrix(expr.data)
es.dif <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
es.dif <- t(es.dif)
es.dif <- data.frame(es.dif)
es.dif$Common_score <- es.dif$X1 - es.dif$X2
z_scores <- es.dif
z_scores$X1 <- NULL
z_scores$X2 <- NULL
QS <- anno
QS <- QS[QS$cell %in% rownames(z_scores),]
QS$z_score <- z_scores$Common_score
QS$product_name <- as.character(QS$product_name)
QS$pathway_level_1 <- as.character(QS$pathway_level_1)
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Results/")
save(QS, file = "GSE139944_QS_A549.RData")

#Plot % of cells that are estimated to be quiescent after the different treatments
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Figures/")
treatments <- unique(as.character(QS$product_name))
Percentage <- NULL
for (i in treatments) {
  
  test <- QS[QS$product_name %in% i,]
  total <- dim(test)[1]
  quiescent <- test[test$z_score > 0,]
  quiescent <- dim(quiescent)[1]
  percentage <- quiescent / total
  percentage <- percentage * 100
  Percentage <- c(Percentage, percentage)
  
}
results <- data.frame(treatments, Percentage)
results$Change <- sapply(results$Percentage, function(x)
  ifelse(x > 46.37858, "Increase","Decrease")) #46.37858 is the % of cells that are quiescent in the vehicle treatment
#Select top 30 and bottom 30 results
results <- results[order(results$Percentage),]
bottom_30 <- results[1:30,]
bottom_30 <- as.character(bottom_30$treatments)
results <- results[order(-results$Percentage),]
top_30 <- results[1:30,]
top_30 <- as.character(top_30$treatments)
selected_treatments <- c(top_30, bottom_30, "Vehicle")
results <- results[results$treatments %in% selected_treatments,]
p <- ggbarplot(results, x = "treatments", y = "Percentage",
               fill = "Change",           # change fill color by mpg_level
               color = "white",            # Set bar border colors to white
               palette = c("#FC4E07","#00AFBB"),           
               sort.val = "desc",          # Sort the value in descending order
               sort.by.groups = FALSE,     # Don't sort inside each group
               x.text.angle = 90,          # Rotate vertically x axis texts
               ylab = "Percentage of Quiescent Cells",
               legend.title = "Quiescent cell population change",
               xlab = "Treatment",
               rotate = TRUE,
)
pdf("GSE139944_A549_treatment_name_percentage_quiescence.pdf",height = 20,width = 10)
p + theme_classic()
dev.off()

#Same analysis at the pathway level
treatments <- unique(as.character(QS$pathway_level_1))
Percentage <- NULL
for (i in treatments) {
  
  test <- QS[QS$pathway_level_1 %in% i,]
  total <- dim(test)[1]
  quiescent <- test[test$z_score > 0,]
  quiescent <- dim(quiescent)[1]
  percentage <- quiescent / total
  percentage <- percentage * 100
  Percentage <- c(Percentage, percentage)
  
}

results <- data.frame(treatments, Percentage)
results$Change <- sapply(results$Percentage, function(x)
  ifelse(x > 46.37858, "Increase","Decrease"))
library(ggpubr)
p <- ggbarplot(results, x = "treatments", y = "Percentage",
               fill = "Change",           # change fill color by mpg_level
               color = "white",            # Set bar border colors to white
               palette = c("#FC4E07","#00AFBB"),           
               sort.val = "desc",          # Sort the value in descending order
               sort.by.groups = FALSE,     # Don't sort inside each group
               x.text.angle = 90,          # Rotate vertically x axis texts
               ylab = "Percentage of Quiescent Cells",
               legend.title = "Quiescent cell population change",
               xlab = "Treatment",
               rotate = TRUE,
)
pdf("GSE139944_A549_treatment_pathway_percentage_quiescence.pdf",height = 5,width = 10)
p + theme_classic()
dev.off()





######################
#MCF7 cells
######################

#################
#Data processing
#Load expr data
combined.data <- readRDS("GSM4150378_sciPlex3_MCF7_24hrs.RDS") #Data can be downloaded from GEO using accession code GSE139944
###This is annotation data
anno <- combined.data@colData
anno <- data.frame(anno)
#This is expression data:
expr.data <- combined.data@assays$data$counts
#Set up Seurat object
lc.data <- CreateSeuratObject(counts = expr.data, min.cells = 3, min.features  = 200, project = "ScreenStudy", assay = "RNA")
mito.genes <- grep(pattern = "^MT-", x = rownames(lc.data@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(lc.data@assays[["RNA"]][mito.genes, ])/Matrix::colSums(lc.data@assays[["RNA"]])
lc.data <- AddMetaData(object = lc.data, metadata = percent.mito, col.name = "percent.mito") 
VlnPlot(object = lc.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
par(mfrow = c(1, 2))
FeatureScatter(object = lc.data, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object = lc.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# Filter out cells that have unique gene counts over 6,000 or mitochondrial
# content over 15%
lc.data <- subset(x = lc.data, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mito >  -Inf & percent.mito < 0.15 )
#############################
#Normalisign the data:
lc.data <- NormalizeData(object = lc.data, normalization.method = "LogNormalize", scale.factor = 10000)
#Select expression of dormancy genes only
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("downregulated_common.RData")
load("upregulated_common.RData")
genes <- c(upregulated_common, downregulated_common)
gene.anno <- read.table("GSM4150378_sciPlex3_A549_MCF7_K562_screen_gene.annotations.txt", header = TRUE, sep = " ")#Data can be downloaded from GEO using accession code GSE139944
gene.anno <- gene.anno[gene.anno$gene_short_name %in% genes,]
genes <- as.character(gene.anno$id)
subset.seurat <- subset(lc.data, features = genes)
expr.data <- as.matrix(GetAssayData(subset.seurat, slot = "data"))


#############################
#Quiescence Scoring
upregulated_common <- gene.anno[gene.anno$gene_short_name %in% upregulated_common,]
upregulated_common <- upregulated_common$id
downregulated_common <- gene.anno[gene.anno$gene_short_name %in% downregulated_common,]
downregulated_common <- downregulated_common$id
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
QS <- anno
QS <- QS[QS$cell %in% rownames(z_scores),]
QS$z_score <- z_scores$Common_score
QS$product_name <- as.character(QS$product_name)
QS$pathway_level_1 <- as.character(QS$pathway_level_1)
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Results/")
save(QS, file = "GSE139944_QS_MCF7.RData")



###Plot percentage of cells that are quiescent across different treatments:
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Figures/")
treatments <- unique(as.character(QS$product_name))
Percentage <- NULL
for (i in treatments) {
  
  test <- QS[QS$product_name %in% i,]
  total <- dim(test)[1]
  quiescent <- test[test$z_score > 0,]
  quiescent <- dim(quiescent)[1]
  percentage <- quiescent / total
  percentage <- percentage * 100
  Percentage <- c(Percentage, percentage)
  
}

results <- data.frame(treatments, Percentage)
results$Change <- sapply(results$Percentage, function(x)
  ifelse(x > 44.81013, "Increase","Decrease")) #44.81013 is the percentage of cells that are quiescent in Vehicle treated samples
#Select top 30 and bottom 30
results <- results[order(results$Percentage),]
bottom_30 <- results[1:30,]
bottom_30 <- as.character(bottom_30$treatments)
results <- results[order(-results$Percentage),]
top_30 <- results[1:30,]
top_30 <- as.character(top_30$treatments)
selected_treatments <- c(top_30, bottom_30, "Vehicle")
results <- results[results$treatments %in% selected_treatments,]
p <- ggbarplot(results, x = "treatments", y = "Percentage",
               fill = "Change",           # change fill color by mpg_level
               color = "white",            # Set bar border colors to white
               palette = c("#FC4E07","#00AFBB"),           
               sort.val = "desc",          # Sort the value in descending order
               sort.by.groups = FALSE,     # Don't sort inside each group
               x.text.angle = 90,          # Rotate vertically x axis texts
               ylab = "Percentage of Quiescent Cells",
               legend.title = "Quiescent cell population change",
               xlab = "Treatment",
               rotate = TRUE,
)
pdf("GSE134499_MCF7_treatment_name_percentage_quiescence.pdf",height = 20,width = 10)
p + theme_classic()
dev.off()

#Repeat the analysis on a pathway level
treatments <- unique(as.character(QS$pathway_level_1))
Percentage <- NULL
for (i in treatments) {
  
  test <- QS[QS$pathway_level_1 %in% i,]
  total <- dim(test)[1]
  quiescent <- test[test$z_score > 0,]
  quiescent <- dim(quiescent)[1]
  percentage <- quiescent / total
  percentage <- percentage * 100
  Percentage <- c(Percentage, percentage)
  
}

results <- data.frame(treatments, Percentage)
results$Change <- sapply(results$Percentage, function(x)
  ifelse(x > 44.81013, "Increase","Decrease"))
p <- ggbarplot(results, x = "treatments", y = "Percentage",
               fill = "Change",           # change fill color by mpg_level
               color = "white",            # Set bar border colors to white
               palette = c("#FC4E07","#00AFBB"),           
               sort.val = "desc",          # Sort the value in descending order
               sort.by.groups = FALSE,     # Don't sort inside each group
               x.text.angle = 90,          # Rotate vertically x axis texts
               ylab = "Percentage of Quiescent Cells",
               legend.title = "Quiescent cell population change",
               xlab = "Treatment",
               rotate = TRUE,
)
pdf("GSE134499_MCF7_treatment_pathway_percentage_quiescence.pdf",height = 5,width = 10)
p + theme_classic()
dev.off()







