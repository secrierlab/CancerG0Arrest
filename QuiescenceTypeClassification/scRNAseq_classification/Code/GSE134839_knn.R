#############################################################################
##GSE134839 - KNN analsysis 
#############################################################################

#Load required R packages:
library(biomaRt)
library(plyr)
library(GeoTcgaData)
library(data.table)
library(sva)
library(class)
library(ggplot2)
library(ggpubr)



##################################
#Load bulk RNA-seq reference data:
count.data <- read.table("GSE122927_ReadCount.txt", header = TRUE,sep = "\t") # (data can be downloaded from GEO using accession code GSE122927)
rownames(count.data) <- count.data$X
#Change gene anntation to HGNC
human = useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "useast.ensembl.org")
count.data$X <- as.character(count.data$X)
genes <- count.data$X
results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=genes, mart=human, filters = "ensembl_gene_id")
count.matrix <- merge(count.data, results,
                      by.x = "X", by.y = "ensembl_gene_id")
count.matrix$X <- NULL
#Set the rownames as HGNC gene names:
count.matrix <- ddply(count.matrix,"hgnc_symbol",numcolwise(sum))
rownames(count.matrix) <- count.matrix$hgnc_symbol
count.matrix$hgnc_symbol <- NULL
#TPM normalisation
count.matrix <- as.matrix(count.matrix)
TPM_expression <- countToTpm_matrix(count.matrix)
reference_data <- TPM_expression
rm(TPM_expression)
##Make sure this is on the log scale:
reference_data <- log2(reference_data + 1)
reference_data <- data.frame(reference_data)
#The following samples are controls:
reference_data$hgnc_symbol <- rownames(reference_data)
reference_data$Sample3 <- NULL
reference_data$Sample6 <- NULL
reference_data$Sample7 <- NULL
reference_data$Sample9 <- NULL
reference_data$Sample10 <- NULL
reference_data$Sample16 <- NULL
reference_data$Sample17 <- NULL
reference_data$Sample18 <- NULL
reference_data$Sample19 <- NULL
reference_data$Sample20 <- NULL
reference_data$Sample21 <- NULL
reference_data$Sample22 <- NULL
reference_data$Sample31 <- NULL
reference_data$Sample32 <- NULL


#######################
#Load single cell data
#expression data (can be downloaded from GEO using the GSE149383 accession code):
D0 <- read.table("GSM3972657_D0.dge.txt", sep = "\t", header = TRUE)
D1 <- read.table("GSM3972658_D1.dge.txt", sep = "\t", header = TRUE)
D2 <- read.table("GSM3972659_D2.dge.txt", sep = "\t", header = TRUE)
D4 <- read.table("GSM3972660_D4.dge.txt", sep = "\t", header = TRUE)
D9 <- read.table("GSM3972661_D9.dge.txt", sep = "\t", header = TRUE)
D11 <- read.table("GSM3972662_D11.dge.txt", sep = "\t", header = TRUE)
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
#SEURAT ANALYSIS
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
expr.data <- data.frame(expr.data)
expr.data$hgnc_symbol <- rownames(expr.data)


##################
##Merge datasets
#################
matrix_for_combat<-merge(reference_data,expr.data,by.x='hgnc_symbol', by.y = 'hgnc_symbol')
batch<-factor(c(rep(1,26),rep(2,2949)))



#####################################
#ComBat batch effect removal analysis
#####################################
matrix_for_combat2<-as.data.frame(setDT(matrix_for_combat)[, lapply(.SD, mean), by = hgnc_symbol])
rownames(matrix_for_combat2)<-matrix_for_combat2[,1]

prepareCombat<-function(dat,batch){        
  
  dat <- as.matrix(dat)
  batch <- as.factor(batch)
  zero.rows.lst <- lapply(levels(batch), function(batch_level) {
    if (sum(batch == batch_level) > 1) {
      return(which(apply(dat[, batch == batch_level], 1,
                         function(x) {
                           var(x) == 0
                         })))
    }else {
      
      return(which(rep(1, 3) == 2))
      
    }       
  })
  zero.rows <- Reduce(union, zero.rows.lst)
  keep.rows <- setdiff(1:nrow(dat), zero.rows)
  if (length(zero.rows) > 0) {
    cat(sprintf("Found %d genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.\n",
                length(zero.rows)))
    dat.orig <- dat
    dat <- dat[keep.rows, ]
  }
}

combat_temp<-prepareCombat(dat=matrix_for_combat2[,-1],batch)
combat_temp2<- ComBat(dat=combat_temp, batch=batch)
combat_temp2[combat_temp2<0]<-0
combat_temp3<-data.frame(genes=rownames(combat_temp2),combat_temp2) #combat treated object with annotation


###########################
#Remove low variable genes
###########################
combat_temp3$genes <- as.character(combat_temp3$genes)
vardata<-apply(combat_temp3[,-1],1,var)
combat_temp3$vardata<-vardata
combat_temp3<-combat_temp3[order(combat_temp3$vardata,decreasing=T),]
combat_filter<- combat_temp3[1:round((nrow(combat_temp3)*30)/100),-ncol(combat_temp3)]





##############################################
#Principal component analysis
###############################################
####Do this on a day by day basis

#Select reference cells and cells from day 0 of the experiment
Samples <- colnames(combat_filter)
Selected_cells <- Samples[1:27]
Samples <- data.frame(Samples)
Samples$Day <- sapply(Samples$Samples, function(x)
  strsplit(x,"_")[[1]][1])
Samples <- Samples[Samples$Day %in% "D0",]
Samples <- as.character(Samples$Samples)
Selected_cells <- c(Selected_cells, Samples)
combat_filter_test <- combat_filter[,colnames(combat_filter) %in% Selected_cells]
irlba_res <- prcomp(t(combat_filter_test[,-1]))
pca.coordinates <- irlba_res$x
rownames(pca.coordinates)
##Add annotation about quiescence type:
QuiescenceType <- c("spontaneous","spontaneous","control","spontaneous","spontaneous","control","control","spontaneous",
                    "control","control","spontaneous","spontaneous","spontaneous","spontaneous","spontaneous",
                    "control","control","control","control","control","control","control","serum","serum",
                    "mek","mek","cdk","cdk","contact","contact","control","control","serum","serum","mek","mek",
                    "cdk","cdk","contact","contact")
QuiescenceType <- QuiescenceType[!(QuiescenceType %in% "control")]
dim(pca.coordinates)
QuiescenceType <- c(QuiescenceType, rep("single cell",1572))
pca.coordinates <- data.frame(pca.coordinates)
pca.coordinates$QuiescenceType <- QuiescenceType
pca.coordinates_D0 <- pca.coordinates


##Select reference cells and cells from day 0 of the experiment
Samples <- colnames(combat_filter)
Selected_cells <- Samples[1:27]
Samples <- data.frame(Samples)
Samples$Day <- sapply(Samples$Samples, function(x)
  strsplit(x,"_")[[1]][1])
Samples <- Samples[Samples$Day %in% "D1",]
Samples <- as.character(Samples$Samples)
Selected_cells <- c(Selected_cells, Samples)
combat_filter_test <- combat_filter[,colnames(combat_filter) %in% Selected_cells]
irlba_res <- prcomp(t(combat_filter_test[,-1]))
pca.coordinates <- irlba_res$x
rownames(pca.coordinates)
##Add annotation about quiescence type:
QuiescenceType <- c("spontaneous","spontaneous","control","spontaneous","spontaneous","control","control","spontaneous",
                    "control","control","spontaneous","spontaneous","spontaneous","spontaneous","spontaneous",
                    "control","control","control","control","control","control","control","serum","serum",
                    "mek","mek","cdk","cdk","contact","contact","control","control","serum","serum","mek","mek",
                    "cdk","cdk","contact","contact")
QuiescenceType <- QuiescenceType[!(QuiescenceType %in% "control")]
dim(pca.coordinates)
QuiescenceType <- c(QuiescenceType, rep("single cell",309))
pca.coordinates <- data.frame(pca.coordinates)
pca.coordinates$QuiescenceType <- QuiescenceType
pca.coordinates_D1 <- pca.coordinates


##Select reference cells and cells from day 2 of the experiment
Samples <- colnames(combat_filter)
Selected_cells <- Samples[1:27]
Samples <- data.frame(Samples)
Samples$Day <- sapply(Samples$Samples, function(x)
  strsplit(x,"_")[[1]][1])
Samples <- Samples[Samples$Day %in% "D2",]
Samples <- as.character(Samples$Samples)
Selected_cells <- c(Selected_cells, Samples)
combat_filter_test <- combat_filter[,colnames(combat_filter) %in% Selected_cells]
irlba_res <- prcomp(t(combat_filter_test[,-1]))
pca.coordinates <- irlba_res$x
rownames(pca.coordinates)
##Add annotation about quiescence type:
QuiescenceType <- c("spontaneous","spontaneous","control","spontaneous","spontaneous","control","control","spontaneous",
                    "control","control","spontaneous","spontaneous","spontaneous","spontaneous","spontaneous",
                    "control","control","control","control","control","control","control","serum","serum",
                    "mek","mek","cdk","cdk","contact","contact","control","control","serum","serum","mek","mek",
                    "cdk","cdk","contact","contact")
QuiescenceType <- QuiescenceType[!(QuiescenceType %in% "control")]
dim(pca.coordinates)
QuiescenceType <- c(QuiescenceType, rep("single cell",226))
pca.coordinates <- data.frame(pca.coordinates)
pca.coordinates$QuiescenceType <- QuiescenceType
pca.coordinates_D2 <- pca.coordinates


##Select reference cells and cells from day 4 of the experiment:
Samples <- colnames(combat_filter)
Selected_cells <- Samples[1:27]
Samples <- data.frame(Samples)
Samples$Day <- sapply(Samples$Samples, function(x)
  strsplit(x,"_")[[1]][1])
Samples <- Samples[Samples$Day %in% "D4",]
Samples <- as.character(Samples$Samples)
Selected_cells <- c(Selected_cells, Samples)
combat_filter_test <- combat_filter[,colnames(combat_filter) %in% Selected_cells]
irlba_res <- prcomp(t(combat_filter_test[,-1]))
pca.coordinates <- irlba_res$x
rownames(pca.coordinates)
##Add annotation about quiescence type:
QuiescenceType <- c("spontaneous","spontaneous","control","spontaneous","spontaneous","control","control","spontaneous",
                    "control","control","spontaneous","spontaneous","spontaneous","spontaneous","spontaneous",
                    "control","control","control","control","control","control","control","serum","serum",
                    "mek","mek","cdk","cdk","contact","contact","control","control","serum","serum","mek","mek",
                    "cdk","cdk","contact","contact")
QuiescenceType <- QuiescenceType[!(QuiescenceType %in% "control")]
dim(pca.coordinates)
QuiescenceType <- c(QuiescenceType, rep("single cell",200))
pca.coordinates <- data.frame(pca.coordinates)
pca.coordinates$QuiescenceType <- QuiescenceType
pca.coordinates_D4 <- pca.coordinates


##Select reference cells and cells from day 9 of the experiment
Samples <- colnames(combat_filter)
Selected_cells <- Samples[1:27]
Samples <- data.frame(Samples)
Samples$Day <- sapply(Samples$Samples, function(x)
  strsplit(x,"_")[[1]][1])
Samples <- Samples[Samples$Day %in% "D9",]
Samples <- as.character(Samples$Samples)
Selected_cells <- c(Selected_cells, Samples)
combat_filter_test <- combat_filter[,colnames(combat_filter) %in% Selected_cells]
irlba_res <- prcomp(t(combat_filter_test[,-1]))
pca.coordinates <- irlba_res$x
rownames(pca.coordinates)
##Add annotation about quiescence type:
QuiescenceType <- c("spontaneous","spontaneous","control","spontaneous","spontaneous","control","control","spontaneous",
                    "control","control","spontaneous","spontaneous","spontaneous","spontaneous","spontaneous",
                    "control","control","control","control","control","control","control","serum","serum",
                    "mek","mek","cdk","cdk","contact","contact","control","control","serum","serum","mek","mek",
                    "cdk","cdk","contact","contact")
QuiescenceType <- QuiescenceType[!(QuiescenceType %in% "control")]
dim(pca.coordinates)
QuiescenceType <- c(QuiescenceType, rep("single cell",389))
pca.coordinates <- data.frame(pca.coordinates)
pca.coordinates$QuiescenceType <- QuiescenceType
pca.coordinates_D9 <- pca.coordinates


##Select reference cells and cells from day 11 of the experiment
Samples <- colnames(combat_filter)
Selected_cells <- Samples[1:27]
Samples <- data.frame(Samples)
Samples$Day <- sapply(Samples$Samples, function(x)
  strsplit(x,"_")[[1]][1])
Samples <- Samples[Samples$Day %in% "D11",]
Samples <- as.character(Samples$Samples)
Selected_cells <- c(Selected_cells, Samples)
combat_filter_test <- combat_filter[,colnames(combat_filter) %in% Selected_cells]
irlba_res <- prcomp(t(combat_filter_test[,-1]))
pca.coordinates <- irlba_res$x
rownames(pca.coordinates)
##Add annotation about quiescence type:
QuiescenceType <- c("spontaneous","spontaneous","control","spontaneous","spontaneous","control","control","spontaneous",
                    "control","control","spontaneous","spontaneous","spontaneous","spontaneous","spontaneous",
                    "control","control","control","control","control","control","control","serum","serum",
                    "mek","mek","cdk","cdk","contact","contact","control","control","serum","serum","mek","mek",
                    "cdk","cdk","contact","contact")
QuiescenceType <- QuiescenceType[!(QuiescenceType %in% "control")]
dim(pca.coordinates)
QuiescenceType <- c(QuiescenceType, rep("single cell",253))
pca.coordinates <- data.frame(pca.coordinates)
pca.coordinates$QuiescenceType <- QuiescenceType
pca.coordinates_D11 <- pca.coordinates





##########################################################################
###Quiescence scoring of scRNA-seq dataset
##########################################################################

#Load quiescence biomarker genes:
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList")
load("upregulated_common.RData")
load("downregulated_common.RData")
quiescence_genes <- c(downregulated_common, upregulated_common)
expr.data$hgnc_symbol <- NULL
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






#######################################
###kNN analysis 
#######################################
#Do this on a day by day basis (only on quiescent cells)
#Select quiescent cells:
QS <- QS[QS$Common_score > 0,]
quiescent_cells <- as.character(QS$SampleID)

#D0
pca.coordinates_D0 <- pca.coordinates_D0[,colnames(pca.coordinates_D0) %in% c("QuiescenceType","PC1","PC2")]
pca_train <- pca.coordinates_D0[1:26,]
pca_train_category <- pca_train$QuiescenceType
pca_train$QuiescenceType <- NULL
dim(pca.coordinates_D0)
pca_test <- pca.coordinates_D0[27:1598,]
pca.proliferating <- pca_test[!rownames(pca_test) %in% quiescent_cells,]
pca.proliferating$QuiescenceType <- "Proliferating"
pca_test <- pca_test[rownames(pca_test) %in% quiescent_cells,]
pca_test$QuiescenceType <- NULL
pr <- knn(pca_train,pca_test,cl=pca_train_category,k=3, l = 2)
results_D0 <- pca_test
results_D0$QuiescenceType <- pr
results_D0$QuiescenceType <- as.character(results_D0$QuiescenceType)
table(results_D0$QuiescenceType)
results_D0$QuiescenceType <- sapply(results_D0$QuiescenceType, function(x)
  ifelse(x %in% "spontaneous","spontaneous",
         ifelse(x %in% "serum","serum",
                ifelse(x %in% "contact","contact",
                       ifelse(x %in% "cdk","cdk",
                              ifelse(x %in% "mek","mek","uncertain"))))))
pca_train$QuiescenceType <- pca_train_category
results_D0 <- rbind(pca_train, results_D0, pca.proliferating)
table(results_D0$QuiescenceType)




#D1
pca.coordinates_D1 <- pca.coordinates_D1[,colnames(pca.coordinates_D1) %in% c("QuiescenceType","PC1","PC2")]
pca_train <- pca.coordinates_D1[1:26,]
pca_train_category <- pca_train$QuiescenceType
pca_train$QuiescenceType <- NULL
dim(pca.coordinates_D1)
pca_test <- pca.coordinates_D1[27:335,]
pca.proliferating <- pca_test[!rownames(pca_test) %in% quiescent_cells,]
pca.proliferating$QuiescenceType <- "Proliferating"
pca_test <- pca_test[rownames(pca_test) %in% quiescent_cells,]
pca_test$QuiescenceType <- NULL
pr <- knn(pca_train,pca_test,cl=pca_train_category,k=3, l = 2)
results_D1 <- pca_test
results_D1$QuiescenceType <- pr
results_D1$QuiescenceType <- as.character(results_D1$QuiescenceType)
table(results_D1$QuiescenceType)
results_D1$QuiescenceType <- sapply(results_D1$QuiescenceType, function(x)
  ifelse(x %in% "spontaneous","spontaneous",
         ifelse(x %in% "serum","serum",
                ifelse(x %in% "contact","contact",
                       ifelse(x %in% "cdk","cdk",
                              ifelse(x %in% "mek","mek","uncertain"))))))
pca_train$QuiescenceType <- pca_train_category
results_D1 <- rbind(pca_train, results_D1, pca.proliferating)
table(results_D1$QuiescenceType)



#D2
pca.coordinates_D2 <- pca.coordinates_D2[,colnames(pca.coordinates_D2) %in% c("QuiescenceType","PC1","PC2")]
pca_train <- pca.coordinates_D2[1:26,]
pca_train_category <- pca_train$QuiescenceType
pca_train$QuiescenceType <- NULL
dim(pca.coordinates_D2)
pca_test <- pca.coordinates_D2[27:252,]
pca.proliferating <- pca_test[!rownames(pca_test) %in% quiescent_cells,]
pca.proliferating$QuiescenceType <- "Proliferating"
pca_test <- pca_test[rownames(pca_test) %in% quiescent_cells,]
pca_test$QuiescenceType <- NULL
pr <- knn(pca_train,pca_test,cl=pca_train_category,k=3, l = 2)
results_D2 <- pca_test
results_D2$QuiescenceType <- pr
results_D2$QuiescenceType <- as.character(results_D2$QuiescenceType)
table(results_D2$QuiescenceType)
results_D2$QuiescenceType <- sapply(results_D2$QuiescenceType, function(x)
  ifelse(x %in% "spontaneous","spontaneous",
         ifelse(x %in% "serum","serum",
                ifelse(x %in% "contact","contact",
                       ifelse(x %in% "cdk","cdk",
                              ifelse(x %in% "mek","mek","uncertain"))))))
pca_train$QuiescenceType <- pca_train_category
results_D2 <- rbind(pca_train, results_D2, pca.proliferating)
table(results_D2$QuiescenceType)




#D4
pca.coordinates_D4 <- pca.coordinates_D4[,colnames(pca.coordinates_D4) %in% c("QuiescenceType","PC1","PC2")]
pca_train <- pca.coordinates_D4[1:26,]
pca_train_category <- pca_train$QuiescenceType
pca_train$QuiescenceType <- NULL
dim(pca.coordinates_D4)
pca_test <- pca.coordinates_D4[27:226,]
pca.proliferating <- pca_test[!rownames(pca_test) %in% quiescent_cells,]
pca.proliferating$QuiescenceType <- "Proliferating"
pca_test <- pca_test[rownames(pca_test) %in% quiescent_cells,]
pca_test$QuiescenceType <- NULL
pr <- knn(pca_train,pca_test,cl=pca_train_category,k=3, l = 2)
results_D4 <- pca_test
results_D4$QuiescenceType <- pr
results_D4$QuiescenceType <- as.character(results_D4$QuiescenceType)
table(results_D4$QuiescenceType)
results_D4$QuiescenceType <- sapply(results_D4$QuiescenceType, function(x)
  ifelse(x %in% "spontaneous","spontaneous",
         ifelse(x %in% "serum","serum",
                ifelse(x %in% "contact","contact",
                       ifelse(x %in% "cdk","cdk",
                              ifelse(x %in% "mek","mek","uncertain"))))))
pca_train$QuiescenceType <- pca_train_category
results_D4 <- rbind(pca_train, results_D4, pca.proliferating)
table(results_D4$QuiescenceType)



#D9
pca.coordinates_D9 <- pca.coordinates_D9[,colnames(pca.coordinates_D9) %in% c("QuiescenceType","PC1","PC2")]
pca_train <- pca.coordinates_D9[1:26,]
pca_train_category <- pca_train$QuiescenceType
pca_train$QuiescenceType <- NULL
dim(pca.coordinates_D9)
pca_test <- pca.coordinates_D9[27:415,]
pca.proliferating <- pca_test[!rownames(pca_test) %in% quiescent_cells,]
pca.proliferating$QuiescenceType <- "Proliferating"
pca_test <- pca_test[rownames(pca_test) %in% quiescent_cells,]
pca_test$QuiescenceType <- NULL
pr <- knn(pca_train,pca_test,cl=pca_train_category,k=3, l = 2)
results_D9 <- pca_test
results_D9$QuiescenceType <- pr
results_D9$QuiescenceType <- as.character(results_D9$QuiescenceType)
table(results_D9$QuiescenceType)
results_D9$QuiescenceType <- sapply(results_D9$QuiescenceType, function(x)
  ifelse(x %in% "spontaneous","spontaneous",
         ifelse(x %in% "serum","serum",
                ifelse(x %in% "contact","contact",
                       ifelse(x %in% "cdk","cdk",
                              ifelse(x %in% "mek","mek","uncertain"))))))
pca_train$QuiescenceType <- pca_train_category
results_D9 <- rbind(pca_train, results_D9, pca.proliferating)
table(results_D9$QuiescenceType)


#D11
pca.coordinates_D11 <- pca.coordinates_D11[,colnames(pca.coordinates_D11) %in% c("QuiescenceType","PC1","PC2")]
pca_train <- pca.coordinates_D11[1:26,]
pca_train_category <- pca_train$QuiescenceType
pca_train$QuiescenceType <- NULL
dim(pca.coordinates_D11)
pca_test <- pca.coordinates_D11[27:279,]
pca.proliferating <- pca_test[!rownames(pca_test) %in% quiescent_cells,]
pca.proliferating$QuiescenceType <- "Proliferating"
pca_test <- pca_test[rownames(pca_test) %in% quiescent_cells,]
pca_test$QuiescenceType <- NULL
pr <- knn(pca_train,pca_test,cl=pca_train_category,k=3, l = 2)
results_D11 <- pca_test
results_D11$QuiescenceType <- pr
results_D11$QuiescenceType <- as.character(results_D11$QuiescenceType)
table(results_D11$QuiescenceType)
results_D11$QuiescenceType <- sapply(results_D11$QuiescenceType, function(x)
  ifelse(x %in% "spontaneous","spontaneous",
         ifelse(x %in% "serum","serum",
                ifelse(x %in% "contact","contact",
                       ifelse(x %in% "cdk","cdk",
                              ifelse(x %in% "mek","mek","uncertain"))))))
pca_train$QuiescenceType <- pca_train_category
results_D11 <- rbind(pca_train, results_D11, pca.proliferating)
table(results_D11$QuiescenceType)




#########################
#Combine and plot results
########################


#Combine results
results_D0$Day <- "D0"
results_D0$Sample <- rownames(results_D0)
results_D1$Day <- "D1"
results_D1$Sample <- rownames(results_D1)
results_D2$Day <- "D2"
results_D2$Sample <- rownames(results_D2)
results_D4$Day <- "D4"
results_D4$Sample <- rownames(results_D4)
results_D9$Day <- "D9"
results_D9$Sample <- rownames(results_D9)
results_D11$Day <- "D11"
results_D11$Sample <- rownames(results_D11)
results <- rbind(results_D0, results_D1, results_D11, results_D2, results_D4, results_D9)
results$Day <- factor(results$Day, levels = c("D0","D1","D2","D4","D9","D11"))
reference <- results$Sample[1:26]
results$SampleType <- sapply(results$Sample, function(x)
  ifelse(x %in% reference, "RNA-seq reference","scRNA-seq"))
results$SampleType <- factor(results$SampleType, levels = c("scRNA-seq","RNA-seq reference"))

#Plot
setwd("~/Documents/GitHub/CancerDormancy/QuiescenceTypeClassification/scRNAseq_classification/Figures/")
plot <- ggplot(results, aes(x= PC1, y=PC2, color = QuiescenceType, shape = SampleType)) + geom_point(aes(shape=SampleType, color=QuiescenceType), alpha = 0.7, size =1)
pdf("GSE134839_knn.pdf", width = 9, height = 3)
plot + scale_colour_manual(values = c("cdk" = "firebrick3", "contact" = "royalblue4", "mek" = "lightskyblue", "serum" = "goldenrod1", "spontaneous" = "mediumorchid4","uncertain" = "grey81", "Proliferating" = "mediumseagreen")) + facet_wrap(~Day,nrow = 1) + theme_classic() + scale_size_manual(values=c(1, 1))
dev.off()




##############################
##Barplot composition summary:
#############################

QuiescenceType <- c("cdk","contact","mek","serum","spontaneous","uncertain","cdk","contact","mek","serum","spontaneous","uncertain","cdk","contact","mek","serum","spontaneous","uncertain","cdk","contact","mek","serum","spontaneous","uncertain","cdk","contact","mek","serum","spontaneous","uncertain","cdk","contact","mek","serum","spontaneous","uncertain")
Days <- c(rep(0,6),rep(1,6),rep(2,6),rep(4,6),rep(9,6),rep(11,6))
N <- c(55,5,304,12,131,119,6,9,78,138,23,41,13,3,55,64,25,21,5,1,76,52,13,18,14,6,78,64,12,59,9,3,63,4,7,18)
Summary <- data.frame(Days, QuiescenceType, N)
Summary$Days <- as.factor(Summary$Days)
pdf("GSE134839_knn_barplot_cell_composition.pdf", height = 5, width = 4)
p <- ggplot(Summary, aes(fill=QuiescenceType, y=N, x=Days, width = 0.75)) + 
  geom_bar(stat = "identity", position = "fill") + theme_classic()
p + rotate_x_text(45) + scale_colour_manual(values = c("cdk" = "firebrick3", "contact" = "royalblue4", "mek" = "lightskyblue", "serum" = "goldenrod1", "spontaneous" = "mediumorchid4","uncertain" = "grey81"))

dev.off()
