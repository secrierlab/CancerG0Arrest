##############################################
#QS vs Proliferation marker performance
##############################################

#Load the required packages:
library(dplyr)
library(ggplot2)
library(plotROC)
library(RColorBrewer)
library(biomaRt)


#######################################################################################
#Load expression data of the 7 validation datasets for the given proliferation markers
#######################################################################################



##############
#GSE131594

#Load expression data
expr.data <- read.table(file = "GSE131594_RNAseq_DormantCells_FPKM.txt", sep = "\t", header = TRUE) #Can be obtained from GEO using accession code GSE131594
expr.data <- expr.data %>% group_by(gene.symbol) %>% mutate_each(funs(mean)) %>% distinct
expr.data <- data.frame(expr.data)
rownames(expr.data) <- expr.data$gene.symbol
expr.data$gene.symbol <- NULL
expr.data$gene_id <- NULL
expr.data[1:10,1:10]
X <- as.matrix(expr.data)
X <- log2(X+1)
X <- data.frame(t(X))
#Add annotation information
Samples <- colnames(expr.data)
binary_groups <- c(1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0)
annotation <- data.frame(Samples, binary_groups)
#Results:
X[1:10,1:10]
rownames(annotation) <- annotation$Samples
annotation$Samples  <- NULL
all(rownames(annotation) == rownames(X))
results <- annotation
#Add information about the genes of interest to the results data.frame
results$CDKN1A <- X$CDKN1A
results$MKI67 <- X$MKI67
results$CCNB1 <- X$CCNB1
results$CDK2 <- X$CDK4
results$CDK4 <- X$CDK4
results$CDK6 <- X$CDK6
results$CCNE1 <- X$CCNE1
results$RBL2 <- X$RBL2
results$CDKN1A <- X$CDKN1A
#Mean expression of DREAM targets
setwd("~/Documents/GitHub/CancerDormancy/Data/OtherGeneLists/")
DREAM_targets <- read.table(file = "DREAM_target_genelist.txt", header = FALSE, sep = "\t")
DREAM_targets <- as.character(DREAM_targets$V1)
all_genes <- colnames(X)
DREAM_targets <- DREAM_targets[DREAM_targets %in% all_genes]
results$"Mean DREAM target expression" <- rowMeans(subset(X, select = DREAM_targets), na.rm = TRUE)
#Replication complex genes (mcm2-7)
replication_complex_genes <- c("MCMC2","MCM3","MCM4","MCM5","MCM6","MCM7")
replication_complex_genes <- replication_complex_genes[replication_complex_genes %in% all_genes]
results$"Mean RC gene expression" <- rowMeans(subset(X, select = replication_complex_genes), na.rm = TRUE)
GSE131594 <- results
GSE131594$Dataset <- "GSE131594"





##############
#GSE114012

#Load expression data
setwd("~/Documents/GitHub/CancerDormancy/Data/NormalisedRNAseqDatasets/")
load("GSE114012_TPM_expression.RData")
load("GSE114012_sample_annotations.RData")
X <- as.matrix(TPM_expression)
X <- log2(X+1)
X <- data.frame(t(X))
all_genes <- as.vector(colnames(X))
#Create a results dataframe
binary_groups <- sapply(coldata$Group, function(x)
  ifelse(x %in% c("Quiescent"),1,0))
results <- data.frame(binary_groups)
#Add information about the genes of interest to the results data.frame
results$CDKN1A <- X$CDKN1A
results$MKI67 <- X$MKI67
results$CCNB1 <- X$CCNB1
results$CDK2 <- X$CDK4
results$CDK4 <- X$CDK4
results$CDK6 <- X$CDK6
results$CCNE1 <- X$CCNE1
results$RBL2 <- X$RBL2
results$CDKN1A <- X$CDKN1A
#Mean expression of DREAM targets
setwd("~/Documents/GitHub/CancerDormancy/Data/OtherGeneLists/")
DREAM_targets <- read.table(file = "DREAM_target_genelist.txt", header = FALSE, sep = "\t")
DREAM_targets <- as.character(DREAM_targets$V1)
all_genes <- colnames(X)
DREAM_targets <- DREAM_targets[DREAM_targets %in% all_genes]
results$"Mean DREAM target expression" <- rowMeans(subset(X, select = DREAM_targets), na.rm = TRUE)
#Replication complex genes (mcm2-7)
replication_complex_genes <- c("MCMC2","MCM3","MCM4","MCM5","MCM6","MCM7")
replication_complex_genes <- replication_complex_genes[replication_complex_genes %in% all_genes]
results$"Mean RC gene expression" <- rowMeans(subset(X, select = replication_complex_genes), na.rm = TRUE)
GSE114012 <- results
GSE114012$Dataset <- "GSE114012"
rm(results)




##############
#GSE137912

#Load expression data
X <- read.csv("GSE137912_logcounts.csv",header = TRUE) #Data cam be obtained from GEO using accession code GSE137912
X <- as.matrix(normalised_imputed_data)
all_genes <- as.character(X$X)
X$X <- NULL
rownames(X) <- all_genes
##Load annotation data
annotation <- read.csv("GSE137912_cell.annotation.csv",header = TRUE) #Data can be obtained from GEO using accession code GSE137912
annotation_H358 <- annotation[annotation$Line %in% "H358" & annotation$Hours %in% c(0,24),]
##Select an expression matrix with the desired samples
samples_H358 <- as.character(annotation_H358$X)
H358_expression <- X[,colnames(X) %in% samples_H358]
#Create a results dataframe
results <- data.frame(samples_H358)
quiescent_samples <- annotation_H358[annotation_H358$Hours %in% 24,]
quiescent_samples <- quiescent_samples$X
quiescent_samples <- as.character(quiescent_samples)
results$binary_groups <- sapply(results$samples_H358, function(x)
  ifelse(x %in% quiescent_samples,1,0))
rownames(results) <- results$samples_H358
results$samples_H358 <- NULL
H358_expression <- data.frame(t(H358_expression))
X <- H358_expression
#Results:
X[1:10,1:10]
all(rownames(X) == rownames(results))
for (i in colnames(X)) {
  X[[i]] <- as.numeric(as.character(X[[i]]))
}
#Add information about the genes of interest to the results data.frame
results$CDKN1A <- X$CDKN1A
results$MKI67 <- X$MKI67
results$CCNB1 <- X$CCNB1
results$CDK2 <- X$CDK4
results$CDK4 <- X$CDK4
results$CDK6 <- X$CDK6
results$CCNE1 <- X$CCNE1
results$RBL2 <- X$RBL2
results$CDKN1A <- X$CDKN1A
#Mean expression of DREAM targets
setwd("~/Documents/GitHub/CancerDormancy/Data/OtherGeneLists/")
DREAM_targets <- read.table(file = "DREAM_target_genelist.txt", header = FALSE, sep = "\t")
DREAM_targets <- as.character(DREAM_targets$V1)
all_genes <- colnames(X)
DREAM_targets <- DREAM_targets[DREAM_targets %in% all_genes]
results$"Mean DREAM target expression" <- rowMeans(subset(X, select = DREAM_targets), na.rm = TRUE)
#Replication complex genes (mcm2-7)
replication_complex_genes <- c("MCMC2","MCM3","MCM4","MCM5","MCM6","MCM7")
replication_complex_genes <- replication_complex_genes[replication_complex_genes %in% all_genes]
results$"Mean RC gene expression" <- rowMeans(subset(X, select = replication_complex_genes), na.rm = TRUE)
GSE137912 <- results
GSE137912$Dataset <- "GSE137912"
rm(results)




##############
#GSE152699

#Load expression data
expr.data <- read.table(file = "Data/GSE152699_Cells_fpkm.txt", sep = "\t", header = TRUE) #Data can be downloaded from GEO using accession code GSE152699
expr.data <- expr.data %>% group_by(Name) %>% mutate_each(funs(mean)) %>% distinct
expr.data <- data.frame(expr.data)
rownames(expr.data) <- expr.data$Name
expr.data$Name <- NULL
expr.data[1:1,1:10]
X <- as.matrix(expr.data)
X <- log2(X+1)
X <- data.frame(t(X))
#Add annotation information
Samples <- colnames(expr.data)
binary_groups <- c(0,0,0,0,0,0,1,1,1,1,1,1)
annotation <- data.frame(Samples, binary_groups)
#Results:
X[1:10,1:10]
rownames(annotation) <- annotation$Samples
annotation$Samples  <- NULL
all(rownames(annotation) == rownames(X))
results <- annotation
#Add information about the genes of interest to the results data.frame
results$CDKN1A <- X$CDKN1A
results$MKI67 <- X$MKI67
results$CCNB1 <- X$CCNB1
results$CDK2 <- X$CDK4
results$CDK4 <- X$CDK4
results$CDK6 <- X$CDK6
results$CCNE1 <- X$CCNE1
results$RBL2 <- X$RBL2
results$CDKN1A <- X$CDKN1A
#Mean expression of DREAM targets
setwd("~/Documents/GitHub/CancerDormancy/Data/OtherGeneLists/")
DREAM_targets <- read.table(file = "DREAM_target_genelist.txt", header = FALSE, sep = "\t")
DREAM_targets <- as.character(DREAM_targets$V1)
all_genes <- colnames(X)
DREAM_targets <- DREAM_targets[DREAM_targets %in% all_genes]
results$"Mean DREAM target expression" <- rowMeans(subset(X, select = DREAM_targets), na.rm = TRUE)
#Replication complex genes (mcm2-7)
replication_complex_genes <- c("MCMC2","MCM3","MCM4","MCM5","MCM6","MCM7")
replication_complex_genes <- replication_complex_genes[replication_complex_genes %in% all_genes]
results$"Mean RC gene expression" <- rowMeans(subset(X, select = replication_complex_genes), na.rm = TRUE)
GSE152699 <- results
GSE152699$Dataset <- "GSE152699"
rm(results)






##############
#GSE75367

#Load expression data
setwd("~/Documents/GitHub/CancerDormancy/Data/NormalisedScRNAseqDatasets/")
load("GSE75367_normalised_counts.RData")
load("GSE75367_annotation.RData")
X <- as.matrix(normalised_counts)
X <- data.frame(t(X))
#Create the results dataframe
results <- anno
results$binary_groups <- sapply(results$group, function(x)
  ifelse(x %in% c("low"),1,0))
results$group <- NULL
results$sample <- rownames(results)
#make sure all the samples in the results dataframe are also in the X dataframe
samples <- rownames(X)
results <- results[rownames(results) %in% samples,]
all(rownames(results) %in% rownames(X))
results$sample <- NULL
#Add information about the genes of interest to the results data.frame
results$CDKN1A <- X$CDKN1A
results$MKI67 <- X$MKI67
results$CCNB1 <- X$CCNB1
results$CDK2 <- X$CDK4
results$CDK4 <- X$CDK4
results$CDK6 <- X$CDK6
results$CCNE1 <- X$CCNE1
results$RBL2 <- X$RBL2
results$CDKN1A <- X$CDKN1A
#Mean expression of DREAM targets
setwd("~/Documents/GitHub/CancerDormancy/Data/OtherGeneLists/")
DREAM_targets <- read.table(file = "DREAM_target_genelist.txt", header = FALSE, sep = "\t")
DREAM_targets <- as.character(DREAM_targets$V1)
all_genes <- colnames(X)
DREAM_targets <- DREAM_targets[DREAM_targets %in% all_genes]
results$"Mean DREAM target expression" <- rowMeans(subset(X, select = DREAM_targets), na.rm = TRUE)
#Replication complex genes (mcm2-7)
replication_complex_genes <- c("MCMC2","MCM3","MCM4","MCM5","MCM6","MCM7")
replication_complex_genes <- replication_complex_genes[replication_complex_genes %in% all_genes]
results$"Mean RC gene expression" <- rowMeans(subset(X, select = replication_complex_genes), na.rm = TRUE)
GSE75367 <- results
GSE75367$Dataset <- "GSE75367"





##############
#GSE83142

#Load expression data
setwd("~/Documents/GitHub/CancerDormancy/Data/NormalisedScRNAseqDatasets/")
load("GSE83142_normalised_counts.RData")
X <- as.matrix(normalised_counts)
#Form the results dataframe
binary_groups <- c(rep(0,34), rep(1,14))
results <- data.frame(binary_groups)
rownames(results) <- colnames(X)
X <- data.frame(t(X))
#Add information about the genes of interest to the results data.frame
results$CDKN1A <- X$ENSG00000124762
results$MKI67 <- X$ENSG00000148773
results$CCNB1 <- X$ENSG00000134057
results$CDK2 <- X$ENSG00000123374
results$CDK4 <- X$ENSG00000135446
results$CDK6 <- X$ENSG00000105810
results$CCNE1 <- X$ENSG00000105173
results$RBL2 <- X$ENSG00000103479
results$CDKN1A <- X$ENSG00000124762
#Mean expression of DREAM targets
setwd("~/Documents/GitHub/CancerDormancy/Data/OtherGeneLists/")
DREAM_targets <- read.table(file = "DREAM_target_genelist.txt", header = FALSE, sep = "\t")
DREAM_targets <- as.character(DREAM_targets$V1)
all_genes <- colnames(X)
#Replication complex genes (mcm2-7)
replication_complex_genes <- c("ENSG00000073111", "ENSG00000112118","ENSG00000104738","ENSG00000100297","ENSG00000076003","ENSG00000166508")
replication_complex_genes <- replication_complex_genes[replication_complex_genes %in% all_genes]
results$"Mean RC gene expression" <- rowMeans(subset(X, select = replication_complex_genes), na.rm = TRUE)
#Convert using biomart
human = useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "useast.ensembl.org")
results_biomart <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=DREAM_targets, mart=human, filters = "hgnc_symbol")
DREAM_targets <- results_biomart$ensembl_gene_id
DREAM_targets <- DREAM_targets[DREAM_targets %in% all_genes]
results$"Mean DREAM target expression" <- rowMeans(subset(X, select = DREAM_targets), na.rm = TRUE)
GSE83142 <- results
GSE83142$Dataset <- "GSE83142"
rm(results)







##############
#GSE93991
expr.data <- read.table(file = "GSE93991_RNA-seq_Expression_Data.txt", sep = "\t", header = FALSE) #data can be downloaded from GEO using accession code GSE93991
rownames(expr.data) <- expr.data$V1
expr.data$V1 <- NULL
#Expression data starts from V10
expr.data$V2 <- NULL
expr.data$V3 <- NULL
expr.data$V4 <- NULL
expr.data$V5 <- NULL
expr.data$V6 <- NULL
expr.data$V7 <- NULL
expr.data$V8 <- NULL
expr.data$V9 <- NULL
#Add colnames:
Sample <- c("Quiescence1", "Quiescence2", "Quiescence3", "Proliferation1", "Quiescence4", "Proliferation2", "Quiescence5", "Proliferation3","Proliferation4", "Proliferation5", "Proliferation6", "proliferation7", "Proliferation8", "Proliferation9", "Quiescence6" )
colnames(expr.data) <- Sample
X <- as.matrix(expr.data)
all_genes <- as.vector(rownames(X))
X <- data.frame(t(X))
#Create a results dataframe:
binary_groups <- c(1,1,1,0,1,0,1,0,0,0,0,0,0,0,1)
results <- as.data.frame(binary_groups)
#Add information about the genes of interest to the results data.frame
results$CDKN1A <- X$ENSG00000124762
results$MKI67 <- X$ENSG00000148773
results$CCNB1 <- X$ENSG00000134057
results$CDK2 <- X$ENSG00000123374
results$CDK4 <- X$ENSG00000135446
results$CDK6 <- X$ENSG00000105810
results$CCNE1 <- X$ENSG00000105173
results$RBL2 <- X$ENSG00000103479
results$CDKN1A <- X$ENSG00000124762
#Mean expression of DREAM targets
setwd("~/Documents/GitHub/CancerDormancy/Data/OtherGeneLists/")
DREAM_targets <- read.table(file = "DREAM_target_genelist.txt", header = FALSE, sep = "\t")
DREAM_targets <- as.character(DREAM_targets$V1)
all_genes <- colnames(X)
#Replication complex genes (mcm2-7)
replication_complex_genes <- c("ENSG00000073111", "ENSG00000112118","ENSG00000104738","ENSG00000100297","ENSG00000076003","ENSG00000166508")
replication_complex_genes <- replication_complex_genes[replication_complex_genes %in% all_genes]
results$"Mean RC gene expression" <- rowMeans(subset(X, select = replication_complex_genes), na.rm = TRUE)
#Convert using biomart
human = useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "useast.ensembl.org")
results_biomart <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=DREAM_targets, mart=human, filters = "hgnc_symbol")
DREAM_targets <- results_biomart$ensembl_gene_id
DREAM_targets <- DREAM_targets[DREAM_targets %in% all_genes]
results$"Mean DREAM target expression" <- rowMeans(subset(X, select = DREAM_targets), na.rm = TRUE)
GSE93991 <- results
GSE93991$Dataset <- "GSE93991"
rm(results)




##################################
#Combine expression data summaries
##################################
markers <- colnames(GSE137912)
GSE114012 <- GSE114012[,colnames(GSE114012) %in% markers]
GSE131594 <- GSE131594[,colnames(GSE131594) %in% markers]
GSE152699 <- GSE152699[,colnames(GSE152699) %in% markers]
GSE75367 <- GSE75367[,colnames(GSE75367) %in% markers]
GSE93991 <- GSE93991[,colnames(GSE93991) %in% markers]
merged_data <- rbind(GSE114012,GSE131594,GSE137912,GSE152699,GSE75367,GSE83142,GSE93991)



####################################
##Check the accuracy of each marker
####################################
colnames(merged_data)
#p21
basicplot <- ggplot(merged_data, aes(d = binary_groups, m = CDKN1A, color = Dataset)) + geom_roc(n.cuts = 10, labels = FALSE)  
basicplot + labs(x = "False positive fraction", y = "True positive fraction")
#AUC calculations
AUC <- data.frame(calc_auc(basicplot))
AUC <- AUC$AUC
Accuracy <- AUC
CDKN1A <- data.frame(Accuracy)
CDKN1A$Marker <- "CDKN1A"

#MKI67
basicplot <- ggplot(merged_data, aes(d = binary_groups, m = MKI67, color = Dataset)) + geom_roc(n.cuts = 10, labels = FALSE) 
basicplot + labs(x = "False positive fraction", y = "True positive fraction")
#AUC calculations
AUC <- data.frame(calc_auc(basicplot))
AUC <- AUC$AUC
1 - AUC
Accuracy <- 1 - AUC
MKI67 <- data.frame(Accuracy)
MKI67$Marker <- "MKI67"

#CCNB1
basicplot <- ggplot(merged_data, aes(d = binary_groups, m = CCNB1, color = Dataset)) + geom_roc(n.cuts = 10, labels = FALSE) 
basicplot + labs(x = "False positive fraction", y = "True positive fraction")
#AUC calculations
AUC <- data.frame(calc_auc(basicplot))
AUC <- AUC$AUC
1 - AUC
#Combine results for the differnt markers
Accuracy <- 1 - AUC
CCNB1 <- data.frame(Accuracy)
CCNB1$Marker <- "CCNB1"

#CDK2
basicplot <- ggplot(merged_data, aes(d = binary_groups, m = CDK2, color = Dataset)) + geom_roc(n.cuts = 10, labels = FALSE) 
basicplot + labs(x = "False positive fraction", y = "True positive fraction")
#AUC calculations
AUC <- data.frame(calc_auc(basicplot))
AUC <- AUC$AUC
1 - AUC
Accuracy <- 1 - AUC
CDK2 <- data.frame(Accuracy)
CDK2$Marker <- "CDK2"

#CDK4
basicplot <- ggplot(merged_data, aes(d = binary_groups, m = CDK4, color = Dataset)) + geom_roc(n.cuts = 10, labels = FALSE)
basicplot + labs(x = "False positive fraction", y = "True positive fraction")
#AUC calculations
AUC <- data.frame(calc_auc(basicplot))
AUC <- AUC$AUC
1 - AUC
Accuracy <- 1 - AUC
CDK4 <- data.frame(Accuracy)
CDK4$Marker <- "CDK4"

#CDK6
basicplot <- ggplot(merged_data, aes(d = binary_groups, m = CDK6, color = Dataset)) + geom_roc(n.cuts = 10, labels = FALSE) 
basicplot + labs(x = "False positive fraction", y = "True positive fraction")
#AUC calculations
AUC <- data.frame(calc_auc(basicplot))
AUC <- AUC$AUC
1 - AUC
Accuracy <- 1 - AUC
CDK6 <- data.frame(Accuracy)
CDK6$Marker <- "CDK6"

#CCNE1
basicplot <- ggplot(merged_data, aes(d = binary_groups, m = CCNE1, color = Dataset)) + geom_roc(n.cuts = 10, labels = FALSE) 
basicplot + labs(x = "False positive fraction", y = "True positive fraction")
#AUC calculations
AUC <- data.frame(calc_auc(basicplot))
AUC <- AUC$AUC
1 - AUC
#Combine results for the differnt markers
#Mean dream expression -DONE
Accuracy <- 1 - AUC
CCNE1 <- data.frame(Accuracy)
CCNE1$Marker <- "CCNE1"

#RBL2
basicplot <- ggplot(merged_data, aes(d = binary_groups, m = RBL2, color = Dataset)) + geom_roc(n.cuts = 10, labels = FALSE) 
basicplot + labs(x = "False positive fraction", y = "True positive fraction")
#AUC calculations
AUC <- data.frame(calc_auc(basicplot))
AUC <- AUC$AUC
AUC
#Combine results for the differnt markers
#Mean dream expression -DONE
Accuracy <- AUC
RBL2 <- data.frame(Accuracy)
RBL2$Marker <- "RBL2"

#DREAM
basicplot <- ggplot(merged_data, aes(d = binary_groups, m = `Mean DREAM target expression`, color = Dataset)) + geom_roc(n.cuts = 10, labels = FALSE) 
basicplot + labs(x = "False positive fraction", y = "True positive fraction")
#AUC calculations
AUC <- data.frame(calc_auc(basicplot))
AUC <- AUC$AUC
1 - AUC
#Combine results for the differnt markers
#Mean dream expression -DONE
Accuracy <- 1 - AUC
DREAM <- data.frame(Accuracy)
DREAM$Marker <- "Mean DREAM target expression"

#RC
basicplot <- ggplot(merged_data, aes(d = binary_groups, m = `Mean RC gene expression`, color = Dataset)) + geom_roc(n.cuts = 10, labels = FALSE)  
basicplot + labs(x = "False positive fraction", y = "True positive fraction")
#AUC calculations
AUC <- data.frame(calc_auc(basicplot))
AUC <- AUC$AUC
1 - AUC
#Combine results for the differnt markers
#Mean dream expression -DONE
Accuracy <- 1 - AUC
RC <- data.frame(Accuracy)
RC$Marker <- "Mean RC gene expression"

#Combined Z-score
Accuracy <- c(0.7656250,0.9652778,0.9861261,1.0000000,0.9275362,0.9222689,0.8148148) #These are determined in the "QuiescenceScoreValidation_ROC.R" script
combined_z_score <- data.frame(Accuracy)
combined_z_score$Marker <- "Combined Z-Score"

#################################
#Combine the dataframes and plot
################################
combined_data <- rbind(combined_z_score,DREAM,RC,CCNB1,MKI67,CDK4, CDK2,CCNE1,CDKN1A,CDK6,RBL2)
class(combined_data$Marker)
order <- unique(combined_data$Marker)
order <- rev(order)
#Then turn it back into a factor with the levels in the correct order
combined_data$Marker <- factor(combined_data$Marker, levels=order)
myCol <- brewer.pal(11, "Spectral")
setwd("~/Documents/GitHub/CancerDormancy/QuiescenceScoreValidation/Figures/")
p <- ggplot(combined_data, aes(x=Marker, y=Accuracy, fill=Marker)) +
  geom_boxplot() +
  scale_fill_manual(values = myCol) +
  geom_jitter(color="black", size=0.4, alpha=0.9) + theme_classic() +
  ggtitle("") +
  xlab("")
pdf("QS_vs_ProliferationMarkers_Performance.pdf", height = 5, width = 6)
p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()





