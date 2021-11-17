###########################################################
#Bulk RNA-seq dataset quiescence score treatment response
###########################################################

#Load required packages
library(GSVA)
library(biomaRt)
library(plyr)
library(dplyr)
library(GeoTcgaData)
library(ggpubr)


###############################
##GGSE99116 Quiescence Scoring
###############################
#Load expression data:
count.matrix <- read.csv("GSE99116_RPKM_acute.csv", header = TRUE) #Data can be downloaded from GEO using GSE99116 accession code
rownames(count.matrix) <- count.matrix$X
count.matrix$X <- NULL
count.matrix <- data.frame(t(count.matrix))
##Add sample annotation:
Samples <- rownames(count.matrix)
Annotation <- data.frame(Samples)
Annotation$Samples <- as.character(Annotation$Samples)
Annotation$Treatment <- sapply(Annotation$Samples, function(x)
  strsplit(x,"_")[[1]][2])
Annotation <- Annotation[Annotation$Treatment %in% c("Control","Palbociclib"),]
Annotation$Treatment <- sapply(Annotation$Treatment, function(x)
  ifelse(x %in% "Control",0,1))
count.matrix <- count.matrix[rownames(count.matrix) %in% Annotation$Samples,]
#Quiescence scores
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("upregulated_common_ENSG.RData")
load("downregulated_common_ENSG.RData")
#Select only genes involved in quiescence to reduce the size of the dataframe:
expr.data <- count.matrix
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated_common_ENSG, upregulated_common_ENSG))]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
gene_lists <- list(upregulated_common_ENSG, downregulated_common_ENSG)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
all(rownames(z_score) == Annotation$Samples)
Annotation$Common <- z_score$z_score
z_score <- Annotation$Common
Treatment <- Annotation$Treatment
Treatment <- sapply(Treatment, function(x)
  ifelse(x %in% 1,"Treated","Control"))
z_score <- data.frame(z_score, Treatment)
rownames(z_score) <- Annotation$Samples
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Results/")
save(z_score, file = "GSE99116_QS.RData")





###############################
##GSE124854 Quiescence Scoring
###############################
#Load expression data:
count.data <- read.csv("GSE124854_combined_counts.csv", header = TRUE) #Data can be downloaded from the GEO using accession code GSE124854
#####TPM transformation:
Genes <- as.character(count.data$id)
#Convert gene names using biomart:
human = useMart("ensembl", dataset="hsapiens_gene_ensembl")
results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=Genes, mart=human, filters = "ensembl_gene_id")
#Merge the results from biomart with the count matrix
count.matrix <- merge(count.data, results,
                      by.x = "id", by.y = "ensembl_gene_id")
##Remove duplicate values:
count.matrix <- count.matrix %>% group_by(hgnc_symbol) %>% mutate_each(funs(mean)) %>% distinct
count.matrix <- data.frame(count.matrix)
rownames(count.matrix) <- count.matrix$hgnc_symbol
count.matrix$id <- NULL
count.matrix$hgnc_symbol <- NULL
#Convert the dataframe to a matrix
count.matrix <- as.matrix(count.matrix)
#Convert raw counts to TPM values:
TPM_expression <- countToTpm_matrix(count.matrix)
TPM_expression <- log2(TPM_expression + 1)
TPM_expression <- data.frame(t(TPM_expression))
##Add sample annotation:
Samples <- rownames(TPM_expression)
Annotation <- data.frame(Samples)
Annotation$Samples <- as.character(Annotation$Samples)
Annotation$Treatment <- sapply(Annotation$Samples, function(x)
  ifelse(x %in% c("X18L","X18R","X20R","X20L","X24R","X24L","X25R","X25L","X26L","X26R","X37R","X37L","X45L","X45R"),0,
         ifelse(x %in% c("X31L","X31R","X36L","X37R","X40L","X40R","X41R","X47L"),1,2)))
Annotation <- Annotation[Annotation$Treatment %in% c(1,0),]
TPM_expression <- TPM_expression[rownames(TPM_expression) %in% Annotation$Samples,]
#Quiescence score calculation
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("upregulated_common.RData")
load("downregulated_common.RData")
#Select only genes involved in quiescence to reduce the size of the dataframe:
expr.data <- TPM_expression
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated_common, upregulated_common))]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
gene_lists <- list(upregulated_common, downregulated_common)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
all(rownames(z_score) == Annotation$Samples)
Annotation$Common <- z_score$z_score
z_score <- Annotation$Common
Treatment <- Annotation$Treatment
Treatment <- sapply(Treatment, function(x)
  ifelse(x %in% 1,"Treated","Control"))
z_score <- data.frame(z_score, Treatment)
rownames(z_score) <- Annotation$Samples
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Results/")
save(z_score, file = "GSE124854_QS.RData")




###############################
##GSE135215 Quiescence Scoring
###############################
#Load expression data:
count.matrix <- read.csv("GSE135215_ALLsamples.normalizedCounts.csv", header = TRUE) #Data can be downloaded from the GEO using accession code: GSE135215 
count.matrix <- count.matrix %>% group_by(Name) %>% mutate_each(funs(mean)) %>% distinct
count.matrix <- data.frame(count.matrix)
rownames(count.matrix) <- count.matrix$Name
count.matrix$Name <- NULL
count.matrix <- as.matrix(count.matrix)
count.matrix <- log2(count.matrix +1)
count.matrix <- data.frame(t(count.matrix))
##Add sample annotation:
Samples <- rownames(count.matrix)
Annotation <- data.frame(Samples)
Annotation$Samples <- as.character(Annotation$Samples)
Annotation$Treatment <- sapply(Annotation$Samples, function(x)
  strsplit(x,"_")[[1]][3])
Annotation <- Annotation[Annotation$Treatment %in% c("DMSO","Palbo"),]
Annotation$Treatment <- sapply(Annotation$Treatment, function(x)
  ifelse(x %in% "DMSO",0,1))
count.matrix <- count.matrix[rownames(count.matrix) %in% Annotation$Samples,]
#Calculate quiescence scores
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("upregulated_common.RData")
load("downregulated_common.RData")
#Select only genes involved in quiescence to reduce the size of the dataframe:
expr.data <- count.matrix
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated_common, upregulated_common))]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
gene_lists <- list(upregulated_common, downregulated_common)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
all(rownames(z_score) == Annotation$Samples)
Annotation$Common <- z_score$z_score
z_score <- Annotation$Common
Treatment <- Annotation$Treatment
Treatment <- sapply(Treatment, function(x)
  ifelse(x %in% 1,"Treated","Control"))
z_score <- data.frame(z_score, Treatment)
rownames(z_score) <- Annotation$Samples
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Results/")
save(z_score, file = "GSE135215_QS.RData")





###############################
##GSE152699 Quiescence Scoring
###############################
#Load expression data
expr.data <- read.table(file = "Data/GSE152699_Cells_fpkm.txt", sep = "\t", header = TRUE) #Data can be downloaded from the GEO using accession code GSE152699
expr.data <- expr.data %>% group_by(Name) %>% mutate_each(funs(mean)) %>% distinct
expr.data <- data.frame(expr.data)
rownames(expr.data) <- expr.data$Name
expr.data$Name <- NULL
expr.data[1:1,1:10]
X <- as.matrix(expr.data)
X <- log2(X+1)
#Add annotation information
Samples <- colnames(expr.data)
binary_groups <- c(0,0,0,0,0,0,1,1,1,1,1,1)
annotation <- data.frame(Samples, binary_groups)
#Quiescence score calculation
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("upregulated_common.RData")
load("downregulated_common.RData")
gene_lists <- list(upregulated_common, downregulated_common)
es.dif <- gsva(X, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
es.dif <- t(es.dif)
es.dif <- data.frame(es.dif)
es.dif$GSVA_score <- es.dif$X1 - es.dif$X2
z_score <- es.dif$GSVA_score
z_score <- data.frame(z_score)
z_score$Treatment <- c(rep("Control",6),rep("Treated",6))
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Results/")
save(z_score, file = "GSE152699_QS.RData")






###############################
##GSE178839 Quiescence Scoring
###############################
#Load expression data:
expr.data <- read.table("GSE178839_Combined_counts_HTSeq.csv", header = TRUE,sep = ";") #Data can be downloaded from the GEO using accession code GSE178839
#####TPM transformation:
Genes <- as.character(expr.data$GeneID)
#Convert gene names using biomart:
human = useMart("ensembl", dataset="hsapiens_gene_ensembl")
results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=Genes, mart=human, filters = "ensembl_gene_id")
#Merge the results from biomart with the count matrix
count.matrix <- merge(expr.data, results,
                      by.x = "GeneID", by.y = "ensembl_gene_id")
##Remove duplicate values:
count.matrix$GeneID <- NULL
count.matrix <- count.matrix %>% group_by(hgnc_symbol) %>% mutate_each(funs(mean)) %>% distinct
count.matrix <- data.frame(count.matrix)
rownames(count.matrix) <- count.matrix$hgnc_symbol
count.matrix$hgnc_symbol <- NULL
#Convert the dataframe to a matrix
count.matrix <- as.matrix(count.matrix)
#Convert raw counts to TPM values:
TPM_expression <- countToTpm_matrix(count.matrix)
TPM_expression <- log2(TPM_expression + 1)
TPM_expression <- data.frame(t(TPM_expression))
#Quiescence Score calculation
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("upregulated_common.RData")
load("downregulated_common.RData")
#Select only genes involved in quiescence to reduce the size of the dataframe:
expr.data <- TPM_expression
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated_common, upregulated_common))]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
gene_lists <- list(upregulated_common, downregulated_common)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
z_score$X1 <- NULL
z_score$X2 <- NULL
#Add treatment information:
z_score$Treatment <- c("Treated","Treated","Control","Treated","Control","Control","Treated","Treated","Treated")
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Results/")
save(z_score, file = "GSE178839_QS.RData")





########################################
#Summary quiescence score comparison
########################################

#Load and combine quiescence scores from the 5 datasets

setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Results/")
load("GSE124854_QS.RData")
z_score$Study <- "GSE124854"
combined.data <- z_score
rm(z_score)
load("GSE135215_QS.RData")
z_score$Study <- "GSE135215"
combined.data <- rbind(combined.data, z_score)
load("GSE152699_QS.RData")
z_score$Study <- "GSE152699"
combined.data <- rbind(combined.data, z_score)
load("GSE178839_QS.RData")
z_score$Study <- "GSE178839"
combined.data <- rbind(combined.data, z_score)
load("GSE99116_QS.RData")
z_score$Study <- "GSE99116"
combined.data <- rbind(combined.data, z_score)


###Plot:
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Figures/")
pdf("BulkRNAseq_TreatmentResponse.pdf",height = 6, width = 7)


combined.data %>%
  ggplot( aes(x=Study, y=z_score, fill=Treatment)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(), size = 0.1) +
  theme_classic() +
  xlab("") + stat_compare_means(aes(group = Treatment), label = "p.signif") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


dev.off()





