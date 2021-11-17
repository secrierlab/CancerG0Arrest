############################################################################
#Bulk RNA-seq dataset quiescence score palbociclib treatment response
###########################################################################

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
#Common score:
#Load the genes involved in quiescence
load("upregulated_common_ENSG.RData")
load("downregulated_common_ENSG.RData")
expr.data <- count.matrix
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated_common_ENSG, upregulated_common_ENSG))]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
#Now work out the combined z_score for the samples
gene_lists <- list(upregulated_common_ENSG, downregulated_common_ENSG)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
all(rownames(z_score) == Annotation$Samples)
Annotation$Common <- z_score$z_score
#CDK inhibition
load("cdk_upregulated_ENSG.RData")
load("cdk_downregulated_ENSG.RData")
expr.data <- count.matrix
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated_ENSG, upregulated_ENSG))]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
#Now work out the combined z_score for the samples
gene_lists <- list(upregulated_ENSG, downregulated_ENSG)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
all(rownames(z_score) == Annotation$Samples)
Annotation$CDK <- z_score$z_score
#Contact inhibition
load("contact_upregulated_ENSG.RData")
load("contact_downregulated_ENSG.RData")
expr.data <- count.matrix
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated_ENSG, upregulated_ENSG))]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
#Now work out the combined z_score for the samples
gene_lists <- list(upregulated_ENSG, downregulated_ENSG)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
all(rownames(z_score) == Annotation$Samples)
Annotation$ContactInhibition <- z_score$z_score
#MEK inhibition
load("mek_downregulated_ENSG.RData")
load("mek_upregulated_ENSG.RData")
expr.data <- count.matrix
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated_ENSG, upregulated_ENSG))]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
#Now work out the combined z_score for the samples
gene_lists <- list(upregulated_ENSG, downregulated_ENSG)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
all(rownames(z_score) == Annotation$Samples)
Annotation$MEKInhibition <- z_score$z_score
#serum starvation
load("serum_downregulated_ENSG.RData")
load("serum_upregulated_ENSG.RData")
expr.data <- count.matrix
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated_ENSG, upregulated_ENSG))]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
#Now work out the combined z_score for the samples
gene_lists <- list(upregulated_ENSG, downregulated_ENSG)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
all(rownames(z_score) == Annotation$Samples)
Annotation$SerumStarvation <- z_score$z_score
#spontaneous quiescence
load("spontaneous_upregulated_ENSG.RData")
load("spontaneous_downregulated_ENSG.RData")
expr.data <- count.matrix
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated_ENSG, upregulated_ENSG))]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
#Now work out the combined z_score for the samples
gene_lists <- list(upregulated_ENSG, downregulated_ENSG)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
all(rownames(z_score) == Annotation$Samples)
Annotation$Spontaneous <- z_score$z_score
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Results/")
save(Annotation, file = "GSE99116_IndividualProgramme_QS.RData")




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
#CDK inhibition
load("cdk_upregulated.RData")
load("cdk_downregulated.RData")
#Select only genes involved in quiescence to reduce the size of the dataframe:
expr.data <- TPM_expression
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated, upregulated))]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
#Now work out the combined z_score for the samples
gene_lists <- list(upregulated, downregulated)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
all(rownames(z_score) == Annotation$Samples)
Annotation$CDK <- z_score$z_score
#Contact inhibition
load("contact_upregulated.RData")
load("contact_downregulated.RData")
#Select only genes involved in quiescence to reduce the size of the dataframe:
expr.data <- TPM_expression
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated, upregulated))]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
#Now work out the combined z_score for the samples
gene_lists <- list(upregulated, downregulated)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
all(rownames(z_score) == Annotation$Samples)
Annotation$ContactInhibition <- z_score$z_score
#MEK inhibition
load("mek_downregulated.RData")
load("mek_upregulated.RData")
#Select only genes involved in quiescence to reduce the size of the dataframe:
expr.data <- TPM_expression
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated, upregulated))]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
#Now work out the combined z_score for the samples
gene_lists <- list(upregulated, downregulated)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
all(rownames(z_score) == Annotation$Samples)
Annotation$MEKInhibition <- z_score$z_score
#serum starvation
load("serum_downregulated.RData")
load("serum_upregulated.RData")
#Select only genes involved in quiescence to reduce the size of the dataframe:
expr.data <- TPM_expression
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated, upregulated))]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
#Now work out the combined z_score for the samples
gene_lists <- list(upregulated, downregulated)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
all(rownames(z_score) == Annotation$Samples)
Annotation$SerumStarvation <- z_score$z_score
#spontaneous quiescence
load("spontaneous_upregulated.RData")
load("spontaneous_downregulated.RData")
#Select only genes involved in quiescence to reduce the size of the dataframe:
expr.data <- TPM_expression
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated, upregulated))]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
#Now work out the combined z_score for the samples
gene_lists <- list(upregulated, downregulated)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
all(rownames(z_score) == Annotation$Samples)
Annotation$Spontaneous <- z_score$z_score
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Results/")
save(Annotation, file = "GSE124854_IndividualProgramme_QS.RData")




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
#CDK inhibition
load("cdk_upregulated.RData")
load("cdk_downregulated.RData")
#Select only genes involved in quiescence to reduce the size of the dataframe:
expr.data <- count.matrix
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated, upregulated))]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
#Now work out the combined z_score for the samples
gene_lists <- list(upregulated, downregulated)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
all(rownames(z_score) == Annotation$Samples)
Annotation$CDK <- z_score$z_score
#Contact inhibition
#Load the genes involved in quiescence
load("contact_upregulated.RData")
load("contact_downregulated.RData")
#Select only genes involved in quiescence to reduce the size of the dataframe:
expr.data <- count.matrix
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated, upregulated))]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
#Now work out the combined z_score for the samples
gene_lists <- list(upregulated, downregulated)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
all(rownames(z_score) == Annotation$Samples)
Annotation$ContactInhibition <- z_score$z_score
#MEK inhibition
load("mek_downregulated.RData")
load("mek_upregulated.RData")
#Select only genes involved in quiescence to reduce the size of the dataframe:
expr.data <- count.matrix
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated, upregulated))]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
#Now work out the combined z_score for the samples
gene_lists <- list(upregulated, downregulated)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
all(rownames(z_score) == Annotation$Samples)
Annotation$MEKInhibition <- z_score$z_score
#serum starvation
load("serum_downregulated.RData")
load("serum_upregulated.RData")
#Select only genes involved in quiescence to reduce the size of the dataframe:
expr.data <- count.matrix
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated, upregulated))]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
#Now work out the combined z_score for the samples
gene_lists <- list(upregulated, downregulated)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
all(rownames(z_score) == Annotation$Samples)
Annotation$SerumStarvation <- z_score$z_score
#spontaneous quiescence
load("spontaneous_upregulated.RData")
load("spontaneous_downregulated.RData")
#Select only genes involved in quiescence to reduce the size of the dataframe:
expr.data <- count.matrix
expr.data <- expr.data[,which(colnames(expr.data) %in% c(downregulated, upregulated))]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
#Now work out the combined z_score for the samples
gene_lists <- list(upregulated, downregulated)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
all(rownames(z_score) == Annotation$Samples)
Annotation$Spontaneous <- z_score$z_score
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Results/")
save(Annotation, file = "GSE135215_IndividualProgramme_QS.RData")




########################################
#Summary quiescence score comparison
########################################

#Load and combine quiescence scores from the 5 datasets



#################################################
###Combine quiescence score estimates from the 3 studies:
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Results/")
load("GSE124854_IndividualProgramme_QS.RData")
combined.data <- Annotation
combined.data$Study <- "GSE124854"
load("GSE135215_IndividualProgramme_QS.RData")
Annotation$Study <- "GSE135215"
combined.data <- rbind(combined.data, Annotation)
load("GSE99116_IndividualProgramme_QS.RData")
Annotation$Study <- "GSE99116"
combined.data <- rbind(combined.data, Annotation)

##Clean up the dataframe:
combined.data$Treatment <- sapply(combined.data$Treatment, function(x)
  ifelse(x %in% 1,"Palbociclib Inhibition","Control"))
library(reshape)
combined.data <- melt(combined.data, id.vars = c("Samples","Study","Treatment"))
colnames(combined.data) <- c("Samples","Study","Treatment","QuiescenceScoreType","QuiescenceScore")


#########################################
#Plot comparisons for all QS
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Figures/")
pdf("BulkRNAseq_PalbociclibTreatmentResponse.pdf",height = 6, width = 7)


combined.data %>%
  ggplot( aes(x=Study, y=QuiescenceScore, fill=Treatment)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(), size = 0.1) +
  theme_classic() +
  xlab("") + stat_compare_means(aes(group = Treatment), label = "p.signif") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_wrap(~QuiescenceScoreType)


dev.off()

