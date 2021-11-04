########################################
#Obtaining dormancy biomarker gene lists 
########################################

#Code to obtain a list of genes associated with spontaneous, contact inhibition, serum starvation, mek inhibition and cdk inhibition quiescence
#As well as a list of genes associated with all 5 forms of quiescence 
#Based on RNA-seq data from GSE122927

#Load required packages:
library(biomaRt)


########################################################################################################################
##Step 1: Save genes that are diffierentially expressed in each form of quiescence (with an adjusted p value of < 0.05)
#########################################################################################################################
#Load the differential expression analysis results and save genes differentially expressed in each quiescence form
DE_data <- read.csv("journal.pbio.3000178.s008.csv", header = TRUE)
#Data can be downloaded from the Min et al 2019,PLOS BIOL

#Select genes that are differentially expressed in all 5 forms of quiescence with an adjusted p_value of 0.05
#Spontaneous quiescence
spontaneous_quiescence <- DE_data[which(DE_data$spontaneous_padj < 0.05),]
spontaneous_quiescence_genes <- as.vector(spontaneous_quiescence$Ensembl.Gene.ID)
#Contact inhibition quiescence
contact_quiescence <- DE_data[which(DE_data$contact_padj < 0.05),]
contact_quiescence_genes <- as.vector(contact_quiescence$Ensembl.Gene.ID)
#Cdk inhibition quiescence
cdk_quiescence <- DE_data[which(DE_data$cdki_padj < 0.05),]
cdk_quiescence_genes <- as.vector(cdk_quiescence$Ensembl.Gene.ID)
#Mek inhibition quiescence
mek_quiescence <- DE_data[which(DE_data$meki_padj < 0.05),]
mek_quiescence_genes <- as.vector(mek_quiescence$Ensembl.Gene.ID)
#Serum starvation
serum_quiescence <- DE_data[which(DE_data$serum_padj < 0.05),]
serum_quiescence_genes <- as.vector(serum_quiescence$Ensembl.Gene.ID)
#Obtaining a list of genes differentially expressed in all 5 forms of quiescence 
intersect_all <- function(a,b,c,d,e){
  Reduce(intersect, list(a,b,c,d,e))
}
common_genes <- intersect_all(cdk_quiescence_genes, contact_quiescence_genes, mek_quiescence_genes, serum_quiescence_genes, spontaneous_quiescence_genes)





##########################################################################################################################################
##Step2: For the common gene list ONLY, make sure that all of the genes show unidirectional logfold changes across all 5 quiescence forms:
##########################################################################################################################################
all_genes <- c(cdk_quiescence_genes, common_genes, contact_quiescence_genes, mek_quiescence_genes, serum_quiescence_genes, spontaneous_quiescence_genes)
all_genes <- unique(all_genes)
test <- DE_data[which(DE_data$Ensembl.Gene.ID %in% all_genes),]
test <- test[,c(1,6,9,12,15,18)]
rownames(test) <- test$Ensembl.Gene.ID
test$Ensembl.Gene.ID <- NULL
test$spontaneous_log2FoldChange <- sign(test$spontaneous_log2FoldChange)
test$serum_log2FoldChange <- sign(test$serum_log2FoldChange)
test$meki_log2FoldChange <- sign(test$meki_log2FoldChange)
test$cdki_log2FoldChange <- sign(test$cdki_log2FoldChange)
test$contact_log2FoldChange <- sign(test$contact_log2FoldChange)
test$direction <- rowMeans(test)
test <- test[test$direction %in% c(1.0, -1.0),]
remaining_genes <- as.vector(rownames(test))
common_genes <- common_genes[common_genes %in% remaining_genes]




############################################################################################################
##Step 3:Make sure that genes for the other quiescence programs are unique to only this quiescence program
###########################################################################################################
###Spontaneous quiescence:(342 genes left)
spontaneous_quiescence_genes <- spontaneous_quiescence_genes[!(spontaneous_quiescence_genes %in% cdk_quiescence_genes)]
spontaneous_quiescence_genes <- spontaneous_quiescence_genes[!(spontaneous_quiescence_genes %in% serum_quiescence_genes)]
spontaneous_quiescence_genes <- spontaneous_quiescence_genes[!(spontaneous_quiescence_genes %in% mek_quiescence_genes)]
spontaneous_quiescence_genes <- spontaneous_quiescence_genes[!(spontaneous_quiescence_genes %in% contact_quiescence_genes)]
#Contact inhibition quiescence: (359 genes left)
contact_quiescence_genes <- contact_quiescence_genes[!(contact_quiescence_genes %in% cdk_quiescence_genes)]
contact_quiescence_genes <- contact_quiescence_genes[!(contact_quiescence_genes %in% serum_quiescence_genes)]
contact_quiescence_genes <- contact_quiescence_genes[!(contact_quiescence_genes %in% mek_quiescence_genes)]
contact_quiescence_genes <- contact_quiescence_genes[!(contact_quiescence_genes %in% spontaneous_quiescence_genes)]
#Mek inhibition: (1117 genes left)
mek_quiescence_genes <- mek_quiescence_genes[!(mek_quiescence_genes %in% cdk_quiescence_genes)]
mek_quiescence_genes <- mek_quiescence_genes[!(mek_quiescence_genes %in% serum_quiescence_genes)]
mek_quiescence_genes <- mek_quiescence_genes[!(mek_quiescence_genes %in% spontaneous_quiescence_genes)]
mek_quiescence_genes <- mek_quiescence_genes[!(mek_quiescence_genes %in% serum_quiescence_genes)]
#serum starvation (5349 genes lef)t:
serum_quiescence_genes <- serum_quiescence_genes[!(serum_quiescence_genes %in% cdk_quiescence_genes)]
serum_quiescence_genes <- serum_quiescence_genes[!(serum_quiescence_genes %in% spontaneous_quiescence_genes)]
serum_quiescence_genes <- serum_quiescence_genes[!(serum_quiescence_genes %in% mek_quiescence_genes)]
serum_quiescence_genes <- serum_quiescence_genes[!(serum_quiescence_genes %in% contact_quiescence_genes)]
#Contact inhibition (128 genes left)
cdk_quiescence_genes <- cdk_quiescence_genes[!(cdk_quiescence_genes %in% contact_quiescence_genes)]
cdk_quiescence_genes <- cdk_quiescence_genes[!(cdk_quiescence_genes %in% mek_quiescence_genes)]
cdk_quiescence_genes <- cdk_quiescence_genes[!(cdk_quiescence_genes %in% serum_quiescence_genes)]
cdk_quiescence_genes <- cdk_quiescence_genes[!(cdk_quiescence_genes %in% spontaneous_quiescence_genes)]



##################################################################################
##Step 4:Select genes that have mean pancancer expression > 0.5 (log2 FPKM data)
##################################################################################

#Load the combined pancancer expression data
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
load("combined_experssion_FPKM.RData")
#This is an example expression dataframe with 100 entries only
#Remove columns denoting cancer types
combined_data$cancer_type <- NULL
#Calculate the mean expression of each gene:
means <- colMeans(combined_data)
#Look at the distribution of the mean expression values
hist(means, breaks = 100)
#Select genes with mean logfold expression values greater than 0.5
genes_with_high_expression <- means[which(means > 0.5)]
genes_with_high_expression <- data.frame(genes_with_high_expression)
genes_high_expression <- as.vector(rownames(genes_with_high_expression))
#Remove genes from each gene list that show low pancancer expression
common_genes <- common_genes[common_genes %in% genes_high_expression] #(191 genes remain)
cdk_quiescence_genes <- cdk_quiescence_genes[cdk_quiescence_genes %in% genes_high_expression] #(1198 genes remain)
contact_quiescence_genes <- contact_quiescence_genes[contact_quiescence_genes %in% genes_high_expression] #(320 genes remain)
mek_quiescence_genes <- mek_quiescence_genes[mek_quiescence_genes %in% genes_high_expression] #(1021 genes remain)
serum_quiescence_genes <- serum_quiescence_genes[serum_quiescence_genes %in% genes_high_expression] #(4863 genes remain)
spontaneous_quiescence_genes <- spontaneous_quiescence_genes[spontaneous_quiescence_genes %in% genes_high_expression] #(291 genes remain)




###################################################################
##Step 5: Remove genes involved in other stages of the cell cycle:
##################################################################

#Firstly load genes involved in different stages of the cell cycle
setwd("~/Documents/GitHub/CancerDormancy/Data/OtherGeneLists/")
cell_cycle_gene_list <- read.table(file = "CellCycle_genelist.txt", sep = "\t", header = TRUE)
cell_cycle_gene_list <- as.vector(cell_cycle_gene_list$KEGG_CELL_CYCLE)
#Convert this gene list into ENSG numbers
human = useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "useast.ensembl.org")
results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=cell_cycle_gene_list, mart=human, filters = "hgnc_symbol")
cell_cycle_gene_list <- results$ensembl_gene_id
#Now check if any of the cell cycle genes are present in the list of genes used to assess quiescence 
common_genes <- common_genes[!(common_genes %in% cell_cycle_gene_list)] #180 genes remain
cdk_quiescence_genes<- cdk_quiescence_genes[!(cdk_quiescence_genes %in% cell_cycle_gene_list)] #1161 genes remain
contact_quiescence_genes<- contact_quiescence_genes[!(contact_quiescence_genes %in% cell_cycle_gene_list)] #320 genes remain
mek_quiescence_genes<- mek_quiescence_genes[!(mek_quiescence_genes %in% cell_cycle_gene_list)] #1016 genes remain
serum_quiescence_genes<- serum_quiescence_genes[!(serum_quiescence_genes %in% cell_cycle_gene_list)] #4811 genes remain
spontaneous_quiescence_genes<- spontaneous_quiescence_genes[!(spontaneous_quiescence_genes %in% cell_cycle_gene_list)] #289 genes remain


#################################################################
##Step 6:Remove genes with sd of pancancer expression less than 1
#################################################################
standard_deviations <- apply(combined_data, 2, sd)
hist(standard_deviations, breaks = 100)
standard_deviations_1.0 <- standard_deviations[which(means > 1)]
standard_deviations_1.0 <- data.frame(standard_deviations_1.0)
standard_deviations_1.0 <- as.vector(rownames(standard_deviations_1.0))
common_genes <- common_genes[common_genes %in% standard_deviations_1.0] #175 genes remain
cdk_quiescence_genes <- cdk_quiescence_genes[cdk_quiescence_genes %in% standard_deviations_1.0] #1088 genes remain
contact_quiescence_genes <- contact_quiescence_genes[contact_quiescence_genes %in% standard_deviations_1.0] #278 genes remain
mek_quiescence_genes <- mek_quiescence_genes[mek_quiescence_genes %in% standard_deviations_1.0] #920 genes remain
serum_quiescence_genes <- serum_quiescence_genes[serum_quiescence_genes %in% standard_deviations_1.0] #4428 genes remain
spontaneous_quiescence_genes <- spontaneous_quiescence_genes[spontaneous_quiescence_genes %in% standard_deviations_1.0] #273 genes remain




###################################################################################
##Step 6: Split each gene list into downregulated and up-regulated genes:
####################################################################################

#Common gene list:
common_logfoldchanges <- DE_data[,c(1,6,9,12,15,18)]
common_logfoldchanges$average <- rowMeans(common_logfoldchanges[,2:6])
logfoldchanges_positive <- common_logfoldchanges[which(common_logfoldchanges$average > 0),]
logfoldchanges_negative <- common_logfoldchanges[which(common_logfoldchanges$average < 0),]
positive <- as.vector(logfoldchanges_positive$Ensembl.Gene.ID)
negative <- as.vector(logfoldchanges_negative$Ensembl.Gene.ID)
upregulated_common <- unique(positive[positive %in% common_genes])
downregulated_common <- unique(negative[negative %in% common_genes])
#CDK inhibition quiescence
logfoldchanges_positive <- DE_data[which(DE_data$cdki_log2FoldChange > 0),]
logfoldchanges_negative <- DE_data[which(DE_data$cdki_log2FoldChange < 0),]
positive <- as.vector(logfoldchanges_positive$Ensembl.Gene.ID)
negative <- as.vector(logfoldchanges_negative$Ensembl.Gene.ID)
upregulated_cdk <- unique(positive[positive %in% cdk_quiescence_genes])
downregulated_cdk <- unique(negative[negative %in% cdk_quiescence_genes])
#Contact inhibition
logfoldchanges_positive <- DE_data[which(DE_data$contact_log2FoldChange > 0),]
logfoldchanges_negative <- DE_data[which(DE_data$contact_log2FoldChange < 0),]
positive <- as.vector(logfoldchanges_positive$Ensembl.Gene.ID)
negative <- as.vector(logfoldchanges_negative$Ensembl.Gene.ID)
upregulated_contact <- unique(positive[positive %in% contact_quiescence_genes])
downregulated_contact <- unique(negative[negative %in% contact_quiescence_genes])
#Mek inhibition
logfoldchanges_positive <- DE_data[which(DE_data$meki_log2FoldChange > 0),]
logfoldchanges_negative <- DE_data[which(DE_data$meki_log2FoldChange < 0),]
positive <- as.vector(logfoldchanges_positive$Ensembl.Gene.ID)
negative <- as.vector(logfoldchanges_negative$Ensembl.Gene.ID)
upregulated_mek <- unique(positive[positive %in% mek_quiescence_genes])
downregulated_mek <- unique(negative[negative %in% mek_quiescence_genes])
#Serum starvation:
logfoldchanges_positive <- DE_data[which(DE_data$serum_log2FoldChange > 0),]
logfoldchanges_negative <- DE_data[which(DE_data$serum_log2FoldChange < 0),]
positive <- as.vector(logfoldchanges_positive$Ensembl.Gene.ID)
negative <- as.vector(logfoldchanges_negative$Ensembl.Gene.ID)
upregulated_serum <- unique(positive[positive %in% serum_quiescence_genes])
downregulated_serum <- unique(negative[negative %in% serum_quiescence_genes])
#Spontaneous quiescence
logfoldchanges_positive <- DE_data[which(DE_data$spontaneous_log2FoldChange > 0),]
logfoldchanges_negative <- DE_data[which(DE_data$spontaneous_log2FoldChange < 0),]
positive <- as.vector(logfoldchanges_positive$Ensembl.Gene.ID)
negative <- as.vector(logfoldchanges_negative$Ensembl.Gene.ID)
upregulated_spontaneous <- unique(positive[positive %in% spontaneous_quiescence_genes])
downregulated_spontaneous <- unique(negative[negative %in% spontaneous_quiescence_genes])



###########################################################################################################
##Step 7: Remove genes that are not correlated with DREAM complex target expression (common gene list only)
##########################################################################################################
#Work out the mean expression of DREAM targets
#Load DREAM complex targets:
setwd("~/Documents/GitHub/CancerDormancy/Data/OtherGeneLists/")
DREAM_target <- read.table(file = "DREAM_target_genelist.txt", header = FALSE, sep = "\t")
HUGO_genes <- as.character(DREAM_target$V1)
#Convert this gene list into ENSG numbers:
human = useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "useast.ensembl.org")
results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=HUGO_genes, mart=human, filters = "hgnc_symbol")
DREAM_target <- results$ensembl_gene_id
#Work out the mean expression of the DREAM targets for each pancancer sample:
all_genes <- colnames(combined_data)
DREAM_target <- DREAM_target[DREAM_target %in% all_genes]
DREAM_target_expression <- rowMeans(subset(combined_data, select = DREAM_target, na.rm = TRUE))
combined_data$DREAM_target_expression <- DREAM_target_expression
#Calculate the correlation of all genes upregulated in quiescence with DREAM target expression:
upregulated_genes <- c(upregulated_cdk, upregulated_common, upregulated_contact, upregulated_mek, upregulated_serum, upregulated_spontaneous)
#Initialise the variables:
correlations <- NULL
for (i in upregulated_genes) {
  print(i)
  test <- cor.test(combined_data$DREAM_target_expression, combined_data[[(paste(i))]])
  correlation <- test$estimate
  correlations <- c(correlations, correlation)
}
correlations <- data.frame(correlations)
correlations$ENSG_genes <- upregulated_genes
#Look at the distribution of the correlation estimates
par(mar = c(4, 4, 4, 4)) # Set the margin on all sides to 2
hist(correlations$correlations, breaks = 20)
#Save the list of genes that are upregulated in quiescence and are negatively correlated with DREAM target expression 
test <- correlations[correlations$correlations <= 0, ]
upregulated_genes <- unique(as.character(test$ENSG_genes))
upregulated_common <- upregulated_common[upregulated_common %in% upregulated_genes] 
rm(test)
rm(correlations)

#Now repeat this for the downregulated gene list:
downregulated_genes <- c(downregulated_cdk, downregulated_common, downregulated_contact, downregulated_mek, downregulated_serum, downregulated_spontaneous)
#Initialise the variables:
correlations <- NULL
for (i in downregulated_genes) {
  print(i)
  test <- cor.test(combined_data$DREAM_target_expression, combined_data[[(paste(i))]])
  correlation <- test$estimate
  correlations <- c(correlations, correlation)
}
correlations <- data.frame(correlations)
correlations$ENSG_genes <- downregulated_genes
#Look at the distribution of the correlation estimates
par(mar = c(4, 4, 4, 4)) # Set the margin on all sides to 2
hist(correlations$correlations, breaks = 20)
#Save the list of genes that are downregulated in quiescence and are positively correlated with DREAM target expression 
test <- correlations[correlations$correlations >= 0, ]
downregulated_genes <- unique(as.character(test$ENSG_genes))
#Update the gene lists:
downregulated_common <- downregulated_common[downregulated_common %in% downregulated_genes] 
rm(test)
rm(correlations)





###########################################################
#Step8: ensure that all gene lists are of comparable size
###########################################################
#There are 112 downregulated and 27 upregulated genes in all 5 forms of quiescence
#Other forms of quiescence have significantly larger gene lists
#Norrow these down to 10 upregulated and 10 downregulated genes
#Report both HGNC and ENSG gene annotation forms


#Common gene list:
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList")
upregulated_common_ENSG <- upregulated_common
save(upregulated_common_ENSG, file = "upregulated_common_ENSG.RData")
downregulated_common_ENSG <- downregulated_common
save(downregulated_common_ENSG, file = "downregulated_common_ENSG.RData")
human = useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "useast.ensembl.org")
results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=upregulated_common, mart=human, filters = "ensembl_gene_id")
upregulated_common <- results$hgnc_symbol
save(upregulated_common, file = "upregulated_common.RData")
results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=downregulated_common, mart=human, filters = "ensembl_gene_id")
downregulated_common <- results$hgnc_symbol
save(downregulated_common, file = "downregulated_common.RData")


#Spontaneous quiescence
DE_data <- read.csv("journal.pbio.3000178.s008.csv", header = TRUE)
spontaneous_quiescence <- DE_data[DE_data$Ensembl.Gene.ID %in% c(upregulated_spontaneous, downregulated_spontaneous),]
spontaneous_quiescence <- spontaneous_quiescence[order(spontaneous_quiescence$spontaneous_log2FoldChange),]
spontaneous_quiescence <- na.omit(spontaneous_quiescence)
downregulated <- spontaneous_quiescence[1:10,]
upregulated <- spontaneous_quiescence[250:259,]
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList")
downregulated_ENSG <- unique(as.character(downregulated$Ensembl.Gene.ID))
upregulated_ENSG <- unique(as.character(upregulated$Ensembl.Gene.ID))
save(upregulated_ENSG, file = "spontaneous_upregulated_ENSG.RData")
save(downregulated_ENSG, file = "spontaneous_downregulated_ENSG.RData")
downregulated <- unique(as.character(downregulated$Gene.Symbol))
upregulated <- unique(as.character(upregulated$Gene.Symbol))
save(upregulated, file = "spontaneous_upregulated.RData")
save(downregulated, file = "spontaneous_downregulated.RData")



#CDK
setwd("~/Documents/Dormancy_PhD_project_data/RNA_seq_data/GSE122927")
DE_data <- read.csv("journal.pbio.3000178.s008.csv", header = TRUE)
cdk_quiescence <- DE_data[DE_data$Ensembl.Gene.ID %in% c(upregulated_cdk, downregulated_cdk),]
cdk_quiescence <- cdk_quiescence[order(cdk_quiescence$cdki_log2FoldChange),]
cdk_quiescence <- na.omit(cdk_quiescence)
downregulated <- cdk_quiescence[1:10,]
upregulated <- cdk_quiescence[1016:1025,]
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList")
downregulated_ENSG <- unique(as.character(downregulated$Ensembl.Gene.ID))
upregulated_ENSG <- unique(as.character(upregulated$Ensembl.Gene.ID))
save(upregulated_ENSG, file = "cdk_upregulated_ENSG.RData")
save(downregulated_ENSG, file = "cdk_downregulated_ENSG.RData")
downregulated <- unique(as.character(downregulated$Gene.Symbol))
upregulated <- unique(as.character(upregulated$Gene.Symbol))
save(upregulated, file = "cdk_upregulated.RData")
save(downregulated, file = "cdk_downregulated.RData")



#serum
setwd("~/Documents/Dormancy_PhD_project_data/RNA_seq_data/GSE122927")
DE_data <- read.csv("journal.pbio.3000178.s008.csv", header = TRUE)
serum_quiescence <- DE_data[DE_data$Ensembl.Gene.ID %in% c(upregulated_serum, downregulated_serum),]
serum_quiescence <- serum_quiescence[order(serum_quiescence$serum_log2FoldChange),]
serum_quiescence <- na.omit(serum_quiescence)
downregulated <- serum_quiescence[1:10,]
upregulated <- serum_quiescence[4135:4144,]
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList")
downregulated_ENSG <- unique(as.character(downregulated$Ensembl.Gene.ID))
upregulated_ENSG <- unique(as.character(upregulated$Ensembl.Gene.ID))
save(upregulated_ENSG, file = "serum_upregulated_ENSG.RData")
save(downregulated_ENSG, file = "serum_downregulated_ENSG.RData")
downregulated <- unique(as.character(downregulated$Gene.Symbol))
upregulated <- unique(as.character(upregulated$Gene.Symbol))
save(upregulated, file = "serum_upregulated.RData")
save(downregulated, file = "serum_downregulated.RData")


#mek
setwd("~/Documents/Dormancy_PhD_project_data/RNA_seq_data/GSE122927")
DE_data <- read.csv("journal.pbio.3000178.s008.csv", header = TRUE)
mek_quiescence <- DE_data[DE_data$Ensembl.Gene.ID %in% c(upregulated_mek, downregulated_mek),]
mek_quiescence <- mek_quiescence[order(mek_quiescence$meki_log2FoldChange),]
mek_quiescence <- na.omit(mek_quiescence)
downregulated <- mek_quiescence[1:10,]
upregulated <- mek_quiescence[815:824,]
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList")
downregulated_ENSG <- unique(as.character(downregulated$Ensembl.Gene.ID))
upregulated_ENSG <- unique(as.character(upregulated$Ensembl.Gene.ID))
save(upregulated_ENSG, file = "mek_upregulated_ENSG.RData")
save(downregulated_ENSG, file = "mek_downregulated_ENSG.RData")
downregulated <- unique(as.character(downregulated$Gene.Symbol))
upregulated <- unique(as.character(upregulated$Gene.Symbol))
save(upregulated, file = "mek_upregulated.RData")
save(downregulated, file = "mek_downregulated.RData")


#contact
setwd("~/Documents/Dormancy_PhD_project_data/RNA_seq_data/GSE122927")
DE_data <- read.csv("journal.pbio.3000178.s008.csv", header = TRUE)
contact_quiescence <- DE_data[DE_data$Ensembl.Gene.ID %in% c(upregulated_contact, downregulated_contact),]
contact_quiescence <- contact_quiescence[order(contact_quiescence$contact_log2FoldChange),]
contact_quiescence <- na.omit(contact_quiescence)
downregulated <- contact_quiescence[1:10,]
upregulated <- contact_quiescence[229:238,]
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList")
downregulated_ENSG <- unique(as.character(downregulated$Ensembl.Gene.ID))
upregulated_ENSG <- unique(as.character(upregulated$Ensembl.Gene.ID))
save(upregulated_ENSG, file = "contact_upregulated_ENSG.RData")
save(downregulated_ENSG, file = "contact_downregulated_ENSG.RData")
downregulated <- unique(as.character(downregulated$Gene.Symbol))
upregulated <- unique(as.character(upregulated$Gene.Symbol))
save(upregulated, file = "contact_upregulated.RData")
save(downregulated, file = "contact_downregulated.RData")



