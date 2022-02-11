######################################################################################
#######Correlation between quiescence scores and other proliferation markers:
######################################################################################

#Load required libraries:
library(biomaRt)
library(ggcorrplot)


##########################################
#Select markers of interest:
#Replication complex, MKI67, DREAM and CDK2
MKI67 <- "ENSG00000148773"
CDK2 <- "ENSG00000123374"
CDK4 <- "ENSG00000135446"
CCNE1 <- "ENSG00000105173"
CDKN1A <- "ENSG00000124762"
CDK6 <- "ENSG00000105810"
RBL2 <- "ENSG00000103479"
CCNB1 <-"ENSG00000134057" 

#DREAM target genes:
setwd("~/Documents/GitHub/CancerDormancy/Data/OtherGeneLists/")
DREAM_targets <- read.table(file = "DREAM_target_genelist.txt", header = FALSE, sep = "\t")
#Downloaded from the Broad insitite (FISCHER DREAM TARGETS)
DREAM_targets <- as.character(DREAM_targets$V1)
human = useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "useast.ensembl.org")
results_biomart <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=DREAM_targets, mart=human, filters = "hgnc_symbol")
DREAM_targets <- results_biomart$ensembl_gene_id
#Replication complex genes:
replication_complex_genes <- c("ENSG00000073111", "ENSG00000112118","ENSG00000104738","ENSG00000100297","ENSG00000076003","ENSG00000166508")

##########################################
#Load TCGA expression data
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
load("combined_scaled_expr_data.RData")
#This is a toy example dataset with 100 entries only

##########################################
#Load TCGA quiescence scores
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
load("TCGA_common_QS_purity_scaled.RData")
z_score$Barcode <- rownames(z_score)

##########################
#combine the two dataframes
combined.scaled <- combined.scaled[rownames(combined.scaled) %in% z_score$Barcode,]
combined.scaled <- combined.scaled[order(rownames(combined.scaled)),]
z_score <- z_score[order(z_score$Barcode),]
all(rownames(combined.scaled) == z_score$Barcode)
combined.scaled$quiescence_score <- z_score$z_score


#rename genes
combined.scaled$MKI67 <- combined.scaled$ENSG00000148773
combined.scaled$CDK2 <- combined.scaled$ENSG00000123374
combined.scaled$CDK4 <- combined.scaled$ENSG00000135446
combined.scaled$CCNE1 <- combined.scaled$ENSG00000105173
combined.scaled$CDKN1A <- combined.scaled$ENSG00000124762
combined.scaled$CDK6 <- combined.scaled$ENSG00000105810
combined.scaled$CCNB1 <- combined.scaled$ENSG00000134057
combined.scaled$RBL2 <- combined.scaled$ENSG00000103479


#Make sure to remove genes which were used in calculating quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("upregulated_common_ENSG.RData")
load("downregulated_common_ENSG.RData")
quiescence_genes <- c(upregulated_common_ENSG, downregulated_common_ENSG)
#Calculate the mean expression of DREAM complex target genes, as well as RC genes
all_genes <- colnames(combined.scaled)
DREAM_targets <- DREAM_targets[DREAM_targets %in% all_genes]
DREAM_targets <- DREAM_targets[!(DREAM_targets %in% quiescence_genes)]
combined.scaled$DREAM <- rowMeans(subset(combined.scaled, select = DREAM_targets), na.rm = TRUE)
replication_complex_genes <- replication_complex_genes[!(replication_complex_genes %in% quiescence_genes)]
combined.scaled$RC <- rowMeans(subset(combined.scaled, select = replication_complex_genes), na.rm = TRUE)
combined.scaled <- combined.scaled[colnames(combined.scaled) %in% c("MKI67","CDK2","CDK4","CDK6","CCNE1","CCNB1","DREAM","RBL2","RC","CDKN1A","quiescence_score")]

#Work out correlations pancancer
genes <- c("MKI67","CDK2","CDK4","CDK6","CCNE1","CCNB1","DREAM","RC","CDKN1A")
Pval <- NULL
Corr <- NULL
for (i in genes) {
  print(i)
  test <- cor.test(combined.scaled[[i]],combined.scaled$quiescence_score)
  pval <- test$p.value
  corr <- test$estimate
  Pval <- c(Pval,pval)
  Corr <- c(Corr,corr)
  
}



#######################################
#Add EXTEND scores and stemness scores:
telomerase.score <- read.csv("41467_2020_20474_MOESM7_ESM.csv", header = TRUE)
#The scores were obtained from Noureen et al 2021, Nature Communications
telomerase.score$SampleID <-  gsub('\\.', '-', telomerase.score$SampleID)
telomerase.score$SampleID <- sapply(telomerase.score$SampleID, function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
combined.scaled$SampleID <- sapply(rownames(z_score), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
merge.data <- merge(combined.scaled, telomerase.score, 
                    by.x = "SampleID", by.y = "SampleID")
test <- cor.test(merge.data$EXTEND.Scores,merge.data$quiescence_score)
pval <- test$p.value
corr <- test$estimate
Pval <- c(Pval,pval)
Corr <- c(Corr,corr)

###########################################
###Add stemness indices:
stemness.rna.score <- read.csv("ExpressionScore.csv", header = TRUE)
#Stemness indicies were obtained from Mlata et al 2018, Cell
stemness.rna.score$TCGAlong.id <- as.character(stemness.rna.score$TCGAlong.id)
stemness.rna.score$SampleID <- sapply(stemness.rna.score$TCGAlong.id, function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
merged.data <- merge(merge.data, stemness.rna.score,
                     by.x = "SampleID", by.y = "SampleID")
test <- cor.test(merged.data$mRNAsi,merged.data$quiescence_score)
pval <- test$p.value
corr <- test$estimate
Pval <- c(Pval,pval)
Corr <- c(Corr,corr)



#Results dataframe:
genes <- c("MKI67","CDK2","CDK4","CDK6","CCNE1","CCNB1","DREAM","RC","CDKN1A","EXTEND","Stemness")
Results_correlation <- data.frame(genes)
Results_pvalue <- data.frame(genes)
Results_correlation$Pancancer <- Corr
Results_pvalue$Pancancer <- Pval
#Clean the dataframes:
rownames(Results_correlation) <- Results_correlation$genes
rownames(Results_pvalue) <- Results_pvalue$genes
Results_correlation$genes <- NULL
Results_pvalue$genes <- NULL
Results_correlation <- as.matrix(t(Results_correlation))
Results_pvalue <- as.matrix(t(Results_pvalue))


#############Plot the results:
setwd("~/Documents/GitHub/CancerDormancy/TCGA_DormancyEvaluation/Figures/")
mypalette  = colorRampPalette(
  c("darkblue", "lightgoldenrodyellow","darkorange3")
)(100)
pdf("TCGA_QS_correlation_with_ProliferationMarkers.pdf", width = 4, height = 5)
p <- ggcorrplot(Results_correlation, p.mat = Results_pvalue, insig = "blank", sig.level = 0.05, lab = TRUE, lab_size = 2.5, method = "circle",show.legend = TRUE,tl.cex = 10,legend.title = "Correlation",ggtheme = ggplot2::theme_classic(),colors = c("darkorange3", "lightgoldenrodyellow", "darkblue")) 
p 
dev.off()
