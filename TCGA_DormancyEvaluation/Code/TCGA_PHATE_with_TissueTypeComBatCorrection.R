##############################################################################
#Phate dimentionality reduction on ComBat scaled data (tissue effects removed)
##############################################################################

#Load the required R packages:
library(phateR)
library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

#Load the quiescence biomarker gene lists
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("downregulated_common.RData")
load("upregulated_common.RData")


#Load the pancancer expression (ComBat scaled):
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
expr.data <- fread("TCGA_combat_tumor_type_correction.txt", sep = "\t")
#This is an example dataframe with 100 entries only

#Select only genes involved in quiescence to reduce the size of the dataframe:
expr.data <- expr.data[which(expr.data$genes %in% c(downregulated_common, upregulated_common)),]
#Set the genes as rownames
all_genes <- expr.data$genes
expr.data$genes <- NULL
expr.data <- as.data.frame(expr.data)
rownames(expr.data) <- all_genes

#Transpose the dataframe
expr.data <- data.frame(t(expr.data))

#Make sure that you select primary tumour samples:
rownames(expr.data) <- gsub('\\.', '-', rownames(expr.data))
expr.data$SampleType <- sapply(rownames(expr.data), function(x)
  strsplit(x,"-")[[1]][5])
table(expr.data$SampleType)
expr.data <- expr.data[which(expr.data$SampleType %in% c("01A","01B","01C")),]
expr.data$SampleType <- NULL

#Make sure that you only select solid tumour samples (not blood cancers)
expr.data$StudyType <- sapply(rownames(expr.data), function(x)
  strsplit(x,"-")[[1]][1])
table(expr.data$StudyType)
CT <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
expr.data <- expr.data[which(expr.data$StudyType %in% CT),]
expr.data$StudyType <- NULL

#The data is already log transformed and so is ready for phate
expr.data <- as.matrix(expr.data)

#Now the data is ready for phate
setwd("~/Documents/GitHub/CancerDormancy/TCGA_DormancyEvaluation/Results")
set.seed(123456)
data_phate <- phate(expr.data)
phate_coordinates <- data.frame(data_phate$embedding)
save(phate_coordinates, file = "PHATE_Coordinates_with_ComBat.RData")






####################################
#Plot PHATE Coordinates
####################################

#Load phate coordinates
setwd("~/Documents/GitHub/CancerDormancy/TCGA_DormancyEvaluation/Results/")
load("PHATE_Coordinates_with_ComBat.RData")
phate_coordinates$CancerType <- sapply(rownames(phate_coordinates), function(x)
  strsplit(x,"-")[[1]][1])
phate_coordinates$Barcode <- sapply(rownames(phate_coordinates), function(x)
  paste(strsplit(x,"-")[[1]][2:8],collapse="-"))

#Load quiescence scores scaled to take into account tumour purity 
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
load("TCGA_common_QS_purity_scaled.RData")


#Merge the two dataframes
phate_coordinates <- phate_coordinates[phate_coordinates$Barcode %in% rownames(z_score),]
phate_coordinates <- phate_coordinates[order(phate_coordinates$Barcode),]
z_score <- z_score[order(rownames(z_score)),]
all(rownames(z_score) == phate_coordinates$Barcode)
phate_coordinates$"Quiescence Score" <- z_score$z_score


#Plot the coordiantes

#Plot by cancer type:
setwd("~/Documents/GitHub/CancerDormancy/TCGA_DormancyEvaluation/Figures/")
plot <- ggplot(phate_coordinates, aes(x= PHATE2, y=PHATE1, color = CancerType)) + geom_point(size = 1.25, alpha = 0.3)
pdf("PHATE_ComBat_coloured_by_CT.pdf", width = 6, height = 6)
plot + theme_classic()
dev.off()
#colour by QS
colfunc <- colorRampPalette(c("#1B7837","grey95","#762A83"))
colfunc(100)
myCols <- colfunc(100)
#Colour by quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/TCGA_DormancyEvaluation/Figures/")
plot <- ggplot(phate_coordinates, aes(x= PHATE2, y=PHATE1, color = `Quiescence Score`)) + geom_point(size = 2, alpha = 0.3) + theme_classic()
pdf("PHATE_ComBat_coloured_by_QS.pdf", width = 6, height = 6)
plot + scale_colour_gradientn(colours = myCols, limits=c(-20, 20))
dev.off()

setwd("~/Documents/GitHub/CancerDormancy/TCGA_DormancyEvaluation/Figures/")
plot <- ggplot(phate_coordinates, aes(x= PHATE2, y=PHATE1, color = `Quiescence Score`)) + geom_point(size = 3, alpha = 0.3) + theme_classic()
pdf("PHATE_ComBat_coloured_by_QS_v2.pdf", width = 6, height = 6)
plot + scale_colour_gradientn(colours = myCols, limits=c(-20, 20))
dev.off()

