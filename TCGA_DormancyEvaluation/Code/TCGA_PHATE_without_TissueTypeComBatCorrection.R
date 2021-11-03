######################################################################################
#Phate dimentionality reduction without ComBat correction (tissue effects NOT removed)
#####################################################################################

#Load the required R packages:
library(phateR)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

#Load the quiescence biomarker gene lists
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("downregulated_common_ENSG.RData")
load("upregulated_common_ENSG.RData")


#Load the pancancer expression data:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
load("combined_experssion_FPKM.RData")
#This is an example dataframe with 100 entries only


#Perform phate dimentionality reduction only on the genes differentially expressed in all 5 forms of quiescence:
phate_expr <- combined_data[,colnames(combined_data) %in% c(upregulated_common_ENSG, downregulated_common_ENSG)]
phate_expr <- as.matrix(phate_expr)
setwd("~/Documents/GitHub/CancerDormancy/TCGA_DormancyEvaluation/Results")
library(phateR)
set.seed(123456)
data_phate <- phate(phate_expr)
phate_coordinates <- data.frame(data_phate$embedding)
save(phate_coordinates, file = "PHATE_Coordinates_without_ComBat.RData")




####################################
#Plot PHATE Coordinates
####################################

#Load phate coordinates
setwd("~/Documents/GitHub/CancerDormancy/TCGA_DormancyEvaluation/Results")
load("PHATE_Coordinates_without_ComBat.RData")
phate_coordinates$Barcode <- rownames(phate_coordinates)

#Load quiescence scores scaled to take into account tumour purity 
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
load("TCGA_common_QS_purity_scaled.RData")


#Merge the two dataframes
phate_coordinates <- phate_coordinates[phate_coordinates$Barcode %in% rownames(z_score),]
phate_coordinates <- phate_coordinates[order(phate_coordinates$Barcode),]
z_score <- z_score[order(rownames(z_score)),]
all(rownames(z_score) == phate_coordinates$Barcode)
phate_coordinates$"Quiescence Score" <- z_score$z_score
phate_coordinates$CancerType <- z_score$CancerType

#Plot the coordiantes
mypalette  = colorRampPalette(
  c("darkorange3", "lightgoldenrodyellow","darkblue")
)(100)


#Colour by quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/TCGA_DormancyEvaluation/Figures/")
plot <- ggplot(phate_coordinates, aes(x= PHATE2, y=PHATE1, color = `Quiescence Score`)) + geom_point(size = 2, alpha = 0.3) + theme_classic()
pdf("PHATE_NoComBat_coloured_by_QS.pdf", width = 6, height = 6)
plot + scale_colour_gradientn(colours = mypalette, limits=c(-20, 20))
dev.off()


#Plot by cancer type:
setwd("~/Documents/GitHub/CancerDormancy/TCGA_DormancyEvaluation/Figures/")
plot <- ggplot(phate_coordinates, aes(x= PHATE2, y=PHATE1, color = CancerType)) + geom_point(size = 1.25, alpha = 0.3)
pdf("PHATE_noComBat_coloured_by_CT.pdf", width = 6, height = 6)
plot + theme_classic()
dev.off()


