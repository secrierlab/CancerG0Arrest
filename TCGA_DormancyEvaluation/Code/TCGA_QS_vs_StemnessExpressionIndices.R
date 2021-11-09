####################################################
##Quiescence score association with stemness markers:
###################################################

#Load required R packages:
library(ggpubr)
library(RColorBrewer)


####################
#Load quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
load("TCGA_common_QS_purity_scaled.RData")



###################
#Load stemness scores:
stemness.rna.score <- read.csv("ExpressionScore.csv", header = TRUE)
#scores can be downloaded from Malta et al 2018  (Cell)
stemness.rna.score$TCGAlong.id <- as.character(stemness.rna.score$TCGAlong.id)
stemness.rna.score$SampleID <- sapply(stemness.rna.score$TCGAlong.id, function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
z_score$SampleID <- sapply(rownames(z_score), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
merged.data <- merge(z_score, stemness.rna.score,
                     by.x = "SampleID", by.y = "SampleID")




####Plot correlation:
setwd("~/Documents/GitHub/CancerDormancy/TCGA_DormancyEvaluation/Figures/")
pdf("TCGA_QS_vs_StemnessExpressionIndex.pdf", width = 5, height = 5)
p <- ggscatter(merged.data, x = "mRNAsi", y = "z_score",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
               cor.coef = TRUE, alpha = 0.1, size = 2)
p + labs(x = "Expression Stemness Index", y = "Quiescence Score")
dev.off()






