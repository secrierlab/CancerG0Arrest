###########################################################
#####TCGA quiescence  association with CA20 scores
###########################################################

#Load required packages:
library(ggplot2)
library(ggpubr)


#Load quiescence scores
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
load("TCGA_common_QS_purity_scaled.RData")


#Load CA20 scores:
CA20 <- read.table("journal.pcbi.1006832.s018.txt", header = TRUE, sep = "\t")
#CA20 scores were obtained from de Almeida et al 2019, Plos Computational Biology


#Merge the two dataframes:
z_score$PatientID <- sapply(rownames(z_score), function(x)
  paste(strsplit(x,"-")[[1]][1:3],collapse="-"))
CA20$Sample.ID <- gsub('\\.', '-', CA20$Sample.ID)
CA20$PatientID <- sapply(CA20$Sample.ID, function(x)
  paste(strsplit(x,"-")[[1]][1:3],collapse="-"))
z_score <- merge(z_score, CA20, 
                 by.x = "PatientID", by.y = "PatientID")



#Plot the correlation association:
setwd("~/Documents/GitHub/CancerDormancy/CEP89_QuiescenceAssociation/Figures/")
pdf("TCGA_QuiescenceScore_vs_CA20.pdf", width = 7, height = 7)
p <- ggscatter(z_score, x = "CA20", y = "z_score", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE,
               conf.int = TRUE, cor.coeff.args = list(method = "pearson", label.x = 7,label.y = 40, label.sep = "\n"), alpha = 0.2)
p + labs(x = "CA20 Score", y = "Quiescence Score") + geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE)
dev.off()





