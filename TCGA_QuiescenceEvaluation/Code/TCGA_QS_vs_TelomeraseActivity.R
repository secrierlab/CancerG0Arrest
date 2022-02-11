#######################################################
##Quiescence score association with telomerase activity
#######################################################


#Load required R packages:
library(ggpubr)
library(RColorBrewer)

####################
#Load quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
load("TCGA_common_QS_purity_scaled.RData")


###################
#Load telomerase acitivity scores:
telomerase.score <- read.csv("41467_2020_20474_MOESM7_ESM.csv", header = TRUE)
#data was downloaded from Noureen et al 2021 (Nature Communications)
telomerase.score$SampleID <-  gsub('\\.', '-', telomerase.score$SampleID)
telomerase.score$SampleID <- sapply(telomerase.score$SampleID, function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))

######################
#Merge the two dataframes:
z_score$SampleID <- sapply(rownames(z_score), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
merge.data <- merge(z_score, telomerase.score, 
                    by.x = "SampleID", by.y = "SampleID")


####################
##Plot correlation;
mypalette  = colorRampPalette(
  c("darkorange3", "lightgoldenrodyellow","darkblue")
)(100)

#Plot quiescence scores vs log transformed division rates:
setwd("~/Documents/GitHub/CancerDormancy/TCGA_DormancyEvaluation/Figures/")
pdf("TCGA_QS_vs_telomerase_activity.pdf", width = 5, height = 5)
p <- ggscatter(merge.data, x = "EXTEND.Scores", y = "z_score",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"),
               cor.coef = TRUE, alpha = 0.1, size = 2)
p + labs(x = "EXTEND Score", y = "Quiescence Score") 
dev.off()
