############################################################################################
#######Variation in quiescence score correlates with stem cell division rate
############################################################################################

#Load required packages:
library(ggpubr)


#Load stem cell division data
setwd("~/Documents/GitHub/CancerDormancy/Data/StemCellDivisionRates/")
stem_cell_division_rate <- read.csv("Stem_cell_division_rates.csv", header = TRUE) #This data was obtained from data reported by Tomasetti & Vogelstein, Science, 2015 
stem_cell_division_rate <- stem_cell_division_rate[1:11,]

#Load quiescence scores (general programme):
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
load("TCGA_common_QS_purity_scaled.RData")


#Get an average estimate QS estimate for ESCA-SC
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_ClinicalData/")
ESCA_clinical <- read.csv("TCGA-ESCA_clinical.csv")
#Mean unscaled quiescence score:
z_score_ESCA <- z_score[z_score$CancerType %in% "ESCA",]
z_score_ESCA$PatientID <- sapply(rownames(z_score_ESCA), function(x)
  paste(strsplit(x,"-")[[1]][1:3],collapse="-"))
#Merge with the clinical data:
z_score_ESCA <- merge(z_score_ESCA, ESCA_clinical,
                      by.x = "PatientID",by.y = "submitter_id")
table(z_score_ESCA$primary_diagnosis)
z_score_ESCA <- z_score_ESCA[z_score_ESCA$primary_diagnosis %in% c("Squamous cell carcinoma, NOS","Basaloid squamous cell carcinoma","Squamous cell carcinoma, keratinizing, NOS"),]
mean_ESCA <- mean(z_score_ESCA$z_score)


#########Work out means for other cancers:
CT <- as.character(stem_cell_division_rate$Cancer.type)
CT <- CT[-1]
Mean_QS <- mean_ESCA

for (a in CT) {
  
  selected.data <- z_score[z_score$CancerType %in% a, ]
  
  mean <- mean(selected.data$z_score)
  
  Mean_QS <- c(Mean_QS, mean)
  rm(mean)
  
}

#Add this information to the stem cell information dataframe
stem_cell_division_rate$Mean_QS <- Mean_QS

###Log transform the number of stem cell divisions per lifetime:
stem_cell_division_rate$Log10.stem.cell.divisions.per.lifetime <- log10(stem_cell_division_rate$Number.of.divisions.of.each.stem.cell.per.lifetime + 1)
stem_cell_division_rate$Log2.stem.cell.divisions.per.lifetime <- log2(stem_cell_division_rate$Number.of.divisions.of.each.stem.cell.per.lifetime + 1)

#Plot quiescence scores vs log transformed division rates:
setwd("~/Documents/GitHub/CancerDormancy/TCGA_DormancyEvaluation/Figures/")
pdf("AvStemCellDivision_vs_MeanTCGA_QS.pdf", width = 5, height = 5)
p <- ggscatter(stem_cell_division_rate, x = "Log2.stem.cell.divisions.per.lifetime", y = "Mean_QS",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE,label = stem_cell_division_rate$Cancer.type,repel = TRUE,
               conf.int = TRUE, cor.coeff.args = list(method = "pearson", label.x = 7,label.y = 11, label.sep = "\n"))
p + labs(x = "Log2 number of stem cell divisions per lifetime", y = "Mean tissue quiescence score")
dev.off()


