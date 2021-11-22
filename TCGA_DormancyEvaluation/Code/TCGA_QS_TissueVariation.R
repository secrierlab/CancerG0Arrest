#######################################################
###Variation in quiescence scores across tissues (TCGA)
######################################################

#Load the required R packages:
library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)



#Load quiescence scores (common programme)
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
load("TCGA_common_QS_purity_scaled.RData")

#Split ESCA cancers into adenocarcinoma and squamous cell carcinomas
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_ClinicalData/")
ESCA_clinical <- read.csv("TCGA-ESCA_clinical.csv")
z_score$PatientID <- sapply(rownames(z_score), function(x)
  paste(strsplit(x,"-")[[1]][1:3],collapse="-"))
#Merge with the clinical data:
z_score_ESCA <- merge(z_score, ESCA_clinical,
                      by.x = "PatientID",by.y = "submitter_id")
table(z_score_ESCA$primary_diagnosis)
SC <- z_score_ESCA[z_score_ESCA$primary_diagnosis %in% c("Squamous cell carcinoma, NOS","Basaloid squamous cell carcinoma","Squamous cell carcinoma, keratinizing, NOS"),]
SC <- as.character(SC$PatientID)
AC <- z_score_ESCA[z_score_ESCA$primary_diagnosis %in% c("Adenocarcinoma, NOS","Tubular adenocarcinoma","Mucinous adenocarcinoma"),]
AC <- as.character(AC$PatientID)
z_score$CancerType[z_score$PatientID %in% AC] <- 'ESCA-AC'
z_score$CancerType[z_score$PatientID %in% SC] <- 'ESCA-SC'
z_score <- z_score[(!z_score$CancerType %in% "ESCA"),]
table(z_score$CancerType)


#Get a mean quiescence score estimate for all cancer types
CT <- unique(z_score$CancerType)
mean <- NULL
for (i in CT) {
  
  selected.data <- z_score[z_score$CancerType %in% i,]
  x <- median(selected.data$z_score)
  mean <- c(mean, x)
}
means <- data.frame(CT, mean)
means <- means[order(means$mean),]
z_score$"Quiescence Score" <- z_score$z_score


##############Plot the pan-cancer variation:
Ordered_cancers <- means$CT
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
myCols <- myPalette(32)
z_score$CancerType <- factor(z_score$CancerType,
                             levels = Ordered_cancers)
setwd("~/Documents/GitHub/CancerDormancy/TCGA_DormancyEvaluation/Figures/")


colfunc <- colorRampPalette(c("#1B7837","grey95","#762A83"))
colfunc(32)
myCols <- colfunc(32)

pdf("TCGA_QuiescenceScore_TissueVariation.pdf", width = 5, height = 6)
p <- ggplot(z_score, aes(x=CancerType, y=z_score, fill=CancerType)) +
  geom_boxplot() +
  scale_fill_manual(values = myCols) + theme_classic() +
  ggtitle("") +
  xlab("Cancer Type") + ylab("Quiescence Score")
p + rotate_x_text(45) + labs(y = "Quiescence Score", x = "Cancer Type") + theme(legend.position = "bottom")
dev.off()