#############################################################
########CDKN1A expression association with cancer quiescence
#############################################################

#Load required packages:
library(ggplot2)

############################
####Load quiescence scores:
############################
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
load("TCGA_common_QS_purity_scaled.RData")
z_score$QuiescenceScore <- z_score$z_score
z_score$SampleID <- sapply(rownames(z_score), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
z_score$Barcode <- rownames(z_score)


######################
###Load RNA seq data:
######################
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
load("combined_scaled_expr_data.RData") #This is an example dataset with 100 entries only 
combined.scaled$Barcode <- rownames(combined.scaled)
combined.scaled <- combined.scaled[,colnames(combined.scaled) %in% c("Barcode","ENSG00000124762")] #CDKN1A = ENSG00000124762

########################
##Combine the datasets:
########################
combined.data <- merge(combined.scaled, z_score,
                       by.x = "Barcode", by.y = "Barcode")
combined.data$CDKN1A <- combined.data$ENSG00000124762


#####################################################
#Split ESCA patients into SC and AC
#####################################################

####Split ESCA SC and AC patients
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_ClinicalData/")
ESCA_clinical <- read.csv("TCGA-ESCA_clinical.csv")
combined.data$PatientID <- sapply(combined.data$Barcode, function(x)
  paste(strsplit(x,"-")[[1]][1:3],collapse="-"))
#Merge with the clinical data:
combined.data_ESCA <- merge(combined.data, ESCA_clinical,
                            by.x = "PatientID",by.y = "submitter_id")
table(combined.data_ESCA$primary_diagnosis)
SC <- combined.data_ESCA[combined.data_ESCA$primary_diagnosis %in% c("Squamous cell carcinoma, NOS","Basaloid squamous cell carcinoma","Squamous cell carcinoma, keratinizing, NOS"),]
SC <- as.character(SC$PatientID)
AC <- combined.data_ESCA[combined.data_ESCA$primary_diagnosis %in% c("Adenocarcinoma, NOS","Tubular adenocarcinoma","Mucinous adenocarcinoma"),]
AC <- as.character(AC$PatientID)
combined.data$CancerType[combined.data$PatientID %in% AC] <- 'ESCA-AC'
combined.data$CancerType[combined.data$PatientID %in% SC] <- 'ESCA-SC'
combined.data <- combined.data[(!combined.data$CancerType %in% "ESCA"),]
table(combined.data$CancerType)

#######Determine UQ and LQ for each cancer type 
CT <- unique(as.character(combined.data$CancerType))
a <- CT[1]
new.df <- combined.data[combined.data$CancerType %in% a,]
LQ <- summary(new.df$CDKN1A)[2]
UQ <- summary(new.df$CDKN1A)[5]
new.df$CDKN1A_expression <- sapply(new.df$CDKN1A, function(x)
  ifelse(x > UQ, "High",
         ifelse( x < LQ,"Low","Other")))
new.df <- new.df[new.df$CDKN1A_expression %in% c("High","Low"),]
merged.df <- new.df
CT <- CT[-1]
for (a in CT) {
  new.df <- combined.data[combined.data$CancerType %in% a,]
  LQ <- summary(new.df$CDKN1A)[2]
  UQ <- summary(new.df$CDKN1A)[5]
  new.df$CDKN1A_expression <- sapply(new.df$CDKN1A, function(x)
    ifelse(x > UQ, "High",
           ifelse( x < LQ,"Low","Other")))
  new.df <- new.df[new.df$CDKN1A_expression %in% c("High","Low"),]
  merged.df <- rbind(merged.df, new.df)
  
}
#Order cancers by quiescence score
CT <- unique(merged.df$CancerType)
Ordered_cancers <- c("TGCT","THYM","UCS","READ","COAD","CESC","BLCA","CHOL","LUSC","STAD","HNSC","BRCA","UCEC","ESCA-SC","SKCM","ESCA-AC","LUAD","MESO","OV","PAAD","PRAD","UVM","SARC","THCA","LIHC","KIRC","PCPG","KICH","LGG","KIRP","GBM","ACC")
Ordered_cancers <- Ordered_cancers[Ordered_cancers %in% CT]
merged.df <- merged.df[merged.df$CancerType %in% Ordered_cancers,]
merged.df$CancerType <- factor(merged.df$CancerType,
                               levels = Ordered_cancers)
setwd("~/Documents/GitHub/CancerDormancy/TCGA_DormancyEvaluation/Figures/")
pdf("TCGA_CDKN1A_expression_vs_QS.pdf",height = 5,width = 15)
my_comparisons <- list( c("High", "Low"))
p <- ggplot(merged.df, aes(x=CancerType, y=QuiescenceScore, fill=CDKN1A_expression)) + 
  geom_boxplot() +theme_classic() + xlab("Cancer Type") + ylab("Quiescence Score")
p + stat_compare_means(aes(group = CDKN1A_expression), label = "p.signif") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
