################################################
##BRCA_quiescence_apobec_sbs2 association_heatmap
################################################

#Load required packages:
library(pheatmap)




####Load mutational signatures 
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_MutationalSignatures/")
load("sigs.defaultnounknown.BRCA.RData")
sigs.defaultnounknown$Barcode <- rownames(sigs.defaultnounknown)
sigs.defaultnounknown$PatientID <- sapply(sigs.defaultnounknown$Barcode, function(x)
  paste(strsplit(x,"-")[[1]][1:3],collapse="-"))
Summary <- sigs.defaultnounknown[,colnames(sigs.defaultnounknown) %in% c("SBS2","SBS13","Barcode","PatientID")]


###Add quiescence Scores:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
load("TCGA_common_QS_purity_scaled.RData")
z_score$PatientID <- sapply(rownames(z_score), function(x)
  paste(strsplit(x,"-")[[1]][1:3],collapse="-"))
Summary <- merge(Summary, z_score,
                 by.x = "PatientID", by.y = "PatientID")


#########################
#Add Subtype information:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_ClinicalData/")
BRCA.clinical <- read.csv("combined_study_clinical_data_brca.csv", header = TRUE)
table(BRCA.clinical$Subtype)
BRCA_Basal <- BRCA.clinical[BRCA.clinical$Subtype %in% "BRCA_Basal",]
BRCA_Basal <- as.character(BRCA_Basal$Patient.ID)
BRCA_Her2 <- BRCA.clinical[BRCA.clinical$Subtype %in% "BRCA_Her2",]
BRCA_Her2 <- as.character(BRCA_Her2$Patient.ID)
BRCA_LumA <- BRCA.clinical[BRCA.clinical$Subtype %in% "BRCA_LumA",]
BRCA_LumA <- as.character(BRCA_LumA$Patient.ID)
BRCA_LumB <- BRCA.clinical[BRCA.clinical$Subtype %in% "BRCA_LumB",]
BRCA_LumB <- as.character(BRCA_LumB$Patient.ID)
Summary$Subtype <- sapply(Summary$PatientID, function(x)
  ifelse(x %in% BRCA_Basal, "BRCA_Basal",
         ifelse(x %in% BRCA_Her2, "BRCA_Her2",
                ifelse(x %in% BRCA_LumA, "BRCA_LumA",
                       ifelse(x %in% BRCA_LumB, "BRCA_LumB","Other")))))
table(Summary$Subtype)
Summary <- Summary[Summary$Subtype %in% c("BRCA_Basal","BRCA_Her2","BRCA_LumA","BRCA_LumB"),]





###Add enrichment information:
load("apobec.group.replic.Rdata")
apobec.group.replic$samples <- as.character(apobec.group.replic$samples)
apobec.group.replic$PatientID <- sapply(apobec.group.replic$samples, function(x)
  paste(strsplit(x,"-")[[1]][1:3],collapse="-"))
APOBEC_enriched <- apobec.group.replic[apobec.group.replic$Apobec.Total > 50,]
APOBEC_enriched <- as.character(APOBEC_enriched$PatientID)
APOBEC_depleted<- apobec.group.replic[apobec.group.replic$Apobec.Total <= 50,]
APOBEC_depleted <- as.character(APOBEC_depleted$PatientID)
Summary$APOBEC_enriched <- sapply(Summary$PatientID, function(x)
  ifelse(x %in% APOBEC_depleted, "NO",
         ifelse(x %in% APOBEC_enriched,"YES","Other")))
table(Summary$APOBEC_enriched)



#################
##PLOT
################

#SORT BY QS
Summary <- Summary[order(Summary$z_score),]
Summary$Sum <- NULL
Summary$PatientID <- NULL
Summary$Barcode <- NULL
Summary$CancerType <- NULL

results <- Summary[,1:2]
anno <- Summary[,3:5]

setwd("~/Documents/GitHub/CancerDormancy/BRCA_Quiescence/Figures/")
p<-pheatmap(t(results), show_colnames = FALSE, show_rownames = TRUE, cluster_rows = FALSE, cluster_cols = FALSE, annotation = anno, legend = FALSE)
pdf("BRCA_quiescence_apobec_sbs2_association_heatmap.pdf", height = 3,width = 6)
p
dev.off()


