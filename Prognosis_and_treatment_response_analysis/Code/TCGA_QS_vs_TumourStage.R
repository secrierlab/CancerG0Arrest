#######################################################################
####This code looks at how quiescence scores vary across tumour stage:
#######################################################################

#Load required R packages:
library(ggpubr)

###########
#Load clinical information
###########
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_ClinicalData/")
clinical.data <- read.csv("combined_study_clinical_data.csv", header = TRUE)
columns_of_interest <- c("Patient.ID","Sample.ID","Sample.Type","Diagnosis.Age","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Cancer.Type",
                         "TCGA.PanCanAtlas.Cancer.Type.Acronym","Disease.Free..Months.","Disease.Free.Status","Months.of.disease.specific.survival",
                         "Disease.specific.Survival.status","Neoadjuvant.Therapy.Type.Administered.Prior.To.Resection.Text","MSI.MANTIS.Score",
                         "MSIsensor.Score","Mutation.Count","New.Neoplasm.Event.Post.Initial.Therapy.Indicator","Overall.Survival..Months.",
                         "Overall.Survival.Status","Progress.Free.Survival..Months.","Progression.Free.Status","Radiation.Therapy",
                         "Sex","Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code","American.Joint.Committee.on.Cancer.Tumor.Stage.Code",
                         "Primary.Lymph.Node.Presentation.Assessment","Aneuploidy.Score")
clinical.data <- clinical.data[,colnames(clinical.data) %in% columns_of_interest]


##########################
##Load quiescence scores
#########################
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
load("TCGA_common_QS_purity_scaled.RData")
z_score$PatientID <- sapply(rownames(z_score), function(x)
  paste(strsplit(x,"-")[[1]][1:3],collapse="-"))


###########################
#Combine the two dataframes:
###########################
clinical.data <- merge(z_score, clinical.data,
                       by.x = "PatientID", by.y = "Patient.ID")



########################################
#Group tumours according to their stage:
########################################
table(clinical.data$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code)
clinical.data$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code <- as.character(clinical.data$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code)
clinical.data$Tumour_stage <- sapply(clinical.data$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code, function(x)
  ifelse(x %in% c("STAGE I","STAGE IA","STAGE IB"),"Stage 1",
         ifelse(x %in% c("STAGE II","STAGE IIA","STAGE IIB","STAGE IIC"),"Stage 2",
                ifelse(x %in% c("STAGE III","STAGE IIIA","STAGE IIIB","STAGE IIIC"),"Stage 3",
                       ifelse(x %in% c("STAGE IV","STAGE IVA","STAGE IVB","STAGE IVC"),"Stage 4","Other")))))
clinical.data <- clinical.data[clinical.data$Tumour_stage %in% c("Stage 1","Stage 2","Stage 3","Stage 4"),]
table(clinical.data$Tumour_stage)




########################################################################
#Plot a boxplot comparison of quiescence scores across tumour stages:
########################################################################
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Figures/")
my_comparisons <- list( c("Stage 1","Stage 2"), c("Stage 1","Stage 3"), c("Stage 1","Stage 4"), c("Stage 2","Stage 3"), c("Stage 2","Stage 4"), c("Stage 3", "Stage 4") )
p <- ggboxplot(clinical.data, x = "Tumour_stage", y = "z_score",
               fill = "Tumour_stage", shape = "Tumour_stage", order = c("Stage 1","Stage 2", "Stage 3","Stage 4"))
pdf("TCGA_QS_vs_TumourStage.pdf",height = 10, width = 5)
p+stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 60)                   # Add global p-value
dev.off()

