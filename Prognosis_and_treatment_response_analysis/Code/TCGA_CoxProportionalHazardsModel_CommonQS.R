##########################################################################
#########TCGA Cox Proportional Hazards Model 
##########################################################################

#Load required packages:
library(gtsummary)
library("survival")
library("knitr")
library("survutils")
library("dplyr")
library(maxstat)
library(survminer)
library(forestmodel)

###########################################################################
#Load clinical information and clean up:
setwd("~/Documents/GitHub/CancerDormancy/Data/ClinicalData")
clinical.data <- read.csv("combined_study_clinical_data.csv", header = TRUE)
#(data was downloaded from  cbioportal)
#Select columns of interest
columns_of_interest <- c("Patient.ID","Sample.ID","Sample.Type","Diagnosis.Age","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Cancer.Type",
                         "TCGA.PanCanAtlas.Cancer.Type.Acronym","Disease.Free..Months.","Disease.Free.Status","Months.of.disease.specific.survival",
                         "Disease.specific.Survival.status","Neoadjuvant.Therapy.Type.Administered.Prior.To.Resection.Text","MSI.MANTIS.Score",
                         "MSIsensor.Score","Mutation.Count","New.Neoplasm.Event.Post.Initial.Therapy.Indicator","Overall.Survival..Months.",
                         "Overall.Survival.Status","Progress.Free.Survival..Months.","Progression.Free.Status","Radiation.Therapy",
                         "Sex","Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code","American.Joint.Committee.on.Cancer.Tumor.Stage.Code",
                         "Primary.Lymph.Node.Presentation.Assessment","Aneuploidy.Score")
clinical.data <- clinical.data[,colnames(clinical.data) %in% columns_of_interest]
#Select primary tumour samples
clinical.data <- clinical.data[clinical.data$Sample.Type %in% "Primary",]
rownames(clinical.data) <- clinical.data$Patient.ID
#Load quiescence scores and mutation rates
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
load("TCGA_common_QS_purity_scaled_and_TMB.RData")
Pancancer_samples$cluster <- NULL
Pancancer_samples$SampleID <- NULL
Pancancer_samples$CancerType <- NULL
Pancancer_samples$Mutation_rate <- NULL
Pancancer_samples$Total_numer_of_somatic_mutations <- NULL
#Merge dataframes
Pancancer_samples$PatientID <- sapply(rownames(Pancancer_samples), function(x)
  paste(strsplit(x,"-")[[1]][1:3],collapse="-"))
clinical.data <- merge(clinical.data, Pancancer_samples,
                       by.x = "Patient.ID", by.y= "PatientID")
#Select cancer types of interest:
CT <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
clinical.data <- clinical.data[clinical.data$TCGA.PanCanAtlas.Cancer.Type.Acronym %in% CT,]
#Group tumours according to their stage:
table(clinical.data$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code)
clinical.data$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code <- as.character(clinical.data$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code)
clinical.data$Tumour_stage <- sapply(clinical.data$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code, function(x)
  ifelse(x %in% c("STAGE I","STAGE IA","STAGE IB"),"Stage 1",
         ifelse(x %in% c("STAGE II","STAGE IIA","STAGE IIB","STAGE IIC"),"Stage 2",
                ifelse(x %in% c("STAGE III","STAGE IIIA","STAGE IIIB","STAGE IIIC"),"Stage 3",
                       ifelse(x %in% c("STAGE IV","STAGE IVA","STAGE IVB","STAGE IVC"),"Stage 4","Other")))))
clinical.data <- clinical.data[clinical.data$Tumour_stage %in% c("Stage 1","Stage 2","Stage 3","Stage 4"),]
table(clinical.data$Tumour_stage)
clinical.data$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code <- NULL
#Convert features to character or numerical variables:
clinical.data$Cancer.Type <- as.character(clinical.data$TCGA.PanCanAtlas.Cancer.Type.Acronym)
clinical.data$TCGA.PanCanAtlas.Cancer.Type.Acronym <- NULL
clinical.data$Neoadjuvant.Therapy.Type.Administered.Prior.To.Resection.Text <- NULL
clinical.data$MSI.MANTIS.Score <- NULL
clinical.data$MSIsensor.Score <- NULL
clinical.data$Sex <- as.character(clinical.data$Sex)
clinical.data$Radiation.Therapy <- NULL
#Select patients for which their age,sex, tumour burden is known:
clinical.data <- clinical.data[complete.cases(clinical.data[ , "Diagnosis.Age"]),]
clinical.data <- clinical.data[complete.cases(clinical.data[ , "Sex"]),]




#####################################################
#Survival analysis - Disease specific survival status:
#####################################################
colnames(clinical.data)
test.data <- clinical.data
table(test.data$Disease.specific.Survival.status)
test.data$Disease.specific.Survival.status <- as.character(test.data$Disease.specific.Survival.status)
test.data$Disease.specific.Survival.status <- sapply(test.data$Disease.specific.Survival.status, function(x)
  ifelse(x %in% c("0:ALIVE OR DEAD TUMOR FREE"),0, 
         ifelse(x %in% "1:DEAD WITH TUMOR",1,2)))
test.data <- test.data[test.data$Disease.specific.Survival.status %in% c(0,1),]
#Categorize the data/find the optimal cutoff
res.cut <- surv_cutpoint(test.data, time = "Months.of.disease.specific.survival", event = "Disease.specific.Survival.status",
                         variables = c("z_score"))
summary(res.cut)
plot(res.cut, "z_score", palette = "npg")
res.cat <- surv_categorize(res.cut)
test.data$quiescence_group <- res.cat$z_score

#Forest plot
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Figures/")
pdf("TCGA_CoxProportionalHazardsModel.pdf", height = 7, width = 10)
print(forest_model(coxph(Surv(Months.of.disease.specific.survival, Disease.specific.Survival.status) ~ quiescence_group + Cancer.Type + Diagnosis.Age + Tumour_stage + Sex + log_Mutation_rate, test.data)))
dev.off()

#Kaplan-meier survival curve
#(Use a cutoff of 180 months)
test.data.1 <- test.data[test.data$Months.of.disease.specific.survival >= 180,]
test.data.2 <- test.data[test.data$Months.of.disease.specific.survival < 180,]
test.data.1$Months.of.disease.specific.survival <- 180
test.data.1$Disease.specific.Survival.status <- 0
test.data <- rbind(test.data.1, test.data.2)
res.cox <- coxph(Surv(Months.of.disease.specific.survival, Disease.specific.Survival.status) ~ strata(quiescence_group) + Cancer.Type + Diagnosis.Age + Tumour_stage + Sex + log_Mutation_rate, test.data)
pdf("TCGA_KaplanMeier_SurvivalCurve.pdf",height = 5,width = 5,onefile = FALSE)
d.cox<-survfit(res.cox)
p1<-ggsurvplot(d.cox,data = test.data,conf.int=FALSE, risk.table=TRUE)
print(p1)
dev.off()

