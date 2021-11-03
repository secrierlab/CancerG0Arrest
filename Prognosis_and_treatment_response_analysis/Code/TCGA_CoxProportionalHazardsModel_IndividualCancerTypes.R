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
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_ClinicalData/")
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
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Figures/")

#Perform the analysis on a cancer-by-cancer manner
CT <- c("ACC","KIRP","MESO","LIHC","KIRC","BRCA","LUSC","COAD","ESCA")


for (a in CT) {
  
  print(a)
  test.data <- clinical.data[clinical.data$CancerType %in% a,]
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
  pdf(paste("TCGA_CoxProportionalHazardsModel_",a,".pdf",sep = ""), height = 7, width = 10)
  print(forest_model(coxph(Surv(Months.of.disease.specific.survival, Disease.specific.Survival.status) ~ quiescence_group + Tumour_stage, test.data)))
  dev.off()
  
  
}



##########Summary table:
CancerType <- c("ACC","BRCA","COAD","ESCA","KIRC","KIRP","LIHC","LUSC","MESO")
HazardRatio <- c(9.76,1.76,0.37,0.36,2.31,7.49,2.51,0.42,4.55)
HazardRatio <- log2(HazardRatio)
lower_conf <- c(3.51,1.02,0.21,0.14,1.06,2.95,1.48,0.22,2.18)
lower_conf <- log2(lower_conf)
upper_conf <- c(27.12,3.02,0.68,0.93,5.00,19.04,4.25,0.78,9.49)
upper_conf <- log2(upper_conf)
Summary <- data.frame(CancerType, HazardRatio, lower_conf, upper_conf)
Summary <- Summary[order(Summary$HazardRatio),]
Summary$yaxis <- 1:9
#Plot the results:
pdf("TCGA_CoxProportionalHazardsModel_CancerTypeSummary.pdf", width = 4, height = 3)
p <-ggplot(Summary, aes(x = HazardRatio, y = yaxis))
p + geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = upper_conf, xmin = lower_conf), size = .5, height = .2, color = "gray50") +
  geom_point(size = 1, color = "darkblue") +
  theme_classic() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = Summary$yaxis, labels = Summary$CancerType) +
  scale_x_continuous(limits = c(-3,5))+
  ylab("") +
  xlab("log Hazard Ratio (95% CI)")
dev.off()
