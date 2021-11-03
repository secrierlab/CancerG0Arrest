##########################################################################
#########TCGA survival analysis - CEP89 status expression stratification
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



####################
###Add CEP89 expression:
###################
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
load("combined_experssion_FPKM.RData")
#This a limited example dataframe with only 100 entries 
#To download the full dataset follow instructions in the "TCGA_DataDownload" folder
combined_data$PatientID <- sapply(rownames(combined_data), function(x)
  paste(strsplit(x,"-")[[1]][1:3],collapse="-"))
combined_data <- combined_data[,colnames(combined_data) %in% c("PatientID","ENSG00000121289")]
clinical.data <- merge(clinical.data, combined_data,
                       by.x = "Patient.ID", by.y = "PatientID")



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
clinical.data$CEP89_expr <- clinical.data$ENSG00000121289






############################################################
#########Survival analysis - Disease specific survival status:
############################################################

#Perfrom analysis on a cancer-by-cancer basis

HR <- NULL
UPPER_LIMIT <- NULL
LOWER_LIMIT <- NULL
PVAL <- NULL

#ACC
test.data <- clinical.data[clinical.data$Cancer.Type %in% "ACC",]
table(test.data$Disease.specific.Survival.status)
test.data$Disease.specific.Survival.status <- as.character(test.data$Disease.specific.Survival.status)
test.data$Disease.specific.Survival.status <- sapply(test.data$Disease.specific.Survival.status, function(x)
  ifelse(x %in% c("0:ALIVE OR DEAD TUMOR FREE"),0, 
         ifelse(x %in% "1:DEAD WITH TUMOR",1,2)))
test.data <- test.data[test.data$Disease.specific.Survival.status %in% c(0,1),]
#Try and chategorise the data/find the optimal cutoff
res.cut <- surv_cutpoint(test.data, time = "Months.of.disease.specific.survival", event = "Disease.specific.Survival.status",
                         variables = c("CEP89_expr"))
summary(res.cut)
plot(res.cut, "CEP89_expr", palette = "npg")
res.cat <- surv_categorize(res.cut)
test.data$CEP89_expr <- res.cat$CEP89_expr
test.data$CEP89_expr <- factor(test.data$CEP89_expr, levels = c("low","high"))
setwd("~/Documents/GitHub/CancerDormancy/CEP89_QuiescenceAssociation/Figures/")
pdf("TCGA_ACC_CEP89_CoxProportionalHazardsModel.pdf",height = 5, width = 7)
print(forest_model(coxph(Surv(Months.of.disease.specific.survival, Disease.specific.Survival.status) ~ CEP89_expr + Tumour_stage, test.data)))
dev.off()

#Record model values:
test <- coxph(Surv(Months.of.disease.specific.survival, Disease.specific.Survival.status) ~ CEP89_expr + Tumour_stage, test.data)
test <- summary(test)
hr <- test$coefficients[1,2]
upper_limit <- test$conf.int[1,4]
lower_limit <- test$conf.int[1,3]
pval <- test[[7]][1,5]
HR <- c(HR,hr)
UPPER_LIMIT <- c(UPPER_LIMIT, upper_limit)
LOWER_LIMIT <- c(LOWER_LIMIT, lower_limit)
PVAL <- c(PVAL, pval)


#HNSC
test.data <- clinical.data[clinical.data$Cancer.Type %in% "HNSC",]
table(test.data$Disease.specific.Survival.status)
test.data$Disease.specific.Survival.status <- as.character(test.data$Disease.specific.Survival.status)
test.data$Disease.specific.Survival.status <- sapply(test.data$Disease.specific.Survival.status, function(x)
  ifelse(x %in% c("0:ALIVE OR DEAD TUMOR FREE"),0, 
         ifelse(x %in% "1:DEAD WITH TUMOR",1,2)))
test.data <- test.data[test.data$Disease.specific.Survival.status %in% c(0,1),]
#Try and chategorise the data/find the optimal cutoff
res.cut <- surv_cutpoint(test.data, time = "Months.of.disease.specific.survival", event = "Disease.specific.Survival.status",
                         variables = c("CEP89_expr"))
summary(res.cut)
plot(res.cut, "CEP89_expr", palette = "npg")
res.cat <- surv_categorize(res.cut)
test.data$CEP89_expr <- res.cat$CEP89_expr
test.data$CEP89_expr <- factor(test.data$CEP89_expr, levels = c("low","high"))
setwd("~/Documents/GitHub/CancerDormancy/CEP89_QuiescenceAssociation/Figures/")
pdf("TCGA_HNSC_CEP89_CoxProportionalHazardsModel.pdf",height = 5, width = 7)
print(forest_model(coxph(Surv(Months.of.disease.specific.survival, Disease.specific.Survival.status) ~ CEP89_expr + Tumour_stage, test.data)))
dev.off()

#Record model values
test <- coxph(Surv(Months.of.disease.specific.survival, Disease.specific.Survival.status) ~ CEP89_expr + Tumour_stage, test.data)
test <- summary(test)
hr <- test$coefficients[1,2]
upper_limit <- test$conf.int[1,4]
lower_limit <- test$conf.int[1,3]
pval <- test[[7]][1,5]
HR <- c(HR,hr)
UPPER_LIMIT <- c(UPPER_LIMIT, upper_limit)
LOWER_LIMIT <- c(LOWER_LIMIT, lower_limit)
PVAL <- c(PVAL, pval)


#KIRC
test.data <- clinical.data[clinical.data$Cancer.Type %in% "KIRC",]
table(test.data$Disease.specific.Survival.status)
test.data$Disease.specific.Survival.status <- as.character(test.data$Disease.specific.Survival.status)
test.data$Disease.specific.Survival.status <- sapply(test.data$Disease.specific.Survival.status, function(x)
  ifelse(x %in% c("0:ALIVE OR DEAD TUMOR FREE"),0, 
         ifelse(x %in% "1:DEAD WITH TUMOR",1,2)))
test.data <- test.data[test.data$Disease.specific.Survival.status %in% c(0,1),]
#Try and chategorise the data/find the optimal cutoff
res.cut <- surv_cutpoint(test.data, time = "Months.of.disease.specific.survival", event = "Disease.specific.Survival.status",
                         variables = c("CEP89_expr"))
summary(res.cut)
plot(res.cut, "CEP89_expr", palette = "npg")
res.cat <- surv_categorize(res.cut)
test.data$CEP89_expr <- res.cat$CEP89_expr
test.data$CEP89_expr <- factor(test.data$CEP89_expr, levels = c("low","high"))
setwd("~/Documents/GitHub/CancerDormancy/CEP89_QuiescenceAssociation/Figures/")
pdf("TCGA_KIRC_CEP89_CoxProportionalHazardsModel.pdf",height = 5, width = 7)
print(forest_model(coxph(Surv(Months.of.disease.specific.survival, Disease.specific.Survival.status) ~ CEP89_expr + Tumour_stage, test.data)))
dev.off()

#Record model values
test <- coxph(Surv(Months.of.disease.specific.survival, Disease.specific.Survival.status) ~ CEP89_expr + Tumour_stage, test.data)
test <- summary(test)
hr <- test$coefficients[1,2]
upper_limit <- test$conf.int[1,4]
lower_limit <- test$conf.int[1,3]
pval <- test[[7]][1,5]
HR <- c(HR,hr)
UPPER_LIMIT <- c(UPPER_LIMIT, upper_limit)
LOWER_LIMIT <- c(LOWER_LIMIT, lower_limit)
PVAL <- c(PVAL, pval)


#KIRP
test.data <- clinical.data[clinical.data$Cancer.Type %in% "KIRP",]
table(test.data$Disease.specific.Survival.status)
test.data$Disease.specific.Survival.status <- as.character(test.data$Disease.specific.Survival.status)
test.data$Disease.specific.Survival.status <- sapply(test.data$Disease.specific.Survival.status, function(x)
  ifelse(x %in% c("0:ALIVE OR DEAD TUMOR FREE"),0, 
         ifelse(x %in% "1:DEAD WITH TUMOR",1,2)))
test.data <- test.data[test.data$Disease.specific.Survival.status %in% c(0,1),]
#Try and chategorise the data/find the optimal cutoff
res.cut <- surv_cutpoint(test.data, time = "Months.of.disease.specific.survival", event = "Disease.specific.Survival.status",
                         variables = c("CEP89_expr"))
summary(res.cut)
plot(res.cut, "CEP89_expr", palette = "npg")
res.cat <- surv_categorize(res.cut)
test.data$CEP89_expr <- res.cat$CEP89_expr
test.data$CEP89_expr <- factor(test.data$CEP89_expr, levels = c("low","high"))
setwd("~/Documents/GitHub/CancerDormancy/CEP89_QuiescenceAssociation/Figures/")
pdf("TCGA_KIRP_CEP89_CoxProportionalHazardsModel.pdf",height = 5, width = 7)
print(forest_model(coxph(Surv(Months.of.disease.specific.survival, Disease.specific.Survival.status) ~ CEP89_expr + Tumour_stage, test.data)))
dev.off()

#Record model values:
test <- coxph(Surv(Months.of.disease.specific.survival, Disease.specific.Survival.status) ~ CEP89_expr + Tumour_stage, test.data)
test <- summary(test)
hr <- test$coefficients[1,2]
upper_limit <- test$conf.int[1,4]
lower_limit <- test$conf.int[1,3]
pval <- test[[7]][1,5]
HR <- c(HR,hr)
UPPER_LIMIT <- c(UPPER_LIMIT, upper_limit)
LOWER_LIMIT <- c(LOWER_LIMIT, lower_limit)
PVAL <- c(PVAL, pval)



#LIHC
test.data <- clinical.data[clinical.data$Cancer.Type %in% "LIHC",]
table(test.data$Disease.specific.Survival.status)
test.data$Disease.specific.Survival.status <- as.character(test.data$Disease.specific.Survival.status)
test.data$Disease.specific.Survival.status <- sapply(test.data$Disease.specific.Survival.status, function(x)
  ifelse(x %in% c("0:ALIVE OR DEAD TUMOR FREE"),0, 
         ifelse(x %in% "1:DEAD WITH TUMOR",1,2)))
test.data <- test.data[test.data$Disease.specific.Survival.status %in% c(0,1),]
#Try and chategorise the data/find the optimal cutoff
res.cut <- surv_cutpoint(test.data, time = "Months.of.disease.specific.survival", event = "Disease.specific.Survival.status",
                         variables = c("CEP89_expr"))
summary(res.cut)
plot(res.cut, "CEP89_expr", palette = "npg")
res.cat <- surv_categorize(res.cut)
test.data$CEP89_expr <- res.cat$CEP89_expr
test.data$CEP89_expr <- factor(test.data$CEP89_expr, levels = c("low","high"))
setwd("~/Documents/GitHub/CancerDormancy/CEP89_QuiescenceAssociation/Figures/")
pdf("TCGA_LIHC_CEP89_CoxProportionalHazardsModel.pdf",height = 5, width = 7)
print(forest_model(coxph(Surv(Months.of.disease.specific.survival, Disease.specific.Survival.status) ~ CEP89_expr + Tumour_stage, test.data)))
dev.off()

#Record model values
test <- coxph(Surv(Months.of.disease.specific.survival, Disease.specific.Survival.status) ~ CEP89_expr + Tumour_stage, test.data)
test <- summary(test)
hr <- test$coefficients[1,2]
upper_limit <- test$conf.int[1,4]
lower_limit <- test$conf.int[1,3]
pval <- test[[7]][1,5]
HR <- c(HR,hr)
UPPER_LIMIT <- c(UPPER_LIMIT, upper_limit)
LOWER_LIMIT <- c(LOWER_LIMIT, lower_limit)
PVAL <- c(PVAL, pval)

#LUSC
test.data <- clinical.data[clinical.data$Cancer.Type %in% "LUSC",]
table(test.data$Disease.specific.Survival.status)
test.data$Disease.specific.Survival.status <- as.character(test.data$Disease.specific.Survival.status)
test.data$Disease.specific.Survival.status <- sapply(test.data$Disease.specific.Survival.status, function(x)
  ifelse(x %in% c("0:ALIVE OR DEAD TUMOR FREE"),0, 
         ifelse(x %in% "1:DEAD WITH TUMOR",1,2)))
test.data <- test.data[test.data$Disease.specific.Survival.status %in% c(0,1),]
#Try and chategorise the data/find the optimal cutoff
res.cut <- surv_cutpoint(test.data, time = "Months.of.disease.specific.survival", event = "Disease.specific.Survival.status",
                         variables = c("CEP89_expr"))
summary(res.cut)
plot(res.cut, "CEP89_expr", palette = "npg")
res.cat <- surv_categorize(res.cut)
test.data$CEP89_expr <- res.cat$CEP89_expr
test.data$CEP89_expr <- factor(test.data$CEP89_expr, levels = c("low","high"))
setwd("~/Documents/GitHub/CancerDormancy/CEP89_QuiescenceAssociation/Figures/")
pdf("TCGA_LUSC_CEP89_CoxProportionalHazardsModel.pdf",height = 5, width = 7)
print(forest_model(coxph(Surv(Months.of.disease.specific.survival, Disease.specific.Survival.status) ~ CEP89_expr + Tumour_stage, test.data)))
dev.off()

#Record model values
test <- coxph(Surv(Months.of.disease.specific.survival, Disease.specific.Survival.status) ~ CEP89_expr + Tumour_stage, test.data)
test <- summary(test)
hr <- test$coefficients[1,2]
upper_limit <- test$conf.int[1,4]
lower_limit <- test$conf.int[1,3]
pval <- test[[7]][1,5]
HR <- c(HR,hr)
UPPER_LIMIT <- c(UPPER_LIMIT, upper_limit)
LOWER_LIMIT <- c(LOWER_LIMIT, lower_limit)
PVAL <- c(PVAL, pval)


#PAAD
test.data <- clinical.data[clinical.data$Cancer.Type %in% "PAAD",]
table(test.data$Disease.specific.Survival.status)
test.data$Disease.specific.Survival.status <- as.character(test.data$Disease.specific.Survival.status)
test.data$Disease.specific.Survival.status <- sapply(test.data$Disease.specific.Survival.status, function(x)
  ifelse(x %in% c("0:ALIVE OR DEAD TUMOR FREE"),0, 
         ifelse(x %in% "1:DEAD WITH TUMOR",1,2)))
test.data <- test.data[test.data$Disease.specific.Survival.status %in% c(0,1),]
#Try and chategorise the data/find the optimal cutoff
res.cut <- surv_cutpoint(test.data, time = "Months.of.disease.specific.survival", event = "Disease.specific.Survival.status",
                         variables = c("CEP89_expr"))
summary(res.cut)
plot(res.cut, "CEP89_expr", palette = "npg")
res.cat <- surv_categorize(res.cut)
test.data$CEP89_expr <- res.cat$CEP89_expr
test.data$CEP89_expr <- factor(test.data$CEP89_expr, levels = c("low","high"))
setwd("~/Documents/GitHub/CancerDormancy/CEP89_QuiescenceAssociation/Figures/")
pdf("TCGA_PAAD_CEP89_CoxProportionalHazardsModel.pdf",height = 5, width = 7)
print(forest_model(coxph(Surv(Months.of.disease.specific.survival, Disease.specific.Survival.status) ~ CEP89_expr + Tumour_stage, test.data)))
dev.off()

#Record model values
test <- coxph(Surv(Months.of.disease.specific.survival, Disease.specific.Survival.status) ~ CEP89_expr + Tumour_stage, test.data)
test <- summary(test)
hr <- test$coefficients[1,2]
upper_limit <- test$conf.int[1,4]
lower_limit <- test$conf.int[1,3]
pval <- test[[7]][1,5]
HR <- c(HR,hr)
UPPER_LIMIT <- c(UPPER_LIMIT, upper_limit)
LOWER_LIMIT <- c(LOWER_LIMIT, lower_limit)
PVAL <- c(PVAL, pval)


#STAD
table(clinical.data$Cancer.Type)
test.data <- clinical.data[clinical.data$Cancer.Type %in% "STAD",]
table(test.data$Disease.specific.Survival.status)
test.data$Disease.specific.Survival.status <- as.character(test.data$Disease.specific.Survival.status)
test.data$Disease.specific.Survival.status <- sapply(test.data$Disease.specific.Survival.status, function(x)
  ifelse(x %in% c("0:ALIVE OR DEAD TUMOR FREE"),0, 
         ifelse(x %in% "1:DEAD WITH TUMOR",1,2)))
test.data <- test.data[test.data$Disease.specific.Survival.status %in% c(0,1),]
#Try and chategorise the data/find the optimal cutoff
res.cut <- surv_cutpoint(test.data, time = "Months.of.disease.specific.survival", event = "Disease.specific.Survival.status",
                         variables = c("CEP89_expr"))
summary(res.cut)
plot(res.cut, "CEP89_expr", palette = "npg")
res.cat <- surv_categorize(res.cut)
test.data$CEP89_expr <- res.cat$CEP89_expr
test.data$CEP89_expr <- factor(test.data$CEP89_expr, levels = c("low","high"))
setwd("~/Documents/GitHub/CancerDormancy/CEP89_QuiescenceAssociation/Figures/")
pdf("TCGA_STAD_CEP89_CoxProportionalHazardsModel.pdf",height = 5, width = 7)
print(forest_model(coxph(Surv(Months.of.disease.specific.survival, Disease.specific.Survival.status) ~ CEP89_expr + Tumour_stage, test.data)))
dev.off()

#Record model values
test <- coxph(Surv(Months.of.disease.specific.survival, Disease.specific.Survival.status) ~ CEP89_expr + Tumour_stage, test.data)
test <- summary(test)
hr <- test$coefficients[1,2]
upper_limit <- test$conf.int[1,4]
lower_limit <- test$conf.int[1,3]
pval <- test[[7]][1,5]
HR <- c(HR,hr)
UPPER_LIMIT <- c(UPPER_LIMIT, upper_limit)
LOWER_LIMIT <- c(LOWER_LIMIT, lower_limit)
PVAL <- c(PVAL, pval)






CancerType <- c("ACC","HNSC","KIRC","KIRP","LIHC","LUSC","PAAD","STAD")
Survival_analysus_results <- data.frame(CancerType, HR, UPPER_LIMIT, LOWER_LIMIT, PVAL)


######################################
#Plot the results (only significant)
#####################################
CancerType <- c("ACC","HNSC","KIRC","KIRP","LIHC","LUSC","PAAD","STAD")
Survival_analysus_results <- data.frame(CancerType, HR, UPPER_LIMIT, LOWER_LIMIT, PVAL)
Survival_analysus_results <- Survival_analysus_results[Survival_analysus_results$PVAL < 0.05,]
Survival_analysus_results <- Survival_analysus_results[order(Survival_analysus_results$HR),]
Survival_analysus_results$yaxis <- 1:8
Survival_analysus_results$HR <- log2(Survival_analysus_results$HR)
Survival_analysus_results$UPPER_LIMIT <- log2(Survival_analysus_results$UPPER_LIMIT)
Survival_analysus_results$LOWER_LIMIT <- log2(Survival_analysus_results$LOWER_LIMIT)
#Plot the results:
pdf("TCGA_CEP89_CoxProportionalHazardsModel_CancerwiseSummary.pdf", width = 4, height = 7)
p <-ggplot(Survival_analysus_results, aes(x = HR, y = yaxis))
p + geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = UPPER_LIMIT, xmin = LOWER_LIMIT), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange") +
  theme_classic() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = Survival_analysus_results$yaxis, labels = Survival_analysus_results$CancerType) +
  scale_x_continuous(limits = c(-5,5))+
  ylab("Cancer Type") +
  xlab("log2 Hazard Ratio")
dev.off()


