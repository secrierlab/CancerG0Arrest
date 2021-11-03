########################################################################
####TCGA scaling expression data by purity
########################################################################


#Load the required packages
library("MOFA")
library("MOFAdata")
library(MultiAssayExperiment)



####################################################################
######Load tumour purity data and identify samples with purity > 0.3
tumour.purity <- read.table("TCGA_mastercalls.abs_tables_JSedit.fixed.txt", header = TRUE, sep = "\t")
#tumour purity data was downloaded from the PanCanAtlas Publications
#https://gdc.cancer.gov/about-data/publications/pancanatlas
tumour.purity <- tumour.purity[tumour.purity$purity >= 0.3,]
tumour.purity$sample <- as.character(tumour.purity$sample)
tumour.purity$SampleID <- sapply(tumour.purity$sample, function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
pure.samples <- as.character(tumour.purity$SampleID)



#########################################################################
####Load expression data and select samples with purity greater than 0.3
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
load("combined_experssion_FPKM.RData")
#This a limited example dataframe with only 100 entries 
#To download the full dataset follow instructions in the "TCGA_DataDownload" folder
combined_data$SampleID <- sapply(rownames(combined_data), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
combined_data$PatientID <- sapply(rownames(combined_data), function(x)
  paste(strsplit(x,"-")[[1]][1:3],collapse="-"))
combined_data <- combined_data[combined_data$SampleID %in% pure.samples,]
combined_data$Barcode <- rownames(combined_data)



######################################################################
######Some patients have two or more samples (randomly select only one)
Barcodes <- as.character(combined_data$Barcode)
Patients <- as.character(combined_data$PatientID)
uniquePatients <- unique(as.character(combined_data$PatientID))
selection.df <- data.frame(Barcodes, Patients)
Barcodes <- NULL
set.seed(123)
for (i in uniquePatients) {
  print(i)
  selected_data <- selection.df[selection.df$Patients %in% i,]
  selected_data <- as.character(selected_data$Barcodes)
  random_sample <- sample(selected_data, size = 1, replace = FALSE)
  Barcodes <- c(Barcodes, random_sample)
  rm(random_sample)
}
combined_data <- combined_data[combined_data$Barcode %in% Barcodes,]




########################################################################
##Scale the expression data by tumour purity on a cancer by cancer basis

CT <- unique(as.character(combined_data$cancer_type))

for (a in CT) {
  
  print(a)
  
  #Select data from the cancer type of interest
  selected.data <- combined_data[combined_data$cancer_type %in% a,]
  
  #Merge data with purity scores
  selected.data <- merge(selected.data, tumour.purity,
                         by.x = "SampleID", by.y = "SampleID")
  
  #Create the colname information
  Purity <- selected.data$purity
  Purity <- data.frame(Purity)
  rownames(Purity) <- selected.data$Barcode
  ColData <- Purity
  
  #Prepare the RNA expression data
  rownames(selected.data) <- selected.data$Barcode
  selected.data$cancer_type <- NULL
  selected.data$SampleID <- NULL
  selected.data$ploidy <- NULL
  selected.data$purity <- NULL
  selected.data$Genome.doublings <- NULL
  selected.data$call.status <- NULL
  selected.data$Cancer.DNA.fraction <- NULL
  selected.data$Subclonal.genome.fraction <- NULL
  selected.data$solution <- NULL
  selected.data$Barcode <- NULL
  selected.data$sample <- NULL
  selected.data$Coverage.for.80..power <- NULL
  selected.data$PatientID <- NULL
  selected.data$array <- NULL
  selected.data <- data.frame(t(selected.data))
  colnames(selected.data) <- gsub('\\.', '-', colnames(selected.data))
  all(rownames(ColData) == colnames(selected.data))
  Data <- list(selected.data)
  names(Data)[1] <- "mRNA"
  
  #Creae a multiassayexperiment object
  
  rna.seq.scale <- MultiAssayExperiment(
    experiments = Data, 
    colData = ColData
  )
  
  MOFAobject <- createMOFAobject(rna.seq.scale)
  MOFAobject <- prepareMOFA(MOFAobject)
  MOFAobject_reg <- regressCovariates(
    object = MOFAobject,
    views = c("mRNA"),
    covariates = InputData(MOFAobject)$Purity
  )
  
  #Extract the scaled RNA seq data:
  
  scaled.data <- MOFAobject_reg@TrainData$mRNA 
  
  #Save the scaled data
  setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
  save(scaled.data, file = paste("scaled_expr_data_",a,".RData",sep = ""))
  
}


#############################################################################
#Combine scaled RNA-seq data from various cancer studies

setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
load("scaled_expr_data_ACC.RData")
combined.scaled <- scaled.data
combined.scaled <- data.frame(t(combined.scaled))
combined.scaled$CancerType <- "ACC"
rm(scaled.data)
CT <- CT[-1]

for (a in CT) {
  
  print(a)
  
  load(paste("scaled_expr_data_",a,".RData",sep = ""))
  
  scaled.data <- data.frame(t(scaled.data))
  scaled.data$CancerType <- a
  
  combined.scaled <- rbind(combined.scaled, scaled.data)
  
}


####Save the combined data:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
save(combined.scaled, file = "combined_scaled_expr_data.RData")



