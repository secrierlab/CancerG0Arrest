#########################################################
#TCGA combining CNV data
#########################################################

#Load the CNV data for ACC
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
load("TCGA-ACCCopy_Number_Variation.RData")
CNV <- as.data.frame(data)
rownames(CNV) <- CNV$`Gene Symbol`
CNV$`Gene Symbol` <- NULL
CNV$`Gene ID` <- NULL
CNV$Cytoband <- NULL
CNV <- CNV[,order(colnames(CNV))]
CNV <- data.frame(t(CNV))

CT <- CT[!(CT %in% "ACC")]
#Load combine the CNV dataframe with CNV data for ther remaining cancer types:
for (a in CT) {
  load(paste("TCGA-",a,"Copy_Number_Variation.RData",sep = ""))
  data <- as.data.frame(data)
  rownames(data) <- data$`Gene Symbol`
  data$`Gene Symbol` <- NULL
  data$`Gene ID` <- NULL
  data$Cytoband <- NULL
  data <- data.frame(t(data))
  data <- data[,order(colnames(data))]
  
  #Merge the dataframes:
  CNV <- rbind(CNV, data)
}

#Remove the numbers after decimal points in the ENSG annotation (these denote the version of the gene and can be ignored)
names <- as.character(colnames(CNV))
names2 <- gsub('\\.', '-', names)
colnames(CNV) <- sapply(names2, function(x)
  strsplit(x,"-")[[1]][1])
save(CNV, file = "combined_CNV.RData")


