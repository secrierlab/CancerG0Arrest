######################################################
##TCGA CEP89 and CA20 association
######################################################

#Load required pacakges
library(ggplot2)
library(ggpubr)

#Select TCGA studies for the analysis:
#
CT <- c( 'ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'ESCA',  'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP',  'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC',  'UCS','UVM' 
)


#######
#Load CA20 scores:
CA20 <- read.table("journal.pcbi.1006832.s018.txt", header = TRUE, sep = "\t")
#CA20 scores were obtained from de Almeida et al 2019, Plos Computational Biology
CA20$Sample.ID <- gsub('\\.', '-', CA20$Sample.ID)
CA20$PatientID <- sapply(CA20$Sample.ID, function(x)
  paste(strsplit(x,"-")[[1]][1:3],collapse="-"))

#########
#Load CNV data:
#Load the CNV data for ACC
setwd("~/Documents/GitHub/CancerDormancy/Data/")
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
#Identify patients with CEP89 amplifications
CNV_samples <- sapply(rownames(CNV), function(x)
  paste(strsplit(x,"-")[[1]][1:3],collapse="-"))
CA20 <- CA20[CA20$PatientID %in% CNV_samples,]
CEP89_samples <- CNV[CNV$ENSG00000121289 %in% c(1,2,3,4,5),]
CEP89_samples <- sapply(rownames(CEP89_samples), function(x)
  paste(strsplit(x,"-")[[1]][1:3],collapse="-"))
CA20$CEP89_status <- sapply(CA20$PatientID, function(x)
  ifelse(x %in% CEP89_samples,"Amp","Other"))
table(CA20$CEP89_status)


#Plot comparison:
setwd("~/Documents/GitHub/CancerDormancy/CEP89_QuiescenceAssociation/Figures/")
pdf("TCGA_CEP89_status_vs_CA20.pdf",height = 5,width = 5)
my_comparisons <- list( c("Amp", "Other"))
p <- ggplot(CA20, aes(x=CEP89_status, y=CA20, fill = CEP89_status)) + 
  geom_boxplot()+
  labs(x="CEP89 Status", y = "CA20 Score")+
  theme_classic()
p +  stat_compare_means(comparisons = my_comparisons, label = "p.signif") # Add pairwise comparisons p-value
dev.off()







