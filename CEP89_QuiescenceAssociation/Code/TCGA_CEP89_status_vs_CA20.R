######################################################
##TCGA CEP89 and CA20 association
######################################################

#Load required pacakges
library(ggplot2)
library(ggpubr)


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
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
load("combined_CNV.RData")
#This a limited example dataframe with only 100 entries 
#To download the full dataset follow instructions in the "TCGA_DataDownload" folder

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







