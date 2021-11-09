###################################################################
#Association between quiescence and microsatelite instability
###################################################################

#Load required packages:
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(grid)


####Load quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
load("TCGA_common_QS_purity_scaled.RData")
z_score$PatientID <- sapply(rownames(z_score), function(x)
  paste(strsplit(x,"-")[[1]][1:3],collapse="-"))
Pancancer_samples <- z_score
Pancancer_samples$quiescence_score <- Pancancer_samples$z_score


#####Load MSI event data:
#NOTE: this data was downloaded from Cotez-Ciriano et al 2017 
MSI.data <- read.csv("41467_2017_BFncomms15180_MOESM259_ESM.csv",header = TRUE)
MSI.data$Log_total_nb_MSI_events <- log2(MSI.data$Total_nb_MSI_events + 1)
summary(MSI.data$Log_total_nb_MSI_events)

#Merge the data:
Merged.data <- merge(Pancancer_samples, MSI.data,
                     by.x = "PatientID", by.y = "Barcode")


#Clean up the dataframe:
rownames(Merged.data) <- Merged.data$PatientID
Merged.data <- Merged.data[,colnames(Merged.data) %in% c("Cancer_type","MSI_category_nb_from_TCGA_consortium","quiescence_score")]

#Select cancer types which have MSI annotation
Merged.data <- Merged.data[Merged.data$MSI_category_nb_from_TCGA_consortium %in% c("msi-h","msi-l","mss"),]
Merged.data$MSI_category_nb_from_TCGA_consortium <- as.character(Merged.data$MSI_category_nb_from_TCGA_consortium)
table(Merged.data$MSI_category_nb_from_TCGA_consortium)
table(Merged.data$Cancer_type)
Merged.data$Cancer_type <- as.character(Merged.data$Cancer_type)
Merged.data <- Merged.data[Merged.data$Cancer_type %in% c("COAD","STAD","UCEC"),]

#Group the MSI-l and MSI-h categories
Merged.data$MSI.status <- sapply(Merged.data$MSI_category_nb_from_TCGA_consortium, function(x)
  ifelse(x %in% c("msi-h","msi-l"),"MSI",
         ifelse(x %in% "mss","MSS","other")))
table(Merged.data$MSI.status)
Merged.data$MSI_category_nb_from_TCGA_consortium <- NULL


#Plot association betweem MSI status and quiescence scores:
table(Merged.data$Cancer_type)
text_high <- textGrob("n = 224", gp=gpar(fontsize=10, fontface="plain"))
text_medium <- textGrob("n = 193", gp=gpar(fontsize=10, fontface="plain"))
text_low <- textGrob("n = 250",gp=gpar(fontsize=10, fontface="plain"))
myCol <- brewer.pal(2, "Set2")
setwd("~/Documents/GitHub/CancerDormancy/PanCancer_ElasticNetRegression/FeatureExploration/Figures/")
pdf("MSI_vs_quiescence.pdf",height = 9,width = 9)
my_comparisons <- list( c("MSS", "MSI"))
p <- ggplot(Merged.data, aes(x=Cancer_type, y=quiescence_score, fill=MSI.status)) + 
  geom_boxplot() +theme_classic() + xlab("\n\nCancer Type") + ylab("Quiescence Score")  + scale_fill_manual(values = myCol) +
  annotation_custom(text_low, xmin = 1, xmax = 1, ymin = -23, ymax = -23) + 
  annotation_custom(text_medium, xmin = 2, xmax = 2, ymin = -23, ymax = -23) + 
  annotation_custom(text_high, xmin = 3, xmax = 3, ymin = -23, ymax = -23)
p + stat_compare_means(label.y = 15, label = "p.signif") + coord_cartesian(clip = "off")    # Add global p-value

dev.off()
