#######################################################################
#Genome doubling events in relation to quiescene scores (purity scaled)
#######################################################################


###Load the required packages:
library(ggplot2)
library(ggpubr)


###########################
#Load the quiescence scores
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
load("TCGA_common_QS_purity_scaled.RData")

######################################################
#Load WGD event annotation and merge with QS dataframe
Absolute_data <- read.table("TCGA_mastercalls.abs_tables_JSedit.fixed.txt", header = TRUE, sep = "\t")
#(data was downloaded from https://gdc.cancer.gov/about-data/publications/pancanatlas)
#Data was obtained from PanCanAtlas publications
Absolute_data$SampleID <- Absolute_data$array
Absolute_data$array <- NULL
Absolute_data$call.status <- NULL
Absolute_data$purity <- NULL
Absolute_data$Coverage.for.80..power <- NULL
Absolute_data$Cancer.DNA.fraction <- NULL
Absolute_data$Subclonal.genome.fraction <- NULL
Absolute_data$solution <- NULL
Absolute_data$sample <- NULL
#Alter the Merged_data SampleID column to match those reported in the PanCanAtlas publication
Absolute_data$SampleID <- as.character(Absolute_data$SampleID)
z_score$SampleID <- sapply(rownames(z_score), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
z_score$SampleID <-  gsub('\\01A', '1', z_score$SampleID)
z_score$SampleID <-  gsub('\\01B', '1', z_score$SampleID)
z_score$SampleID <-  gsub('\\01C', '1', z_score$SampleID)


#Merge  the dataframes :
Pancancer_samples <- merge(z_score, Absolute_data,
                           by.x = "SampleID", by.y = "SampleID")
colnames(Pancancer_samples)
table(Pancancer_samples$Genome.doublings)
Pancancer_samples <- Pancancer_samples[Pancancer_samples$Genome.doublings %in% c("0","1","2"),]



####################
#Boxplot comparison
setwd("~/Documents/GitHub/CancerDormancy/TCGA_DormancyEvaluation/Figures/")
Pancancer_samples$Genome.doublings <- as.factor(Pancancer_samples$Genome.doublings)
pdf("TCGA_QS_vs_WGD.pdf",height = 5,width = 5)
my_comparisons <- list( c("0", "1"), c("0", "2"), c("1", "2") )
p <- ggplot(Pancancer_samples, aes(x=Genome.doublings, y=z_score, fill = Genome.doublings)) + 
  geom_boxplot()+
  labs(x="Genome doublings", y = "Quiescence score")+
  theme_classic()
p +  stat_compare_means(comparisons = my_comparisons, label = "p.signif") # Add pairwise comparisons p-value
dev.off()

