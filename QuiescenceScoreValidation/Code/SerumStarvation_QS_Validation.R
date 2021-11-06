################################################
##serum starvation quiescence programme validaton
#Measure correlation between different quiescence scores and serum starvation gene activity

#Load required packages:
library(GSVA)
library(biomaRt)
library(corrplot)
library(reshape)
library(RColorBrewer)
library(ggpubr)


####################################
##Load pan-cancer quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
load("TCGA_common_QS_purity_scaled.RData")
QS <- z_score
QS$CommonQuiescenceScore <- QS$z_score
QS$z_score <- NULL
load("TCGA_contact_inhibition_QS_purity_scaled.RData")
QS$ContactInhibitionScore <- z_score$z_score
load("TCGA_CDK_inhibition_QS_purity_scaled.RData")
QS$CDKInhibitionScore <- z_score$z_score
load("TCGA_MEK_inhibition_QS_purity_scaled.RData")
QS$MEKInhibitionScore <- z_score$z_score
load("TCGA_serum_starvation_QS_purity_scaled.RData")
QS$SerumStarvationScore <- z_score$z_score
load("TCGA_spotaneous_quiescence_QS_purity_scaled.RData")
QS$SpontaneousQuiescenceScore <- z_score$z_score

######################################
#Calculate serum starvation scores:

setwd("~/Documents/GitHub/CancerDormancy/Data/OtherGeneLists/")
Genes <- read.table("GO_term_summary_20210204_100045.csv", header = TRUE, sep = ",")
Genes <- as.character(Genes$Symbol)
Genes <- toupper(Genes)
human = useMart("ensembl", dataset="hsapiens_gene_ensembl")
biomart_conversion <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=Genes, mart=human, filters = "hgnc_symbol")
ENSG_genes <- biomart_conversion$ensembl_gene_id
#Select only genes of interest to reduce the size of the dataframe:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
load("combined_scaled_expr_data.RData") #This is an example dataset with 100 entries only
expr.data <- combined.scaled
expr.data <- expr.data[,which(colnames(expr.data) %in% ENSG_genes)]
#The data is already log transformed and so is ready for z_score calculation:
expr.data <- as.matrix(t(expr.data))
#Now work out the combined z_score for the samples
gene_lists <- list(ENSG_genes)
z_score <- gsva(expr.data, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
all(rownames(z_score) == rownames(QS))
QS$Score <- z_score$z_score



################################################
###work out the correlations between serum starvation scores and the 5 quiescence scores (pancancer and in individual cancer types)
Scores <- c("CommonQuiescenceScore","ContactInhibitionScore","CDKInhibitionScore","MEKInhibitionScore","SerumStarvationScore","SpontaneousQuiescenceScore")
P_val <- NULL
Corr <- NULL
for (i in Scores) {
  test <- cor.test(QS$Score, QS[[i]])
  p_val <- test$p.value
  corr <- test$estimate
  P_val <- c(P_val, p_val)
  Corr <- c(Corr, corr)
}
Correlations.df <- data.frame(Scores)
Correlations.df$Pancancer <- Corr
Pvalue.df <- data.frame(Scores)
Pvalue.df$Pancancer <- P_val
CT <- unique(as.character(QS$CancerType))
CT <- CT[!(CT %in% c("CHOL","READ","TGCT"))]
for (a in CT) {
  P_val <- NULL
  Corr <- NULL
  test.data <- QS[QS$CancerType %in% a,]
  for (i in Scores) {
    test <- cor.test(test.data$Score, test.data[[i]])
    p_val <- test$p.value
    corr <- test$estimate
    P_val <- c(P_val, p_val)
    Corr <- c(Corr, corr)
  }
  Correlations.df[[a]] <- Corr
  Pvalue.df[[a]] <- P_val
  
}
##Clean up the dataframes:
rownames(Correlations.df) <- Correlations.df$Scores
Correlations.df$Scores <- NULL
rownames(Pvalue.df) <- Pvalue.df$Scores
Pvalue.df$Scores <- NULL
Correlations.df <- as.matrix(Correlations.df)
Pvalue.df <- as.matrix(Pvalue.df)



##################
setwd("~/Documents/GitHub/CancerDormancy/QuiescenceScoreValidation/Figures/")
pdf("SerumStarvation_QS_ValidationHeatplot.pdf",height = 3,width = 16)
corrplot(Correlations.df,
         method = "color",
         insig = "label_sig",
         p.mat = Pvalue.df,
         cl.pos = "b",
         pch.cex = 1,
         is.corr = F,
         tl.col = "black",
         cl.length = 3,
         cl.cex = 0.8)
dev.off()





#####Also plot these as boxplots:
Pvalue.df[Pvalue.df > 0.05] <- NA
Correlations.df <- Correlations.df + Pvalue.df
Correlations.df <- Correlations.df - Pvalue.df
Correlations.df <- data.frame(Correlations.df)
Correlations.df$Pancancer <- NULL
Correlations.df <- data.frame(t(Correlations.df))
Scores <- colnames(Correlations.df)
mean <- NULL
for (i in Scores) {
  test <- Correlations.df[[i]]
  test <- na.omit(test)
  m <- mean(test)
  mean <- c(mean, m)
  
  
}
Order <- data.frame(Scores, mean)
Order <- Order[order(Order$mean),]
Order <- as.character(Order$Scores)
Correlations.df <- melt(Correlations.df)
colnames(Correlations.df) <- c("Score","SerumStarvation_markers")
Correlations.df$Score <- factor(Correlations.df$Score, levels = Order)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
myCols <- myPalette(6)
pdf("SerumStarvation_QS_ValidationBoxplot.pdf",height = 7, width = 7)
p <- ggplot(Correlations.df, aes(x=Score, y=SerumStarvation_markers, fill=Score)) +
  geom_boxplot() +
  scale_fill_manual(values = myCols) + theme_classic() +
  ggtitle("") +
  xlab("Quiescence Score") + ylab("Cancerwise correlation with SerumStarvation markers")
p + rotate_x_text(45) + theme(legend.position = "bottom")
dev.off()



