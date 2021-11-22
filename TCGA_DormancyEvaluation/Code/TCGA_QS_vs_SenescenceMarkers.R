#####################################################
#Quiescence score correlation with senescence markers
#####################################################


#Load required R packages:
library(biomaRt)
library(ggpubr)


#Load quiescence scores
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
load("TCGA_common_QS_purity_scaled.RData")

#Load expression data
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
load("combined_scaled_expr_data.RData") #This is an example dataset with only 100 entries 
all(rownames(combined.scaled) == rownames(z_score))

#Load SASP genes:
SASP <- read.csv("journal.pbio.3000599.s009.csv")
#This data can be downloaded from Basisty et al 2020 (PLOS Biology)
SASP <- as.character(SASP$Genes)
human = useMart("ensembl", dataset="hsapiens_gene_ensembl")
biomart_conversion <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=SASP, mart=human, filters = "hgnc_symbol")
SASP <- biomart_conversion$ensembl_gene_id
SASP <- SASP[SASP %in% colnames(combined.scaled)]
combined.scaled$SASP <- rowMeans(subset(combined.scaled, select = SASP), na.rm = TRUE)


################
#Plot correlations:
combined.scaled$QuiescenceScore <- z_score$z_score
setwd("~/Documents/GitHub/CancerDormancy/TCGA_DormancyEvaluation/Figures/")
pdf("TCGA_QS_vs_SASP.pdf", width = 5, height = 5)
p <- ggscatter(combined.scaled, x = "SASP", y = "QuiescenceScore",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE,
               conf.int = TRUE, alpha = 0.1, size = 2)
p + labs(x = "SASP", y = "Quiescence Score")
dev.off()


setwd("~/Documents/GitHub/CancerDormancy/TCGA_DormancyEvaluation/Figures/")
pdf("TCGA_QS_vs_GLB1.pdf", width = 5, height = 5)
p <- ggscatter(combined.scaled, x = "ENSG00000170266", y = "QuiescenceScore",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE,
               conf.int = TRUE, alpha = 0.1, size = 2)
p + labs(x = "GLB1", y = "Quiescence Score")
dev.off()

