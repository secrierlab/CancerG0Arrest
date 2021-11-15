#########################################################
####GSE137912 treatment response analysis
#########################################################

#Load required packages:
library(GSVA)
library(umap)
library(ggplot2)
library(ggpubr)


##Load expression data
X <- read.csv("GSE137912_logcounts.csv",header = TRUE) #Data can be obtained from GEO using accession code GSE137912
X <- as.matrix(normalised_imputed_data)
all_genes <- as.character(X$X)
X$X <- NULL
rownames(X) <- all_genes
##Load annotation data
annotation <- read.csv("GSE137912_cell.annotation.csv",header = TRUE)#Data can be obtained from GEO using accession code GSE137912
annotation_H358 <- annotation[annotation$Line %in% "H358",]
##Select an expression matrix with the desired samples
samples_H358 <- as.character(annotation_H358$X)
H358_expression <- X[,colnames(X) %in% samples_H358]

##Load genes differentially expressed in all 5 forms of quiescence 
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("upregulated_common.RData")
load("downregulated_common.RData")
quiescence_genes <- c(downregulated_common, upregulated_common)


##Calculating quiescence scores 
gene_lists <- list(upregulated_common, downregulated_common)
H358_expression <- as.matrix(H358_expression)
es.dif <- gsva(H358_expression, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
es.dif <- t(es.dif)
es.dif <- data.frame(es.dif)
es.dif$GSVA_score <- es.dif$X1 - es.dif$X2
all(rownames(es.dif) == annotation_H358$X)
es.dif$Hours <- annotation_H358$Hours
es.dif$Common_score <- es.dif$GSVA_score
es.dif$GSVA_score <- NULL
z_scores <- es.dif
QS <- z_scores
QS$X1 <- NULL
QS$X2 <- NULL
#Save quiescence score annotation:
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Results/")
save(z_scores, file = "GSE137912_QS_H348.RData")



##########
#Summarise the data:
QS$CellStatus <- sapply(QS$Common_score, function(x)
  ifelse(x < 0, "Proliferating","Quiescent"))
table(QS$Hours)
Hours <- c(0,0,4,4,24,24,72,72)
CellStatus <- c("Proliferating","Quiescent","Proliferating","Quiescent","Proliferating","Quiescent","Proliferating","Quiescent")
N <- NULL
for (i in c(0,4,24,72)) {
  
  print(i)
  test <- QS[QS$Hours %in% i,]
  test <- table(test$CellStatus)
  n <- test[1]
  N <- c(N,n)
  n <- test[2]
  N <- c(N,n)
  
}
Summary <- data.frame(Hours, CellStatus, N)
Summary$Hours <- as.factor(Summary$Hours)
setwd("~/Documents/GitHub/CancerDormancy/Prognosis_and_treatment_response_analysis/Figures/")
pdf("GSE137912_barplot_cell_composition.pdf", height = 5, width = 4)
p <- ggplot(Summary, aes(fill=CellStatus, y=N, x=Hours, width = 0.75)) + 
  geom_bar(stat = "identity", position = "fill") + theme_classic()
p + rotate_x_text(45) + scale_fill_manual(values = c("#666666", "#D95F02"))
dev.off()




#####################################################################
####run umap with 0,4,24,72 hr time points
timepoint0 <- QS[QS$Hours %in% 0,]
timepoint0 <- as.character(rownames(timepoint0))
H358_expression <- as.data.frame(H358_expression)
umap.expr <- H358_expression[,colnames(H358_expression) %in% timepoint0]
umap.expr <- umap.expr[rownames(umap.expr) %in% c(upregulated_common,downregulated_common),]
umap.expr <- as.matrix(t(umap.expr))

#Running UMAP - hr 0:
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
##This is the matrix with cooridnates
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_0 <- UMAP_coordinates
UMAP_coordinates_0$Hour <- "0hr"


###############################
timepoint4 <- QS[QS$Hours %in% 4,]
timepoint4 <- as.character(rownames(timepoint4))
umap.expr <- H358_expression[,colnames(H358_expression) %in% timepoint4]
umap.expr <- umap.expr[rownames(umap.expr) %in% c(upregulated_common,downregulated_common),]
umap.expr <- as.matrix(t(umap.expr))

#Running UMAP - 4hr:
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
##This is the matrix with cooridnates
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_4 <- UMAP_coordinates
UMAP_coordinates_4$Hour <- "4hr"









###############################
timepoint24 <- QS[QS$Hours %in% 24,]
timepoint24 <- as.character(rownames(timepoint24))
umap.expr <- H358_expression[,colnames(H358_expression) %in% timepoint24]
umap.expr <- umap.expr[rownames(umap.expr) %in% c(upregulated_common,downregulated_common),]
umap.expr <- as.matrix(t(umap.expr))
#Running UMAP - 24hr:
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
##This is the matrix with cooridnates
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_24 <- UMAP_coordinates
UMAP_coordinates_24$Hour <- "24hr"





###############################
timepoint72 <- QS[QS$Hours %in% 72,]
timepoint72 <- as.character(rownames(timepoint72))
umap.expr <- H358_expression[,colnames(H358_expression) %in% timepoint72]
umap.expr <- umap.expr[rownames(umap.expr) %in% c(upregulated_common,downregulated_common),]
umap.expr <- as.matrix(t(umap.expr))
#Running UMAP - 72hr:
common.umap = umap(umap.expr, random_state=123)
common.umap$layout 
##This is the matrix with cooridnates
UMAP_coordinates <- data.frame(common.umap$layout)
colnames(UMAP_coordinates) <- c("UMAP1","UMAP2")
UMAP_coordinates$Sample <- rownames(UMAP_coordinates)
UMAP_coordinates_72 <- UMAP_coordinates
UMAP_coordinates_72$Hour <- "72hr"



###Plot the UMAP plots together:
UMAP_coordinates <- rbind(UMAP_coordinates_0, UMAP_coordinates_4, UMAP_coordinates_24, UMAP_coordinates_72)
QS$SampleID <- rownames(QS)
UMAP_coordinates <- merge(UMAP_coordinates, QS,
                          by.x = "Sample", by.y = "SampleID")
UMAP_coordinates$QuiescenceScore <- UMAP_coordinates$Common_score
pdf("GSE137912_UMAP.pdf", height = 6, width = 8)
ggplot(UMAP_coordinates, aes(x=UMAP1, y=UMAP2, colour = QuiescenceScore)) +
  geom_point() +
  scale_color_gradient2(low = "#59ac53", midpoint = 0,mid = "grey95", high = "#8b5aa8") + theme_classic() + facet_wrap(~Hour,dir = "v")
dev.off()



