################################################
#MCF7 PTEN mutation vs doubling time
################################################

#Load required R packages:
library(GSVA)
library(pheatmap)


###Load RNA-seq data:
rna.seq <- read.csv("Supplementary_Table_12_MCF7_Gene_expression_profiles_BING.csv", header = TRUE)
#Data was obtained from Ben-David et al 2018 (Nature)
rownames(rna.seq) <- rna.seq$Gene_Symbol
rna.seq$Gene_Symbol <- NULL
rna.seq$Probeset_ID <- NULL


#Load quiescence biomarker genes
setwd("~/Documents/GitHub/CancerDormancy/Data/DormancyGeneList/")
load("downregulated_common.RData")
load("upregulated_common.RData")



#Select only genes involved in quiescence to reduce the size of the dataframe:
rna.seq <- data.frame(t(rna.seq))
rna.seq <- rna.seq[,which(colnames(rna.seq) %in% c(downregulated_common, upregulated_common))]
#The data is already log transformed and so is ready for z_score calculation:
rna.seq <- as.matrix(t(rna.seq))
#Now work out the combined z_score for the samples
gene_lists <- list(upregulated_common, downregulated_common)
z_score <- gsva(rna.seq, gene_lists, method = "zscore", verbose = FALSE, parallel.sz=1)
z_score <- t(z_score)
z_score <- data.frame(z_score)
z_score$z_score <- z_score$X1 - z_score$X2
z_score$X1 <- NULL
z_score$X2 <- NULL


############
#Add information about the doubling times:
low <- c("MCF7.U","MCF7.Z","MCF7.T","MCF7.X","MCF7.Y","MCF7.W","MCF7.Q","MCF7.V","MCF7.A")
medium <- c("MCF7.F","MCF7.S","MCF7.R","MCF7.H","MCF7.N","MCF7.B","MCF7.D","MCF7.M","MCF7.E")
high <- c("MCF7.K","MCF7.I","MCF7.J","MCF7.AA","MCF7.L","MCF7.G","MCF7.P","MCF7.O","MCF7.C")
z_score$DoublingTime <- sapply(rownames(z_score), function(x)
  ifelse(x %in% low, "Low",
         ifelse(x %in% medium,"Medium",
                ifelse(x %in% high, "High","Other"))))

table(z_score$DoublingTime)



###Load information about mutations:
CellLines <- rownames(z_score)
mutation.data <- read.csv("Supplementary_Table_5_MCF7_SNV-Indel-Phased.csv", header = TRUE)
#Data was obtained from Ben-David et al 2018 (Nature)
mutation.data$tumor_sample_name <- as.character(mutation.data$tumor_sample_name)
mutation.data$tumor_sample_name <- sapply(mutation.data$tumor_sample_name, function(x)
  paste(strsplit(x,"-")[[1]][1:2],collapse="-"))
mutation.data$tumor_sample_name <- gsub('\\-', '.', mutation.data$tumor_sample_name)
#Remove snvs that are likey to not be deletorious
table(mutation.data$Canonical_Variant_Classification)
mutation.data <- mutation.data[mutation.data$Canonical_Variant_Classification %in% c("Frameshift","Inframe_Del","Initiator_Codon","Missense","Nonsense"),]


##Create a results table 
results <- data.frame(CellLines)
##LOF genes
LOF <- c("TP53","PTEN")
for (i in LOF) {
  
  print(i)
  selected_snv <- mutation.data[mutation.data$Canonical_Hugo_Symbol %in% i,]
  selected.snv <- unique(as.character(selected_snv$tumor_sample_name))
  
  results[[paste(i," MUT",sep = "")]] <- sapply(results$CellLines, function(x)
    ifelse(x %in% c(selected.snv),1,0))
  
  
  
}



##Order cell libes according to their quiescence score:
z_score <- z_score[order(z_score$z_score),]
z_score$QuiescenceScore <- z_score$z_score
z_score$z_score <- NULL
sample_order <- rownames(z_score)
rownames(results) <- results$CellLines
results <- results[order(match(rownames(results), sample_order)), , drop = FALSE]
results$DoublingTime <- z_score$DoublingTime
z_score$DoublingTime <- NULL
results$CellLines <- NULL
z_score$DoublingTime <- NULL
results$`PTEN MUT` <- as.character(results$`PTEN MUT`)

colfunc <- colorRampPalette(c("#00441B","#F7F7F7","#40004B"))
colfunc(27)
colours <- colfunc(27)
#Plot:
p<-pheatmap(t(z_score), show_colnames = FALSE, show_rownames = TRUE, cluster_rows = FALSE, cluster_cols = FALSE, annotation = results, legend = TRUE, color = colours)
setwd("~/Documents/GitHub/CancerDormancy/PanCancer_ElasticNetRegression/FeatureExploration/Figures/")
pdf("MCF7_heatmap.pdf",height = 4,width = 8)
p 
dev.off()