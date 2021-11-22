###########################################################################
####Validation ballonplots (TCGA and cBioPortal studies)
###########################################################################

#Load required R packages:
library(TCGAbiolinks)
library(biomaRt)
library(reshape)
library(ggplot2)
library(ggpubr)

###########################################################################
####TCGA analysis
###########################################################################

########
#Load quiescence scores
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
load("TCGA_common_QS_purity_scaled.RData")
z_score$QuiescenceScore <- z_score$z_score
z_score$SampleID <- sapply(rownames(z_score), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))

######
#Load CNV data:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
load("combined_CNV.RData") #This is a toy example dataset with 100 entires only
CNV$SampleID <- sapply(rownames(CNV), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
all.cnv.samples <- unique(as.character(CNV$SampleID))
z_score <- z_score[z_score$SampleID %in% all.cnv.samples,]

##############
###Load SNV data:

####Download maf data:
solid_tumours <- unique(z_score$CancerType)
for (tumour in solid_tumours){
  print(paste0(tumour, ": This is ", which(tumour == solid_tumours), " out of ", length(solid_tumours)))    #Show progress
  df.maf <- GDCquery_Maf(tumour, pipelines = "mutect2")
  df.maf.short <- df.maf[,c("Tumor_Sample_Barcode","Variant_Classification","Hugo_Symbol")]
  df.maf.short$Cancer <- as.character(tumour)
  
  if (tumour == solid_tumours[1]){
    maf_all <- df.maf.short
  } else {
    maf_all <- rbind(maf_all, df.maf.short)
  }
}
pancancer_snv <- maf_all
pancancer_snv$SampleID <- sapply(pancancer_snv$Tumor_Sample_Barcode, function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
#Select snvs likely to be deleterious
pancancer_snv <- pancancer_snv[which(pancancer_snv$Variant_Classification
                                     %in% c("Frame_Shift_Del","Frame_Shift_Ins ","In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation")),]
all.snv.samples <- unique(as.character(pancancer_snv$SampleID))
z_score <- z_score[z_score$SampleID %in% all.snv.samples,]

##############################
##Create a data frame to summarise genes that need to be validated:
Gene <- c("TP53","MYC","RB1","CDKN2A","PTEN","RAD21","LMNA","CEP89","PTEN","LRP1B","ZMYM2","FHIT","PRRX1","ZFHX3","SDHA","FSTL3","ETV6","ATF1","PLAG1","CYLD","CCNB1IP1","SMAD4","SIRPA","RECQL4","ATF1","STAG1","SDHB","FUS","RB1","KLF6","GATA3")
Feature_type <- c("SNV","AMP","DEL","DEL","DEL","AMP","AMP","AMP","SNV","DEL","DEL","DEL","AMP","DEL","AMP","DEL","AMP","DEL","AMP","DEL","DEL","DEL","DEL","DEL","AMP","AMP","AMP","AMP","AMP","DEL","DEL")
Gene_annotation <- data.frame(Gene, Feature_type)
##Add ensg annotation:
human = useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "useast.ensembl.org")
results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=Gene, mart=human, filters = "hgnc_symbol")
results <- results[results$ensembl_gene_id %in% colnames(CNV),]

############################
#Set up results dataframes  
Difference_matrix <- data.frame(solid_tumours)
Pvalue_matrix <- data.frame(solid_tumours)

########################################
#######Fill information for SNV features
SNV_genes <- Gene_annotation[Gene_annotation$Feature_type %in% "SNV",]
SNV_genes <- as.character(SNV_genes$Gene)

for (a in SNV_genes) {
  
  print(a)

  ##Select samples with mutations
  snvs_selected <- pancancer_snv[pancancer_snv$Hugo_Symbol %in% a,]
  snvs_selected <- unique(as.character(snvs_selected$SampleID))
  z_score$Mutational.status <- sapply(z_score$SampleID, function(x)
    ifelse(x %in% snvs_selected, "MT","WT"))
  
  ####For each cancer type record change in mean between MT and WT, as well as significance
  mean_change <- NULL
  pvalue <- NULL
  for (c in solid_tumours) {
    
    test.data <- z_score[z_score$CancerType %in% c,]
    print(table(test.data$Mutational.status))
    WT <- test.data[test.data$Mutational.status %in% "WT",]
    WT <- mean(WT$QuiescenceScore)
    MT <- test.data[test.data$Mutational.status %in% "MT",]
    MT <- mean(MT$QuiescenceScore)
    change <- MT - WT
    #Record change in quiescence score
    mean_change <- c(mean_change, change)
    #record p value:
    test.data <- z_score[z_score$CancerType %in% c,]
    WT <- test.data[test.data$Mutational.status %in% "WT",]
    MT <- test.data[test.data$Mutational.status %in% "MT",]
    p <- 1
    
    tryCatch({
      res <- wilcox.test(WT$z_score, MT$z_score)
    }, error=function(e){})
    
    
    tryCatch({
      p <- res$p.value
    }, error=function(e){})
    
    
    pvalue <- c(pvalue,p)
    rm(res)
    rm(p)
    
    
  }
  
  Difference_matrix[[paste(a,"_SNV",sep = "")]] <- mean_change
  Pvalue_matrix[[paste(a,"_SNV",sep = "")]] <- pvalue
  
}


###################################################
#######Fill the results table with CNV-loss events
DEL_genes <- Gene_annotation[Gene_annotation$Feature_type %in% "DEL",]
DEL_genes <- as.character(DEL_genes$Gene)


for (a in DEL_genes) {
  
  print(a)
  ensg <- results[results$hgnc_symbol %in% a,]
  ensg <- ensg$ensembl_gene_id
  
  #####select samples with the alteration of interest:
  CNV_selected <- CNV[CNV[[ensg]] %in% -1,]
  CNV_selected <- as.character(CNV_selected$SampleID)
  z_score$Mutational.status <- sapply(z_score$SampleID, function(x)
    ifelse(x %in% CNV_selected, "MT","WT"))
  
  ####For each cancer type record change in mean between MT and WT, as well as significance
  
  mean_change <- NULL
  pvalue <- NULL
  
  for (c in solid_tumours) {
    
    test.data <- z_score[z_score$CancerType %in% c,]
    print(table(test.data$Mutational.status))
    WT <- test.data[test.data$Mutational.status %in% "WT",]
    WT <- mean(WT$QuiescenceScore)
    MT <- test.data[test.data$Mutational.status %in% "MT",]
    MT <- mean(MT$QuiescenceScore)
    change <- MT - WT
    #Record change in quiescence score
    mean_change <- c(mean_change, change)
    #record p value:
    test.data <- z_score[z_score$CancerType %in% c,]
    WT <- test.data[test.data$Mutational.status %in% "WT",]
    MT <- test.data[test.data$Mutational.status %in% "MT",]
    p <- 1
    
    tryCatch({
      res <- wilcox.test(WT$z_score, MT$z_score)
    }, error=function(e){})
    
    
    tryCatch({
      p <- res$p.value
    }, error=function(e){})
    
    
    pvalue <- c(pvalue,p)
    rm(res)
    rm(p)
    
    
  }
  
  Difference_matrix[[paste(a,"_DEL",sep = "")]] <- mean_change
  Pvalue_matrix[[paste(a,"_DEL",sep = "")]] <- pvalue
  
}

###################################################
#######Fill the results table with CNV-loss events
AMP_genes <- Gene_annotation[Gene_annotation$Feature_type %in% "AMP",]
AMP_genes <- as.character(AMP_genes$Gene)


for (a in AMP_genes) {
  
  print(a)
  ensg <- results[results$hgnc_symbol %in% a,]
  ensg <- ensg$ensembl_gene_id
  
  ##Select samples with the alteration of interest:
  CNV_selected <- CNV[CNV[[ensg]] %in% c(1,2,3,4,5),]
  CNV_selected <- as.character(CNV_selected$SampleID)
  
  
  #####Determine the mutational status of each sample
  z_score$Mutational.status <- sapply(z_score$SampleID, function(x)
    ifelse(x %in% CNV_selected, "MT","WT"))
  
  ####For each cancer type record change in mean between MT and WT, as well as significance
  
  mean_change <- NULL
  pvalue <- NULL
  
  for (c in solid_tumours) {
    
    test.data <- z_score[z_score$CancerType %in% c,]
    print(table(test.data$Mutational.status))
    WT <- test.data[test.data$Mutational.status %in% "WT",]
    WT <- mean(WT$QuiescenceScore)
    MT <- test.data[test.data$Mutational.status %in% "MT",]
    MT <- mean(MT$QuiescenceScore)
    change <- MT - WT
    #Record change in quiescence score
    mean_change <- c(mean_change, change)
    #record p value:
    test.data <- z_score[z_score$CancerType %in% c,]
    WT <- test.data[test.data$Mutational.status %in% "WT",]
    MT <- test.data[test.data$Mutational.status %in% "MT",]
    p <- 1
    
    tryCatch({
      res <- wilcox.test(WT$z_score, MT$z_score)
    }, error=function(e){})
    
    
    tryCatch({
      p <- res$p.value
    }, error=function(e){})
    
    
    pvalue <- c(pvalue,p)
    rm(res)
    rm(p)
    
    
  }
  
  Difference_matrix[[paste(a,"_AMP",sep = "")]] <- mean_change
  Pvalue_matrix[[paste(a,"_AMP",sep = "")]] <- pvalue
  
}


######save the final dataframe:
setwd("~/Documents/GitHub/CancerDormancy/PanCancer_ElasticNetRegression/ValidationBallonPlot/Results/")
save(Difference_matrix, file = "TCGA_difference_matrix.RData")
save(Pvalue_matrix, file = "TCGA_pvalue_marix.RData")












###########################################################################
####cBioPortal Validation Dataset analysis
###########################################################################


###List genes which should be validated:
Gene <- c("TP53","MYC","RB1","CDKN2A","PTEN","RAD21","LMNA","CEP89","PTEN","LRP1B","ZMYM2","FHIT","PRRX1","ZFHX3","SDHA","FSTL3","ETV6","ATF1","PLAG1","CYLD","CCNB1IP1","SMAD4","SIRPA","RECQL4","ATF1","STAG1","SDHB","FUS","RB1","KLF6","GATA3")
Feature_type <- c("SNV","AMP","DEL","DEL","DEL","AMP","AMP","AMP","SNV","DEL","DEL","DEL","AMP","DEL","AMP","DEL","AMP","DEL","AMP","DEL","DEL","DEL","DEL","DEL","AMP","AMP","AMP","AMP","AMP","DEL","DEL")
Gene_annotation <- data.frame(Gene, Feature_type)
Gene_annotation <- Gene_annotation[rev(order(Gene_annotation$Feature_type)),]
Features <- paste(Gene_annotation$Gene, Gene_annotation$Feature_type, sep = "_")
##Add ensg annotation:
human = useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "useast.ensembl.org")
results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=Gene, mart=human, filters = "hgnc_symbol")
##Two ensg numbers were reported for PTEN (use ENSG00000171862 not  ENSG00000284792)
results <- results[!(results$ensembl_gene_id %in% "ENSG00000284792"),]
#Set up results dataframes  
Difference_matrix <- data.frame(Features)
Pvalue_matrix <- data.frame(Features)


##############################################
####Dataset 1: BRCA METABRIC (both maf and snv)
#Load the quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/Data/cBioPortal_Studies/QuiescenceScores/")
load("metabrick_quiescence_score.RData")
z_score$sample <- rownames(z_score)
#Load maf file:
maf <- read.csv("data_mutations_mskcc.csv",header = TRUE) #Data can be obtained from cbiportal
maf$Variant_Classification <- as.character(maf$Variant_Classification)
maf <- maf[maf$Variant_Classification %in% c("Frame_Shift_Ins","In_Frame_Ins","Missense_Mutation","Frame_Shift_Del","Nonsense_Mutation","In_Frame_Del","Nonstop_Mutation"),]
##CNV data
CNV.data <- read.table("data_CNA.txt",sep = "\t",header = TRUE) #Data can be obtained from cbioportal
CNV.data <- data.frame(CNV.data)
rownames(CNV.data) <- CNV.data$Hugo_Symbol
CNV.data$Hugo_Symbol <- NULL
CNV.data$Entrez_Gene_Id <- NULL
CNV.data <- data.frame(t(CNV.data))
rownames(CNV.data) <- gsub("\\.", "-", rownames(CNV.data))
CNV.data$Sample <- rownames(CNV.data)
####Look at samples which have reported mutations
snv.samples <- unique(as.character(maf$Tumor_Sample_Barcode))
cnv.samples <- unique(as.character(CNV.data$Sample))
z_score <- z_score[z_score$sample %in% c(snv.samples, cnv.samples),]
############
#Check which genes can be checked:
test.snv <- Gene_annotation[Gene_annotation$Feature_type %in% "SNV",]
test.snv <- as.character(test.snv$Gene)
all.genes.snv <- unique(as.character(maf$Hugo_Symbol))
test.snv <- test.snv[test.snv %in% all.genes.snv]
test.amp  <- Gene_annotation[Gene_annotation$Feature_type %in% "AMP",]
test.amp <- as.character(test.amp$Gene)
all.genes.cnv <- unique(as.character(colnames(CNV.data)))
test.amp <- test.amp[test.amp %in% all.genes.cnv]
test.del  <- Gene_annotation[Gene_annotation$Feature_type %in% "DEL",]
test.del <- as.character(test.del$Gene)
all.genes.del <- unique(as.character(colnames(CNV.data)))
test.del <- test.del[test.del %in% all.genes.cnv]
Features <- c(paste(test.snv,"_SNV",sep = ""),paste(test.del,"_DEL",sep = ""),paste(test.amp,"_AMP",sep = ""))
test.df.change <- data.frame(Features)
test.df.pval <- data.frame(Features)

mean_change <- NULL
pvalue <- NULL
################SNVs:
for (a in test.snv) {
  
  print(a)
  selected_maf <- maf[maf$Hugo_Symbol %in% a,]
  selected_maf <- selected_maf[selected_maf$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation"),]
  mutated_samples <- as.character(selected_maf$Tumor_Sample_Barcode)
  mutated_samples <- unique(mutated_samples)
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% mutated_samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}
################DELs:
for (a in test.del) {
  
  print(a)
  cnv.samples <- CNV.data[CNV.data[[a]] %in% c(-1,-2,-3,-4,-5),]
  cnv.samples <- unique(as.character(cnv.samples$Sample))
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% cnv.samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}


################AMPs:
for (a in test.amp) {
  
  print(a)
  cnv.samples <- CNV.data[CNV.data[[a]] %in% c(1,2,3,4,5),]
  cnv.samples <- unique(as.character(cnv.samples$Sample))
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% cnv.samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}
test.df.change$BRCA_metabric <- mean_change
test.df.pval$BRCA_metabric <- pvalue
#Merge with the original dataframe:
Difference_matrix <- merge(Difference_matrix, test.df.change,
                           by.x = "Features", by.y = "Features",all = TRUE)
Pvalue_matrix <- merge(Pvalue_matrix, test.df.pval,
                       by.x = "Features", by.y = "Features",all = TRUE)





###################################################################
####Dataset 2: BRCA smc (only maf)
setwd("~/Documents/GitHub/CancerDormancy/Data/cBioPortal_Studies/QuiescenceScores/")
load("brca_smc_2018_quiescence_score.RData")
z_score$sample <- rownames(z_score)
#Load maf file:
maf <- read.table("data_mutations_mskcc.txt",sep = "\t", header = TRUE) #Data can be obtained from cbioportal
####Look at samples which have reported mutations
snv.samples <- unique(as.character(maf$Tumor_Sample_Barcode))
z_score <- z_score[z_score$sample %in% snv.samples,]
############
#Check which genes can be checked:
test.snv <- Gene_annotation[Gene_annotation$Feature_type %in% "SNV",]
test.snv <- as.character(test.snv$Gene)
all.genes.snv <- unique(as.character(maf$Hugo_Symbol))
test.snv <- test.snv[test.snv %in% all.genes.snv]
Features <- c(paste(test.snv,"_SNV",sep = ""))
test.df.change <- data.frame(Features)
test.df.pval <- data.frame(Features)
mean_change <- NULL
pvalue <- NULL
################SNVs:
for (a in test.snv) {
  
  print(a)
  selected_maf <- maf[maf$Hugo_Symbol %in% a,]
  selected_maf <- selected_maf[selected_maf$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation"),]
  mutated_samples <- as.character(selected_maf$Tumor_Sample_Barcode)
  mutated_samples <- unique(mutated_samples)
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% mutated_samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}
test.df.change$BRCA_smc <- mean_change
test.df.pval$BRCA_smc <- pvalue
#Merge with the original dataframe:
Difference_matrix <- merge(Difference_matrix, test.df.change,
                           by.x = "Features", by.y = "Features",all = TRUE)
Pvalue_matrix <- merge(Pvalue_matrix, test.df.pval,
                       by.x = "Features", by.y = "Features",all = TRUE)






###################################################################
####Dataset 3: GIS021 (maf and cnv)
#Load the quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/Data/cBioPortal_Studies/QuiescenceScores/")
load("GIS031_quiescence_score.RData")
z_score$sample <- rownames(z_score)
#Load maf file:
maf <- read.table("snv_indel.maf",sep = "\t", header = TRUE) #Data can be obtained from cbioportal
snv.samples <- unique(as.character(maf$Tumor_Sample_Barcode))
maf <- maf[maf$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation"),]
#####Load CNV data:
CNV.data <- read.table("cnv.tsv",sep = "\t",header = TRUE) #Data can be obtained from cbioportal
CNV.data <- CNV.data[CNV.data$Hugo_Symbol %in% Gene,]
rownames(CNV.data) <- CNV.data$Hugo_Symbol
CNV.data$Hugo_Symbol <- NULL
CNV.data$Entrez_Gene_Id <- NULL
CNV.data <- data.frame(t(CNV.data))
rownames(CNV.data) <- gsub("\\.", "-", rownames(CNV.data))
CNV.data$Sample <- rownames(CNV.data)
cnv.samples <- unique(as.character(CNV.data$Sample))
z_score <- z_score[z_score$sample %in% c(snv.samples, cnv.samples),]
############
#Check which genes can be checked:
test.snv <- Gene_annotation[Gene_annotation$Feature_type %in% "SNV",]
test.snv <- as.character(test.snv$Gene)
all.genes.snv <- unique(as.character(maf$Hugo_Symbol))
test.snv <- test.snv[test.snv %in% all.genes.snv]
test.amp  <- Gene_annotation[Gene_annotation$Feature_type %in% "AMP",]
test.amp <- as.character(test.amp$Gene)
all.genes.cnv <- unique(as.character(colnames(CNV.data)))
test.amp <- test.amp[test.amp %in% all.genes.cnv]
test.del  <- Gene_annotation[Gene_annotation$Feature_type %in% "DEL",]
test.del <- as.character(test.del$Gene)
all.genes.del <- unique(as.character(colnames(CNV.data)))
test.del <- test.del[test.del %in% all.genes.cnv]
Features <- c(paste(test.del,"_DEL",sep = ""),paste(test.amp,"_AMP",sep = ""))
test.df.change <- data.frame(Features)
test.df.pval <- data.frame(Features)
mean_change <- NULL
pvalue <- NULL
################SNVs:
for (a in test.snv) {
  
  print(a)
  selected_maf <- maf[maf$Hugo_Symbol %in% a,]
  mutated_samples <- as.character(selected_maf$Tumor_Sample_Barcode)
  mutated_samples <- unique(mutated_samples)
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% mutated_samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}
################DELs:
for (a in test.del) {
  
  print(a)
  cnv.samples <- CNV.data[CNV.data[[a]] %in% c(-1,-2,-3,-4,-5),]
  cnv.samples <- unique(as.character(cnv.samples$Sample))
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% cnv.samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}
################AMPs:
for (a in test.amp) {
  
  print(a)
  cnv.samples <- CNV.data[CNV.data[[a]] %in% c(1,2,3,4,5),]
  cnv.samples <- unique(as.character(cnv.samples$Sample))
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% cnv.samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}
test.df.change$GIS031 <- mean_change
test.df.pval$GIS031 <- pvalue
#Merge with the original dataframe:
Difference_matrix <- merge(Difference_matrix, test.df.change,
                           by.x = "Features", by.y = "Features",all = TRUE)
Pvalue_matrix <- merge(Pvalue_matrix, test.df.pval,
                       by.x = "Features", by.y = "Features",all = TRUE)






###################################################################
####Dataset 4: mskcc_solit_2012 (maf only)
##Load the quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/Data/cBioPortal_Studies/QuiescenceScores/")
load("blca_mskcc_solit_2012_quiescence_score.RData")
z_score$sample <- rownames(z_score)
#Load maf file:
maf <- read.table("data_mutations_mskcc.txt",sep = "\t", header = TRUE) #Data can be obtained from cbioportal
####Look at samples which have reported mutations
snv.samples <- unique(as.character(maf$Tumor_Sample_Barcode))
z_score <- z_score[z_score$sample %in% snv.samples,]
#Check which genes can be checked:
test.snv <- Gene_annotation[Gene_annotation$Feature_type %in% "SNV",]
test.snv <- as.character(test.snv$Gene)
all.genes.snv <- unique(as.character(maf$Hugo_Symbol))
test.snv <- test.snv[test.snv %in% all.genes.snv]
Features <- c(paste(test.snv,"_SNV",sep = ""))
test.df.change <- data.frame(Features)
test.df.pval <- data.frame(Features)
mean_change <- NULL
pvalue <- NULL
################SNVs:
for (a in test.snv) {
  
  print(a)
  selected_maf <- maf[maf$Hugo_Symbol %in% a,]
  selected_maf <- selected_maf[selected_maf$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation"),]
  mutated_samples <- as.character(selected_maf$Tumor_Sample_Barcode)
  mutated_samples <- unique(mutated_samples)
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% mutated_samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}
test.df.change$mskcc_solit_2012 <- mean_change
test.df.pval$mskcc_solit_2012 <- pvalue
#Merge with the original dataframe:
Difference_matrix <- merge(Difference_matrix, test.df.change,
                           by.x = "Features", by.y = "Features",all = TRUE)
Pvalue_matrix <- merge(Pvalue_matrix, test.df.pval,
                       by.x = "Features", by.y = "Features",all = TRUE)






###################################################################
####Dataset 5: paad_qcmg_uq (maf only)
setwd("~/Documents/GitHub/CancerDormancy/Data/cBioPortal_Studies/QuiescenceScores/")
load("paad_qcmg_uq_2016_quiescence_score.RData")
z_score$sample <- rownames(z_score)
#Load maf file:
maf <- read.table("data_mutations_mskcc.txt",sep = "\t", header = TRUE) #Data can be obtained from cbioportal
####Look at samples which have reported mutations
snv.samples <- unique(as.character(maf$Tumor_Sample_Barcode))
z_score <- z_score[z_score$sample %in% snv.samples,]
#Check which genes can be checked:
test.snv <- Gene_annotation[Gene_annotation$Feature_type %in% "SNV",]
test.snv <- as.character(test.snv$Gene)
all.genes.snv <- unique(as.character(maf$Hugo_Symbol))
test.snv <- test.snv[test.snv %in% all.genes.snv]
Features <- c(paste(test.snv,"_SNV",sep = ""))
test.df.change <- data.frame(Features)
test.df.pval <- data.frame(Features)
mean_change <- NULL
pvalue <- NULL
################SNVs:
for (a in test.snv) {
  
  print(a)
  selected_maf <- maf[maf$Hugo_Symbol %in% a,]
  selected_maf <- selected_maf[selected_maf$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation"),]
  mutated_samples <- as.character(selected_maf$Tumor_Sample_Barcode)
  mutated_samples <- unique(mutated_samples)
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% mutated_samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}
test.df.change$paad_qcmg_uq <- mean_change
test.df.pval$paad_qcmg_uq <- pvalue
#Merge with the original dataframe:
Difference_matrix <- merge(Difference_matrix, test.df.change,
                           by.x = "Features", by.y = "Features",all = TRUE)
Pvalue_matrix <- merge(Pvalue_matrix, test.df.pval,
                       by.x = "Features", by.y = "Features",all = TRUE)






###################################################################
####Dataset 6: prad bread (maf and cnv)
#Load the quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/Data/cBioPortal_Studies/QuiescenceScores/")
load("prad_broad_quiescence_score.RData")
z_score$sample <- rownames(z_score)
#Load maf file
maf <- read.table("data_mutations_mskcc.txt",sep = "\t", header = TRUE) #Data can be obtained from cbioportal
#Load cnv data
CNV.data <- read.table("data_CNA.txt",sep = "\t",header = TRUE)
CNV.data$Hugo_Symbol <- as.character(CNV.data$Hugo_Symbol)
duplicated.genes <- CNV.data$Hugo_Symbol[duplicated(CNV.data$Hugo_Symbol)]
CNV.data <- CNV.data[!(CNV.data$Hugo_Symbol %in% duplicated.genes),]
rownames(CNV.data) <- CNV.data$Hugo_Symbol
CNV.data$Hugo_Symbol <- NULL
CNV.data$Entrez_Gene_Id <- NULL
CNV.data <- data.frame(t(CNV.data))
rownames(CNV.data) <- gsub("\\.", "-", rownames(CNV.data))
CNV.data$Sample <- rownames(CNV.data)
####Look at samples which have reported mutations
snv.samples <- unique(as.character(maf$Tumor_Sample_Barcode))
cnv.samples <- unique(as.character(CNV.data$Sample))
z_score <- z_score[z_score$sample %in% c(snv.samples, cnv.samples),]
############
#Check which genes can be checked:
test.snv <- Gene_annotation[Gene_annotation$Feature_type %in% "SNV",]
test.snv <- as.character(test.snv$Gene)
all.genes.snv <- unique(as.character(maf$Hugo_Symbol))
test.snv <- test.snv[test.snv %in% all.genes.snv]
test.amp  <- Gene_annotation[Gene_annotation$Feature_type %in% "AMP",]
test.amp <- as.character(test.amp$Gene)
all.genes.cnv <- unique(as.character(colnames(CNV.data)))
test.amp <- test.amp[test.amp %in% all.genes.cnv]
test.del  <- Gene_annotation[Gene_annotation$Feature_type %in% "DEL",]
test.del <- as.character(test.del$Gene)
all.genes.del <- unique(as.character(colnames(CNV.data)))
test.del <- test.del[test.del %in% all.genes.cnv]
Features <- c(paste(test.snv,"_SNV",sep = ""),paste(test.del,"_DEL",sep = ""),paste(test.amp,"_AMP",sep = ""))
test.df.change <- data.frame(Features)
test.df.pval <- data.frame(Features)
mean_change <- NULL
pvalue <- NULL
################SNVs:
for (a in test.snv) {
  
  print(a)
  
  selected_maf <- maf[maf$Hugo_Symbol %in% a,]
  selected_maf <- selected_maf[selected_maf$Variant_Classification %in% c("Frame_Shift_Ins","In_Frame_Ins","Missense_Mutation","Frame_Shift_Del","Nonsense_Mutation","In_Frame_Del","Nonstop_Mutation"),]
  mutated_samples <- as.character(selected_maf$Tumor_Sample_Barcode)
  mutated_samples <- unique(mutated_samples)
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% mutated_samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}
################DELs:
for (a in test.del) {
  
  print(a)
  cnv.samples <- CNV.data[CNV.data[[a]] %in% c(-1,-2,-3,-4,-5),]
  cnv.samples <- unique(as.character(cnv.samples$Sample))
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% cnv.samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}
################AMPs:
for (a in test.amp) {
  
  print(a)
  cnv.samples <- CNV.data[CNV.data[[a]] %in% c(1,2,3,4,5),]
  cnv.samples <- unique(as.character(cnv.samples$Sample))
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% cnv.samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}
test.df.change$prad_broad <- mean_change
test.df.pval$prad_broad <- pvalue
#Merge with the original dataframe:
Difference_matrix <- merge(Difference_matrix, test.df.change,
                           by.x = "Features", by.y = "Features",all = TRUE)
Pvalue_matrix <- merge(Pvalue_matrix, test.df.pval,
                       by.x = "Features", by.y = "Features",all = TRUE)









###################################################################
####Dataset 7: prad mskcc (maf and cnv)
#Load the quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/Data/cBioPortal_Studies/QuiescenceScores/")
load("prad_mskcc_quiescence_score.RData")
z_score$sample <- rownames(z_score)
#Load maf file:
maf <- read.table("data_mutations_mskcc.txt",sep = "\t", header = TRUE) #Data can be obtained from cbioportal
#Load CNV data
CNV.data <- read.table("data_CNA.txt",sep = "\t",header = TRUE) #Data can be obtained from cbioportal
CNV.data$Hugo_Symbol <- as.character(CNV.data$Hugo_Symbol)
duplicated.genes <- CNV.data$Hugo_Symbol[duplicated(CNV.data$Hugo_Symbol)]
CNV.data <- CNV.data[!(CNV.data$Hugo_Symbol %in% duplicated.genes),]
rownames(CNV.data) <- CNV.data$Hugo_Symbol
CNV.data$Hugo_Symbol <- NULL
CNV.data$Entrez_Gene_Id <- NULL
CNV.data <- data.frame(t(CNV.data))
rownames(CNV.data) <- gsub("\\.", "-", rownames(CNV.data))
CNV.data$Sample <- rownames(CNV.data)
####Look at samples which have reported mutations
snv.samples <- unique(as.character(maf$Tumor_Sample_Barcode))
cnv.samples <- unique(as.character(CNV.data$Sample))
z_score <- z_score[z_score$sample %in% c(snv.samples, cnv.samples),]
############
#Check which genes can be checked:
test.snv <- Gene_annotation[Gene_annotation$Feature_type %in% "SNV",]
test.snv <- as.character(test.snv$Gene)
all.genes.snv <- unique(as.character(maf$Hugo_Symbol))
test.snv <- test.snv[test.snv %in% all.genes.snv]
test.amp  <- Gene_annotation[Gene_annotation$Feature_type %in% "AMP",]
test.amp <- as.character(test.amp$Gene)
all.genes.cnv <- unique(as.character(colnames(CNV.data)))
test.amp <- test.amp[test.amp %in% all.genes.cnv]
test.del  <- Gene_annotation[Gene_annotation$Feature_type %in% "DEL",]
test.del <- as.character(test.del$Gene)
all.genes.del <- unique(as.character(colnames(CNV.data)))
test.del <- test.del[test.del %in% all.genes.cnv]
Features <- c(paste(test.snv,"_SNV",sep = ""),paste(test.del,"_DEL",sep = ""),paste(test.amp,"_AMP",sep = ""))
test.df.change <- data.frame(Features)
test.df.pval <- data.frame(Features)
mean_change <- NULL
pvalue <- NULL
################SNVs:
for (a in test.snv) {
  
  print(a)
  
  selected_maf <- maf[maf$Hugo_Symbol %in% a,]
  selected_maf <- selected_maf[selected_maf$Variant_Classification %in% c("Frame_Shift_Ins","In_Frame_Ins","Missense_Mutation","Frame_Shift_Del","Nonsense_Mutation","In_Frame_Del","Nonstop_Mutation"),]
  mutated_samples <- as.character(selected_maf$Tumor_Sample_Barcode)
  mutated_samples <- unique(mutated_samples)
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% mutated_samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}
################DELs:
for (a in test.del) {
  
  print(a)
  cnv.samples <- CNV.data[CNV.data[[a]] %in% c(-1,-2,-3,-4,-5),]
  cnv.samples <- unique(as.character(cnv.samples$Sample))
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% cnv.samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}
################AMPs:
for (a in test.amp) {
  
  print(a)
  cnv.samples <- CNV.data[CNV.data[[a]] %in% c(1,2,3,4,5),]
  cnv.samples <- unique(as.character(cnv.samples$Sample))
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% cnv.samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}
test.df.change$prad_mskcc <- mean_change
test.df.pval$prad_mskcc <- pvalue
#Merge with the original dataframe:
Difference_matrix <- merge(Difference_matrix, test.df.change,
                           by.x = "Features", by.y = "Features",all = TRUE)
Pvalue_matrix <- merge(Pvalue_matrix, test.df.pval,
                       by.x = "Features", by.y = "Features",all = TRUE)






###################################################################
####Dataset 8: prostate dkfz (maf only)
##Load the quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/Data/cBioPortal_Studies/QuiescenceScores/")
load("prostate_dkfz_2018_quiescence_score.RData")
z_score$sample <- rownames(z_score)
#Load maf file:
maf <- read.table("data_mutations_mskcc.txt",sep = "\t", header = TRUE) #Data can be obtained from cbioportal
####Look at samples which have reported mutations
snv.samples <- unique(as.character(maf$Tumor_Sample_Barcode))
z_score <- z_score[z_score$sample %in% snv.samples,]
#Check which genes can be checked:
test.snv <- Gene_annotation[Gene_annotation$Feature_type %in% "SNV",]
test.snv <- as.character(test.snv$Gene)
all.genes.snv <- unique(as.character(maf$Hugo_Symbol))
test.snv <- test.snv[test.snv %in% all.genes.snv]
Features <- c(paste(test.snv,"_SNV",sep = ""))
test.df.change <- data.frame(Features)
test.df.pval <- data.frame(Features)
mean_change <- NULL
pvalue <- NULL
################SNVs:
for (a in test.snv) {
  
  print(a)
  selected_maf <- maf[maf$Hugo_Symbol %in% a,]
  selected_maf <- selected_maf[selected_maf$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation"),]
  mutated_samples <- as.character(selected_maf$Tumor_Sample_Barcode)
  mutated_samples <- unique(mutated_samples)
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% mutated_samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}
test.df.change$prostate_dkfz_2018 <- mean_change
test.df.pval$prostate_dkfz_2018 <- pvalue
#Merge with the original dataframe:
Difference_matrix <- merge(Difference_matrix, test.df.change,
                           by.x = "Features", by.y = "Features",all = TRUE)
Pvalue_matrix <- merge(Pvalue_matrix, test.df.pval,
                       by.x = "Features", by.y = "Features",all = TRUE)










###################################################################
####Dataset 9: sarc mskcc (maf and cnv)
setwd("~/Documents/GitHub/CancerDormancy/Data/cBioPortal_Studies/QuiescenceScores/")
load("sarc_mskcc_quiescence_score.RData")
z_score$sample <- rownames(z_score)
#Load maf file:
maf <- read.table("data_mutations_mskcc.txt",sep = "\t", header = TRUE) #Data can be obtained from cbioportal
##CNV data
CNV.data <- read.table("data_CNA.txt",sep = "\t",header = TRUE) #Data can be obtained from cbioportal
rownames(CNV.data) <- CNV.data$Entrez_Gene_Id
CNV.data$Entrez_Gene_Id <- NULL
CNV.data <- data.frame(t(CNV.data))
rownames(CNV.data) <- gsub("\\.", "-", rownames(CNV.data))
CNV.data$Sample <- rownames(CNV.data)
#Change names of genes from entez to hugo
CNV.data$TP53 <- CNV.data$X7157
CNV.data$RB1 <- CNV.data$X5925
CNV.data$CDKN2A <- CNV.data$X1029
CNV.data$PTEN <- CNV.data$X5728
CNV.data$LRP1B <- CNV.data$X53353
CNV.data$ZMYM2 <- CNV.data$X7750
CNV.data$FHIT <- CNV.data$X2272
CNV.data$ZFHX3 <- CNV.data$X463
CNV.data$FSTL3 <- CNV.data$X10272
CNV.data$ATF1 <- CNV.data$X466
CNV.data$CYLD <- CNV.data$X1540
CNV.data$CCNB1IP1 <- CNV.data$X57820
CNV.data$SMAD4 <- CNV.data$X4089
CNV.data$SIRPA <- CNV.data$X140885
CNV.data$RECQL4 <- CNV.data$X9401
CNV.data$BCL3 <- CNV.data$X602
CNV.data$NCOR1 <- CNV.data$X9611
CNV.data$KLF6 <- CNV.data$X1316
CNV.data$MYC <- CNV.data$X4609
CNV.data$RAD21 <- CNV.data$X5885
CNV.data$LMNA <- CNV.data$X4000
CNV.data$CEP89 <- CNV.data$X84902
CNV.data$PRRX1 <- CNV.data$X5396
CNV.data$ETV6 <- CNV.data$X2120
CNV.data$PLAG1 <- CNV.data$X5324
CNV.data$ATF1 <- CNV.data$X466
CNV.data$SDHB <- CNV.data$X6390
####Look at samples which have reported mutations
snv.samples <- unique(as.character(maf$Tumor_Sample_Barcode))
cnv.samples <- unique(as.character(CNV.data$Sample))
z_score <- z_score[z_score$sample %in% c(snv.samples, cnv.samples),]
############
#Check which genes can be checked:
test.snv <- Gene_annotation[Gene_annotation$Feature_type %in% "SNV",]
test.snv <- as.character(test.snv$Gene)
all.genes.snv <- unique(as.character(maf$Hugo_Symbol))
test.snv <- test.snv[test.snv %in% all.genes.snv]
test.amp  <- Gene_annotation[Gene_annotation$Feature_type %in% "AMP",]
test.amp <- as.character(test.amp$Gene)
all.genes.cnv <- unique(as.character(colnames(CNV.data)))
test.amp <- test.amp[test.amp %in% all.genes.cnv]
test.del  <- Gene_annotation[Gene_annotation$Feature_type %in% "DEL",]
test.del <- as.character(test.del$Gene)
all.genes.del <- unique(as.character(colnames(CNV.data)))
test.del <- test.del[test.del %in% all.genes.cnv]
Features <- c(paste(test.snv,"_SNV",sep = ""),paste(test.del,"_DEL",sep = ""),paste(test.amp,"_AMP",sep = ""))
test.df.change <- data.frame(Features)
test.df.pval <- data.frame(Features)
mean_change <- NULL
pvalue <- NULL
################SNVs:
for (a in test.snv) {
  
  print(a)
  
  selected_maf <- maf[maf$Hugo_Symbol %in% a,]
  selected_maf <- selected_maf[selected_maf$Variant_Classification %in% c("Frame_Shift_Ins","In_Frame_Ins","Missense_Mutation","Frame_Shift_Del","Nonsense_Mutation","In_Frame_Del","Nonstop_Mutation"),]
  mutated_samples <- as.character(selected_maf$Tumor_Sample_Barcode)
  mutated_samples <- unique(mutated_samples)
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% mutated_samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}
################DELs:
for (a in test.del) {
  
  print(a)
  cnv.samples <- CNV.data[CNV.data[[a]] %in% c(-1,-2,-3,-4,-5),]
  cnv.samples <- unique(as.character(cnv.samples$Sample))
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% cnv.samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}
################AMPs:
for (a in test.amp) {
  
  print(a)
  cnv.samples <- CNV.data[CNV.data[[a]] %in% c(1,2,3,4,5),]
  cnv.samples <- unique(as.character(cnv.samples$Sample))
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% cnv.samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}
test.df.change$sarc_mskcc <- mean_change
test.df.pval$sarc_mskcc <- pvalue
#Merge with the original dataframe:
Difference_matrix <- merge(Difference_matrix, test.df.change,
                           by.x = "Features", by.y = "Features",all = TRUE)
Pvalue_matrix <- merge(Pvalue_matrix, test.df.pval,
                       by.x = "Features", by.y = "Features",all = TRUE)




###################################################################
####Dataset 10: Wt target 2018
#Load the quiescence scores:
setwd("~/Documents/GitHub/CancerDormancy/Data/cBioPortal_Studies/QuiescenceScores/")
load("wt_target_2018_pub_quiescence_score.RData")
z_score$sample <- rownames(z_score)
#Load maf file:
maf <- read.table("data_mutations_mskcc.txt",sep = "\t", header = TRUE) #Data can be downloaded from cbioportal
#Load cnv data
CNV.data <- read.table("data_CNA.txt",sep = "\t",header = TRUE) #Data can be downloaded from cbioportal
rownames(CNV.data) <- CNV.data$Hugo_Symbol
CNV.data$Hugo_Symbol <- NULL
CNV.data$Entrez_Gene_Id <- NULL
CNV.data <- data.frame(t(CNV.data))
rownames(CNV.data) <- gsub("\\.", "-", rownames(CNV.data))
CNV.data$Sample <- rownames(CNV.data)
####Look at samples which have reported mutations
snv.samples <- unique(as.character(maf$Tumor_Sample_Barcode))
cnv.samples <- unique(as.character(CNV.data$Sample))
z_score <- z_score[z_score$sample %in% c(snv.samples, cnv.samples),]
############
#Check which genes can be checked:
test.snv <- Gene_annotation[Gene_annotation$Feature_type %in% "SNV",]
test.snv <- as.character(test.snv$Gene)
all.genes.snv <- unique(as.character(maf$Hugo_Symbol))
test.snv <- test.snv[test.snv %in% all.genes.snv]
test.amp  <- Gene_annotation[Gene_annotation$Feature_type %in% "AMP",]
test.amp <- as.character(test.amp$Gene)
all.genes.cnv <- unique(as.character(colnames(CNV.data)))
test.amp <- test.amp[test.amp %in% all.genes.cnv]
test.del  <- Gene_annotation[Gene_annotation$Feature_type %in% "DEL",]
test.del <- as.character(test.del$Gene)
all.genes.del <- unique(as.character(colnames(CNV.data)))
test.del <- test.del[test.del %in% all.genes.cnv]
Features <- c(paste(test.snv,"_SNV",sep = ""),paste(test.del,"_DEL",sep = ""),paste(test.amp,"_AMP",sep = ""))
test.df.change <- data.frame(Features)
test.df.pval <- data.frame(Features)
mean_change <- NULL
pvalue <- NULL
################SNVs:
for (a in test.snv) {
  
  print(a)
  
  selected_maf <- maf[maf$Hugo_Symbol %in% a,]
  selected_maf <- selected_maf[selected_maf$Variant_Classification %in% c("Frame_Shift_Ins","In_Frame_Ins","Missense_Mutation","Frame_Shift_Del","Nonsense_Mutation","In_Frame_Del","Nonstop_Mutation"),]
  mutated_samples <- as.character(selected_maf$Tumor_Sample_Barcode)
  mutated_samples <- unique(mutated_samples)
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% mutated_samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}
################DELs:
for (a in test.del) {
  
  print(a)
  cnv.samples <- CNV.data[CNV.data[[a]] %in% c(-1,-2,-3,-4,-5),]
  cnv.samples <- unique(as.character(cnv.samples$Sample))
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% cnv.samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}
################AMPs:
for (a in test.amp) {
  
  print(a)
  cnv.samples <- CNV.data[CNV.data[[a]] %in% c(1,2,3,4,5),]
  cnv.samples <- unique(as.character(cnv.samples$Sample))
  z_score$Mutation.status <- sapply(z_score$sample, function(x)
    ifelse(x %in% cnv.samples,"MT","WT"))
  
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  WT <- mean(WT$z_score)
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  MT <- mean(MT$z_score)
  change <- MT - WT
  #Record change in quiescence score
  mean_change <- c(mean_change, change)
  #record p value:
  WT <- z_score[z_score$Mutation.status %in% "WT",]
  MT <- z_score[z_score$Mutation.status %in% "MT",]
  p <- 1
  tryCatch({
    res <- wilcox.test(WT$z_score, MT$z_score)
  }, error=function(e){})
  tryCatch({
    p <- res$p.value
  }, error=function(e){})
  pvalue <- c(pvalue,p)
  rm(res)
  rm(p)
}
test.df.change$wt_target_2018 <- mean_change
test.df.pval$wt_target_2018 <- pvalue
#Merge with the original dataframe:
Difference_matrix <- merge(Difference_matrix, test.df.change,
                           by.x = "Features", by.y = "Features",all = TRUE)
Pvalue_matrix <- merge(Pvalue_matrix, test.df.pval,
                       by.x = "Features", by.y = "Features",all = TRUE)

#############
##Save the annotation
setwd("~/Documents/GitHub/CancerDormancy/PanCancer_ElasticNetRegression/ValidationBallonPlot/Results/")
save(Difference_matrix, file = "cBioPortal_difference_matrix.RData")
save(Pvalue_matrix, file = "cBioPortal_pvalue_matrix.RData")






###########################################################################
#### Ballonplot visualisation
###########################################################################


############
##Load TCGA summary:
setwd("~/Documents/GitHub/CancerDormancy/PanCancer_ElasticNetRegression/ValidationBallonPlot/Results/")
load("TCGA_difference_matrix.RData")
load("TCGA_pvalue_marix.RData")
Difference_matrix_TCGA <- Difference_matrix
Pvalue_matrix_TCGA <- Pvalue_matrix

#########################
#Load external validation:
load("cBioPortal_difference_matrix.RData")
load("cBioPortal_pvalue_matrix.RData")
Difference_matrix_val <- Difference_matrix
Pvalue_matrix_val <- Pvalue_matrix
rm(Difference_matrix)
rm(Pvalue_matrix)

##############
#Combine the two:
rownames(Difference_matrix_TCGA) <- Difference_matrix_TCGA$solid_tumours
Difference_matrix_TCGA$solid_tumours <- NULL
rownames(Pvalue_matrix_TCGA) <- Pvalue_matrix_TCGA$solid_tumours
Pvalue_matrix_TCGA$solid_tumours <- NULL
rownames(Difference_matrix_val) <- Difference_matrix_val$Features
Difference_matrix_val$Features <- NULL
Difference_matrix_val <- data.frame(t(Difference_matrix_val))
rownames(Pvalue_matrix_val) <- Pvalue_matrix_val$Features
Pvalue_matrix_val$Features <- NULL
Pvalue_matrix_val <- data.frame(t(Pvalue_matrix_val))
#Make sure features are in the same order
Difference_matrix_TCGA <- Difference_matrix_TCGA[,order(colnames(Difference_matrix_TCGA))]
Difference_matrix_val <- Difference_matrix_val[,order(colnames(Difference_matrix_val))]
Pvalue_matrix_TCGA <- Pvalue_matrix_TCGA[,order(colnames(Pvalue_matrix_TCGA))]
Pvalue_matrix_val <- Pvalue_matrix_val[,order(colnames(Pvalue_matrix_val))]
#Merge datafranes
Difference_matrix <- rbind(Difference_matrix_TCGA, Difference_matrix_val)
Pvalue_matrix <- rbind(Pvalue_matrix_TCGA, Pvalue_matrix_val)
Difference_matrix$CT <- rownames(Difference_matrix)
Pvalue_matrix$CT <- rownames(Pvalue_matrix)
#Remove nas
Difference_matrix[is.na(Difference_matrix)] <- 0

##########
#Melt dataframes:
Difference_matrix <- melt(Difference_matrix, id.vars = "CT")
Pvalue_matrix <- melt(Pvalue_matrix, id.vars = "CT")
colnames(Difference_matrix) <- c("CancerType", "Gene", "QuiescenceScoreChange")
Difference_matrix$pvalue  <- Pvalue_matrix$value
Difference_matrix$pvalue <- as.numeric(as.character(Difference_matrix$pvalue))
Difference_matrix$pvalue[Difference_matrix$pvalue > 0.05] <- NA
Difference_matrix$minuslog10pval <- -log10(Difference_matrix$pvalue)

###Exlude the folloing cancer types (no significant hits) and order alphabetically
exclude <- c("THYM","prostate_dkfz_2018","prad_broad","paad_qcmg_uq","PAAD","KICH","ESCA","CHOL")
CT <- unique(Difference_matrix$CancerType)
CT <- CT[!(CT %in% exclude)]
Difference_matrix <- Difference_matrix[Difference_matrix$CancerType %in% CT,]
Levels <- c("ACC","BLCA","BRCA","CESC","COAD","GBM","HNSC","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","UCEC","UCS","UVM","BRCA_metabric","BRCA_smc","GIS031","mskcc_solit_2012","prad_mskcc","sarc_mskcc","wt_target_2018")
Levels <- rev(Levels)
Difference_matrix$CancerType <- factor(Difference_matrix$CancerType, levels = Levels)

###Plot ballonplots:
setwd("~/Documents/GitHub/CancerDormancy/PanCancer_ElasticNetRegression/ValidationBallonPlot/Figure")
pdf("TCGA_cBioPortal_Ballonplot_v2.pdf", height = 7,width = 14)
p <- ggballoonplot(Difference_matrix, x = "Gene", y = "CancerType",
                   size = "minuslog10pval", fill = "QuiescenceScoreChange", size.range = c(3,10)) 
p +  scale_fill_gradient2(low = "#00441B", high = "#40004B", mid = "#F7F7F7")

dev.off()




p +  scale_fill_gradient2(low = "#40004B", high = "#00441B", mid = "#F7F7F7")












