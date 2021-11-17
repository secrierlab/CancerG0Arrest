#######################################################################################
#Fisher's exact test to check for depletion of DDR pathway gene mutations in quiescence
######################################################################################


#Load required packages:
library(biomaRt)
library(ggplot2)
library(TCGAbiolinks)



##Load quiescence score information along with kmeans clusters (obtained using by clusterin on ComBat treated data)
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
load("TCGA_common_QS_and_kmeans_clusters.RData")
z_score$quiescence_group <- z_score$cluster
solid_tumours <- unique(as.character(z_score$CancerType))


####Download maf data:
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


#Remove samples which have no mutational data
pancancer_snv <- maf_all
pancancer_snv$SampleID <- sapply(pancancer_snv$Tumor_Sample_Barcode, function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
z_score$SampleID <- sapply(rownames(z_score), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
snv.samples <- unique(as.character(pancancer_snv$SampleID))
z_score <- z_score[z_score$SampleID %in% snv.samples,]
Pancancer_samples <- z_score


##Select only mutations of interest which are likely to be deleterious
snvs_selected <- pancancer_snv[which(pancancer_snv$Variant_Classification
                                     %in% c("Frame_Shift_Del","Frame_Shift_Ins ","In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation")),]
#The following mutation types are ommited are it is more challenging to predict if they are deleterious:
#3'Flank, 3'UTR, 5'Flank, 5'UTR, IGR, Intron, RNA, Silent, Splice_Region, Splice_Site, Translation_Start_Site 



##Load CNV information:
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
load("combined_CNV.RData") #This is a toy example data with only 100 entries


#Load genes of interest for core DNA repair pathways
DDR <- read.csv("DDRpathways.csv",header = TRUE) #This data was obtained from Pearl et al 2015 (Nature Reviews Cancer)
#Select genes which are involved in repair pathways
DDR <- DDR[DDR$DDR %in% "Repair pathway",]
table(DDR$Pathway.1)
BER <- DDR[DDR$Pathway.1 %in% "BER",]
BER <- as.character(BER$Gene.ID)
Alt_NHEJ <- DDR[DDR$Pathway.1 %in% "Alt-NHEJ",]
Alt_NHEJ <- as.character(Alt_NHEJ$Gene.ID)
FA <- DDR[DDR$Pathway.1 %in% "FA (Fanconi anemia pathway)",]
FA <- as.character(FA$Gene.ID)
HR <- DDR[DDR$Pathway.1 %in% "HR (Homologous Recombination)",]
HR <- as.character(HR$Gene.ID)
MMR <- DDR[DDR$Pathway.1 %in% "MMR",]
MMR <- as.character(MMR$Gene.ID)
NER <- DDR[DDR$Pathway.1 %in% "NER",]
NER <- as.character(NER$Gene.ID)
NHEJ <- DDR[DDR$Pathway.1 %in% "NHEJ",]
NHEJ <- as.character(NHEJ$Gene.ID)
TLS <- DDR[DDR$Pathway.1 %in% "TLS",]
TLS <- as.character(TLS$Gene.ID)
DR <- DDR[DDR$Pathway.1 %in% "Direct Repair",]
DR <- as.character(DR$Gene.ID)



#Load genes of interest for associated pathways
DDR <- read.csv("DDRpathways.csv",header = TRUE) #This data was obtained from Pearl et al 2015 (Nature Reviews Cancer)
DDR <- DDR[DDR$DDR %in% "Associated process",]
table(DDR$Process)
CR <- DDR[DDR$Process %in% "Chromatin remodelling",]
CR <- as.character(CR$Gene.ID)
TM <- DDR[DDR$Process %in% "Telomere maintenance",]
TM <- as.character(TM$Gene.ID)
CPF <- DDR[DDR$Process %in% "Checkpoint factors",]
CPF <- as.character(CPF$Gene.ID)
UR <- DDR[DDR$Process %in% "Ubiquitin response",]
UR <- as.character(UR$Gene.ID)
p53 <- DDR[DDR$Process %in% "p53 pathway",]
p53 <- as.character(p53$Gene.ID)
CS <- DDR[DDR$Process %in% "Chromosome segregation",]
CS <- as.character(CS$Gene.ID)


#Convert gene lists to ENSG:
human = useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "useast.ensembl.org")

results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=p53, mart=human, filters = "hgnc_symbol")
p53_ENSG <- results$ensembl_gene_id

results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=Alt_NHEJ, mart=human, filters = "hgnc_symbol")
Alt_NHEJ_ENSG <- results$ensembl_gene_id

results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=BER, mart=human, filters = "hgnc_symbol")
BER_ENSG <- results$ensembl_gene_id

results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=FA, mart=human, filters = "hgnc_symbol")
FA_ENSG <- results$ensembl_gene_id

results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=HR, mart=human, filters = "hgnc_symbol")
HR_ENSG <- results$ensembl_gene_id

results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=FA, mart=human, filters = "hgnc_symbol")
FA_ENSG <- results$ensembl_gene_id

results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=MMR, mart=human, filters = "hgnc_symbol")
MMR_ENSG <- results$ensembl_gene_id

results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=NER, mart=human, filters = "hgnc_symbol")
NER_ENSG <- results$ensembl_gene_id

results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=NHEJ, mart=human, filters = "hgnc_symbol")
NHEJ_ENSG <- results$ensembl_gene_id

results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=TLS, mart=human, filters = "hgnc_symbol")
TLS_ENSG <- results$ensembl_gene_id

results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=DR, mart=human, filters = "hgnc_symbol")
DR_ENSG <- results$ensembl_gene_id

results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=CR, mart=human, filters = "hgnc_symbol")
CR_ENSG <- results$ensembl_gene_id

results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=TM, mart=human, filters = "hgnc_symbol")
TM_ENSG <- results$ensembl_gene_id

results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=CPF, mart=human, filters = "hgnc_symbol")
CPF_ENSG <- results$ensembl_gene_id

results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=UR, mart=human, filters = "hgnc_symbol")
UR_ENSG <- results$ensembl_gene_id

results<- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=CS, mart=human, filters = "hgnc_symbol")
CS_ENSG <- results$ensembl_gene_id


########Add columns to the Pancancer Samples dataframe indicating whether there are CN losses or deleterious SNVs in any of the DDR pathways
#p53
revelant_mutations <- snvs_selected[snvs_selected$Hugo_Symbol %in% p53,]
mutated_samples <- revelant_mutations$SampleID
mutated_samples <- unique(mutated_samples)

CN_loss_samples <- CNV[,colnames(CNV) %in% p53_ENSG]
for (i in colnames(CN_loss_samples)) {
  CN_loss_samples[[i]] <- sapply(CN_loss_samples[[i]], function(x)
    ifelse(x <= -1, 1,0))
}
CN_loss_samples$Total <- rowSums(CN_loss_samples)
CN_loss_samples <- CN_loss_samples[CN_loss_samples$Total >= 1,]
CN_loss_samples <- sapply(rownames(CN_loss_samples), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
Pancancer_samples$p53 <- sapply(Pancancer_samples$SampleID, function(x)
  ifelse(x %in% c(CN_loss_samples,mutated_samples), "MT","WT"))
#CPF
revelant_mutations <- snvs_selected[snvs_selected$Hugo_Symbol %in% CPF,]
mutated_samples <- revelant_mutations$SampleID
mutated_samples <- unique(mutated_samples)

CN_loss_samples <- CNV[,colnames(CNV) %in% CPF_ENSG]
for (i in colnames(CN_loss_samples)) {
  CN_loss_samples[[i]] <- sapply(CN_loss_samples[[i]], function(x)
    ifelse(x <= -1, 1,0))
}
CN_loss_samples$Total <- rowSums(CN_loss_samples)
CN_loss_samples <- CN_loss_samples[CN_loss_samples$Total >= 1,]
CN_loss_samples <- sapply(rownames(CN_loss_samples), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
Pancancer_samples$CPF <- sapply(Pancancer_samples$SampleID, function(x)
  ifelse(x %in% c(CN_loss_samples,mutated_samples), "MT","WT"))
#UR
revelant_mutations <- snvs_selected[snvs_selected$Hugo_Symbol %in% UR,]
mutated_samples <- revelant_mutations$SampleID
mutated_samples <- unique(mutated_samples)
CN_loss_samples <- CNV[,colnames(CNV) %in% UR_ENSG]
for (i in colnames(CN_loss_samples)) {
  CN_loss_samples[[i]] <- sapply(CN_loss_samples[[i]], function(x)
    ifelse(x <= -1, 1,0))
}
CN_loss_samples$Total <- rowSums(CN_loss_samples)
CN_loss_samples <- CN_loss_samples[CN_loss_samples$Total >= 1,]
CN_loss_samples <- sapply(rownames(CN_loss_samples), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
Pancancer_samples$UR <- sapply(Pancancer_samples$SampleID, function(x)
  ifelse(x %in% c(CN_loss_samples,mutated_samples), "MT","WT"))
#NER
revelant_mutations <- snvs_selected[snvs_selected$Hugo_Symbol %in% NER,]
mutated_samples <- revelant_mutations$SampleID
mutated_samples <- unique(mutated_samples)
CN_loss_samples <- CNV[,colnames(CNV) %in% NER_ENSG]
for (i in colnames(CN_loss_samples)) {
  CN_loss_samples[[i]] <- sapply(CN_loss_samples[[i]], function(x)
    ifelse(x <= -1, 1,0))
}
CN_loss_samples$Total <- rowSums(CN_loss_samples)
CN_loss_samples <- CN_loss_samples[CN_loss_samples$Total >= 1,]
CN_loss_samples <- sapply(rownames(CN_loss_samples), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
Pancancer_samples$NER <- sapply(Pancancer_samples$SampleID, function(x)
  ifelse(x %in% c(CN_loss_samples,mutated_samples), "MT","WT"))
#MMR
revelant_mutations <- snvs_selected[snvs_selected$Hugo_Symbol %in% MMR,]
mutated_samples <- revelant_mutations$SampleID
mutated_samples <- unique(mutated_samples)

CN_loss_samples <- CNV[,colnames(CNV) %in% MMR_ENSG]
for (i in colnames(CN_loss_samples)) {
  CN_loss_samples[[i]] <- sapply(CN_loss_samples[[i]], function(x)
    ifelse(x <= -1, 1,0))
}
CN_loss_samples$Total <- rowSums(CN_loss_samples)
CN_loss_samples <- CN_loss_samples[CN_loss_samples$Total >= 1,]
CN_loss_samples <- sapply(rownames(CN_loss_samples), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
Pancancer_samples$MMR <- sapply(Pancancer_samples$SampleID, function(x)
  ifelse(x %in% c(CN_loss_samples,mutated_samples), "MT","WT"))
#BER
revelant_mutations <- snvs_selected[snvs_selected$Hugo_Symbol %in% BER,]
mutated_samples <- revelant_mutations$SampleID
mutated_samples <- unique(mutated_samples)

CN_loss_samples <- CNV[,colnames(CNV) %in% BER_ENSG]
for (i in colnames(CN_loss_samples)) {
  CN_loss_samples[[i]] <- sapply(CN_loss_samples[[i]], function(x)
    ifelse(x <= -1, 1,0))
}
CN_loss_samples$Total <- rowSums(CN_loss_samples)
CN_loss_samples <- CN_loss_samples[CN_loss_samples$Total >= 1,]
CN_loss_samples <- sapply(rownames(CN_loss_samples), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
Pancancer_samples$BER <- sapply(Pancancer_samples$SampleID, function(x)
  ifelse(x %in% c(CN_loss_samples,mutated_samples), "MT","WT"))
#DR
revelant_mutations <- snvs_selected[snvs_selected$Hugo_Symbol %in% DR,]
mutated_samples <- revelant_mutations$SampleID
mutated_samples <- unique(mutated_samples)

CN_loss_samples <- CNV[,colnames(CNV) %in% DR_ENSG]
for (i in colnames(CN_loss_samples)) {
  CN_loss_samples[[i]] <- sapply(CN_loss_samples[[i]], function(x)
    ifelse(x <= -1, 1,0))
}
CN_loss_samples$Total <- rowSums(CN_loss_samples)
CN_loss_samples <- CN_loss_samples[CN_loss_samples$Total >= 1,]
CN_loss_samples <- sapply(rownames(CN_loss_samples), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
Pancancer_samples$DR <- sapply(Pancancer_samples$SampleID, function(x)
  ifelse(x %in% c(CN_loss_samples,mutated_samples), "MT","WT"))
#TM
revelant_mutations <- snvs_selected[snvs_selected$Hugo_Symbol %in% TM,]
mutated_samples <- revelant_mutations$SampleID
mutated_samples <- unique(mutated_samples)

CN_loss_samples <- CNV[,colnames(CNV) %in% TM_ENSG]
for (i in colnames(CN_loss_samples)) {
  CN_loss_samples[[i]] <- sapply(CN_loss_samples[[i]], function(x)
    ifelse(x <= -1, 1,0))
}
CN_loss_samples$Total <- rowSums(CN_loss_samples)
CN_loss_samples <- CN_loss_samples[CN_loss_samples$Total >= 1,]
CN_loss_samples <- sapply(rownames(CN_loss_samples), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
Pancancer_samples$TM <- sapply(Pancancer_samples$SampleID, function(x)
  ifelse(x %in% c(CN_loss_samples,mutated_samples), "MT","WT"))
#CS
revelant_mutations <- snvs_selected[snvs_selected$Hugo_Symbol %in% CS,]
mutated_samples <- revelant_mutations$SampleID
mutated_samples <- unique(mutated_samples)

CN_loss_samples <- CNV[,colnames(CNV) %in% CS_ENSG]
for (i in colnames(CN_loss_samples)) {
  CN_loss_samples[[i]] <- sapply(CN_loss_samples[[i]], function(x)
    ifelse(x <= -1, 1,0))
}
CN_loss_samples$Total <- rowSums(CN_loss_samples)
CN_loss_samples <- CN_loss_samples[CN_loss_samples$Total >= 1,]
CN_loss_samples <- sapply(rownames(CN_loss_samples), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
Pancancer_samples$CS <- sapply(Pancancer_samples$SampleID, function(x)
  ifelse(x %in% c(CN_loss_samples,mutated_samples), "MT","WT"))
#CR
revelant_mutations <- snvs_selected[snvs_selected$Hugo_Symbol %in% CR,]
mutated_samples <- revelant_mutations$SampleID
mutated_samples <- unique(mutated_samples)

CN_loss_samples <- CNV[,colnames(CNV) %in% CR_ENSG]
for (i in colnames(CN_loss_samples)) {
  CN_loss_samples[[i]] <- sapply(CN_loss_samples[[i]], function(x)
    ifelse(x <= -1, 1,0))
}
CN_loss_samples$Total <- rowSums(CN_loss_samples)
CN_loss_samples <- CN_loss_samples[CN_loss_samples$Total >= 1,]
CN_loss_samples <- sapply(rownames(CN_loss_samples), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
Pancancer_samples$CR <- sapply(Pancancer_samples$SampleID, function(x)
  ifelse(x %in% c(CN_loss_samples,mutated_samples), "MT","WT"))
#TLS
revelant_mutations <- snvs_selected[snvs_selected$Hugo_Symbol %in% TLS,]
mutated_samples <- revelant_mutations$SampleID
mutated_samples <- unique(mutated_samples)

CN_loss_samples <- CNV[,colnames(CNV) %in% TLS_ENSG]
for (i in colnames(CN_loss_samples)) {
  CN_loss_samples[[i]] <- sapply(CN_loss_samples[[i]], function(x)
    ifelse(x <= -1, 1,0))
}
CN_loss_samples$Total <- rowSums(CN_loss_samples)
CN_loss_samples <- CN_loss_samples[CN_loss_samples$Total >= 1,]
CN_loss_samples <- sapply(rownames(CN_loss_samples), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
Pancancer_samples$TLS <- sapply(Pancancer_samples$SampleID, function(x)
  ifelse(x %in% c(CN_loss_samples,mutated_samples), "MT","WT"))
#NHEJ
revelant_mutations <- snvs_selected[snvs_selected$Hugo_Symbol %in% NHEJ,]
mutated_samples <- revelant_mutations$SampleID
mutated_samples <- unique(mutated_samples)

CN_loss_samples <- CNV[,colnames(CNV) %in% NHEJ_ENSG]
for (i in colnames(CN_loss_samples)) {
  CN_loss_samples[[i]] <- sapply(CN_loss_samples[[i]], function(x)
    ifelse(x <= -1, 1,0))
}
CN_loss_samples$Total <- rowSums(CN_loss_samples)
CN_loss_samples <- CN_loss_samples[CN_loss_samples$Total >= 1,]
CN_loss_samples <- sapply(rownames(CN_loss_samples), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
Pancancer_samples$NHEJ <- sapply(Pancancer_samples$SampleID, function(x)
  ifelse(x %in% c(CN_loss_samples,mutated_samples), "MT","WT"))
#Alt-NHEJ
revelant_mutations <- snvs_selected[snvs_selected$Hugo_Symbol %in% Alt_NHEJ,]
mutated_samples <- revelant_mutations$SampleID
mutated_samples <- unique(mutated_samples)

CN_loss_samples <- CNV[,colnames(CNV) %in% Alt_NHEJ_ENSG]
for (i in colnames(CN_loss_samples)) {
  CN_loss_samples[[i]] <- sapply(CN_loss_samples[[i]], function(x)
    ifelse(x <= -1, 1,0))
}
CN_loss_samples$Total <- rowSums(CN_loss_samples)
CN_loss_samples <- CN_loss_samples[CN_loss_samples$Total >= 1,]
CN_loss_samples <- sapply(rownames(CN_loss_samples), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
Pancancer_samples$Alt_NHEJ <- sapply(Pancancer_samples$SampleID, function(x)
  ifelse(x %in% c(CN_loss_samples,mutated_samples), "MT","WT"))
#FA
revelant_mutations <- snvs_selected[snvs_selected$Hugo_Symbol %in% FA,]
mutated_samples <- revelant_mutations$SampleID
mutated_samples <- unique(mutated_samples)

CN_loss_samples <- CNV[,colnames(CNV) %in% FA_ENSG]
for (i in colnames(CN_loss_samples)) {
  CN_loss_samples[[i]] <- sapply(CN_loss_samples[[i]], function(x)
    ifelse(x <= -1, 1,0))
}
CN_loss_samples$Total <- rowSums(CN_loss_samples)
CN_loss_samples <- CN_loss_samples[CN_loss_samples$Total >= 1,]
CN_loss_samples <- sapply(rownames(CN_loss_samples), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
Pancancer_samples$FA <- sapply(Pancancer_samples$SampleID, function(x)
  ifelse(x %in% c(CN_loss_samples,mutated_samples), "MT","WT"))
#HR
revelant_mutations <- snvs_selected[snvs_selected$Hugo_Symbol %in% HR,]
mutated_samples <- revelant_mutations$SampleID
mutated_samples <- unique(mutated_samples)

CN_loss_samples <- CNV[,colnames(CNV) %in% HR]
for (i in colnames(CN_loss_samples)) {
  CN_loss_samples[[i]] <- sapply(CN_loss_samples[[i]], function(x)
    ifelse(x <= -1, 1,0))
}
CN_loss_samples$Total <- rowSums(CN_loss_samples)
CN_loss_samples <- CN_loss_samples[CN_loss_samples$Total >= 1,]
CN_loss_samples <- sapply(rownames(CN_loss_samples), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
Pancancer_samples$HR <- sapply(Pancancer_samples$SampleID, function(x)
  ifelse(x %in% c(CN_loss_samples,mutated_samples), "MT","WT"))





#####################
##FISHER'S EXACT TEST:
#####################

#initiatlise the variables:
p_value <- NULL
OR <- NULL
Upper_confidence_interval <- NULL
Lower_confidence_interval <- NULL

Pathways <- c("Alt_NHEJ", "BER", "FA","HR","MMR","NER","NHEJ","TLS", "DR","CR","TM","CPF","UR","p53","CS")

###Fisher's exact test
for (i in Pathways) {
  test_table <- table(Pancancer_samples[[i]], Pancancer_samples$quiescence_group)
  print(test_table)
  test <- fisher.test(test_table)
  p <- test$p.value
  p_value <- c(p_value,p)
  or <- test$estimate
  OR <- c(OR, or)
  down <- test$conf.int[1]
  Lower_confidence_interval <- c(Lower_confidence_interval, down)
  up <- test$conf.int[2]
  Upper_confidence_interval <- c(Upper_confidence_interval, up)
  
}

Adjusted_p_value <- p.adjust(p_value, method = "BH")
setwd("~/Documents/GitHub/CancerDormancy/PanCancer_ElasticNetRegression/FeatureExploration/Results/")
Fishers_test_results <- data.frame(Pathways, OR, p_value, Adjusted_p_value, Lower_confidence_interval, Upper_confidence_interval)
save(Fishers_test_results, file = "FishersExactTest_DDRmutation_vs_quiescence.RData")



#########################
##Visualising the results

#Filter the results to include significant pathways only:
Fishers_test_results <- Fishers_test_results[Fishers_test_results$Adjusted_p_value < 0.05,]
Fishers_test_results$yaxis <- 1:12


#Plot the results:
setwd("~/Documents/GitHub/CancerDormancy/PanCancer_ElasticNetRegression/FeatureExploration/Figures/")
pdf("FishersExactTest_DDRmutation_vs_quiescence.pdf", width = 3, height = 3)
p <-ggplot(Fishers_test_results, aes(x = OR, y = yaxis))
p + geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = Upper_confidence_interval, xmin = Lower_confidence_interval), size = .5, height = .2, color = "gray50") +
  geom_point(size = 1, color = "darkorange") +
  theme_classic() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = Fishers_test_results$yaxis, labels = Fishers_test_results$Pathways) +
  scale_x_continuous(limits = c(0.5,1.5))+
  ylab("") +
  xlab("Odds ratio")
dev.off()





