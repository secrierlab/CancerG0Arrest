########################################################################
##BRCA analysis (BRCA subtypes association with APOBEC mutagenesis)
########################################################################

#Load required packages:
library('TCGAbiolinks')
library('BSgenome.Hsapiens.UCSC.hg38')
library(deconstructSigs)
library(ggpubr)


##############################
#Load mutational data
a <- "BRCA"
mutations <- GDCquery_Maf(a, pipelines = "mutect2")
# select columns of interest:
snvs <- data.frame(mutations[,c("Tumor_Sample_Barcode","Hugo_Symbol",
                                "Chromosome","Start_Position","End_Position",
                                "Variant_Classification","Variant_Type",
                                "Reference_Allele","Tumor_Seq_Allele1",
                                "Tumor_Seq_Allele2")])
# only interested in point mutations, not small insertions/deletions:
snvs <- snvs[which((snvs$Tumor_Seq_Allele2 %in% 
                      c("A","C","G","T"))&
                     (snvs$Reference_Allele %in%
                        c("A","C","G","T"))),]
# Use annotations to remove samples that are not primary tumours
snvs$Sample_Type_Code <- sapply(snvs$Tumor_Sample_Barcode,
                                function(tsb) substr(strsplit(tsb,'-')[[1]][4],1,2))
snvs$Sample_Type <- sapply(snvs$Sample_Type_Code,
                           function(x) ifelse(x=='01','PrimaryTumour',
                                              ifelse(x =='03', 'PrimaryBloodDerivedCancer',
                                                     ifelse(x %in% c('10','11'),'Normal',
                                                            ifelse(x=='06','Metastatic','Other')))))
snvs <- snvs[snvs$Sample_Type == 'PrimaryTumour',]
# Convert to deconstructSigs input:
sigs.input <- mut.to.sigs.input(mut.ref = snvs, 
                                sample.id = "Tumor_Sample_Barcode", 
                                chr = "Chromosome", 
                                pos = "Start_Position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele2",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)
sigs.input$Total <- rowSums(sigs.input)
sigs.input$Barcode <- rownames(sigs.input)



##############################
####Load mutational signatures 
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_MutationalSignatures/")
load("sigs.defaultnounknown.BRCA.RData")
sigs.defaultnounknown$Barcode <- rownames(sigs.defaultnounknown)
sigs.defaultnounknown <- merge(sigs.defaultnounknown, sigs.input,
                               by.x = "Barcode", by.y = "Barcode")
sigs.defaultnounknown$PatientID <- sapply(sigs.defaultnounknown$Barcode, function(x)
  paste(strsplit(x,"-")[[1]][1:3],collapse="-"))
sigs.defaultnounknown$SBS2_mutations <- sigs.defaultnounknown$Total * sigs.defaultnounknown$SBS2
sigs.defaultnounknown$log_SBS2_mutations <- log2(sigs.defaultnounknown$SBS2_mutations + 1)

#############################
####Split by cancer subtyper
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_ClinicalData/")
BRCA.clinical <- read.csv("combined_study_clinical_data_brca.csv", header = TRUE)
table(BRCA.clinical$Subtype)
BRCA_Basal <- BRCA.clinical[BRCA.clinical$Subtype %in% "BRCA_Basal",]
BRCA_Basal <- as.character(BRCA_Basal$Patient.ID)
BRCA_Her2 <- BRCA.clinical[BRCA.clinical$Subtype %in% "BRCA_Her2",]
BRCA_Her2 <- as.character(BRCA_Her2$Patient.ID)
BRCA_LumA <- BRCA.clinical[BRCA.clinical$Subtype %in% "BRCA_LumA",]
BRCA_LumA <- as.character(BRCA_LumA$Patient.ID)
BRCA_LumB <- BRCA.clinical[BRCA.clinical$Subtype %in% "BRCA_LumB",]
BRCA_LumB <- as.character(BRCA_LumB$Patient.ID)
sigs.defaultnounknown$Subtype <- sapply(sigs.defaultnounknown$PatientID, function(x)
  ifelse(x %in% BRCA_Basal, "BRCA_Basal",
         ifelse(x %in% BRCA_Her2, "BRCA_Her2",
                ifelse(x %in% BRCA_LumA, "BRCA_LumA",
                       ifelse(x %in% BRCA_LumB, "BRCA_LumB","Other")))))
table(sigs.defaultnounknown$Subtype)
sigs.defaultnounknown <- sigs.defaultnounknown[sigs.defaultnounknown$Subtype %in% c("BRCA_Basal","BRCA_Her2","BRCA_LumA","BRCA_LumB"),]



###Plot contribution of SBS2 in the 4 categories
library(ggpubr)
my_comparisons <- list( c("BRCA_Basal","BRCA_Her2"), c("BRCA_Basal","BRCA_LumA"), c("BRCA_Basal","BRCA_LumB"), c("BRCA_Her2","BRCA_LumA"), c("BRCA_Her2","BRCA_LumB"), c("BRCA_LumA","BRCA_LumB") )
setwd("~/Documents/GitHub/CancerDormancy/BRCA_Quiescence/Figures/")
p <- ggboxplot(sigs.defaultnounknown, x = "Subtype", y = "log_SBS2_mutations",
               color = "Subtype", palette =c("#44355B", "#274029", "#ECA72C","#EE5622"),
               add = "jitter")
pdf("BRCA_subtype_SBS2_prevelance.pdf",height = 7, width = 7)
p+stat_compare_means(comparisons = my_comparisons,label = "p.signif") # Add pairwise comparisons p-value
dev.off()
