###################################################################
##Combine RNA-seq data for all cancer types into a single dataframe
###################################################################


#Load packages:
library(SummarizedExperiment)
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")


#Load RData files with expression data for all solid tumours
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/ACC")
load("TCGA-ACCTranscriptome_ProfilingThu_Oct_17_12:45:41_2019.RData")
ACC <- as.data.frame(t(assay(data)))
ACC$SampleType <- sapply(rownames(ACC), function(x)
  strsplit(x,"-")[[1]][4])
table(ACC$SampleType)
ACC <- ACC[which(ACC$SampleType %in% c("01A","01B","01C")),]
ACC$SampleType <- NULL
ACC <- as.matrix(t(ACC))
ACC <- log2(ACC+1)
ACC <- data.frame(t(ACC))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/BLCA/")
load("TCGA-BLCATranscriptome_ProfilingThu_Oct_17_12:55:20_2019.RData")
BLCA <- as.data.frame(t(assay(data)))
BLCA$SampleType <- sapply(rownames(BLCA), function(x)
  strsplit(x,"-")[[1]][4])
table(BLCA$SampleType)
BLCA <- BLCA[which(BLCA$SampleType %in% c("01A","01B","01C")),]
BLCA$SampleType <- NULL
BLCA <- as.matrix(t(BLCA))
BLCA <- log2(BLCA+1)
BLCA <- data.frame(t(BLCA))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/BRCA/")
load("TCGA-BRCATranscriptome_ProfilingThu_Oct_17_14:32:08_2019.RData")
BRCA <- as.data.frame(t(assay(data)))
BRCA$SampleType <- sapply(rownames(BRCA), function(x)
  strsplit(x,"-")[[1]][4])
table(BRCA$SampleType)
BRCA <- BRCA[which(BRCA$SampleType %in% c("01A","01B","01C")),]
BRCA$SampleType <- NULL
BRCA <- as.matrix(t(BRCA))
BRCA <- log2(BRCA+1)
BRCA <- data.frame(t(BRCA))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/CESC/")
load("TCGA-CESCTranscriptome_ProfilingThu_Oct_17_14:54:55_2019.RData")
CESC <- as.data.frame(t(assay(data)))
CESC$SampleType <- sapply(rownames(CESC), function(x)
  strsplit(x,"-")[[1]][4])
table(CESC$SampleType)
CESC <- CESC[which(CESC$SampleType %in% c("01A","01B","01C")),]
CESC$SampleType <- NULL
CESC <- as.matrix(t(CESC))
CESC <- log2(CESC+1)
CESC <- data.frame(t(CESC))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/CHOL/")
load("TCGA-CHOLTranscriptome_ProfilingThu_Oct_17_15:00:20_2019.RData")
CHOL <- as.data.frame(t(assay(data)))
CHOL$SampleType <- sapply(rownames(CHOL), function(x)
  strsplit(x,"-")[[1]][4])
table(CHOL$SampleType)
CHOL <- CHOL[which(CHOL$SampleType %in% c("01A","01B","01C")),]
CHOL$SampleType <- NULL
CHOL <- as.matrix(t(CHOL))
CHOL <- log2(CHOL+1)
CHOL <- data.frame(t(CHOL))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/COAD/")
load("TCGA-COADTranscriptome_ProfilingThu_Oct_17_15:11:13_2019.RData")
COAD <- as.data.frame(t(assay(data)))
COAD$SampleType <- sapply(rownames(COAD), function(x)
  strsplit(x,"-")[[1]][4])
table(COAD$SampleType)
COAD <- COAD[which(COAD$SampleType %in% c("01A","01B","01C")),]
COAD$SampleType <- NULL
COAD <- as.matrix(t(COAD))
COAD <- log2(COAD+1)
COAD <- data.frame(t(COAD))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/ESCA/")
load("TCGA-ESCATranscriptome_ProfilingThu_Oct_17_15:19:42_2019.RData")
ESCA <- as.data.frame(t(assay(data)))
ESCA$SampleType <- sapply(rownames(ESCA), function(x)
  strsplit(x,"-")[[1]][4])
table(ESCA$SampleType)
ESCA <- ESCA[which(ESCA$SampleType %in% c("01A","01B","01C")),]
ESCA$SampleType <- NULL
ESCA <- as.matrix(t(ESCA))
ESCA <- log2(ESCA+1)
ESCA <- data.frame(t(ESCA))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/GBM/")
load("TCGA-GBMTranscriptome_ProfilingThu_Oct_17_15:29:03_2019.RData")
GBM <- as.data.frame(t(assay(data)))
GBM$SampleType <- sapply(rownames(GBM), function(x)
  strsplit(x,"-")[[1]][4])
table(GBM$SampleType)
GBM <- GBM[which(GBM$SampleType %in% c("01A","01B","01C")),]
GBM$SampleType <- NULL
GBM <- as.matrix(t(GBM))
GBM <- log2(GBM+1)
GBM <- data.frame(t(GBM))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/HNSC/")
load("TCGA-HNSCTranscriptome_ProfilingThu_Oct_17_15:44:26_2019.RData")
HNSC <- as.data.frame(t(assay(data)))
HNSC$SampleType <- sapply(rownames(HNSC), function(x)
  strsplit(x,"-")[[1]][4])
table(HNSC$SampleType)
HNSC <- HNSC[which(HNSC$SampleType %in% c("01A","01B","01C")),]
HNSC$SampleType <- NULL
HNSC <- as.matrix(t(HNSC))
HNSC <- log2(HNSC+1)
HNSC <- data.frame(t(HNSC))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/KICH/")
load("TCGA-KICHTranscriptome_ProfilingThu_Oct_17_15:56:59_2019.RData")
KICH <- as.data.frame(t(assay(data)))
KICH$SampleType <- sapply(rownames(KICH), function(x)
  strsplit(x,"-")[[1]][4])
table(KICH$SampleType)
KICH <- KICH[which(KICH$SampleType %in% c("01A","01B","01C")),]
KICH$SampleType <- NULL
KICH <- as.matrix(t(KICH))
KICH <- log2(KICH+1)
KICH <- data.frame(t(KICH))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/KIRC/")
load("TCGA-KIRCTranscriptome_ProfilingThu_Oct_17_16:14:04_2019.RData")
KIRC <- as.data.frame(t(assay(data)))
KIRC$SampleType <- sapply(rownames(KIRC), function(x)
  strsplit(x,"-")[[1]][4])
table(KIRC$SampleType)
KIRC <- KIRC[which(KIRC$SampleType %in% c("01A","01B","01C")),]
KIRC$SampleType <- NULL
KIRC <- as.matrix(t(KIRC))
KIRC <- log2(KIRC+1)
KIRC <- data.frame(t(KIRC))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/KIRP/")
load("TCGA-KIRPTranscriptome_ProfilingThu_Oct_17_16:35:58_2019.RData")
KIRP <- as.data.frame(t(assay(data)))
KIRP$SampleType <- sapply(rownames(KIRP), function(x)
  strsplit(x,"-")[[1]][4])
table(KIRP$SampleType)
KIRP <- KIRP[which(KIRP$SampleType %in% c("01A","01B","01C")),]
KIRP$SampleType <- NULL
KIRP <- as.matrix(t(KIRP))
KIRP <- log2(KIRP+1)
KIRP <- data.frame(t(KIRP))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/LGG/")
load("TCGA-LGGTranscriptome_ProfilingThu_Oct_17_13:59:44_2019.RData")
LGG <- as.data.frame(t(assay(data)))
LGG$SampleType <- sapply(rownames(LGG), function(x)
  strsplit(x,"-")[[1]][4])
table(LGG$SampleType)
LGG <- LGG[which(LGG$SampleType %in% c("01A","01B","01C")),]
LGG$SampleType <- NULL
LGG <- as.matrix(t(LGG))
LGG <- log2(LGG+1)
LGG <- data.frame(t(LGG))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/LIHC/")
load("TCGA-LIHCTranscriptome_ProfilingFri_Oct_18_10:13:28_2019.RData")
LIHC <- as.data.frame(t(assay(data)))
LIHC$SampleType <- sapply(rownames(LIHC), function(x)
  strsplit(x,"-")[[1]][4])
table(LIHC$SampleType)
LIHC <- LIHC[which(LIHC$SampleType %in% c("01A","01B","01C")),]
LIHC$SampleType <- NULL
LIHC <- as.matrix(t(LIHC))
LIHC <- log2(LIHC+1)
LIHC <- data.frame(t(LIHC))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/LUAD/")
load("TCGA-LUADTranscriptome_ProfilingFri_Oct_18_10:28:26_2019.RData")
LUAD <- as.data.frame(t(assay(data)))
LUAD$SampleType <- sapply(rownames(LUAD), function(x)
  strsplit(x,"-")[[1]][4])
table(LUAD$SampleType)
LUAD <- LUAD[which(LUAD$SampleType %in% c("01A","01B","01C")),]
LUAD$SampleType <- NULL
LUAD <- as.matrix(t(LUAD))
LUAD <- log2(LUAD+1)
LUAD <- data.frame(t(LUAD))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/LUSC/")
load("TCGA-LUSCTranscriptome_ProfilingFri_Oct_18_10:50:12_2019.RData")
LUSC <- as.data.frame(t(assay(data)))
LUSC$SampleType <- sapply(rownames(LUSC), function(x)
  strsplit(x,"-")[[1]][4])
table(LUSC$SampleType)
LUSC <- LUSC[which(LUSC$SampleType %in% c("01A","01B","01C")),]
LUSC$SampleType <- NULL
LUSC <- as.matrix(t(LUSC))
LUSC <- log2(LUSC+1)
LUSC <- data.frame(t(LUSC))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/MESO/")
load("TCGA-MESOTranscriptome_ProfilingFri_Oct_18_11:45:50_2019.RData")
MESO <- as.data.frame(t(assay(data)))
MESO$SampleType <- sapply(rownames(MESO), function(x)
  strsplit(x,"-")[[1]][4])
table(MESO$SampleType)
MESO <- MESO[which(MESO$SampleType %in% c("01A","01B","01C")),]
MESO$SampleType <- NULL
MESO <- as.matrix(t(MESO))
MESO <- log2(MESO+1)
MESO <- data.frame(t(MESO))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/OV/")
load("TCGA-OVTranscriptome_ProfilingFri_Oct_18_16:12:06_2019.RData")
OV <- as.data.frame(t(assay(data)))
OV$SampleType <- sapply(rownames(OV), function(x)
  strsplit(x,"-")[[1]][4])
table(OV$SampleType)
OV <- OV[which(OV$SampleType %in% c("01A","01B","01C")),]
OV$SampleType <- NULL
OV <- as.matrix(t(OV))
OV <- log2(OV+1)
OV <- data.frame(t(OV))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/PAAD/")
load("TCGA-PAADTranscriptome_ProfilingFri_Oct_18_16:39:53_2019.RData")
PAAD <- as.data.frame(t(assay(data)))
PAAD$SampleType <- sapply(rownames(PAAD), function(x)
  strsplit(x,"-")[[1]][4])
table(PAAD$SampleType)
PAAD <- PAAD[which(PAAD$SampleType %in% c("01A","01B","01C")),]
PAAD$SampleType <- NULL
PAAD <- as.matrix(t(PAAD))
PAAD <- log2(PAAD+1)
PAAD <- data.frame(t(PAAD))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/PCPG/")
load("TCGA-PCPGTranscriptome_ProfilingSun_Oct_20_14:19:25_2019.RData")
PCPG <- as.data.frame(t(assay(data)))
PCPG$SampleType <- sapply(rownames(PCPG), function(x)
  strsplit(x,"-")[[1]][4])
table(PCPG$SampleType)
PCPG <- PCPG[which(PCPG$SampleType %in% c("01A","01B","01C")),]
PCPG$SampleType <- NULL
PCPG <- as.matrix(t(PCPG))
PCPG <- log2(PCPG+1)
PCPG <- data.frame(t(PCPG))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/PRAD/")
load("TCGA-PRADTranscriptome_ProfilingSun_Oct_20_14:41:06_2019.RData")
PRAD <- as.data.frame(t(assay(data)))
PRAD$SampleType <- sapply(rownames(PRAD), function(x)
  strsplit(x,"-")[[1]][4])
table(PRAD$SampleType)
PRAD <- PRAD[which(PRAD$SampleType %in% c("01A","01B","01C")),]
PRAD$SampleType <- NULL
PRAD <- as.matrix(t(PRAD))
PRAD <- log2(PRAD+1)
PRAD <- data.frame(t(PRAD))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/READ/")
load("TCGA-READTranscriptome_ProfilingSun_Oct_20_16:43:45_2019.RData")
READ <- as.data.frame(t(assay(data)))
READ$SampleType <- sapply(rownames(READ), function(x)
  strsplit(x,"-")[[1]][4])
table(READ$SampleType)
READ <- READ[which(READ$SampleType %in% c("01A","01B","01C")),]
READ$SampleType <- NULL
READ <- as.matrix(t(READ))
READ <- log2(READ+1)
READ <- data.frame(t(READ))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/SARC/")
load("TCGA-SARCTranscriptome_ProfilingSun_Oct_20_18:23:40_2019.RData")
SARC <- as.data.frame(t(assay(data)))
SARC$SampleType <- sapply(rownames(SARC), function(x)
  strsplit(x,"-")[[1]][4])
table(SARC$SampleType)
SARC <- SARC[which(SARC$SampleType %in% c("01A","01B","01C")),]
SARC$SampleType <- NULL
SARC <- as.matrix(t(SARC))
SARC <- log2(SARC+1)
SARC <- data.frame(t(SARC))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/SKCM/")
load("TCGA-SKCMTranscriptome_ProfilingSun_Oct_20_20:32:33_2019.RData")
SKCM <- as.data.frame(t(assay(data)))
SKCM$SampleType <- sapply(rownames(SKCM), function(x)
  strsplit(x,"-")[[1]][4])
table(SKCM$SampleType)
SKCM <- SKCM[which(SKCM$SampleType %in% c("01A","01B","01C")),]
SKCM$SampleType <- NULL
SKCM <- as.matrix(t(SKCM))
SKCM <- log2(SKCM+1)
SKCM <- data.frame(t(SKCM))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/STAD/")
load("TCGA-STADTranscriptome_ProfilingSun_Oct_20_23:32:03_2019.RData")
STAD <- as.data.frame(t(assay(data)))
STAD$SampleType <- sapply(rownames(STAD), function(x)
  strsplit(x,"-")[[1]][4])
table(STAD$SampleType)
STAD <- STAD[which(STAD$SampleType %in% c("01A","01B","01C")),]
STAD$SampleType <- NULL
STAD <- as.matrix(t(STAD))
STAD <- log2(STAD+1)
STAD <- data.frame(t(STAD))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/TGCT/")
load("TCGA-TGCTTranscriptome_ProfilingSun_Oct_20_23:44:13_2019.RData")
TGCT <- as.data.frame(t(assay(data)))
TGCT$SampleType <- sapply(rownames(TGCT), function(x)
  strsplit(x,"-")[[1]][4])
table(TGCT$SampleType)
TGCT <- TGCT[which(TGCT$SampleType %in% c("01A","01B","01C")),]
TGCT$SampleType <- NULL
TGCT <- as.matrix(t(TGCT))
TGCT <- log2(TGCT+1)
TGCT <- data.frame(t(TGCT))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/THCA/")
load("TCGA-THCATranscriptome_ProfilingMon_Oct_21_08:04:11_2019.RData")
THCA <- as.data.frame(t(assay(data)))
THCA$SampleType <- sapply(rownames(THCA), function(x)
  strsplit(x,"-")[[1]][4])
table(THCA$SampleType)
THCA <- THCA[which(THCA$SampleType %in% c("01A","01B","01C")),]
THCA$SampleType <- NULL
THCA <- as.matrix(t(THCA))
THCA <- log2(THCA+1)
THCA <- data.frame(t(THCA))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/THYM/")
load("TCGA-THYMTranscriptome_ProfilingSun_Oct_20_23:52:16_2019.RData")
THYM <- as.data.frame(t(assay(data)))
THYM$SampleType <- sapply(rownames(THYM), function(x)
  strsplit(x,"-")[[1]][4])
table(THYM$SampleType)
THYM <- THYM[which(THYM$SampleType %in% c("01A","01B","01C")),]
THYM$SampleType <- NULL
THYM <- as.matrix(t(THYM))
THYM <- log2(THYM+1)
THYM <- data.frame(t(THYM))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/UCEC/")
load("TCGA-UCECTranscriptome_ProfilingMon_Oct_21_09:54:58_2019.RData")
UCEC <- as.data.frame(t(assay(data)))
UCEC$SampleType <- sapply(rownames(UCEC), function(x)
  strsplit(x,"-")[[1]][4])
table(UCEC$SampleType)
UCEC <- UCEC[which(UCEC$SampleType %in% c("01A","01B","01C")),]
UCEC$SampleType <- NULL
UCEC <- as.matrix(t(UCEC))
UCEC <- log2(UCEC+1)
UCEC <- data.frame(t(UCEC))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/UCS/")
load("TCGA-UCSTranscriptome_ProfilingMon_Oct_21_08:21:24_2019.RData")
UCS <- as.data.frame(t(assay(data)))
UCS$SampleType <- sapply(rownames(UCS), function(x)
  strsplit(x,"-")[[1]][4])
table(UCS$SampleType)
UCS <- UCS[which(UCS$SampleType %in% c("01A","01B","01C")),]
UCS$SampleType <- NULL
UCS <- as.matrix(t(UCS))
UCS <- log2(UCS+1)
UCS <- data.frame(t(UCS))
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/UVM")
load("TCGA-UVMTranscriptome_ProfilingMon_Oct_21_10:12:39_2019.RData")
UVM <- as.data.frame(t(assay(data)))
UVM$SampleType <- sapply(rownames(UVM), function(x)
  strsplit(x,"-")[[1]][4])
table(UVM$SampleType)
UVM <- UVM[which(UVM$SampleType %in% c("01A","01B","01C")),]
UVM$SampleType <- NULL
UVM <- as.matrix(t(UVM))
UVM <- log2(UVM+1)
UVM <- data.frame(t(UVM))

setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/Combined_TCGA_expression_data")
##Make sure that columns of each data frame are the same and are in the same order
all(colnames(ACC) == colnames(THCA))
combined_data <- rbind(ACC, BLCA, BRCA, CESC, CHOL, COAD, ESCA, GBM, HNSC, KICH, KIRC, KIRP, LGG, LIHC, LUAD, LUSC, MESO, OV, PAAD, PCPG, PRAD, READ, SARC, SKCM, STAD, TGCT,THCA, THYM, UCEC, UCS, UVM)
cancer_type <- c(rep("ACC", 79), rep("BLCA",414), rep("BRCA",1102), rep("CESC",304), rep("CHOL",36), rep("COAD",478), rep("ESCA",161), rep("GBM",156), rep("HNSC",500), rep("KICH",65), rep("KIRC",538), rep("KIRP",288), rep("LGG",511), rep("LIHC",371), rep("LUAD",533), rep("LUSC",502), rep("MESO",86), rep("OV",374), rep("PAAD",177), rep("PCPG",178), rep("PRAD",498), rep("READ",166), rep("SARC",259), rep("SKCM",103), rep("STAD",375), rep("TGCT",150), rep("THCA",502), rep("THYM",119), rep("UCEC",551), rep("UCS",56), rep("UVM",80))
combined_data$cancer_type <- cancer_type
save(combined_data, file = "combined_experssion_FPKM.RData")

