##########################################
#TCGA expression and CNV data download
##########################################


#Load required packages:
library(TCGAbiolinks)
library(SummarizedExperiment)
library(RTCGAToolbox)
#This code allows you to download Harmonized TCGA Hg38 data

#Set the directory to which you want to download the data (a separate folder will be created for each cancer/study type)
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")

#Type in here the list of cancer types/studies for which you want to download the data. The for loop will subseqeuntly download data for each study.
CT <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

current_dir<-getwd()

for (i in CT) {
  
  print(i)
  
  dir.create(i) 
  
  setwd(paste(current_dir,'/',i,sep=''))
  
  cancer_type<-paste("TCGA-",i,sep="")
  
  ##
  ## Copy Number Variation
  ##
  
  #discrete copy number variation
  query1h<- GDCquery(project = cancer_type,data.category="Copy Number Variation",data.type ="Gene Level Copy Number Scores",legacy = F )
  GDCdownload(query1h)
  dataCNV <- as.data.frame(GDCprepare(query1h, remove.files.prepared=T,save=T))
  write.table(file=paste(paste(i,"CNV_gene_level_copy.harmonized",sep='_'),'.txt',sep=''),dataCNV,sep='\t',row.names=T,quote=F)
  
  #segments copy number variation
  query2h<- GDCquery(project = cancer_type,data.category="Copy Number Variation",data.type ="Copy Number Segment",legacy = F)
  GDCdownload(query2h)
  dataCNV2 <- as.data.frame(GDCprepare(query2h, remove.files.prepared=T,save=T))
  write.table(file=paste(paste(i,"CNV_segments.harmonized",sep='_'),'.txt',sep=''),dataCNV2,sep='\t',row.names=T,quote=F)
  
  ##
  ## Transcriptome Profiling
  ##
  
  #FKPM
  query4h<- GDCquery(project = cancer_type,
                     data.category="Transcriptome Profiling",
                     data.type ="Gene Expression Quantification",
                     workflow.type = "HTSeq - FPKM",legacy = FALSE
  )
  
  GDCdownload(query4h)
  data4 <-GDCprepare(query4h, remove.files.prepared=T,save= T)
  matrixRNAseq<-assay(data4)
  
  write.table(file=paste(paste(i,"RNAseqdata_HTSEQ_FKPM.harmonized",sep='_'),'.txt',sep=''),matrixRNAseq,sep='\t',row.names=T,quote=F)
  
  #COUNTS
  query5h<- GDCquery(project = cancer_type,
                     data.category="Transcriptome Profiling",
                     data.type ="Gene Expression Quantification",
                     workflow.type = "HTSeq - Counts",legacy = FALSE
  )
  
  GDCdownload(query5h)
  data5 <-GDCprepare(query5h, remove.files.prepared=T,save=T)
  matrixRNAseq2<-assay(data5)
  
  write.table(file=paste(paste(i,"RNAseqdata_HTSEQ_COUNTS.harmonized",sep='_'),'.txt',sep=''),matrixRNAseq2,sep='\t',row.names=T,quote=F)
  
  clinical<-GDCquery_clinic(project = cancer_type,type="Clinical")
  
  write.table(file=paste(paste(i,"Clinical_harmonized",sep='_'),'.txt',sep=''),clinical,sep='\t',row.names=T,quote=F)
  
  
  setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_Expr_CNV_Data/")
  
  system('rm *.RData')
  
}
