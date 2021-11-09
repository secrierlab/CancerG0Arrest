######################################
##TCGA Quiescence Type Classification
######################################

#Load required packages:
library(ggpubr)



#Load and combine the different quiescence scores 
setwd("~/Documents/GitHub/CancerDormancy/Data/TCGA_QuiescenceScores/")
load("TCGA_common_QS_purity_scaled.RData")
combined.scores <- z_score
combined.scores$common <- combined.scores$z_score
combined.scores$z_score <- NULL
load("TCGA_CDK_inhibition_QS_purity_scaled.RData")
combined.scores$cdk <- z_score$z_score
load("TCGA_contact_inhibition_QS_purity_scaled.RData")
combined.scores$contact <- z_score$z_score
load("TCGA_MEK_inhibition_QS_purity_scaled.RData")
combined.scores$mek <- z_score$z_score
load("TCGA_serum_starvation_QS_purity_scaled.RData")
combined.scores$serum <- z_score$z_score
load("TCGA_spotaneous_quiescence_QS_purity_scaled.RData")
combined.scores$spontaneous <- z_score$z_score

#Select only quiescent samples to classify 
combined.scores <- combined.scores[combined.scores$common > 0,]

#Save cancer type annotation:
CT <- combined.scores$CancerType
Barcodes <- rownames(combined.scores)
CT<- data.frame(CT,Barcodes)
combined.scores$CancerType <- NULL
combined.scores$common <- NULL

##Initially identify quiescence type with the highest quiescence score
#This will be the preliminary classification prior to determining uncertain samples
combined.scores$QuiescneceType <- colnames(combined.scores)[apply(combined.scores,1,which.max)]
table(combined.scores$QuiescneceType)

#######################################################################################################################################
#For each type of quiescecne prediction check if the corresponding score is significantly higher than that for other quiescence types 
#Also mark samples with the corresponding quiescence score < 1 as uncertain
#######################################################################################################################################


#CDK4/6 inhibition quiescence 

cdk <- combined.scores[combined.scores$QuiescneceType %in% "cdk",]
Uncertainty <- NULL
for (i in rownames(cdk)) {
  
  test.data <- combined.scores[rownames(combined.scores) %in% i,]
  mean <- test.data$cdk
  values <- c(test.data$contact, test.data$mek, test.data$serum, test.data$spontaneous)
  test <- t.test(values, mu = mean,
                 alternative = "less")
  pval <- test$p.value
  Uncertainty <- c(Uncertainty, pval)
  
}
Uncertainty <- data.frame(Uncertainty)
rownames(Uncertainty) <- rownames(cdk)
Uncertainty <- Uncertainty[Uncertainty$Uncertainty < 0.05,]
Uncertainty <- as.character(rownames(Uncertainty))
cdk <- combined.scores[combined.scores$QuiescneceType %in% "cdk",]
Uncertainty2 <- cdk[cdk$cdk < 1,]
Uncertainty2 <- as.character(rownames(Uncertainty2))
cdk$QuiescneceType <- sapply(rownames(cdk), function(x)
  ifelse(x %in% c(Uncertainty,Uncertainty2),"Uncertain","cdk"))
table(cdk$QuiescneceType)


#Contact inhibition
contact <- combined.scores[combined.scores$QuiescneceType %in% "contact",]
Uncertainty <- NULL
for (i in rownames(contact)) {
  
  test.data <- combined.scores[rownames(combined.scores) %in% i,]
  mean <- test.data$contact
  values <- c(test.data$cdk, test.data$mek, test.data$serum, test.data$spontaneous)
  test <- t.test(values, mu = mean,
                 alternative = "less")
  pval <- test$p.value
  Uncertainty <- c(Uncertainty, pval)
  
}
Uncertainty <- data.frame(Uncertainty)
rownames(Uncertainty) <- rownames(contact)
Uncertainty <- Uncertainty[Uncertainty$Uncertainty < 0.05,]
Uncertainty <- as.character(rownames(Uncertainty))
contact <- combined.scores[combined.scores$QuiescneceType %in% "contact",]
Uncertainty2 <- contact[contact$contact < 1,]
Uncertainty2 <- as.character(rownames(Uncertainty2))
contact$QuiescneceType <- sapply(rownames(contact), function(x)
  ifelse(x %in% c(Uncertainty,Uncertainty2),"Uncertain","contact"))
table(contact$QuiescneceType)


##MEK inhibition 
mek <- combined.scores[combined.scores$QuiescneceType %in% "mek",]
Uncertainty <- NULL
for (i in rownames(mek)) {
  
  test.data <- combined.scores[rownames(combined.scores) %in% i,]
  mean <- test.data$mek
  values <- c(test.data$cdk, test.data$contact, test.data$serum, test.data$spontaneous)
  test <- t.test(values, mu = mean,
                 alternative = "less")
  pval <- test$p.value
  Uncertainty <- c(Uncertainty, pval)
  
}
Uncertainty <- data.frame(Uncertainty)
rownames(Uncertainty) <- rownames(mek)
Uncertainty <- Uncertainty[Uncertainty$Uncertainty < 0.05,]
Uncertainty <- as.character(rownames(Uncertainty))
mek <- combined.scores[combined.scores$QuiescneceType %in% "mek",]
Uncertainty2 <- mek[mek$mek < 1,]
Uncertainty2 <- as.character(rownames(Uncertainty2))
mek$QuiescneceType <- sapply(rownames(mek), function(x)
  ifelse(x %in% c(Uncertainty,Uncertainty2),"Uncertain","mek"))
table(mek$QuiescneceType)


#Serum starvation
serum <- combined.scores[combined.scores$QuiescneceType %in% "serum",]
Uncertainty <- NULL
for (i in rownames(serum)) {
  
  test.data <- combined.scores[rownames(combined.scores) %in% i,]
  mean <- test.data$serum
  values <- c(test.data$cdk, test.data$contact, test.data$mek, test.data$spontaneous)
  test <- t.test(values, mu = mean,
                 alternative = "less")
  pval <- test$p.value
  Uncertainty <- c(Uncertainty, pval)
  
}
Uncertainty <- data.frame(Uncertainty)
rownames(Uncertainty) <- rownames(serum)
Uncertainty <- Uncertainty[Uncertainty$Uncertainty < 0.05,]
Uncertainty <- as.character(rownames(Uncertainty))
serum <- combined.scores[combined.scores$QuiescneceType %in% "serum",]
Uncertainty2 <- serum[serum$serum < 1,]
Uncertainty2 <- as.character(rownames(Uncertainty2))
serum$QuiescneceType <- sapply(rownames(serum), function(x)
  ifelse(x %in% c(Uncertainty,Uncertainty2),"Uncertain","serum"))
table(serum$QuiescneceType)


#Spontaneous quiescence
spontaneous <- combined.scores[combined.scores$QuiescneceType %in% "spontaneous",]
Uncertainty <- NULL
for (i in rownames(spontaneous)) {
  
  test.data <- combined.scores[rownames(combined.scores) %in% i,]
  mean <- test.data$spontaneous
  values <- c(test.data$cdk, test.data$contact, test.data$mek, test.data$serum)
  test <- t.test(values, mu = mean,
                 alternative = "less")
  pval <- test$p.value
  Uncertainty <- c(Uncertainty, pval)
  
}
Uncertainty <- data.frame(Uncertainty)
rownames(Uncertainty) <- rownames(spontaneous)
Uncertainty <- Uncertainty[Uncertainty$Uncertainty < 0.05,]
Uncertainty <- as.character(rownames(Uncertainty))
spontaneous <- combined.scores[combined.scores$QuiescneceType %in% "spontaneous",]
Uncertainty2 <- spontaneous[spontaneous$spontaneous < 1,]
Uncertainty2 <- as.character(rownames(Uncertainty2))
spontaneous$QuiescneceType <- sapply(rownames(spontaneous), function(x)
  ifelse(x %in% c(Uncertainty,Uncertainty2),"Uncertain","spontaneous"))
table(spontaneous$QuiescneceType)


#Combine the results
combined.scores <- rbind(contact, cdk, mek, serum, spontaneous)


#######Summarise the information (for each cancer type)
combined.scores$Barcodes <- rownames(combined.scores)
combined.scores <- merge(combined.scores, CT, 
                         by.x = "Barcodes", by.y = "Barcodes")
table(combined.scores$QuiescneceType)
combined.scores$CancerType <- combined.scores$CT
CT <- unique(as.character(combined.scores$CancerType))
a <- "GBM"
selected.cancer <- combined.scores[combined.scores$CancerType %in% a,]
QuiescenceTypes <- c("cdk","contact","mek","serum","spontaneous","Uncertain")
N <- NULL
for (i in QuiescenceTypes) {
  print(i)
  selected.type <- selected.cancer[selected.cancer$QuiescneceType %in% i,]
  n <- dim(selected.type)[1]
  N <- c(N,n)
  
  
}
CancerType <- rep(a, 6)
Summary <- data.frame(CancerType, QuiescenceTypes, N)

CT <- CT[-1]
for (a in CT) {
  print(a)
  selected.cancer <- combined.scores[combined.scores$CancerType %in% a,]
  QuiescenceTypes <- c("cdk","contact","mek","serum","spontaneous","Uncertain")
  N <- NULL
  for (i in QuiescenceTypes) {
    print(i)
    selected.type <- selected.cancer[selected.cancer$QuiescneceType %in% i,]
    n <- dim(selected.type)[1]
    N <- c(N,n)
    
    
  }
  CancerType <- rep(a, 6)
  Summary_new <- data.frame(CancerType, QuiescenceTypes, N)
  Summary <- rbind(Summary, Summary_new)
  
}
Summary$QuiescenceTypes <- factor(Summary$QuiescenceTypes, levels = c("Uncertain","cdk","contact","mek","serum","spontaneous"))




###Order cancer types according to % uncertainty 
CT <- unique(as.character(Summary$CancerType))
Uncertain <- NULL
Percentage <- NULL
for (a in CT) {
  print(a)
  test <- Summary[Summary$CancerType %in% a,]
  total <- sum(test$N)
  test <- test[test$QuiescenceTypes %in% "Uncertain",]
  Uncertain <- test$N
  percentage <- Uncertain / total
  print(percentage)
  Percentage <- c(Percentage, percentage)
  
}
Percentages <- data.frame(CT, Percentage)
Percentages <- Percentages[order(Percentages$Percentage),]
CT <- as.character(Percentages$CT)
Summary$CancerType <- factor(Summary$CancerType,levels = CT)



#Plot the classification summary:
setwd("~/Documents/GitHub/CancerDormancy/QuiescenceTypeClassification/TCGA_classification/Figures")
pdf("TCGA_QuiescenceTypeClassification.pdf", height = 5, width = 10)
p <- ggplot(Summary, aes(fill=QuiescenceTypes, y=N, x=CancerType, width = 0.7)) + 
  geom_bar(stat = "identity", position = "fill") + theme_classic()
p + rotate_x_text(45) + scale_fill_manual(values = c("#666666", "#D95F02", "#7570B3", "#E7298A", "#E6AB02", "#1B9E77"))
dev.off()

