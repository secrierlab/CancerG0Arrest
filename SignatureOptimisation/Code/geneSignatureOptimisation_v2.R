#### DEPENDENCIES ####
library(randomForest)
library(GSVA)
library(dplyr)
library(ggplot2)

#### DATA LOADING ####
# combined rnaseq exprsesion data and quiescence annotations, for 3 separate laboratory lines
setwd("~/Documents/GitHub/CancerDormancy/Data/SignatureOptimisation")
load("S1.Rdata")
load("S2.Rdata")
load("S3.Rdata")

# list of known DEGs associated with quiescence
load("downregulated_common.RData")
load("upregulated_common.RData")
known_genes <- c(downregulated_common, upregulated_common)

#Check which genes are reported in the 3 datasets:
S1_expr <- S1_combined$expr_data
S1_genes <- rownames(S1_expr)
S2_expr <- S2_combined$expr_data
S2_genes <- rownames(S2_expr)
S3_expr <- S3_combined$expr_data
S3_genes <- rownames(S3_expr)
known_genes <- known_genes[known_genes %in% S1_genes]
known_genes <- known_genes[known_genes %in% S2_genes]
known_genes <- known_genes[known_genes %in% S3_genes]
downregulated_common <- downregulated_common[downregulated_common %in% known_genes]
upregulated_common <- upregulated_common[upregulated_common %in% known_genes]

#### FUNCTIONS #####
# function to calculate a vector that is normalized from 0 to 1
normalize_vector <- function(input_vec){
  maxval <- max(input_vec, na.rm = TRUE)
  return_vec <- input_vec/maxval
  return(return_vec)
}

#format input data for better random forest model compatibility
dataModelPrep <- function(expr_ann_data){
  expr_data <- as.data.frame(expr_ann_data[1])
  ann_data <- as.data.frame(expr_ann_data[2])
  PQ <- rep(NA, ncol(expr_data))
  FormDat <- rbind(PQ, expr_data); rownames(FormDat)[1] <- "Q"
  # replace text based description of quiescence state with true/false version
  for (i in 1:ncol(expr_data)){
    if (ann_data[i,] == "Proliferating"){
      FormDat[1,i] <- FALSE
    }
    else if (ann_data[i,] == "Quiescent"){
      FormDat[1,i] <- TRUE
    }
  }
  MDat <- as.data.frame(t(FormDat))
  MDat$Q <- as.logical(MDat$Q)
  return(MDat)
}

# function to calculate the relative importance of genes using their gini coefficient in the rf model
giniIndexCalc <- function(expr_ann_data,nIter=1000){
  # preset seeds
  seeds <- seq(1, 1+nIter-1, by = 1)
  # format data to work with random forest
  rf_data <- dataModelPrep(expr_ann_data)
  # filter for the 120 genes
  rf_data <- rf_data[ ,c(1, which(colnames(rf_data) %in% c(downregulated_common, upregulated_common))) ]
  rf_data$Q <- as.factor(rf_data$Q)
  colnames(rf_data) <- make.names(colnames(rf_data), unique = TRUE)
  # set up dataframe to average 5 results
  giniTable <- data.frame(matrix(nrow = 120, ncol = nIter)); rownames(giniTable) <- known_genes
  # do 5 iterations to get a larger spread of results
  for (i in 1:nIter){
    # run random forest model, and extract importance of each variable (gene)
    set.seed(seeds[i])
    model_rf <- randomForest(Q ~ ., data=rf_data, proximity=TRUE, ntree=500)
    gini_values <- importance(model_rf)
    # write importance(gini index) back into the dataframe
    for (gene in rownames(gini_values)){
      giniTable[match(gene, known_genes), i] <- gini_values[gene,]
    }
    cat(i); cat("\n")
  }
  # calculate and normalize mean values
  meanGini <- normalize_vector(rowMeans(giniTable, na.rm = TRUE))
  return(as.numeric(meanGini))
}

# function gradually removes genes from the z-score model, based on the gini index
gene_score_filter <- function(val_data, rank_dataset){
  # set up matrix to hold results
  results_df <- as.data.frame(matrix(nrow = nrow(rank_dataset)-2, ncol = 8))
  colnames(results_df) <- c( "corr_EuD", "corr_pRB", "p_EuD", "p_pRB", "avg_corr", "avg_p" , "gene_number", "gene_names")
  for (i in 1:(nrow(rank_dataset)-2)){
    gene_num <- i+2
    # select the top n genes, and filter the validation dataset to only indclude these
    filtered_df <- rank_dataset[rownames(rank_dataset) %in% rownames(slice_max(GeneRank_Gini_full, order_by = max, n = gene_num)), ]
    included_genes <- rownames(filtered_df)
    upreg_filtered <- intersect(included_genes, upregulated_common)
    downreg_filtered <- intersect(included_genes, downregulated_common)
    gene_list_filtered <- list(upreg_filtered, downreg_filtered)
    # calculate quiescence score using z-score approach with data from only the genes selected above
    z_score <- gsva(refined_expr, gene_list_filtered, method = "zscore", verbose = FALSE, parallel.sz=8)
    z_score <- t(z_score); z_score <- data.frame(z_score)
    if (ncol(z_score) == 1){
      val_data$z_score <- z_score$z_score
    } else {
      z_score$z_score <- z_score$X1 - z_score$X2
      val_data$z_score <- z_score$X1 - z_score$X2
    }
    # calculate correlations, p-values
    m2 <- as.numeric(cor.test(val_data$Percentage_G0_EuD, val_data$z_score)[c(3,4)])
    m3 <- as.numeric(cor.test(val_data$Percentage_G0_phosphoRB, val_data$z_score)[c(3,4)])
    # collect results in a vector
    corr_vec <- c(m2[2], m3[2], m2[1], m3[1], mean(c(m2[2],m3[2])) , mean(c(m2[1],m3[1])), gene_num)
      # use absolute values to account for negative correlations, and round values to 4 digits
      corr_vec <- abs(corr_vec); corr_vec <- round(corr_vec, digits = 4)
      # add the names of all genes present in the combination
      corr_vec <- c(corr_vec, toString(included_genes))
    results_df[i,] <- corr_vec
  }
  for (i in 1:7){
    results_df[,i] <- as.numeric(results_df[,i])
  }
  return(results_df)
}

#### DATA LOADING ####
# combined rnaseq exprsesion data and quiescence annotations, for 3 separate laboratory lines
setwd("~/Documents/GitHub/CancerDormancy/Data/SignatureOptimisation")
load("S1.Rdata")
load("S2.Rdata")
load("S3.Rdata")


# dataset used for validation
  # gene expression values
load("expr_data_deENSMBLd.Rdata")
  # annotated results from laboratory scoring methods
load("annotation.RData")

#### GENE COMBINATION SCORING ASSESMENT ####
# calculate how informative a given gene is in a given dataset, and summarize results into dataframe
GeneRank_Gini_full <- data.frame(matrix(nrow = 120, ncol = 3))
rownames(GeneRank_Gini_full) <- known_genes; colnames(GeneRank_Gini_full) <- c("dataset_1","dataset_2","dataset_3")
GeneRank_Gini_full$dataset_1 <- giniIndexCalc(S1_combined,nIter=1000)
GeneRank_Gini_full$dataset_2 <- giniIndexCalc(S2_combined,nIter=1000)
GeneRank_Gini_full$dataset_3 <- giniIndexCalc(S3_combined,nIter=1000)
GeneRank_Gini_full<- na.omit(GeneRank_Gini_full)
GeneRank_Gini_full[, "max"] <- apply(GeneRank_Gini_full[, 1:3], 1, max)

setwd("~/Documents/GitHub/CancerDormancy/SignatureOptimisation/Results")
save(GeneRank_Gini_full, file = "GeneRank_Gini_full_v2.RData")
# calculate average correlation and p values at different thresholds (will take approx. 1-2 minutes)
geneCombinationScore <- gene_score_filter(Validation_dataset, GeneRank_Gini_full)
save(geneCombinationScore, file = "geneCombinationScore_v2.RData")



# plot of correlation change
corplot <- ggplot(geneCombinationScore, aes(x = gene_number, y = avg_corr)) +
  geom_line(aes(color = "average"), lwd = 1.5, alpha = 0.65) +
  geom_point(aes(y = corr_EuD, color = "eud")) +
  geom_point(aes(y = corr_pRB, color = "prb")) +
  ylim(0.6, 0.85) +
  labs(x = "gene_number",
       y = "correlation",
       color = "Legend") +
  scale_color_manual(name = "Method",
                     values = c(
                       "eud" = "royalblue3",
                       "prb" = "green4",
                       "average" = "black")) +
  theme(legend.position = "bottom") + geom_vline(xintercept=35,linetype="dashed", color = "red") + geom_vline(xintercept=32,linetype="dashed", color = "red")
corplot


pval_plot <- ggplot(geneCombinationScore, aes(x = gene_number, y = avg_p)) +
  geom_line(aes(color = "average"), lwd = 1.5, alpha = 0.65) +
  geom_point(aes(y = p_EuD, color = "eud")) +
  geom_point(aes(y = p_pRB, color = "prb")) +
  ylim(0, 0.08) +
  labs(x = "gene_number",
       y = "p_value",
       color = "Legend") +
  scale_color_manual(name = "Method",
                     values = c(
                       "eud" = "royalblue3",
                       "prb" = "green4",
                       "average" = "black")) +
  theme(legend.position = "bottom") + geom_hline(yintercept=0.05,linetype="dashed", color = "black")+ geom_vline(xintercept=35,linetype="dashed", color = "red") + geom_vline(xintercept=32,linetype="dashed", color = "red")
pval_plot

GeneRank_Gini_full <- GeneRank_Gini_full[order(-GeneRank_Gini_full$max),]
GeneRank_Gini_full$Number_of_genes <- 1:120
gini_threshold <- ggplot(GeneRank_Gini_full, aes(x = Number_of_genes, y = max)) +
  geom_line(aes(color = "average"), lwd = 1.5, alpha = 0.65) +
  labs(x = "gene_number",
       y = "Gini Index Threshold",
       color = "Legend") +
  scale_color_manual(name = "Method",
                     values = c(
                       "eud" = "royalblue3",
                       "prb" = "green4",
                       "average" = "black")) +
  theme(legend.position = "bottom") + geom_vline(xintercept=35,linetype="dashed", color = "red") + geom_vline(xintercept=32,linetype="dashed", color = "red")
gini_threshold

setwd("~/Documents/GitHub/CancerDormancy/SignatureOptimisation/Figures")
pdf("SignatureOptimisationSummary_v2.pdf",height = 8,width = 5)
ggarrange(corplot, pval_plot, gini_threshold ,
          ncol = 1, nrow = 3)
dev.off()
