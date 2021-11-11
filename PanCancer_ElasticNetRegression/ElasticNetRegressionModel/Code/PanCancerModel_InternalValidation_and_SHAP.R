#####################################################################################################################
#####SHAP analysis based on ensemble elastic net regression model (where quiescence scores have been scaled by purity)
####################################################################################################################

#Load the required R packages:
library(ggpubr)
require('ggforce') # for `geom_sina`
library(fastshap)


#Load the train dataset:
setwd("~/Documents/GitHub/CancerDormancy/Data/PanCancer_ElasticNetRegression_Input/")
load("interal_traindataset.RData")
#Select only the features of interest:
setwd("~/Documents/GitHub/CancerDormancy/Data/PanCancer_ElasticNetRegression_Output/Summary/")
load("coefficients.RData")
features <- as.character(coefficients$features)
features <- c(features,"quiescence_score")
trainData <- trainData[,colnames(trainData) %in% features]
load("intercept.RData")


# Fit a projection pursuit regression model 
fit <- glm(quiescence_score ~ ., data = trainData)
new_coefficients <- data.frame(fit$coefficients)
#Change the intercept:
fit$coefficients[1] <- mean_intercept
new_coefficients$features <- rownames(new_coefficients)
new_coefficients <- new_coefficients[-1,]
feature_order <- new_coefficients$features
coefficients <- coefficients[match(feature_order, coefficients$features),]
#Check that everything is in the same order
all(new_coefficients$features == coefficients$features)
coefficients <- coefficients$coefficients
coefficients <- c(mean_intercept,coefficients)
#Change all of the coefficients of the model
fit$coefficients[1:61] <- coefficients[1:61]



#Load test dataset:
setwd("~/Documents/GitHub/CancerDormancy/Data/PanCancer_ElasticNetRegression_Input/")
load("interal_testdataset.RData")
# Predict on testData
predicted <- predict(fit, testData)
head(predicted)


#Plot the results:
results_plot <- data.frame(testData$quiescence_score, predicted)
results_plot$observed <- results_plot$testData.quiescence_score
setwd("~/Documents/GitHub/CancerDormancy/PanCancer_ElasticNetRegression/ElasticNetRegressionModel/Figures/")
pdf("LinearRegressionModelPerformance.pdf",height = 8.5,width = 4)
ggscatter(results_plot, x = "observed", y="predicted", cor.coef = TRUE, add = "reg.line") + ylab("Predicted Quiescence Score") + xlab("Observed Quiescence Score")
dev.off()


##Calculate SHAP scores
# Compute approximate Shapley values using 10 Monte Carlo simulations set.seed(101) # for reproducibility
shap_score <- explain(fit, X = subset(trainData, select = -quiescence_score), nsim = 10,
                      pred_wrapper = predict, adjust = TRUE)
shap_score
#Melt the dataset:
features <- colnames(shap_score)
i <- features[1]
features <- features[-1]
feature_data <- shap_score[[i]]
results <- data.frame(feature_data)
results$variable <- i
results$feature_data <- NULL
results$value <- feature_data
results$rfvalue <- 0
results$stdfvalue <- 0
feature_data <- abs(feature_data)
mean <- mean(feature_data)
results$mean_value <- mean
feature_value<- trainData[[i]]
max <- max(feature_value)
min <- min(feature_value)
feature_value <- scale(feature_value, center = min, scale = max - min)
results$feature_value <- feature_value
for (i in features) {
  print(i)
  
  feature_data <- shap_score[[i]]
  results_new <- data.frame(feature_data)
  results_new$variable <- i
  results_new$feature_data <- NULL
  results_new$value <- feature_data
  results_new$rfvalue <- 0
  results_new$stdfvalue <- 0
  feature_data <- abs(feature_data)
  mean <- mean(feature_data)
  results_new$mean_value <- mean
  feature_value<- trainData[[i]]
  max <- max(feature_value)
  min <- min(feature_value)
  feature_value <- scale(feature_value, center = min, scale = max - min)
  results_new$feature_value <- feature_value
  results <- rbind(results, results_new)
}
results$feature_value <- as.numeric(results$feature_value)
#order the results so that most significant features are at the top:
results <- results[order(-results$mean),]
results$variable <- as.character(results$variable)
ordered_features <- unique(results$variable)
results$variable <- as.factor(results$variable)
results$variable <- factor(results$variable, levels = rev(ordered_features))

## Plot shap overall metrics
plot1 <- ggplot(data = results)+
  coord_flip() + 
  # sina plot: 
  geom_sina(aes(x = variable, y = value, color = feature_value)) + scale_color_gradient(low="#FFCC33", high="#6600CC", 
                                                                                        breaks=c(0,1), labels=c("Low","High")) + labs(y = "SHAP value (impact on model output)", x = "", color = "Feature value") + theme_classic()
setwd("~/Documents/GitHub/CancerDormancy/PanCancer_ElasticNetRegression/ElasticNetRegressionModel/Figures/")
pdf("SHAP_plot_all_features.pdf",height = 10,width = 10)
plot1
dev.off()



#####Plot only noncancer type features:
results$variable <- as.character(results$variable)
removed.features <- c("CancerType.COAD","CancerType.CESC","CancerType.BLCA","CancerType.KIRP","CancerType.LIHC","CancerType.GBM","CancerType.BRCA","CancerType.OV","CancerType.READ","CancerType.SARC","CancerType.STAD","CancerType.KIRC","CancerType.LGG","CancerType.ACC","CancerType.UCS","CancerType.LUSC","CancerType.ESCA_adernocarcinoma","CancerType.CHOL")
filtered.results <- results[(!results$variable %in% removed.features),]
#order the results so that most significant features are at the top:
filtered.results <- filtered.results[order(-filtered.results$mean),]
filtered.results$variable <- as.character(filtered.results$variable)
ordered_features <- unique(filtered.results$variable)
filtered.results$variable <- as.factor(filtered.results$variable)
filtered.results$variable <- factor(filtered.results$variable, levels = rev(ordered_features))
plot1 <- ggplot(data = filtered.results)+
  coord_flip() + 
  # sina plot: 
  geom_sina(aes(x = variable, y = value, color = feature_value)) + scale_color_gradient(low="#FFCC33", high="#6600CC", 
                                                                                        breaks=c(0,1), labels=c("Low","High")) + labs(y = "SHAP value (impact on model output)", x = "", color = "Feature value") + theme_classic()
setwd("~/Documents/GitHub/CancerDormancy/PanCancer_ElasticNetRegression/ElasticNetRegressionModel/Figures/")
pdf("SHAP_plot_non_cancer_type_features.pdf",height = 10,width = 10)
plot1
dev.off()




