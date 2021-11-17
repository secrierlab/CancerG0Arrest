#################################################################################################
##Pan-cancer elastic net regression model to identify genomic features associated with quiescence
#################################################################################################

#Load required R packages:
library(caret)
library("glmnet")



#Load cleaned input data:
setwd("~/Documents/GitHub/CancerDormancy/Data/PanCancer_ElasticNetRegression_Input")
load("input_data_cleaned.RData")


####Split the data into interal test and train datasets
set.seed(123)
trainRowNumbers <- createDataPartition(Merged_data$quiescence_score, p=0.8, list=FALSE)
trainData <- Merged_data[trainRowNumbers,]
testData <- Merged_data[-trainRowNumbers,]
#Save these datasets:
save(trainData, file = "interal_traindataset.RData")
save(testData, file = "interal_testdataset.RData")


###################################################################
#Using the traindataset build a linear regression model 1000 times
##################################################################

#Load the train data:
setwd("~/Documents/GitHub/CancerDormancy/Data/PanCancer_ElasticNetRegression_Input")
load("interal_traindataset.RData")
#Rename this data:
Merged_data <- trainData
Merged_data$SampleID <- NULL


for (i in 1:1000) {
  
  #Print the iteration number
  print(i)
  
  ###Split the data into training and test datasets:
  # Step 1: Get row numbers for the training data
  trainRowNumbers <- createDataPartition(Merged_data$quiescence_score, p=0.9, list=FALSE)
  # Step 2: Create the training  dataset
  trainData <- Merged_data[trainRowNumbers,]
  # Step 3: Create the test dataset
  testData <- Merged_data[-trainRowNumbers,]
  # Store X and Y for later use.
  x = trainData
  x$quiescence_score <- NULL
  y = trainData$quiescence_score
  ########Build an elastic net regression model:
  grid <- expand.grid(.alpha = seq(0,1,by = 0.2),
                      .lambda = seq(0.01, 0.03, by = 0.002))
  control <- caret::trainControl(method = "cv", number = 5)
  model1 <- caret::train(quiescence_score ~., data = trainData, method = "glmnet", trControl = control, tuneGrid = grid)
  fitted <- predict(model1)
  model1$bestTune
  #Coefficients of the final model. You need to sepcify the best lambda:
  coefficients <- coef(model1$finalModel, model1$bestTune$lambda)
  coefficients <- as.matrix(coefficients)
  coefficients <- data.frame(coefficients)
  coefficients$variable <- rownames(coefficients)
  coefficients <- coefficients[coefficients$X1 > 0 | coefficients$X1 < 0,]
  #Save the coefficients for each iteration 
  setwd("~/Documents/GitHub/CancerDormancy/Data/PanCancer_ElasticNetRegression_Output/IterationData/")
  save(coefficients, file = paste("model_coefficients_interation",i,".RData",sep = ""))
  write.csv(coefficients, file = paste("model_coefficients_interation",i,".csv",sep = ""))
  ########PREPARE the test data set and predict:
  # Predict on testData
  predicted <- predict(model1, testData)
  #Model performance metrics:
  performance <- data.frame(
    RMSE = RMSE(predicted, testData$quiescence_score),
    Rsquare = R2(predicted, testData$quiescence_score)
  )
  #save model performance information:
  setwd("~/Documents/GitHub/CancerDormancy/Data/PanCancer_ElasticNetRegression_Output/IterationData/")
  save(performance, file = paste("model_performance_interation_",i,".RData",sep = ""))
}





################################################################################################################
#######################Check whith features occur in all iterations:
################################################################################################################

#Load the train data:
setwd("~/Documents/GitHub/CancerDormancy/Data/PanCancer_ElasticNetRegression_Input/")
load("interal_traindataset.RData")
#Rename this data:
Merged_data <- trainData
Merged_data$SampleID <- NULL
Merged_data$quiescence_score <- NULL
#Record all of the genomic features used to build the model
Features <- colnames(Merged_data)
Features <- data.frame(Features)
###Load and combine results from each iteration
for (i in 1:1000) {
  print(i)
  setwd("~/Documents/GitHub/CancerDormancy/Data/PanCancer_ElasticNetRegression_Output/IterationData/")
  results <- read.csv(paste("model_coefficients_interation",i,".csv",sep = ""))
  results <- results[results$X1 > 0 | results$X1 < 0, ]
  #Merge with Features dataframe 
  Features <- merge(Features, results,
                    by.x = "Features", by.y = "variable")
  
}
####60 features in total
#Save the results
setwd("~/Documents/GitHub/CancerDormancy/Data/PanCancer_ElasticNetRegression_Output/Summary/")
Features <- as.character(Features$Features)
fileConn<-file("results_overview.txt")
writeLines(Features, fileConn)
close(fileConn)



################################################################################################################
#######################Find out the average coefficients for each of the variables:
################################################################################################################
setwd("~/Documents/GitHub/CancerDormancy/Data/PanCancer_ElasticNetRegression_Output/Summary/")
#Load the features:
features <- read.table("results_overview.txt", sep = "\t", header = FALSE)
features <- as.character(features$V1)
coefficients <- NULL
setwd("~/Documents/GitHub/CancerDormancy/Data/PanCancer_ElasticNetRegression_Output/IterationData/")
for (a in features) {
  
  print(a)
  mean_coefficients <- NULL
  
  for (i in 1:1000) {
    
    results <- read.csv(paste("model_coefficients_interation",i,".csv",sep = ""))
    results <- results[results$variable %in% a,]
    coeff <- as.numeric(results$X1)
    mean_coefficients <- c(mean_coefficients, coeff)
    
  }
  
  mean_coefficients <- mean(mean_coefficients)
  coefficients <- c(coefficients, mean_coefficients)
}

#Save the result:
coefficients <- data.frame(features, coefficients)
setwd("~/Documents/GitHub/CancerDormancy/Data/PanCancer_ElasticNetRegression_Output/Summary/")
save(coefficients, file = "coefficients.RData")





################################################################################################################
#######################Find out the mean intercept value:
#################################################################################################################

setwd("~/Documents/GitHub/CancerDormancy/Data/PanCancer_ElasticNetRegression_Output/IterationData/")

mean_intercept <- NULL
for (i in 1:1000) {
  
  results <- read.csv(paste("model_coefficients_interation",i,".csv",sep = ""))
  results <- results[results$variable %in% "(Intercept)",]
  coeff <- as.numeric(results$X1)
  mean_intercept <- c(mean_intercept, coeff)
  
}

mean_intercept <- mean(mean_intercept)
setwd("~/Documents/GitHub/CancerDormancy/Data/PanCancer_ElasticNetRegression_Output/Summary/")
save(mean_intercept, file = "intercept.RData")












