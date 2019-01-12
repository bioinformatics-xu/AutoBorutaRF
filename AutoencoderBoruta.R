library(Boruta)
library(caret)
library(randomForest)
library(ROCR)
library(dplyr)
library(h2o)
h2o.init()


######  training data prepare   ######
dataSource <- "GDSClung"
impRate <- 0.5
correlationRate <- 0.95
seeds <- 777

############   class data  ############
class.file <- paste(dataSource,"_class.csv",sep = "")
cellLine.class <- read.csv(class.file,row.names = 1,check=F)
cellLine.class <- cellLine.class[order(rownames(cellLine.class),decreasing = F),]

############## expression data ###########
rna.file <- paste(dataSource,"_rna.csv",sep = "")
rna.data <- read.csv(rna.file,row.names = 1,check=F)
colnames(rna.data) <- paste(colnames(rna.data),"rna",sep = "_")
rna.data <- rna.data[order(rownames(rna.data),decreasing = F),]

features <- colnames(rna.data)
train.data.h2o <- as.h2o(rna.data)
anomaly_model <- h2o.deeplearning(x = features, training_frame = train.data.h2o, 
                                  activation = "Tanh",   autoencoder = TRUE,   hidden = 100,   
                                  epochs = 100,variable_importances=TRUE,seed=seeds,standardize=TRUE)
rna.import.matrix <- anomaly_model@model$variable_importances
import.feature <- rna.import.matrix$variable[1:(nrow(rna.import.matrix)*impRate)]
rna.data <- rna.data[,which(colnames(rna.data) %in% import.feature)]

omic.name <- "rna"
variable.import <- anomaly_model@model$variable_importances
out.variable.file <- paste(dataSource,"_",omic.name,"_autoVarImport.txt",sep = "")
write.table(variable.import,file = out.variable.file,row.names = F,col.names = T,sep = "\t",quote = F)

# let's take the third hidden layer
train_features <- as.data.frame(h2o.deepfeatures(anomaly_model, train.data.h2o, layer = 1)) 
feature.file.name <- paste(dataSource,"_",omic.name,"_autoLayer1.txt",sep = "")
write.table(train_features,file = feature.file.name,row.names = F,col.names = T,sep = "\t",quote = F)



###############  cnv  data   ##############
cnv.file <- paste(dataSource,"_cnv.csv",sep = "")
cnv.data <- read.csv(cnv.file,row.names = 1,check=F)
colnames(cnv.data) <- paste(colnames(cnv.data),"cnv",sep = "_")
cnv.data <- cnv.data[order(rownames(cnv.data),decreasing = F),]
features <- colnames(cnv.data)
train.data.h2o <- as.h2o(cnv.data)
anomaly_model <- h2o.deeplearning(x = features,   training_frame = train.data.h2o, 
                                  activation = "Tanh",   autoencoder = TRUE,   hidden = 100,   
                                  epochs = 100,variable_importances=TRUE,seed=seeds,standardize=TRUE)
cnv.import.matrix <- anomaly_model@model$variable_importances
import.feature <- cnv.import.matrix$variable[1:(nrow(cnv.import.matrix)*impRate)]
cnv.data <- cnv.data[,which(colnames(cnv.data) %in% import.feature)]

omic.name <- "cnv"
variable.import <- anomaly_model@model$variable_importances
out.variable.file <- paste(dataSource,"_",omic.name,"_autoVarImport.txt",sep = "")
write.table(variable.import,file = out.variable.file,row.names = F,col.names = T,sep = "\t",quote = F)

# let's take the third hidden layer
train_features <- as.data.frame(h2o.deepfeatures(anomaly_model, train.data.h2o, layer = 1)) 
feature.file.name <- paste(dataSource,"_",omic.name,"_autoLayer1.txt",sep = "")
write.table(train_features,file = feature.file.name,row.names = F,col.names = T,sep = "\t",quote = F)


##############   mutation   data   ##################
mut.file <- paste(dataSource,"_mut.csv",sep = "")
mut.data <- read.csv(mut.file,row.names = 1,check=F)
colnames(mut.data) <- paste(colnames(mut.data),"mut",sep = "_")
mut.data <- mut.data[order(rownames(mut.data),decreasing = F),]

#############  correlation detection   ########
combine.feature <- cbind(rna.data,cnv.data)
cor.matrix <- cor(combine.feature)
highlyCorrelated <- findCorrelation(cor.matrix, cutoff=correlationRate)

if(length(highlyCorrelated) > 0){
  combine.feature <- combine.feature[,-highlyCorrelated]
}
combine.feature <- cbind(combine.feature,mut.data)
combine.feature <- combine.feature[order(rownames(combine.feature),decreasing = F),]

predictors.result <- list()

set.seed(seeds)

#####  Boruta method to detect features #####

for (d in 1:ncol(cellLine.class)) {
  print(d)
  combine.feature$class <- cellLine.class[,d]
  train.data <- combine.feature
  
  ######## check if the class is unbalanced ###########
  ensemble.class.flag <- (1/2 < length(which(train.data$class == 0))/length(which(train.data$class == 1))) & (2 > length(which(train.data$class == 0))/length(which(train.data$class == 1)))
  if(ensemble.class.flag){
    ########### not unbalanced #########
    train.data$class <- as.factor(as.character(train.data$class))
    set.seed(seeds)
    boruta.train <- Boruta(class~.,data=train.data,pValue=0.01,maxRuns = 200,mcAdj=FALSE)
    feature.boruta <- gsub("\`","",getSelectedAttributes(boruta.train, withTentative = T))
    predictors.result[[d]] <- feature.boruta
    t <- 1
    while ((length(feature.boruta) <= 3)&t<=20) {
      print("few predictors left, run again")
      boruta.train <- Boruta(class~.,data=train.data,pValue=0.01,maxRuns = 200,mcAdj=FALSE)
      feature.boruta <- gsub("\`","",getSelectedAttributes(boruta.train, withTentative = T))
      predictors.result[[d]] <- feature.boruta
      t <- t+1
    }
    
  }else{
    ######### unbalanced #######
    train.data$class <- as.factor(as.character(train.data$class))
    
    flag <- which.max(c(length(which(train.data$class == 0)),length(which(train.data$class == 1))))
    
    minority.class.data <- train.data[which(train.data[,ncol(train.data)] != (flag-1)),]
    majority.class.data <- train.data[which(train.data[,ncol(train.data)] == (flag-1)),]
    class.times <- floor(nrow(majority.class.data)/nrow(minority.class.data))
    set.seed(seeds)
    folds.ensemble <- createFolds(y=majority.class.data[,1],k = class.times)
    ensemble.i <- 1
    predictors.ensemble <- list()
    
    for (ensemble.i in 1:class.times) {
      print(ensemble.i)
      set.seed(seeds)
      ensemble.i.data <- rbind(minority.class.data,majority.class.data[folds.ensemble[[ensemble.i]],])
      boruta.train <- Boruta(class~.,data=ensemble.i.data,pValue=0.01,maxRuns = 200,mcAdj=FALSE)
      feature.boruta <- gsub("\`","",getSelectedAttributes(boruta.train, withTentative = T))
      predictors.ensemble[[ensemble.i]] <- feature.boruta
      t <- 1
      while ((length(feature.boruta) <= 1)&t<=20) {
        print("few predictors left, run again")
        boruta.train <- Boruta(class~.,data=ensemble.i.data,pValue=0.01,maxRuns = 200,mcAdj=FALSE)
        feature.boruta <- gsub("\`","",getSelectedAttributes(boruta.train, withTentative = T))
        predictors.ensemble[[ensemble.i]] <- feature.boruta
        t <- t+1
      }
    }
    
    predictors.result[[d]] <- unique(unlist(predictors.ensemble))
    
  }
  
}

######### take the predictors as a matrix #######
max.feature <- max(lengths(predictors.result))
predictor.result.matrix <- matrix(data = NA, nrow = length(predictors.result), ncol = max.feature)

for(i in 1:length(predictors.result)){
  result.i.length <- length(predictors.result[[i]])
  if(max.feature - result.i.length != 0){
    predictor.result.matrix[i,] <- c(unlist(predictors.result[[i]]),rep(NA,times=max.feature-result.i.length))
  }else{
    predictor.result.matrix[i,] <- unlist(predictors.result[[i]])
  }
}

############ output the predictor matrix ############
rownames(predictor.result.matrix) <- colnames(cellLine.class)
write.table(predictor.result.matrix,
            file = paste(dataSource,"_predictor_",impRate,"_",correlationRate,".csv",sep = ""),
            row.names=TRUE, na="",col.names=FALSE, sep=",")
