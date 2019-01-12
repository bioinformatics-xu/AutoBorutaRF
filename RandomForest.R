#  drug response prediction classifier training and testing ##

######  Classifier  #####

rm(list=ls())

library(caret)
library(randomForest)
library(ROCR)
library(dplyr)


######  train data prepare   ######
dataSource <- "GDSClung"
impRate <- 0.5
correlationRate <- 0.95
seeds <- 777
fold.num <- 10


############   class data  ############
class.file <- paste(dataSource,"_class.csv",sep = "")
cellLine.class <- read.csv(class.file,row.names = 1,check=F)
cellLine.class <- cellLine.class[order(rownames(cellLine.class),decreasing = F),]

############## expression data preprocess ###########
rna.file <- paste(dataSource,"_rna.csv",sep = "")
rna.data <- read.csv(rna.file,row.names = 1,check=F)
colnames(rna.data) <- paste(colnames(rna.data),"rna",sep = "_")
rna.data <- rna.data[order(rownames(rna.data),decreasing = F),]

###############  cnv  data   ##############
cnv.file <- paste(dataSource,"_cnv.csv",sep = "")
cnv.data <- read.csv(cnv.file,row.names = 1,check=F)
colnames(cnv.data) <- paste(colnames(cnv.data),"cnv",sep = "_")
cnv.data <- cnv.data[order(rownames(cnv.data),decreasing = F),]

##############   mutation   data   ##################
mut.file <- paste(dataSource,"_mut.csv",sep = "")
mut.data <- read.csv(mut.file,row.names = 1,check=F)
colnames(mut.data) <- paste(colnames(mut.data),"mut",sep = "_")
mut.data <- mut.data[order(rownames(mut.data),decreasing = F),]

combine.feature <- cbind(rna.data,cnv.data,mut.data)
combine.feature <- combine.feature[order(rownames(combine.feature),decreasing = F),]


#  read the feature matrix selected by AutoencoderBoruta  
predictor.result.matrix <- as.data.frame(read.csv(paste(dataSource,"_predictor_",impRate,"_",correlationRate,".csv",sep = ""),header = F,row.names = 1,check=F))

all.feature <- unique(unlist(predictor.result.matrix))
all.feature <- all.feature[-which(all.feature == "")]

result.table <- as.data.frame(matrix(data = 0, nrow = ncol(cellLine.class), ncol = 21))
confusion.matrix <- as.data.frame(matrix(data = 0, nrow = ncol(cellLine.class), ncol = 5))
colnames(confusion.matrix) <- c("drugName","TP","FP","TN","FN")

set.seed(seeds)


########  start training and testing    ########
for (d in 1:ncol(cellLine.class)) {
    print(d)
    predictors <- predictor.result.matrix[d,]
    predictors.filter <- unlist(predictors[-which((predictors == ""))])
    
    if(is.null(predictors.filter)){
        input.data <- combine.feature[,which(colnames(combine.feature) %in% all.feature)]
    }else{
        input.data <- combine.feature[,which(colnames(combine.feature) %in% predictors.filter)]
    }
    
    if(is.numeric(input.data)){
        input.data <- cbind(input.data,as.factor(as.character(cellLine.class[,d])))
        colnames(input.data) <- c(predictors.filter,"class")
        input.data <- as.data.frame(input.data)
        input.data$class <- as.factor(as.character(cellLine.class[,d]))
    }else{
        input.data$class <- as.factor(as.character(cellLine.class[,d]))
    }
    
    ensemble.class.flag <- (1/2 < length(which(input.data$class == 0))/length(which(input.data$class == 1))) & (2 > length(which(input.data$class == 0))/length(which(input.data$class == 1)))
    if(ensemble.class.flag){
        
        set.seed(seeds)
        folds<-createFolds(y=input.data$class,k = fold.num) 
        result.overall.cv <- as.data.frame(matrix(data = NA, nrow = fold.num, ncol = 7))
        result.classes.cv <- as.data.frame(matrix(data = NA, nrow = fold.num, ncol = 11))
        confusion.matrix.cv <- as.data.frame(matrix(data = NA, nrow = fold.num, ncol = 4))
        auc.cv <- numeric()
        mcc.cv <- numeric()
        
        colnames(input.data) <- c(paste("V",c(1:(ncol(input.data)-1)),sep = ""),"class")
        
        for(fold.i in 1:fold.num){
            
            traindata <- input.data[-folds[[fold.i]],]
            testdata <- input.data[folds[[fold.i]],]
            testClass <- input.data[folds[[fold.i]],ncol(input.data)]
            set.seed(seeds)
            rf.reslut = randomForest(class ~ .,data = traindata,mtry=1, importance = F,ntree = 1000)
            
            pred = predict(rf.reslut,testdata)
            prob <- predict(rf.reslut,testdata,type = "prob")
            
            pred.confusion.input <- factor(as.vector(pred),levels=c("0","1"))
            testClass.confusion.input <- factor(as.vector(testClass),levels=c("0","1"))
            results <- confusionMatrix(pred.confusion.input,testClass.confusion.input, positive = "1")
            
            result.overall.cv[fold.i,] <- t(as.matrix(results,what = "overall"))
            result.classes.cv[fold.i,] <- t(as.matrix(results,what = "classes"))
            confusion.matrix.cv[fold.i,] <- c(results$table[4],results$table[2],results$table[1],results$table[3])
            
            mcc.numerator.cv <- (confusion.matrix.cv[fold.i,1]*confusion.matrix.cv[fold.i,3]-confusion.matrix.cv[fold.i,2]*confusion.matrix.cv[fold.i,4])
            mcc.denominator.cv <- sqrt((confusion.matrix.cv[fold.i,1]+confusion.matrix.cv[fold.i,2])*(confusion.matrix.cv[fold.i,1]+confusion.matrix.cv[fold.i,4])*
                                           (confusion.matrix.cv[fold.i,3]+confusion.matrix.cv[fold.i,2])*(confusion.matrix.cv[fold.i,3]+confusion.matrix.cv[fold.i,4]))
            mcc.cv[fold.i] <- mcc.numerator.cv/mcc.denominator.cv
            
            if(length(unique(testClass)) == 1){
                auc.cv <- auc.cv
            }else{
                pred.roc <- prediction(prob[,2],as.vector(as.character(testClass)))
                auc.cv[fold.i] <-unlist(performance(pred.roc, "auc")@y.values)
            }
            
        }
        
        result.table[d,2:ncol(result.table)] <- c(mean(auc.cv,na.rm=TRUE),mean(mcc.cv,na.rm = TRUE),apply(result.overall.cv,2,mean,na.rm=TRUE),apply(result.classes.cv,2,mean,na.rm=TRUE)) 
        confusion.matrix[d,2:ncol(confusion.matrix)] <- apply(confusion.matrix.cv,2,mean,na.rm=TRUE)
        
        colnames(result.table) <- c("drugName","AUC","MCC",rownames(as.matrix(results,what = "overall")),rownames(as.matrix(results,what = "classes")))
        result.table$drugName <- colnames(cellLine.class)
        confusion.matrix$drugName <- colnames(cellLine.class)
        
    }else{
        
        set.seed(seeds)
        folds <- createFolds(y = input.data$class, k=fold.num)
        auc.cv <- numeric()
        mcc.cv <- numeric()
        max.flag <- which.max(c(length(which(input.data$class == 0)),length(which(input.data$class == 1))))
        result.overall.cv <- as.data.frame(matrix(data = NA, nrow = fold.num, ncol = 7))
        result.classes.cv <- as.data.frame(matrix(data = NA, nrow = fold.num, ncol = 11))
        confusion.matrix.cv <- as.data.frame(matrix(data = NA, nrow = fold.num, ncol = 4))
        
        colnames(input.data) <- c(paste("V",c(1:(ncol(input.data)-1)),sep = ""),"class")
        for (fold.i in 1:fold.num) {
            
            traindata <- input.data[-folds[[fold.i]],]
            testdata <- input.data[folds[[fold.i]],]
            testClass <- input.data[folds[[fold.i]],ncol(input.data)]
            
            minority.class.data <- traindata[which(traindata[,ncol(traindata)] != (max.flag-1)),]
            majority.class.data <- traindata[which(traindata[,ncol(traindata)] == (max.flag-1)),]
            class.times <- floor(nrow(majority.class.data)/nrow(minority.class.data))
            set.seed(seeds)
            folds.ensemble <- createFolds(y=majority.class.data[,1],k = class.times)
            prob.ensemble <- matrix(data = 0,nrow = nrow(testdata),ncol = 2)
            
            for (ensemble.i in 1:class.times) {
                print(ensemble.i)
                ensemble.data.i <- rbind(minority.class.data,majority.class.data[folds.ensemble[[ensemble.i]],])
                set.seed(seeds)
                ensemble.classifier.i = randomForest(class ~ .,
                                                     data = ensemble.data.i,
                                                     importance = F,
                                                     mtry=1,
                                                     ntree=1000)
                prob.i <- predict(ensemble.classifier.i,testdata,type = "prob")
                prob.ensemble <- prob.ensemble + prob.i
                
            }
            
            prob.ensemble.final <- prob.ensemble
            prob.ensemble.final[,2] <- prob.ensemble.final[,2]/(prob.ensemble.final[,1]+prob.ensemble.final[,2])
            prob.ensemble.final[,1] <- 1 - prob.ensemble.final[,2] 
            pred.class <- vector(length = nrow(testdata))
            pred.class[which(prob.ensemble.final[,2] >= 0.5)] <- 1
            pred.class[which(prob.ensemble.final[,2] < 0.5)] <- 0
            
            pred.confusion.input <- factor(as.vector(as.character(pred.class)),levels=c("0","1"))
            testClass.confusion.input <- factor(as.vector(as.character(testClass)),levels=c("0","1"))
            results <- confusionMatrix(pred.confusion.input,testClass.confusion.input, positive = "1")
            
            result.overall.cv[fold.i,] <- t(as.matrix(results,what = "overall"))
            result.classes.cv[fold.i,] <- t(as.matrix(results,what = "classes"))
            confusion.matrix.cv[fold.i,] <- c(results$table[4],results$table[2],results$table[1],results$table[3])
            
            mcc.numerator.cv <- (confusion.matrix.cv[fold.i,1]*confusion.matrix.cv[fold.i,3]-confusion.matrix.cv[fold.i,2]*confusion.matrix.cv[fold.i,4])
            mcc.denominator.cv <- sqrt((confusion.matrix.cv[fold.i,1]+confusion.matrix.cv[fold.i,2])*(confusion.matrix.cv[fold.i,1]+confusion.matrix.cv[fold.i,4])*
                                           (confusion.matrix.cv[fold.i,3]+confusion.matrix.cv[fold.i,2])*(confusion.matrix.cv[fold.i,3]+confusion.matrix.cv[fold.i,4]))
            mcc.cv[fold.i] <- mcc.numerator.cv/mcc.denominator.cv
            
            if(length(unique(testClass)) == 1){
                auc.cv <- auc.cv
            }else{
                pred.roc <- prediction(prob.ensemble.final[,2],as.vector(as.character(testClass)))
                auc.cv[fold.i] <- unlist(performance(pred.roc, "auc")@y.values)
            }
            
        }
        
        result.table[d,2:ncol(result.table)] <- c(mean(auc.cv,na.rm=TRUE),mean(mcc.cv,na.rm = TRUE),apply(result.overall.cv,2,mean,na.rm=TRUE),apply(result.classes.cv,2,mean,na.rm=TRUE)) 
        confusion.matrix[d,2:ncol(confusion.matrix)] <- apply(confusion.matrix.cv,2,mean,na.rm=TRUE)
        
        colnames(result.table) <- c("drugName","AUC","MCC",rownames(as.matrix(results,what = "overall")),rownames(as.matrix(results,what = "classes")))
        result.table$drugName <- colnames(cellLine.class)
        confusion.matrix$drugName <- colnames(cellLine.class)
        
    }
}

write.table(result.table,file = paste(dataSource,"_resultTable_",impRate,"_",correlationRate,".csv",sep = "") ,row.names = F,col.names = T,quote = F,sep = ",")
write.table(confusion.matrix,file = paste(dataSource,"_confusionMatrix_",impRate,"_",correlationRate,".csv",sep = ""),row.names = F,col.names = T,quote = F,sep = ",")





