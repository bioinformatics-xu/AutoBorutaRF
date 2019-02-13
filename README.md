# AutoBorutaRF R source code
AutoBorutaRF is a project by using Autoencoder network and Boruta algorithm for feature selection and random forest for classification to build drug response prediction model.

## Computational procedure
The computational procedure implemented in AutoBorutaRF is aimed to identify the genetic features responsible for drug response classification and predict the drug response. The procedure includes two steps:
1. Selecting the genetic features
2. Predicting the drug response

These two steps use AutoencoderBoruta.R and RandomForest.R functions respectively.

### AutoencoderBoruta.R
AutoencoderBoruta.R is a two-stepwise feature selection function using the autoencoder network and Boruta algorithm to identify the genetic features responsible for drug response.
 
#### Input files:   
class.csv 
(class.csv is the category label file containing a binary-valued matrix, with rows representing samples and columns representing drugs. The element in row i and column j is "1" if sample i is "sensitive" to drug j, "0" if "non-sensitive".) 

rna.csv
(rna.csv contains the gene expression feature matrix, with rows representing samples and columns representing genes. The elements in the gene expression feature matrix are real-valued.)

cnv.csv
(cnv.csv contains the gene copy number alteration feature matrix, with rows representing samples and columns representing genes. The elements in the gene copy number alteration feature matrix are real-valued.)

mut.csv
(mut.csv contains the gene mutation feature matrix, with rows representing samples and columns representing genes. The elements in the gene mutation feature matrix are binary-valued, "1" for mutation and "0" for wild type.)

#### Notice: 
Class.csv is the label file, and it is mandatory for input; rna.csv, cnv.csv, and mut.csv are the feature files, and they are optional, users could input one, two, or all of them.

#### Output files:
predictor.csv
(predictor.csv contains the selected features, with rows representing drugs. The first element in every row is the drug name, the rest are selected features of the drug.)

### RandomForest.R
RandomForest.R is a classification function to predict drug response using the genetic features selected by AutoencoderBoruta.R.

#### Input file: 
class.csv (the same class.csv as AutoencoderBoruta.R)

rna.csv (the same rna.csv as AutoencoderBoruta.R)

cnv.csv (the same cnv.csv as AutoencoderBoruta.R)

mut.csv (the same mut.csv as AutoencoderBoruta.R)

predictor.csv (predictor.csv file contains the genetic features selected by AutoBorutaRF.R)

#### Output files: 
resultTable.csv
(resultTable.csv contains the classification results of ten-fold cross-validation including the Accuracy, AUC, MCC, Recall, Specificity, and so on.)

confusionMatrix.csv
(confusionMatrix.csv contains the confusion matrices of ten-fold cross-validation of drugs.)

### Example
An example run is:

First, set "dataSource <- 'GDSClung'" in AutoencoderBoruta.R function to use "GDSClung dataset" to identify the genetic features.

Second, set "dataSource <- 'GDSClung'" in RandomForest.R function to use the selected features of "GDSClung dataset" to predict drug response.

### Result
In the resultTable.csv file, we show the mean of some important evaluation criteria obtained by 10-fold cross validation for all drugs. All the evaluation criteria except Matthews correlation coefficient (MCC) are obtained by "confusionMatrix" function (from R package "caret"). Each row corresponds to the results of each drug, for example, drug thapsigargin showing AUC 0.954, MCC 0.7436, accuracy 0.8556 and so on.

In the confusionMatrix.csv file, we show the mean of true positive (TP), false positive(FP), ture negtive(TN), and false negtive(FN) obtained by 10-fold cross validation for all drugs. Each row corresponds to the results of each drug, for example, drug thapsigargin showing TP 4.2, FP	0.9, TN	3.9, and TN	0.5.

### Memory
Because of the training of autoencoder with H2O package, we normally suggest 100 times the size of the dataset for the amount of memory required.


### GDSC and CCLE datasets are available in the following URL.  
https://1drv.ms/u/s!AuqTdFjFP1Q8arSHSFFCco2nwRw

Xu xiaolu 12/01/2019
