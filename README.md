## AutoBorutaRF R source code
AutoBorutaRF is a program by using Autoencoder network and Boruta algorithm for feature selection and random forest for classifier to build drug response prediction model

### AutoencoderBoruta.R is a feature selection function 

Arguments:  
input file       
Processed GDSC lung cell line data set  
GDSClung_class.csv  
GDSClung_rna.csv  
GDSClung_cnv.csv  
GDSClung_mut.csv  
Output file:
selected feature (csv file)  

### RandomForest.R is a classifier using features selected by AutoencoderBoruta.R  

Arguments:
input file     
selected feature file (output by AutoBorutaRF.R)  
Output         
classification results of ten-fold cross validation  

### GDSC and CCLE datasets are available in following URL.  
https://pan.baidu.com/s/1Ey3kmorA7P_cHuWt-qr66g

Xu xiaolu 12/01/2019
