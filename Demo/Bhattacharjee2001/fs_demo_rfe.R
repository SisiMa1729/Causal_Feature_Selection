#################################################################################
# Sisi Ma, 11/02/2018 for AMIA Causal Feature Selection Workshop
# sisima@umn.edu
#################################################################################
# using lung cancer datasets from Bhattacharjee2
# SVM RFE as feature selection and SVM as classifier
rm(list=ls())
print("Predict Lung Cancer vs Normal (SVM RFE)")
setwd('C://Projects/Causal_feature_selection/AMIA2018/Demo/Bhattacharjee2001')
library(e1071)
library(caret)
library(pROC)
source('svm_rfe.R')
set.seed(1)
# load data
data<-read.csv('lung_cancer.csv',header=F)
print(paste('dataset has',nrow(data),'rows, ',ncol(data), 'columns.'))
data[,1]<-as.factor(data[,1])# 0: Adenocarcinoma; 1:Normal
print("# Samples per Class")
print(summary(data[,1]))
levels(data[,1])<-list(c0="0",c1="1")
x<-data[,-1] 
y<-data[,1]
# cross validation
folds<-createFolds(data[,1], k = 5)
mod<-list()
res<-list()
perf<-numeric()
features<-list()
for (f in 1:5){
  test_idx<-folds[[f]]
  train_idx<-setdiff(unlist(folds),test_idx)
  # feature selection
  f_tt<-(f)%%length(folds)+1
  train_test_idx<-folds[[f_tt]]
  train_train_idx<-setdiff(unlist(folds),c(train_test_idx,test_idx))
  features[[f]] <-svm_rfe(x,y, train_train_idx,train_test_idx,2)
  # train
 mod <- svm(x[train_idx,features[[f]]],y[train_idx], 
        kernel = "linear",
        cost=1,
        scale=T,
        probability=T)
  # test
  pred<-predict(mod,newdata=x[test_idx,features[[f]]],probability=T,decision.values=T)
  perf[f]<-roc(y[test_idx],as.vector(attr(pred,"decision.values")),levels=c("c0","c1"),direction=">",plot=F)$auc[[1]]
}


 print(paste('mean AUC = ', mean(perf),sep=''))
 selected_features<-unlist(features)
 freq_table<-as.data.frame(table(selected_features))
 freq_table<-freq_table[with(freq_table, order(-Freq)), ]
 print(freq_table)
 barplot(as.vector(freq_table[,2]),names.arg=freq_table[,1],horiz=TRUE,las=1)


 