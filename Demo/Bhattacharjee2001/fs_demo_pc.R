#################################################################################
# Sisi Ma, 11/02/2018 for AMIA Causal Feature Selection Workshop
# sisima@umn.edu
#################################################################################
# using lung cancer datasets from Bhattacharjee2(a subset: Adenocarcinoma vs Normal)
# pcsimple as feature selection algorithm and SVM as classifier

print("Predict Lung Cancer vs Normal (PC-simple as feature selection)")
rm(list=ls())
library(e1071)
library(caret)
source('fs_pcsimple.R')
setwd('C://Projects/Causal_feature_selection/AMIA2018/Demo/Bhattacharjee2001')

# load data
set.seed(1)
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
features<-list()
perf<-numeric()
for (f in 1:5){
  test_idx<-folds[[f]]
  train_idx<-setdiff(unlist(folds),test_idx)
  #feature selection
  features[[f]] <-fs_pcsimple(x,y, train_idx)
  #train
  mod <- svm(x[train_idx,features[[f]]],y[train_idx], 
        kernel = "linear",
        cost=1,
        scale=T,
        probability=T)
  #test  
  pred<-predict(mod,newdata=x[test_idx,features[[f]]],probability=T,decision.values=T)
  perf[f]<-roc(y[test_idx],as.vector(attr(pred,"decision.values")),levels=c("c0","c1"),direction=">",plot=F)$auc[[1]]
}

print(paste('mean AUC = ', mean(perf),sep=''))
selected_features<-unlist(features)
freq_table<-as.data.frame(table(selected_features))
freq_table<-freq_table[with(freq_table, order(-Freq)), ]
print(freq_table)
barplot(as.vector(freq_table[,2]),names.arg=freq_table[,1],horiz=TRUE,las=1)