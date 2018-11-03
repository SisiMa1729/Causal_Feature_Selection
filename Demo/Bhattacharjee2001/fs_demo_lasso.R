#################################################################################
# Sisi Ma, 11/02/2018 for AMIA Causal Feature Selection Workshop
# sisima@umn.edu
#################################################################################
# using lung cancer datasets from Bhattacharjee2(a subset: Adenocarcinoma vs Normal)
# read data
print("Predict Lung Cancer vs Normal (penalized logistic regression)")
rm(list=ls())
setwd('C://Projects/Causal_feature_selection/AMIA2018/Demo/Bhattacharjee2001')
library(pROC)
library(caret)
library(glmnet)
set.seed(1)
data<-read.csv('lung_cancer.csv',header=F)
print(paste('dataset has',nrow(data),'rows, ',ncol(data), 'columns.'))
data[,1]<-as.factor(data[,1]) # 0: Adenocarcinoma; 1:Normal
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
  # training (features selection and model building is done in one step)
  cv_mod <-cv.glmnet(as.matrix(x[train_idx,]), y[train_idx], family="binomial",nfolds=4)           
  features[[f]]<-setdiff(which(coef(cv_mod)!=0),1)-1 # the first coefficient correspond to intercept
  # testing
  pred<-predict(cv_mod,newx=as.matrix(x[test_idx,]),type="response")# alternative setting s="lambda.min"
  perf[f]<-roc(y[test_idx],as.vector(pred),levels=c("c0","c1"),direction="<",plot=F)$auc[[1]]
}
 
print(paste('mean AUC = ', mean(perf),sep=''))
selected_features<-unlist(features)
freq_table<-as.data.frame(table(selected_features))
freq_table<-freq_table[with(freq_table, order(-Freq)), ]
print(freq_table)
barplot(as.vector(freq_table[,2]),names.arg=freq_table[,1],horiz=TRUE,las=1)


