#################################################################################
# Sisi Ma, 11/02/2018 for AMIA Causal Feature Selection Workshop
# sisima@umn.edu
#################################################################################
# feature selection for survival outcome with lasso
rm(list=ls())
library(survival)
library(caret)
library(glmnet)
library(cancerdata) # data from VIJVER2002
setwd('C://Projects/Causal_feature_selection/AMIA2018/Demo/Vijver2002/')

set.seed(1)
# extrect the 70 genes that are predictive of the prognosis as full feature set
gene<-read.csv("./70genelist.csv",header=FALSE)
gene<-as.character(gene[,1])
data('VIJVER')
idx<-which(is.element(rownames(VIJVER@assayData$exprs),gene))
x<-as.data.frame(t(VIJVER@assayData$exprs[idx,]))
# construct survival outcome
y<-Surv(VIJVER$Follow_up_time_or_metastasis,VIJVER$event_metastasis)
folds<-createFolds(VIJVER$event_metastasis, k = 10)
mod<-list()
features<-list()
perf<-numeric()
# cross validation
for (f in 1:10){
  test_idx<-folds[[f]]
  train_idx<-setdiff(unlist(folds),test_idx)
  # train (feature selection and training classifier is done in one step)
  cv_mod <-cv.glmnet(as.matrix(x[train_idx,]), y[train_idx], family="cox",nfolds=9)
  features[[f]]<-which(coef(cv_mod)!=0)
  # test
  pred<-predict(cv_mod,newx=as.matrix(x[test_idx,,drop=FALSE]))# alternative setting s="lambda.min"
  perf[f]<-survConcordance(y[test_idx]~pred)$concordance
}
print(paste('c-idx = ', mean(perf),sep=''))
selected_features<-unlist(features)
freq_table<-as.data.frame(table(selected_features))
freq_table<-freq_table[with(freq_table, order(-Freq)), ]
print(freq_table)
barplot(as.vector(freq_table[,2]),names.arg=freq_table[,1],horiz=TRUE,las=1)

