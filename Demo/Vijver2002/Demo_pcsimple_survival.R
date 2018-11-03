#################################################################################
# Sisi Ma, 11/02/2018 for AMIA Causal Feature Selection Workshop
# sisima@umn.edu
#################################################################################
# causal feature selection (modified pc simple) for survival 
rm(list=ls())
library(survival)
library(caret)
library(pcalg)
library(cancerdata) # data from VIJVER2002
setwd('C://Projects/Causal_feature_selection/AMIA2018/Demo/Vijver2002/')
source('pcsimple_survival.R')

set.seed(1)
# extrect the 70 genes that are predictive of the prognosis as full feature set
data('VIJVER')
save(VIJVER,file='./VIJVER.rda')
gene<-read.csv("./70genelist.csv",header=FALSE)
gene<-as.character(gene[,1])
idx<-which(is.element(rownames(VIJVER@assayData$exprs),gene))
x<-as.data.frame(t(VIJVER@assayData$exprs[idx,]))
# construct survival outcome
y<-Surv(VIJVER$Follow_up_time_or_metastasis,VIJVER$event_metastasis) # time to metastasis
# cross validation
folds<-createFolds(VIJVER$event_metastasis, k = 10)
features<-list()
perf<-numeric()
for (f in 1:10){
  test_idx<-folds[[f]]
  train_idx<-setdiff(unlist(folds),test_idx)
  # fs
  features[[f]]<-which(pcsimple_survival(y[train_idx,], x[train_idx,], 0.05)$G)
  # train
  cdata<-cbind(y,x[,features[[f]],drop=FALSE])[train_idx,]
  cfit<-coxph(y~.,data=cdata)
  # test
  pred<-predict(cfit,newdata=x[test_idx,features[[f]],drop=FALSE])
  perf[f]<-survConcordance(y[test_idx]~pred)$concordance
}
print(paste('c-idx = ', mean(perf),sep=''))
selected_features<-unlist(features)
freq_table<-as.data.frame(table(selected_features))
freq_table<-freq_table[with(freq_table, order(-Freq)), ]
print(freq_table)
barplot(as.vector(freq_table[,2]),names.arg=freq_table[,1],horiz=TRUE,las=1)

