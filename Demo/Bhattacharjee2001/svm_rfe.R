# svm rfe feature selection
#################################################################################
# Sisi Ma, 11/02/2018 for AMIA Causal Feature Selection Workshop
# sisima@umn.edu
#################################################################################
svm_rfe <- function(x,y, train_train_idx,train_test_idx,reduction_coeff){
# input: x feature matrix, columns are variables, rows are observation
#        y outcome, binary vector
#        train_train_idx, indices for train train set
#        train_test_idx, indices for train test set
#        reduction_coeff, scalar 1/reduction_coeff specifies proporion of variables to drop per interation.
# output: selected features
##################################################################################
# this is a bare bone implementation for illustrative purpose and can be improved, see
# comments in code
#################################################################################
  library(e1071)
# generate array defining number of features
N_f_set<-unique(round(exp(seq(from=log(ncol(x)),to=log(1),by=-log(reduction_coeff)))));
tr_x<-x[train_train_idx,]
tr_y<-y[train_train_idx]
tst_x<-x[train_test_idx,]
tst_y<-y[train_test_idx]
# loop over sets of features
cnt<-0
features<-list()
perf<-double()
for (N_f in N_f_set){
cnt<-cnt+1
#print(cnt)
  if (cnt==1){
  features[[1]]<-1:ncol(x)
  } else {
    st<-sort(abs(w),decreasing=T,index.return=T)
    sorted_previous_features<-features[[cnt-1]][st$ix]
    features[[cnt]]<-sorted_previous_features[1:N_f]
  }
  # train model
   mod<- svm(tr_x[,features[[cnt]]],tr_y, 
	kernel = "linear",
        cost=1,
	      scale=T,
        probability=T)
  pred<-predict(mod,newdata=x[train_test_idx,features[[cnt]]],probability=T,decision.values=T)
   perf[cnt]<-roc(tst_y,as.vector(attr(pred,"decision.values")),levels=c("c0","c1"),direction=">",plot=F)$auc[[1]]
  w <- t(mod$coefs) %*%as.matrix(tr_x[,features[[cnt]],drop=F][mod$index,]) 
}
mperf<-max(perf)
set_idx<-max(which(perf==mperf)) 
# can add statisticall test here to pick a set that have performance that is statistically not
# significantly different from the set with max performance, but with the least number of variables
#print(set_idx)
final_features<-features[[set_idx]]
return(final_features)

}

