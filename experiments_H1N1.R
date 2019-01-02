library(e1071)
library(DMwR)
library(outliers)
library(ggplot2)
library("ROCR")
library(caret)
library(Boruta)
library(igraph)
library(ape)
library(phangorn)


source("/Users/maryam/Desktop/Research/ImperialCollege/Prediction/MyTree/H1N1/code/traintest_H1N1.R")
setwd("/Users/maryam/Desktop/Research/ImperialCollege/Prediction/MyTree/H1N1/pnd_H1N1/RemoveLongBranches")
df=read.csv("df_H1N1.csv",sep= ",",header=T,stringsAsFactors=FALSE)
df=df[,2:ncol(df)]
load("flutreeH1N1_RLB.Rdata")
tree=flutree
plot(tree,show.tip.label=FALSE)

allHeights=node.depth.edgelength(tree)
allD=allDescendants(tree) 

ind=which(df$Labels <= 1.1)
length(ind)

df$Labels[ind]=0
df$Labels[-ind]=1
round(length(which(df$Labels==1))/nrow(df),1)
#remove the clades which are in the last 3.4 years of the tree
res=numeric()
for(i in 1:nrow(df)){
  if(df$Labels[i]==0){
    print(i)
    root=df$Clade[i]
    if( max(allHeights)-allHeights[root]<=3.4){res=c(res,i)}
  }
}
length(res)

df=df[-res,]
df=df[,3:ncol(df)]
df=df[,-29]

#=========================================================================================
#remove the very big clades 
outliers=order(df$numberTipsTrimmed,decreasing = T)[1:11]
df$numberTipsTrimmed[outliers]
df=df[-outliers,]
#########################train and test on H1N1###################
#change the labels
ii=which(df$Labels==1)
df$Labels[ii]=0
df$Labels[-ii]=1
#=========================================================================================
set.seed(123)
DataTT=PrepareData(df)
train=DataTT[[1]]
test=DataTT[[2]]

ind=which(train$Labels==1)
train$Labels[ind]=0
train$Labels[-ind]=1
ind=which(test$Labels==1)
test$Labels[ind]=0
test$Labels[-ind]=1

A=t(as.matrix(TuneParametersR(train,test,"linear")))
A

svm.fit = svm(data = train, Labels ~ .,
              kernel ="linear", degree = 3, gamma =   0.03125 , 
              coef0 = 0, cost =32, nu = 0.5,class.weigth=c("0"=0.50,"1"=0.50))


svm.prob <- predict(svm.fit, newdata = test)
summary(svm.prob)

svmmodel.predict<-predict(svm.fit, newdata = test,decision.values=TRUE)
svmmodel.probs<-attr(svmmodel.predict,"decision.values")
svmmodel.class<-predict(svm.fit,test,type="class")
svmmodel.labels<-test$Labels

#roc analysis for test data
svmmodel.prediction<-prediction(svmmodel.probs,svmmodel.labels)
svmmodel.performance<-performance(svmmodel.prediction,"tpr","fpr")
svmmodel.auc<-performance(svmmodel.prediction,"auc")@y.values[[1]]
plot(svmmodel.performance,type="l", col="green")
round(svmmodel.auc,2)
par(new=TRUE)
##########################################################################################
ii=which(df$Labels==1)
df$Labels[ii]=0
df$Labels[-ii]=1
test=df
#=========================================================================================
df=read.csv("/Users/maryam/Desktop/Research/ImperialCollege/Prediction/MyTree/NewEperiments/ALT0/df_2018-5.csv",sep= ",",header=T,stringsAsFactors=FALSE)
df=df[,2:ncol(df)]

names(df)=c("Clade","numberTipsClade","numberTipsTrimmed","sackin",
            "colless","Variance","I2","B1","B2","avgLadder","ILnumber","pitchforks",
            "maxHeight","MaxWidth","DelW","Stairs1","Stairs2","Cherries","BS","descinm","getstattest","skewness","kurtosis","MeanPairwiseDist","MaxPairwiseDist",
            "diameter", "WienerIndex", "betweenness", "closeness", "eigenvector","MeadianEp","MaxEp","MeanEp","numberTipsTrimmed_3.4","Labels")

ind=which(df$Labels <= 1.1)
length(ind)

df$Labels[ind]=0
df$Labels[-ind]=1
round(length(which(df$Labels==1))/nrow(df),1)

load("/Users/maryam/Desktop/Research/ImperialCollege/Prediction/MyTree/NewEperiments/Trees/flutree2018-5.Rdata")
tree=flutree
allHeights=node.depth.edgelength(tree)
allD=allDescendants(tree) 
#remove recent clades which we are not sure about the labesl
res=numeric()
for(i in 1:nrow(df)){
  if(df$Labels[i]==0){
    print(i)
    root=df$Clade[i]
    if( max(allHeights)-allHeights[root]<=3.4){res=c(res,i)}
  }
}
length(res)
df=df[-res,]
#remove unimportant columns
df=df[,3:ncol(df)]
df=df[,-c(29,30,31,32)]
#=========================================================================================
#remove very big clades 
outliers=order(df$numberTipsTrimmed)[(nrow(df)-3):nrow(df)]
df=df[-outliers,]
train=df
ind=nrow(train)
#==========================================================================================
#######################Merge H1N1 and H3N2#################################################
df=rbind(train,test)
#######################train and test on mereged dataset######################################## 
set.seed(123)
DataTT=PrepareDataRemoveOutliers(df)
train=DataTT[[1]]
test=DataTT[[2]]

A=t(as.matrix(TuneParametersR(train,test,"linear")))
A
svm.fit = svm(data = train, Labels ~ .,
              kernel ="linear", degree = 3, gamma =   0.03125 , 
              coef0 = 0, cost =32, nu = 0.5,class.weigth=c("0"=0.50,"1"=0.50))


svm.prob <- predict(svm.fit, newdata = test)
agreement <- svm.prob == test$Labels
acc=length(which(svm.prob == test$Labels))/length(test$Labels)
acc

svmmodel.predict<-predict(svm.fit, newdata = test,decision.values=TRUE)
svmmodel.probs<-attr(svmmodel.predict,"decision.values")
svmmodel.class<-predict(svm.fit,test,type="class")
svmmodel.labels<-test$Labels

#roc analysis for test data
svmmodel.prediction<-prediction(svmmodel.probs,svmmodel.labels)
svmmodel.performance<-performance(svmmodel.prediction,"tpr","fpr")
svmmodel.auc<-performance(svmmodel.prediction,"auc")@y.values[[1]]
plot(svmmodel.performance,type="l", col="red")
round(svmmodel.auc,2)
par(new=TRUE)
#=========================================================================================
#########################train on H3N2 and test on H1N1###################################
df$numberTipsTrimmed <- as.numeric(df$numberTipsTrimmed)
df$sackin<- as.numeric(df$sackin)
df$colless <- as.numeric(df$colless)
df$Variance <- as.numeric(df$Variance)
df$I2 <- as.numeric(df$I2)
df$B1 <- as.numeric(df$B1)
df$B2<- as.numeric(df$B2)
df$avgLadder <- as.numeric(df$avgLadder)
df$ILnumber  <- as.numeric(df$ILnumber)
df$pitchforks <- as.numeric(df$pitchforks)
df$maxHeight <- as.numeric(df$maxHeight)
df$MaxWidth <-as.numeric(df$MaxWidth)
df$DelW <- as.numeric(df$DelW)
df$Stairs1 <- as.numeric(df$Stairs1)
df$Stairs2 <- as.numeric(df$Stairs2)
df$Cherries <- as.numeric(df$Cherries)
df$BS <- as.numeric(df$BS )
df$descinm <- as.numeric(df$descinm)
df$getstattest <- as.numeric(df$getstattest)
df$skewness<- as.numeric(df$skewness)
df$kurtosis <- as.numeric(df$kurtosis )
df$MeanPairwiseDist<- as.numeric(df$MeanPairwiseDist)
df$MaxPairwiseDist<- as.numeric(df$MaxPairwiseDist)
df$diameter<- as.numeric(df$diameter)
df$WienerIndex<- as.numeric(df$WienerIndex)
df$betweenness<- as.numeric(df$betweenness)
df$closeness<- as.numeric(df$closeness)
df$eigenvector <-as.numeric(df$eigenvector)

df=scaleData(df)
round(length(which(df$Labels==1))/nrow(df),1)
train_ind=c(1:ind)
train=df[train_ind,]
test=df[-train_ind,]

set.seed(321)
kernel="linear"
tune.out=tune(svm, Labels ~ .,
              data = train,kernel =kernel, ranges=list(gamma = 2^c(-5:5),cost=2^c(-5:5)),
              coef0 =0,degree =3,nu = 0.5,class.weigth=c("0"=0.5,"1"=0.5))


gamma=tune.out$best.parameters$gamma
cost=tune.out$best.parameters$cost
degree = 3
coef0 = 0
svm.fit = svm(data = train, Labels ~ .,
              kernel =kernel, degree = 3, gamma = gamma,
              coef0 = 0, cost =cost, nu = 0.5, class.weigth=c("0"=0.5,"1"=0.5))


svm.prob <- predict(svm.fit, newdata = test)
agreement <- svm.prob == test$Labels
acc=length(which(svm.prob == test$Labels))/length(test$Labels)
acc
svmmodel.predict<-predict(svm.fit, newdata = test,decision.values=TRUE)
svmmodel.probs<-attr(svmmodel.predict,"decision.values")
svmmodel.class<-predict(svm.fit,test,type="class")
svmmodel.labels<-test$Labels

#roc analysis for test data
svmmodel.prediction<-prediction(svmmodel.probs,svmmodel.labels)
svmmodel.performance<-performance(svmmodel.prediction,"tpr","fpr")
svmmodel.auc<-performance(svmmodel.prediction,"auc")@y.values[[1]]
plot(svmmodel.performance,type="l", col="blue")
round(svmmodel.auc,2)
par(new=TRUE)

#================================================================================
legend("bottomright", legend=c("Influenza H1N1 tree-AUC:0.86","merge H1N1 and H3N2 clades-AUC:0.85","train on H3N2 and test on H1N1-AUC:0.75"), fill=c("green","red","blue"), cex=0.7)
#================================================================================
