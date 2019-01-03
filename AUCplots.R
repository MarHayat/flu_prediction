#plot the AUC of three models ALT0,ALT1 and ALT2
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
#ALT0
setwd("/Users/maryam/Desktop/Research/ImperialCollege/Prediction/MyTree/NewEperiments/ALT0")
df=read.csv("df_2018-5.csv",sep= ",",header=T,stringsAsFactors=FALSE)
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

load("/Users/maryam/Desktop/Research/ImperialCollege/Prediction/MyTree/NewEperiments/trees/flutree2018-5.Rdata")
tree=flutree
allHeights=node.depth.edgelength(tree)
#find the clades which their root is less than 3.4 years from the height of the tree and their labels are 0
res=numeric()
for(i in 1:nrow(df)){
  if(df$Labels[i]==0){
    print(i)
    root=df$Clade[i]
    if( max(allHeights)-allHeights[root]<=3.4){res=c(res,i)}
  }
}
length(res)


ind=which(df$Labels==1)
df$Labels[ind]=0
df$Labels[-ind]=1

df=df[-res,]
df=df[,3:ncol(df)]
df=df[,-32]

set.seed(123)
DataTT=PrepareDataRemoveOutliers(df)
train=DataTT[[1]]
test=DataTT[[2]]

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
plot(svmmodel.performance,type="l", col="red")
round(svmmodel.auc,2)
par(new=TRUE)

#ALT1
setwd("/Users/maryam/Desktop/Research/ImperialCollege/Prediction/MyTree/NewEperiments/ALT1")
df=read.csv("df_2018-5",sep= ",",header=T,stringsAsFactors=FALSE)
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

load("/Users/maryam/Desktop/Research/ImperialCollege/Prediction/MyTree/NewEperiments/trees/flutree2018-5.Rdata")
tree=flutree
allHeights=node.depth.edgelength(tree)
#find the clades which their root is less than 3.4 years from the height of the tree and their labels are 0
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
df=df[,-32]

set.seed(123)
DataTT=PrepareDataRemoveOutliers(df)
train=DataTT[[1]]
test=DataTT[[2]]

svm.fit = svm(data = train, Labels ~ .,
              kernel ="linear", degree = 3, gamma =   0.03125 , 
              coef0 = 0, cost =8, nu = 0.5,class.weigth=c("0"=0.50,"1"=0.50))


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
plot(svmmodel.performance,type="l", col="blue")
round(svmmodel.auc,2)
par(new=TRUE)

#ALT2
setwd("/Users/maryam/Desktop/Research/ImperialCollege/Prediction/MyTree/NewEperiments/ALT2")
df=read.csv("df_2018-5",sep= ",",header=T,stringsAsFactors=FALSE)
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

load("/Users/maryam/Desktop/Research/ImperialCollege/Prediction/MyTree/NewEperiments/trees/flutree2018-5.Rdata")
tree=flutree
allHeights=node.depth.edgelength(tree)
#find the clades which their root is less than 3.4 years from the height of the tree and their labels are 0
res=numeric()
for(i in 1:nrow(df)){
  if(df$Labels[i]==0){
    print(i)
    root=df$Clade[i]
    if( max(allHeights)-allHeights[root]<=3){res=c(res,i)}
  }
}
length(res)

df=df[-res,]
df=df[,3:ncol(df)]
df=df[,-32]

ind=which(df$Labels==1)
df$Labels[ind]=0
df$Labels[-ind]=1


set.seed(123)
DataTT=PrepareDataRemoveOutliers(df)
train=DataTT[[1]]
test=DataTT[[2]]

svm.fit = svm(data = train, Labels ~ .,
              kernel ="linear", degree = 3, gamma =   0.03125 , 
              coef0 = 0, cost =2, nu = 0.5,class.weigth=c("0"=0.50,"1"=0.50))


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

legend("bottomright", legend=c("ALT0-AUC:0.89","ALT1-AUC:0.78","ALT2-AUC:0.80"), fill=c("red","blue","green"), cex=0.7)

