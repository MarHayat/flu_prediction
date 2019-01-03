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
source("~/traintest_H3N2")
setwd("~/ALT0")
#load the tree
load("~/flutree2018-5.Rdata")
tree=flutree
#compute the height of each node
allHeights=node.depth.edgelength(tree)
allD=allDescendants(tree) 
#read the dataset(each clade with its features)
df=read.csv("df_2018-5.csv",sep= ",",header=T,stringsAsFactors=FALSE)
df=df[,2:ncol(df)]

names(df)=c("Clade","numberTipsClade","numberTipsTrimmed","sackin","colless","Variance","I2","B1","B2",
            "avgLadder","ILnumber","pitchforks","maxHeight","MaxWidth","DelW","Stairs1","Stairs2","Cherries",
            "BS","descinm","getstattest","skewness","kurtosis","MeanPairwiseDist","MaxPairwiseDist", "diameter", 
            "WienerIndex", "betweenness", "closeness", "eigenvector","MeadianEp","MaxEp","MeanEp",
            "numberTipsTrimmed_3.4","Labels")

ind=which(df$Labels <= 1.1)
length(ind)
#assign the labels
df$Labels[ind]=0
df$Labels[-ind]=1
round(length(which(df$Labels==1))/nrow(df),1)

#find the clades which their root is in the last 3.4 years of the tree and their labels are 0
#we are not sure about the labels of these clades since they do not have enough time to growth
res=numeric()
for(i in 1:nrow(df)){
  if(df$Labels[i]==0){
    print(i)
    root=df$Clade[i]
    if( max(allHeights)-allHeights[root]<=3.4){res=c(res,i)}
  }
}
length(res)
#extract the clades which we are not sure about their labels
DD=df[res,c(1,2,3,34,35)]
#remove the unimportant columns
df=df[,3:ncol(df)]
df=df[,-32]
#we then use our model to predict the labels of the recent clades which we removed from the dataset
df_predict=df[res,]
#remove the uncertain clades from dataset (we do not want to include them in train or test)
df=df[-res,]
#change the labels if neccessary (otherwise AUC is less than 0.50)
ii=which(df$Labels==1)
df$Labels[ii]=0
df$Labels[-ii]=1

set.seed(123)
DataTT=PrepareDataRemoveOutliers(df,1)
train=DataTT[[1]]
test=DataTT[[2]]
#find the best hyperparametersfor the model
B=TuneParametersR(train,test,"linear")
#==================================================================================================
#use the best hyperparameters to train the model
svm.fit = svm(data = train, Labels ~ .,
              kernel ="linear",degree = 3, gamma = 0.03125,
              coef0 = 0, cost =32, nu = 0.5,class.weigth=c("0"=0.5,"1"=0.5))


svm.prob <- predict(svm.fit, newdata = test)
agreement <- svm.prob == test$Labels
acc=length(which(svm.prob == test$Labels))/length(test$Labels)
round(acc,2)
svmmodel.predict<-predict(svm.fit, newdata = test,decision.values=TRUE)
svmmodel.probs<-attr(svmmodel.predict,"decision.values")
svmmodel.class<-predict(svm.fit,test,type="class")
svmmodel.labels<-test$Labels

#roc analysis for test data
svmmodel.prediction<-prediction(svmmodel.probs,svmmodel.labels)
svmmodel.performance<-performance(svmmodel.prediction,"tpr","fpr")
svmmodel.auc<-performance(svmmodel.prediction,"auc")@y.values[[1]]
round(svmmodel.auc,2)
plot(svmmodel.performance)
#==================================================================================================
####################################predicton on recent clades#####################################
df_predict=scaleData(df_predict)
svm.prob <- predict(svm.fit, newdata = df_predict)
table(svm.prob)
DD=cbind(DD,svm.prob)
#because we changed the labels of the data before training
ii=which(DD$svm.prob==1)
DD$svm.prob[ii]=0
DD$svm.prob[-ii]=1
write.csv(DD,"Prediction_Recent.csv")
#==================================================================================================
####################################predicton on clades after 2016#################################
df_predict=scaleData(df_predict)
df=read.csv("df_2018-5.csv",sep= ",",header=T,stringsAsFactors=FALSE)
df=df[,2:ncol(df)]
names(df)=c("Clade","numberTipsClade","numberTipsTrimmed","sackin","colless","Variance","I2","B1","B2",
            "avgLadder","ILnumber","pitchforks","maxHeight","MaxWidth","DelW","Stairs1","Stairs2","Cherries",
            "BS","descinm","getstattest","skewness","kurtosis","MeanPairwiseDist","MaxPairwiseDist", "diameter", 
            "WienerIndex", "betweenness", "closeness", "eigenvector","MeadianEp","MaxEp","MeanEp",
            "numberTipsTrimmed_3.4","Labels")
ind=which(df$Labels <= 1.1)
length(ind)

df$Labels[ind]=0
df$Labels[-ind]=1
round(length(which(df$Labels==1))/nrow(df),1)
#find the clades which their tips are after 2016
res_2016=numeric()
myclades=getClades2(flutree, MinTotalSize = 8, MinTrimSize = 8, TimeFrame = 1.4)
Final_trimmedClades=myclades$trimclades[myclades$rejected==0]
Final_trimmedClades_root=as.numeric(names(which(myclades$rejected==0)))
Auxdata=read.csv("~/RNAP2018.csv",sep= ",",header=T,stringsAsFactors=FALSE)
Auxdata=Auxdata[,2:4]
allHeights=node.depth.edgelength(flutree); max(allHeights)  
hdata=data.frame(tiplab=flutree$tip.label, height=allHeights[1:length(flutree$tip.label)]) 

res_2016=numeric()
for(i in 1:length(df$Clade)){
  print(i)
  tr1=extract.clade(flutree,df$Clade[i])
  tr=drop.tip(tr1,setdiff(tr1$tip.label,hdata$tiplab[Final_trimmedClades[[i]]]), trim.internal = TRUE)
  ind=match(tr$tip.label,Auxdata[,1])
  date=as.Date(Auxdata[ind,3])
  if(min(date)>"2016-01-01"){res_2016=c(res_2016,i)}
}

df$Clade[res_2016]
DD=df[res_2016,c(1,2,3,34,35)]
df_predict=df[res_2016,]
#predicton on clades after 2016
df_predict=scaleData(df_predict)
svm.prob <- predict(svm.fit, newdata = df_predict)
table(svm.prob)
DD=cbind(DD,svm.prob )
#because we changed the labels of the data before training
ii=which(DD$svm.prob==1)
DD$svm.prob[ii]=0
DD$svm.prob[-ii]=1
write.csv(DD,"Prediction_2016.csv")

PredLab=rep(NA, length(res_2016)) 
ind=match(df$Clade[res_2016],dd$Clade)
ii=which(is.na(ind))
X=dd[ind[-ii],]
PredLab
ind=match(X$Clade,DD$Clade)
PredLab[ind]=X$svm.prob
DD=cbind(DD,PredLab)
write.csv(DD,"prediction_2016")
dd=read.csv("prediction_2016")
#==================================================================================================
#############################clades that includes just 2016 strains################################
res_2016=numeric()
for(i in 1:length(df$Clade)){
  print(i)
  tr1=extract.clade(flutree,df$Clade[i])
  tr=drop.tip(tr1,setdiff(tr1$tip.label,hdata$tiplab[Final_trimmedClades[[i]]]), trim.internal = TRUE)
  ind=match(tr$tip.label,Auxdata[,1])
  date=as.Date(Auxdata[ind,3])
  if(min(date)>="2016-01-01"&&max(date) <= "2016-12-31"){res_2016=c(res_2016,i)}
}
#==================================================================================================
##############################train on the past and test in recent clades##########################
load("~/flutree2018-5.Rdata")
tree=flutree
allHeights=node.depth.edgelength(tree)
allD=allDescendants(tree) 

setwd("~/ALT0")
df=read.csv("df_2018-5.csv",sep= ",",header=T,stringsAsFactors=FALSE)
df=df[,2:ncol(df)]
names(df)=c("Clade","numberTipsClade","numberTipsTrimmed","sackin","colless","Variance","I2","B1","B2",
            "avgLadder","ILnumber","pitchforks","maxHeight","MaxWidth","DelW","Stairs1","Stairs2","Cherries",
            "BS","descinm","getstattest","skewness","kurtosis","MeanPairwiseDist","MaxPairwiseDist", "diameter", 
            "WienerIndex", "betweenness", "closeness", "eigenvector","MeadianEp","MaxEp","MeanEp",
            "numberTipsTrimmed_3.4","Labels")
ind=which(df$Labels <= 1.1)
length(ind)
df$Labels[ind]=0
df$Labels[-ind]=1
round(length(which(df$Labels==1))/nrow(df),1)
#find the clades which their root is in the last 3.4 years of the tree and their labels are 0
#we are not sure about the labels of these clades since they do not have enough time to growth
res=numeric()
for(i in 1:nrow(df)){
  if(df$Labels[i]==0){
    print(i)
    root=df$Clade[i]
    if( max(allHeights)-allHeights[root]<=3.4){res=c(res,i)}
  }
}
length(res)

resCl=df$Clade[res]
res_2015=numeric()
for(i in 1:length(df$Clade)){
  print(i)
  tr1=extract.clade(flutree,df$Clade[i])
  tr=drop.tip(tr1,setdiff(tr1$tip.label,hdata$tiplab[Final_trimmedClades[[i]]]), trim.internal = TRUE)
  ind=match(tr$tip.label,Auxdata[,1])
  date=as.Date(Auxdata[ind,3])
  if(min(date)>"2015-01-01"){res_2015=c(res_2015,i)}
}
res_2015Cl=df$Clade[res_2015]

test=df[res_2015,]
test=test[which(is.na(match(test$Clade,df$Clade[res]))),]
dim(test)
train=df[-res_2015,]
train=train[which(is.na(match(train$Clade,df$Clade[res]))),]
dim(train)
df=rbind(test,train)
df=df[,3:ncol(df)]
df=df[,-32]
test=df[1:dim(test)[1],]
train=df[(dim(test)[1]+1):nrow(df),]

B=TuneParametersR(train,test,"linear")

