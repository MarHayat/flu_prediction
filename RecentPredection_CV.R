#prediction using cross validation
df=read.csv("df_2018-5.csv",sep= ",",header=T,stringsAsFactors=FALSE)
df=df[,2:ncol(df)]
names(df)=c("Clade","numberTipsClade","numberTipsTrimmed","sackin","colless","Variance","I2","B1","B2",
            "avgLadder","ILnumber","pitchforks","maxHeight","MaxWidth","DelW","Stairs1","Stairs2","Cherries",
            "BS","descinm","getstattest","skewness","kurtosis","MeanPairwiseDist","MaxPairwiseDist", "diameter", 
            "WienerIndex", "betweenness", "closeness", "eigenvector","MeadianEp","MaxEp","MeanEp",
            "numberTipsTrimmed_3.4","Labels")

load("~/flutree2018-5.Rdata")
tree=flutree
allHeights=node.depth.edgelength(tree)
allD=allDescendants(tree) 

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

RecentData=df[res,]
DD=RecentData[,c(1,2,3,34,35)]
df=df[-res,]
df=df[,3:ncol(df)]
df=df[,-32]

kernel="linear"
r=numeric()
set.seed(123)
data=Preparedf(df)
folds <- cut(seq(1,nrow(data)),breaks=10,labels=FALSE)
result_i=numeric()
labels=numeric()
for(j in 1:10){
  print(j)
  testIndexes <- which(folds==j,arr.ind=TRUE)
  test <- data[testIndexes,2:ncol(data) ]
  train <- data[-testIndexes,2:ncol(data) ]
  gamma=0.03125
  cost=64
  svm.fit = svm(data = train, Labels ~ .,
                kernel =kernel, degree = 3, gamma = gamma,
                coef0 = 0, cost =cost, nu = 0.5, class.weigth=c("0"=0.5,"1"=0.5))
  
  
  svm.prob <- predict(svm.fit, newdata = test)
  agreement <- svm.prob == test$Labels
  acc=length(which(svm.prob == test$Labels))/length(test$Labels)
  svmmodel.predict<-predict(svm.fit, newdata = test,decision.values=TRUE)
  svmmodel.probs<-attr(svmmodel.predict,"decision.values")
  svmmodel.class<-predict(svm.fit,test,type="class")
  svmmodel.labels<-test$Labels
  
  #roc analysis for test data
  svmmodel.prediction<-prediction(svmmodel.probs,svmmodel.labels)
  svmmodel.performance<-performance(svmmodel.prediction,"tpr","fpr")
  svmmodel.auc<-performance(svmmodel.prediction,"auc")@y.values[[1]]
  svmmodel.auc
  res=c(acc,svmmodel.auc)
  result_i=rbind(result_i,res)
  labels=c(labels,svm.prob )
  
  #prediction on recent clades
  testData=RecentData[,3:ncol(RecentData)]
  testData=testData[,-32]
  testData=scaleData(testData)
  svm.prob <- predict(svm.fit, newdata = testData)
  DD=cbind(DD,svm.prob)
}

write.csv(DD,"RecentPredictionCV")
#==================================================================================================
####################################predicton on clades after 2016#################################
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
Auxdata=read.csv("/Users/maryam/Desktop/Research/ImperialCollege/Prediction/MyTree/NewEperiments/Trees/RNAP2018.csv",sep= ",",header=T,stringsAsFactors=FALSE)
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

RecentData=df[res_2016,]
DD=RecentData[,c(1,2,3,34,35)]
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

df=df[-res,]
df=df[,3:ncol(df)]
df=df[,-32]

kernel="linear"
r=numeric()
set.seed(123)
data=Preparedf(df)
folds <- cut(seq(1,nrow(data)),breaks=10,labels=FALSE)
result_i=numeric()
labels=numeric()
for(j in 1:10){
  print(j)
  testIndexes <- which(folds==j,arr.ind=TRUE)
  test <- data[testIndexes,2:ncol(data) ]
  train <- data[-testIndexes,2:ncol(data) ]
  gamma=0.03125
  cost=64
  svm.fit = svm(data = train, Labels ~ .,
                kernel =kernel, degree = 3, gamma = gamma,
                coef0 = 0, cost =cost, nu = 0.5, class.weigth=c("0"=0.5,"1"=0.5))
  
  
  svm.prob <- predict(svm.fit, newdata = test)
  agreement <- svm.prob == test$Labels
  acc=length(which(svm.prob == test$Labels))/length(test$Labels)
  svmmodel.predict<-predict(svm.fit, newdata = test,decision.values=TRUE)
  svmmodel.probs<-attr(svmmodel.predict,"decision.values")
  svmmodel.class<-predict(svm.fit,test,type="class")
  svmmodel.labels<-test$Labels
  
  #roc analysis for test data
  svmmodel.prediction<-prediction(svmmodel.probs,svmmodel.labels)
  svmmodel.performance<-performance(svmmodel.prediction,"tpr","fpr")
  svmmodel.auc<-performance(svmmodel.prediction,"auc")@y.values[[1]]
  svmmodel.auc
  res=c(acc,svmmodel.auc)
  result_i=rbind(result_i,res)
  labels=c(labels,svm.prob )
  
  #prediction on recent clades
  testData=RecentData[,3:ncol(RecentData)]
  testData=testData[,-32]
  testData=scaleData(testData)
  svm.prob <- predict(svm.fit, newdata = testData)
  DD=cbind(DD,svm.prob)
}

write.csv(DD,"2016PredictionCV")
