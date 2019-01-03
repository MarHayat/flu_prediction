source("~/CrossValidation.R")
load("~/flutree2018-5.Rdata")
tree=flutree
allHeights=node.depth.edgelength(tree)
allD=allDescendants(tree) 

df=read.csv("~/df_2018-5.csv",sep= ",",header=T,stringsAsFactors=FALSE)
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

df=df[-res,]
Clades=df$Clade
df=df[,3:ncol(df)]
df=df[,-32]
ind=which(df$Labels==1)
df$Labels[ind]=0
df$Labels[-ind]=1
#================================================================================================
#================================================================================================
kernel="linear"
r=numeric()
  data=Preparedf(df)
  folds <- cut(seq(1,nrow(data)),breaks=10,labels=FALSE)
  for(cost in c(2^c(-7:7))){ 
    for(gamma in c(2^c(-7:7))){
      result_i=numeric()
      for(j in 1:10){
        print(j)
        testIndexes <- which(folds==j,arr.ind=TRUE)
        test <- data[testIndexes,]
        train <- data[-testIndexes, ]
        print(gamma) 
        print(cost)
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
      }
      r=rbind(r,c(mean(result_i[,1]),mean(result_i[,2]),gamma,cost))
    }
  }
#================================================================================================
#compute the prediction label for each clade in cross validation
kernel="linear"
r=numeric()
set.seed(123)
data=Preparedf(df)
data_clade=data[,1]
folds <- cut(seq(1,nrow(data)),breaks=10,labels=FALSE)
col=c("red","blue","green","yellow","violet","tomato","orange3","magenta","gray64","goldenrod4")
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
  plot(svmmodel.performance,type="l", col=col[j])
  par(new=TRUE)
  res=c(acc,svmmodel.auc)
  result_i=rbind(result_i,res)
  labels=c(labels,svm.prob )
  
}
legend("bottomright", legend=c("fold1:0.73","fold2: 0.78","fold3: 0.89","fold4: 0.84","fold5: 0.90","fold6:0.80","fold7:0.90","fold8:0.80",
                               "fold9: 0.74","fold10: 0.77"), fill=c("red","blue","green","yellow","violet","tomato","orange3","magenta","gray64","goldenrod4"), cex=0.5)
