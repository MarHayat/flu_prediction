#function to randomley choose train and test datasets
PrepareDataRemoveOutliers=function(df){
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
  df$MeadianEp <-as.numeric(df$MeadianEp)
  df$MaxEp <-as.numeric(df$MaxEp)
  df$MeanEp <-as.numeric(df$MeanEp)
  #================================================================================================
  #remove the outliers
  df_outlier=df[,1:(ncol(df)-1)]
  outlier.scores <- lofactor(df_outlier, k=5)
  #plot(density(outlier.scores))
  # pick top 5 as outliers
  outliers <- order(outlier.scores, decreasing=T)[1:5]
  # who are outliers
  print(outliers)
  df=df[-outliers,]
  #================================================================================================
  df=scaleData(df)
  #================================================================================================
  #shuffle the data
  set.seed(17)
  df_cur=df[sample(seq_len(nrow(df)), size = nrow(df)),]
  round(length(which(df$Labels==1))/nrow(df),1)
  #================================================================================================
  #randomley choose 75% of data for train and the rest for test
    smp_size <- floor(0.75* nrow(df_cur))
    set.seed(123)
    train_ind <- sample(seq_len(nrow(df_cur)), size = smp_size)
    train<- df_cur[train_ind, ]
    test <- df_cur[-train_ind, ]
    LTT=list(train=train,test=test)
    return(LTT)
}
#================================================================================================
TuneParametersR=function(train,test,kernel){
  if(kernel=="linear" | kernel=="radial"){ 
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
    svmmodel.predict<-predict(svm.fit, newdata = test,decision.values=TRUE)
    svmmodel.probs<-attr(svmmodel.predict,"decision.values")
    svmmodel.class<-predict(svm.fit,test,type="class")
    svmmodel.labels<-test$Labels
    
    #roc analysis for test data
    svmmodel.prediction<-prediction(svmmodel.probs,svmmodel.labels)
    svmmodel.performance<-performance(svmmodel.prediction,"tpr","fpr")
    svmmodel.auc<-performance(svmmodel.prediction,"auc")@y.values[[1]]
    svmmodel.auc
    
    table(test$Labels,svm.prob)
    agreement <- svm.prob == test$Labels
    table(agreement)
    prop.table(table(agreement))
    
    #computing the true possitive rate
    TPR=as.numeric(table(test$Labels,svm.prob)[2,2]/table(test$Labels)[2])
    #computing the false negative rate
    FNR=1-TPR 
    #computing the true negative rate
    TNR=as.numeric(table(test$Labels,svm.prob)[1,1]/table(test$Labels)[1])
    #computing the false positive rate
    FPR=1-TNR
    r=c(prop.table(table(agreement))[[2]],TNR,TPR,svmmodel.auc,gamma,cost,degree,coef0)
    #results_i=round(c(prop.table(table(agreement))[[2]],TPR,FNR,TNR,FPR,svmmodel.auc),3)
    
  }
  else{ 
    tune.out=tune(svm, Labels ~ .,
                  data = train,kernel =kernel, ranges=list(gamma = 2^c(-4:4),cost=2^c(-4:4),degree=c(3,4,5),coef0=c(0,1)),
                  nu = 0.5,class.weigth=c("0"=0.5,"1"=0.5))
    
    
    gamma=tune.out$best.parameters$gamma
    cost=tune.out$best.parameters$cost
    degree=tune.out$best.parameters$degree
    coef0=tune.out$best.parameters$coef0
    
    svm.fit = svm(data = train, Labels ~ .,
                  kernel =kernel, degree = degree, gamma = gamma,
                  coef0 = coef0, cost =cost, nu = 0.5, class.weigth=c("0"=0.5,"1"=0.5))
    
    
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
    
    table(test$Labels,svm.prob)
    agreement <- svm.prob == test$Labels
    table(agreement)
    prop.table(table(agreement))
    
    #computing the true possitive rate
    TPR=as.numeric(table(test$Labels,svm.prob)[2,2]/table(test$Labels)[2])
    #computing the false negative rate
    FNR=1-TPR 
    #computing the true negative rate
    TNR=as.numeric(table(test$Labels,svm.prob)[1,1]/table(test$Labels)[1])
    #computing the false positive rate
    FPR=1-TNR
    r=c(prop.table(table(agreement))[[2]],TNR,TPR,svmmodel.auc,gamma,cost,degree,coef0)
    #results_i=round(c(prop.table(table(agreement))[[2]],TPR,FNR,TNR,FPR,svmmodel.auc),3)
    
  }
  return(r)
}

scaleData=function(data){
  scaled.data <- scale(data[,1:(ncol(data)-1)])
  scaled.data=as.data.frame(scaled.data)
  data$Labels=as.factor(data$Labels)
  data=cbind(scaled.data ,data$Labels)
  colnames(data)=c("numberTipsTrimmed","sackin",
                   "colless","Variance","I2","B1","B2","avgLadder","ILnumber","pitchforks",
                   "maxHeight","MaxWidth","DelW","Stairs1","Stairs2","Cherries","BS","descinm","getstattest","skewness","kurtosis","MeanPairwiseDist","MaxPairwiseDist",
                   "diameter", "WienerIndex", "betweenness", "closeness", "eigenvector","MeadianEp","MaxEp","MeanEp","Labels")
  return(data)
}


RandomFr=function(train,test){
  set.seed(17)
  rf_model<-train(Labels~.,data=train,method="rf",
                  trControl=trainControl(method="repeatedcv", number=10, repeats=3),
                  prox=TRUE,allowParallel=TRUE)
  print(rf_model)
  print(rf_model$finalModel)
  rf.fit=rf_model$finalModel
  rf.prob <- predict(rf.fit, newdata = test)
  rf.prob.res = summary(rf.prob)
  
  table(test$Labels,rf.prob)
  agreement <- rf.prob == test$Labels
  table(agreement)
  prop.table(table(agreement))
  
  #computing the true possitive rate
  TPR=as.numeric(table(test$Labels,rf.prob)[2,2]/table(test$Labels)[2])
  #computing the false negative rate
  FNR=1-TPR
  #computing the true negative rate
  TNR=as.numeric(table(test$Labels,rf.prob)[1,1]/table(test$Labels)[1])
  #computing the false positive rate
  FPR=1-TNR
  r=c(prop.table(table(agreement))[[2]],TNR,TPR)
  return(r)
}

