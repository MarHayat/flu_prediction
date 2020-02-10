# Predicting the short-term success of human influenza virus variants with machine learning
This repository contains the material for the paper:

  . Maryam Hayati, Priscila Biller, and Caroline Colijn. Predicting the short-term successof human influenza a variants with machine learning.bioRxivdoi:10.1101/609248,2019.
  
  In this paper we use longitudinally sampled phylogenetic trees based on hemagglutinin sequences from human influenza viruses, together with counts of epitope site polymorphisms in hemagglutinin,  to predict which influenza virus strains are likely to be successful.  We extract small groups of taxa (subtrees) and use a suite of features of these subtrees as key inputs to the machine learning tools. Using a range of training and testing strategies, including training on H3N2 and testing on H1N1, we find that successful prediction of future expansion of small subtrees is possible from these data, with accuracies of 0.71-0.85 and a classifier 'area under the curve' (AUC) 0.75-0.9.
  
  
  flutreeH1N1_RLB.Rdata, flutree2018-5.Rdata, InfluenzaB.Rdata are the trees reconstructed from H1N1 strains, H3N2 strains and B strains. 
  
  df_H1N1.csv, df_2018-5.csv and df_B.csv contain the features of H1N1 tree, H3N2 tree and B tree respectively.
  
  Here, we explain the procedure to classify the subtrees of H3N2 tree and it is the same for the other types:
  
  `data=read.csv("~/df_2018-5.csv",sep= ",",header=T,stringsAsFactors=FALSE)
  load("~/flutree2018-5.Rdata")`
  
  "readData" function read the dataset and do some preproseeing on the data:
  
  `df=readData(data,timeFrame=3.4,alpha=1.1,tree,changeLabels=TRUE)`
    Here alpha is the growth ratio, timeFrame is the time from the root of the tree that we prune the future tips and we set changeLabels equal to TRUE in some cases when the AUC is less than 0.50.
    
   #choose the train and test data
  `set.seed(123)
  DataTT=PrepareData(df,RemoveOutliers=TRUE,TT=TRUE)
  train=DataTT[[1]]
  test=DataTT[[2]]`
  
  #train the model using the best hyperparameters
  `svm.fit = svm(data = train, Labels ~ .,
              kernel ="linear", degree = 3, gamma =   0.03125 , 
              coef0 = 0, cost =32, nu = 0.5,class.weigth=c("0"=0.50,"1"=0.50))`
              
  `svm.prob <- predict(svm.fit, newdata = test)
  summary(svm.prob)
  agreement <- svm.prob == test$Labels
  acc=length(which(svm.prob == test$Labels))/length(test$Labels)
  #roc analysis for test data
  svmmodel.predict<-predict(svm.fit, newdata = test,decision.values=TRUE)
  svmmodel.probs<-attr(svmmodel.predict,"decision.values")
  svmmodel.class<-predict(svm.fit,test,type="class")
  svmmodel.labels<-test$Labels
  svmmodel.prediction<-prediction(svmmodel.probs,svmmodel.labels)
  svmmodel.performance<-performance(svmmodel.prediction,"tpr","fpr")
  svmmodel.auc<-performance(svmmodel.prediction,"auc")@y.values[[1]]
  round(svmmodel.auc,2)`
  
