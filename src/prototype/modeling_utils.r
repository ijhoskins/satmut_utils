# Modeling utils


#' Creates a 2x2 contigency table
#'
#' @param truth logical or character vector
#' @param predicted logical or character vector
#' @return table 2x2 truth by predicted
create_cont_table<- function(truth, predicted){
  res<- table(truth,predicted, dnn=c("Truth", "Predicted"))
  return(res)
}

#' Gets the model accuracy from a 2x2 table
#'
#' @param cont_table table of counts
#' @return numeric prediction accuracy
get_accuracy<- function(cont_table){
  res<- sum(diag(cont_table))/sum(cont_table)
  return(res)
}

#' Gets the model sensitivity from a 2x2 table
#'
#' @param cont_table table of counts
#' @return numeric sensitivity (recall)
get_sensitivity<- function(cont_table){
  res<- cont_table[2,2]/sum(cont_table[2,])
  return(res)
}

#' Gets the model specificity from a 2x2 table
#'
#' @param cont_table table of counts
#' @return numeric specificity (TN/FN+TN)
get_specificity<- function(cont_table){
  res<- cont_table[1,1]/sum(cont_table[1,])
  return(res)
}

#' Gets the model precision from a 2x2 table
#'
#' @param cont_table table of counts
#' @return numeric precision (TP/TP+FP)
get_precision<- function(cont_table){
  res<- cont_table[2,2]/sum(cont_table[,2])
  return(res)
}


#' Gets the F score from precision and recall
#'
#' @param prec precision value
#' @param rec recall value
#' @return numeric F score
get_fscore<- function(prec, rec){
  res<- 2 * prec*rec / (prec+rec)
  return(res)
}

#' Gets the F score from a 2x2 table
#'
#' @param cont_table table of counts
#' @return numeric F score
get_fscore_from_table<- function(cont_table){
  prec<- get_precision(cont_table)
  rec<- get_sensitivity(cont_table)
  res<- get_fscore(prec, rec)
  return(res)
}

#' Creates a PR curve used test and logistic model
#'
#' @param test_dt data.table of test data and selected features
#' @param logistic_model glm lm fitted object
#' @param probs numeric vector of prob cutoffs
#' @return data.table raw values for plotting
get_logistic_pr_curve_dt<- function(test_dt, logistic_model, probs=seq(0.0001,0.9999,0.0001)){
  
  logistic_predictions<- predict(
    logistic_model, test_dt, type="response")
  
  probs_len<- length(probs)
  
  res_dt<- data.table(probs=probs, 
                      precision=vector("numeric", probs_len),
                      recall=vector("numeric", probs_len))
  
  for(i in 1:probs_len){
    
    logistic_prediction_labels<- ifelse(
      logistic_predictions>=probs[i], TRUE, FALSE)
    
    logistic_contigency<- table(
      test_dt[,Truth],
      logistic_prediction_labels,
      dnn=c("Truth", "Predicted"))
    
    precision_res<- get_precision(logistic_contigency)
    recall_res<- get_sensitivity(logistic_contigency)
    
    res_dt[i,precision:=precision_res]
    res_dt[i,recall:=recall_res]
  }
  return(res_dt)
}


#' Runs kNN-CV to determine best k
#' 
#' @param in_dt data.table table of modeling features; can include factors
#' @return train caret::knn.cv result
knn_cv<- function(in_dt){
  
  copy_dt<- copy(in_dt)
  
  copy_dt[,Truth:=as.factor(as.character(Truth))]
  
  train.control <- trainControl(method="repeatedcv", 
                                number=10, repeats=5)
  
  knn_cv_res<- train(form=Truth~., data=copy_dt, method="knn",
                  trControl=train.control)
  
  return(knn_cv_res)
}

#' Creates PR curves for a kNN model
#' 
#' @param test_dt data.table of test data and selected features
#' @param knn_model knn model for predictions
#' @param probs numeric vector of prob cutoffs
#' @return data.table raw values for plotting
get_knn_pr_curve_dt<- function(test_dt, knn_model, probs=seq(0.0001,0.9999,0.0001)){
  
  probs_len<- length(probs)
  
  res_dt<- data.table(probs=probs, 
                      precision=vector("numeric", probs_len),
                      recall=vector("numeric", probs_len))
  
  test_dt_no_y<- test_dt[,.SD,.SDcols=!c("Truth")]
    
  # This performs knn.cv
  knn_predictions<- predict(
    knn_model, test_dt_no_y, type="prob")
  
  for(j in 1:probs_len){
    
    knn_prediction_labels<- knn_predictions[,2]>=probs[j]
    
    knn_contigency<- table(
      test_dt[,Truth],
      knn_prediction_labels,
      dnn=c("Truth", "Predicted"))
    
    precision_res<- get_precision(knn_contigency)
    recall_res<- get_sensitivity(knn_contigency)
    
    res_dt[j,precision:=precision_res]
    res_dt[j,recall:=recall_res]
  }
  return(res_dt)
}


#' Creates PR curves for a random forest
#' 
#' @param test_dt data.table of test data and selected features
#' @param rf_model randomForest model for predictions
#' @param probs numeric vector of prob cutoffs
#' @return data.table raw values for plotting
get_rf_pr_curve_dt<- function(test_dt, rf_model, probs=seq(0.0001,0.9999,0.0001)){
  
  probs_len<- length(probs)
  
  res_dt<- data.table(probs=probs, 
                      precision=vector("numeric", probs_len),
                      recall=vector("numeric", probs_len))
  
  rf_predictions<- predict(
    rf_model, test_dt, type="prob")
  
  for(j in 1:probs_len){
    
    rf_prediction_labels<- rf_predictions[,2]>=probs[j]
    
    rf_contigency<- table(
      test_dt[,Truth],
      rf_prediction_labels,
      dnn=c("Truth", "Predicted"))
    
    precision_res<- get_precision(rf_contigency)
    recall_res<- get_sensitivity(rf_contigency)
    
    res_dt[j,precision:=precision_res]
    res_dt[j,recall:=recall_res]
  }
  
  res_dt<- rbind(res_dt, data.table(
    probs=c(0,1), precision=c(0,1), recall=c(1,0)))

  return(res_dt)
}


#' Creates PR curves for a boosted model
#' 
#' @param test_dt data.table of test data and selected features
#' @param gbm_model gbm model for predictions
#' @param probs numeric vector of prob cutoffs
#' @return data.table raw values for plotting
get_gbm_pr_curve_dt<- function(test_dt, gbm_model, ntrees=500, probs=seq(0.0001,0.9999,0.0001)){
  
  probs_len<- length(probs)
  
  res_dt<- data.table(probs=probs, 
                      precision=vector("numeric", probs_len),
                      recall=vector("numeric", probs_len))
  
  gbm_predictions<- predict(
    gbm_model, test_dt, n.trees=ntrees, type="response")
  
  for(j in 1:probs_len){
    
    gbm_prediction_labels<- gbm_predictions>=probs[j]
    
    gbm_contigency<- table(
      test_dt[,Truth],
      gbm_prediction_labels,
      dnn=c("Truth", "Predicted"))
    
    precision_res<- get_precision(gbm_contigency)
    recall_res<- get_sensitivity(gbm_contigency)
    
    res_dt[j,precision:=precision_res]
    res_dt[j,recall:=recall_res]
  }
  return(res_dt)
}

#' Creates k-folds of shuffled the training dt
#' 
#' @param train_dt data.table training data
#' @param kfolds integer number folds to use, ideally 5 or 10
#' @return list of indices for each fold
create_kfolds<- function(train_dt, kfolds=10){
  
  dt_nrows<- nrow(train_dt)
  fold_size<- floor(dt_nrows / kfolds)
  
  set.seed(9)
  shuffled_indices<- sample(1:dt_nrows, dt_nrows, replace=FALSE)
  
  fold_indices<- vector(mode="list", length=kfolds)
  for(i in 1:kfolds){
    if(i==1){
      fold_indices[[i]]<- shuffled_indices[1:fold_size]
    } else {
      fold_indices[[i]]<- shuffled_indices[
        ((i-1)*fold_size + 1):(i*fold_size)]
    }
  }
  return(fold_indices)
}


#' Runs k-fold CV to tune mtry of a RF model
#' 
#' @param train_dt data.table training data
#' @param mtrys integer vector of mtry values to test
#' @param kfolds integer number folds to use, ideally 5 or 10
#' @param samp_sizes integer vector of sample size to draw from each group
#' @return data.table results of the CV for each value of mtry
tune_rf_mtry<- function(train_dt, mtrys=seq(2,5,1), kfolds=10, samp_sizes=c(2300,2300)){
  
  kfold_indices<- create_kfolds(train_dt, kfolds)
  n_mtrys<- length(mtrys)
  res_dt<- data.table(mtry=rep(mtrys, each=kfolds), 
                      kfold=rep(1:kfolds, times=n_mtrys),
                      accuracy=vector(mode="numeric", length=kfolds*n_mtrys),
                      precision=vector(mode="numeric", length=kfolds*n_mtrys),
                      recall=vector(mode="numeric", length=kfolds*n_mtrys),
                      specificity=vector(mode="numeric", length=kfolds*n_mtrys))

  counter<- 1
  for(i in 1:n_mtrys){
    for(j in 1:kfolds){
      
      test_indices<- kfold_indices[[j]]
      rf_train_dt<- train_dt[!test_indices,]
      rf_test_dt<- train_dt[test_indices,]
      
      rf_model<- randomForest(
        Truth~., rf_train_dt, mtry=mtrys[i], 
        importance=TRUE, sampsize=samp_sizes)
      
      rf_predictions<- predict(
        rf_model, rf_test_dt, 
        type="prob")
      
      rf_test_dt[,Truth_pred:=ifelse(
        rf_predictions[,2]>=0.5, TRUE, FALSE)]
      
      rf_test_contigency<- create_cont_table(
        rf_test_dt[,Truth], 
        rf_test_dt[,Truth_pred])
      
      res_dt[counter, `:=`(accuracy=get_accuracy(rf_test_contigency),
                           precision=get_precision(rf_test_contigency),
                           recall=get_sensitivity(rf_test_contigency),
                           specificity=get_specificity(rf_test_contigency))]
      
      counter<- counter+1
    }
  }
  res_dt[,mtry:=as.factor(mtry)]
  return(res_dt)
}


