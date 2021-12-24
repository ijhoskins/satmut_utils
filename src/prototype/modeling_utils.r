# Modeling utils

library(leaps)
library(caret)
library(e1071)
library(class)
library(randomForest)
library(gbm)
library(glmnet)
library(mltools)


# Factors in R should have these levels for the following functions to work
cont_table_factor_levels =c("FALSE", "TRUE")

#' Creates a 2x2 contigency table
#'
#' @param truth logical or character vector
#' @param predicted logical or character vector
#' @return table 2x2 truth by predicted
create_cont_table<- function(truth, predicted){
  res<- table(truth, predicted, dnn=c("Truth", "Predicted"))
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

get_recall<- get_sensitivity

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


#' Prepares a training dataset for modeling by coercing types
#' 
#' @param train_dt data.table training data
#' @param truth_dt data.table containing truth variants
#' @return data.table 
#'
prep_train_dt<- function(train_dt, truth_dt){
  
  copy_train_dt<- copy(train_dt)
  copy_truth_dt<- copy(truth_dt)
  
  copy_truth_dt[,CAO:=as.integer(sapply(strsplit(
    sapply(strsplit(INFO, ";"), "[", 1), "=", fixed=T), "[", 2))]
  
  copy_truth_dt[,CAF:=as.numeric(sapply(strsplit(
    sapply(strsplit(INFO, ";"), "[", 2), "=", fixed=T), "[", 2))]
  
  copy_truth_dt[,VAR_ID:=paste(`#CHROM`, POS, REF, ALT, sep=":")]
  
  copy_train_dt[,VAR_ID:=paste(`#CHROM`, POS, REF, ALT, sep=":")]
  copy_train_dt[,log10_CAF:=log10(CAF)]
  copy_train_dt[,TYPE:=get_var_type(REF, ALT), by=seq_len(nrow(copy_train_dt))]
  copy_train_dt[,TYPE:=factor(TYPE, levels=c("SNP", "di_nt_MNP", "tri_nt_MNP"))]
  copy_train_dt[,UP_REF_NT:=as.factor(UP_REF_NT)]
  copy_train_dt[,DOWN_REF_NT:=as.factor(DOWN_REF_NT)]
  
  # Generate the substitution feature
  copy_train_dt[,SUBST:=paste(REF_NT, ALT_NT, sep=":")]
  copy_train_dt[SUBST%in%c("A:C", "T:G"), SUBST:="A.C_T.G"]
  copy_train_dt[SUBST%in%c("A:G", "T:C"), SUBST:="A.G_T.C"]
  copy_train_dt[SUBST%in%c("A:T", "T:A"), SUBST:="A.T_T.A"]
  copy_train_dt[SUBST%in%c("C:A", "G:T"), SUBST:="C.A_G.T"]
  copy_train_dt[SUBST%in%c("C:G", "G:C"), SUBST:="C.G_G.C"]
  copy_train_dt[SUBST%in%c("C:T", "G:A"), SUBST:="C.T_G.A"]
  copy_train_dt[,SUBST:=as.factor(SUBST)]
  
  # For variants spanning codons, there will be more than one bool
  # For simplicity and to reduce the number of factor levels which
  # will be one-hot encoded in feature selection/modeling, just set
  # to either True or False using "or" operation
  copy_train_dt[grepl("False", MATCHES_MUT_SIG), MATCHES_MUT_SIG:="False"]
  copy_train_dt[grepl("True", MATCHES_MUT_SIG), MATCHES_MUT_SIG:="True"]
  copy_train_dt[,MATCHES_MUT_SIG:=as.factor(MATCHES_MUT_SIG)]
  
  # Some variants may have not been edited; remove them from final truth dataset
  copy_truth_filt<- copy_truth_dt[CAO>0]
  
  copy_train_dt[,Truth:=FALSE]
  copy_train_dt[VAR_ID%in%copy_truth_filt[,VAR_ID], Truth:=TRUE]
  copy_train_dt[,Truth:=factor(Truth, levels=c("FALSE", "TRUE"))]
  
  # Remove any outlier variants with higher frequency (background or systematic errors)
  copy_train_filt_dt<- copy_train_dt[CAF < 0.3]
  
  # Standardize the numeric data; don't worry about this for now as relevant models 
  # typically do this internally
  #copy_train_filt_std_dt<- scale_numeric_features(copy_train_filt_dt, standardize=TRUE)
  
  return(copy_train_filt_dt)
}


#' Hyperparameter tuning for several models using resampling
#' 
#' @param in_dt data.table input data (should contain both train, test) with Truth column
#' @param model character one of c("elasticnet", "knn", "rf", "gbm", "svm"). Default rf.
#' @param nvmin integer min number of variables to tune for knn model. Default 3.
#' @param nvmax integer max number of variables to tune for knn model Default 10.
#' @param prob_cutoff numeric probability cutoff to use for binary classification. Default 0.5.
#' @return train object or character for method="knn", where first element is train object, second is acc for nvmin:nvmax, and third is chosen nvar
hyperparam_tune<- function(in_dt, model="rf", nvmin=3, nvmax=10, prob_cutoff=0.5){
  
  if(model=="elasticnet"){
    
    elasticnet_cv_res<- train(Truth~., in_dt, method="glmnet", family="binomial")
    return(elasticnet_cv_res)
    
  } else if(model=="knn"){
    
    knn_cv_res<- train(Truth~., in_dt, method="knn", tuneGrid=expand.grid(
      k=seq(3,5,7,9)))
    
    # Need to also run manual CV to tune number of features
    knn_kfolds<- create_kfolds(in_dt, kfolds=5)
    knn_train_dt<- in_dt[!knn_kfolds[[1]]]
    knn_test_dt<- in_dt[knn_kfolds[[1]]]
    
    # Use best subset selection to pick features
    model_accs<- vector(mode="numeric", length=nvmax-nvmin+1)
    nvars<- nvmin:nvmax 
    names(model_accs)<- as.character(nvars)
      
    for(i in 1:length(nvars)){
      
      # Use N features in best subset selection
      best_subset_res<- regsubsets(Truth~., knn_train_dt, nvmax=nvars[i])
      best_subset_res_summary<- summary(best_subset_res)
      select_features_counts<- apply(best_subset_res_summary$which[,-1], 2, sum)
      select_features<- names(select_features_counts)[select_features_counts>0]
      select_features<- c(select_features, "Truth")
      train_select_dt<- knn_train_dt[,select_features, with=FALSE]
      
      knn_model<- knn3(Truth~., train_select_dt, k=knn_cv_res$bestTune$k)
      model_predict<- predict(knn_model, knn_test_dt)
      prediction_labels<- model_predict[,2]>=prob_cutoff
      contigency_table<- table(
        knn_test_dt[,Truth], prediction_labels, dnn=c("Truth", "Predicted"))
      model_acc<- get_accuracy(contigency_table)
      model_accs[i]<- model_acc
    }
    
    # Once we have determined accuracy for the selected model for each value
    # of nvars between nvmin and nvmax, select a reasonably low nvars value
    prev_acc<- 0
    for(i in 1:length(model_accs)){
      if(!(model_accs[i] - prev_acc) > 0.01){
        select_nvars<- nvars[i]
        break
      } else {
        select_nvars<- nvars[i]
        prev_acc<- model_accs[i]
      }
    }
    return(list(knn_cv_res, model_accs, select_nvars))
    
  } else if(model=="rf"){
    
    rf_cv_res<- train(Truth~., in_dt, method="rf", tuneGrid=expand.grid(
      mtry=seq(nvmin, nvmax), ntree=c(250,500,1000)))
    
    return(rf_cv_res)
    
  } else if(model=="gbm"){
    
    # ntrees should not be too high for GBM as this can lead to over-fitting
    # also low values for interaction.depth (stumps) typically work well
    gbm_cv_res<- train(Truth~., in_dt, method="gbm", tuneGrid=expand.grid(
      interaction.depth=seq(1,5), ntrees=c(50,100,200)))
    
    return(gbm_cv_res)
    
  } else if(model=="svm"){
    
    # Tune different costs
    svm_cv_res<- train(Truth~., in_dt, method="svmLinear", 
                       tuneGrid=expand.grid(C=seq(0, 2, length=10)))
    return(svm_cv_res)
    
  }
}

#' Hyperparameter tuning of nvars/nfeatures for several models
#' 
#' @param in_dt data.table input data (should contain both train, test) with Truth column
#' @param model character one of c("elasticnet", "knn", "rf", "gbm", "svm). Default rf.
#' @param prob_cutoff numeric probability cutoff to use for binary classification. Default 0.5.
#' @param kfolds integer number of folds for CV. Default 10. Uses one fold training set for tuning.
#' @return character vector of features to use, selected with best subset for optimal nvars/nfeatures
run_nested_cv<- function(in_dt,  model="rf", prob_cutoff=0.5, kfolds=10){
  
  # Train hyper-parameters
  set.seed(9)
  hyper_indices<- sample(1:nrow(in_dt), round(nrow(in_dt)*0.05))
  hyperparam_train_dt<- in_dt[hyper_indices,]
  
  if(model=="elasticnet"){
    elasticnet_cv_res<- hyperparam_tune(hyperparam_train_dt, model="elasticnet")
  } else if(model=="knn"){
    knn_cv_res<- hyperparam_tune(hyperparam_train_dt, model="knn")
  } else if(model=="rf"){
    rf_cv_res<- hyperparam_tune(hyperparam_train_dt, model="rf")
  } else if(model=="gbm"){
    gbm_cv_res<- hyperparam_tune(hyperparam_train_dt, model="gbm")
  } else if(model=="svm"){
    svm_cv_res<- hyperparam_tune(hyperparam_train_dt, model="svm")
  } else {
    stop("Unsupported model.")
  }
  
  cv_dt<- in_dt[!hyper_indices,]
  kfolds_list<- create_kfolds(cv_dt, kfolds)
  
  # Now test the models
  # Omit the first fold because we use it for hyper-parameter
  cv_perf<- matrix(nrow=kfolds, ncol=3)
  for(i in 1:kfolds){
    
    train_dt<- cv_dt[!kfolds_list[[i]],]
    test_dt<- cv_dt[kfolds_list[[i]],]
    
    if(model=="elasticnet"){
      
      elasticnet_model<- glmnet(x=as.matrix(as.data.frame(train_dt[,-c("Truth")])), 
                           y=train_dt[,Truth], family="binomial",
                           alpha=elasticnet_cv_res$bestTune$alpha, 
                           lambda=elasticnet_cv_res$bestTune$lambda)
      
      model_predict<- predict(elasticnet_model, as.matrix(as.data.frame(test_dt[,-c("Truth")])))
      prediction_labels<- factor(model_predict>=prob_cutoff, levels=c("FALSE", "TRUE"))
      
    } else if(model=="knn"){
      
      # After choosing nvars, select the best features
      best_subset_res<- regsubsets(Truth~., train_dt, nvmax=knn_cv_res[[3]])
      best_subset_res_summary<- summary(best_subset_res)
      select_features_counts<- apply(best_subset_res_summary$which[,-1], 2, sum)
      select_features<- names(select_features_counts)[select_features_counts>0]
      select_features<- c(select_features, "Truth")
      train_select_dt<- train_dt[,..select_features]
      
      knn_model<- knn3(Truth~., train_select_dt, k=knn_cv_res[[1]]$bestTune$k)
      model_predict<- predict(knn_model, test_dt, type="prob")
      prediction_labels<- factor(model_predict[,2]>=prob_cutoff, levels=c("FALSE", "TRUE"))
      
    } else if(model=="rf"){
      
      # Consider tuning mtrees parameter
      rf_model<- randomForest(Truth~., train_dt, mtry=rf_cv_res$bestTune$mtry)
      model_predict<- predict(rf_model, test_dt, type="prob")
      prediction_labels<- factor(model_predict[,2]>=prob_cutoff, levels=c("FALSE", "TRUE"))
      
    } else if(model=="gbm"){
      
      # The outcome needs to be 0 or 1
      train_dt[,Truth_int:=as.integer(as.logical(Truth))]
      test_dt[,Truth_int:=as.integer(as.logical(Truth))]
      
      gbm_model<- gbm(Truth_int~., train_dt, distribution="bernoulli", 
                      n.trees=gbm_cv_res$bestTune$n.trees, 
                      interaction.depth=gbm_cv_res$bestTune$interaction.depth, 
                      shrinkage=gbm_cv_res$bestTune$shrinkage,
                      n.minobsinnode=gbm_cv_res$bestTune$n.minobsinnode)
      
      model_predict<- predict(gbm_model, test_dt, type="response", 
                              n.trees=gbm_cv_res$bestTune$n.trees)
      print(head(model_predict))
      prediction_labels<- factor(model_predict>=prob_cutoff, levels=c("FALSE", "TRUE"))
      
    } else if(model=="svm"){
      
      # This is a "support vector classifier" with linear kernel, but use SVM naming
      svm_model<- svm(Truth~., train_dt, kernel="linear", cost=svm_cv_res$bestTune$C)
      
      # For some reason couldn't get probability=TRUE to give me probs, even with type="prob"
      model_predict<- predict(svm_model, test_dt)
      prediction_labels<- factor(model_predict, levels=c("FALSE", "TRUE"))
      
    } else {
      stop("Unsupported model.")
    }
    
    print(head(prediction_labels))
    print(test_dt[,head(Truth)])
    
    contigency_table<- table(test_dt[,Truth], prediction_labels, dnn=c("Truth", "Predicted"))
    cv_perf[i,1]<- get_accuracy(contigency_table)
    cv_perf[i,2]<- get_precision(contigency_table)
    cv_perf[i,3]<- get_recall(contigency_table)
  }
  colnames(cv_perf)<- c("Accuracy", "Precision", "Recall")
  return(cv_perf)
}


#' Hyperparameter tuning of nvars/nfeatures for several models
#' 
#' @param in_dt data.table input data (should contain both train, test) with Truth column
#' @param model character one of c("elasticnet", "knn", "rf", "gbm", "svm). Default rf.
#' @param prob_cutoff numeric probability cutoff to use for binary classification. Default 0.5.
#' @param kfolds integer number of folds for CV. Default 10. Uses one fold training set for tuning.
#' @return character vector of features to use, selected with best subset for optimal nvars/nfeatures
run_nested_cv_complete<- function(in_dt,  model="rf", prob_cutoff=0.5, kfolds=10){
  
  set.seed(9)
  cv_dt<- copy(in_dt)
  kfolds_list<- create_kfolds(cv_dt, kfolds)
  
  # Test models
  cv_perf<- matrix(nrow=kfolds, ncol=3)
  for(i in 1:kfolds){
    
    train_dt<- cv_dt[!kfolds_list[[i]],]
    test_dt<- cv_dt[kfolds_list[[i]],]
    
    if(model=="elasticnet"){
      
      elasticnet_cv_res<- hyperparam_tune(train_dt, model="elasticnet")
      print(elasticnet_cv_res$bestTune)
      
      elasticnet_model<- glmnet(x=as.matrix(as.data.frame(train_dt[,-c("Truth")])), 
                                y=train_dt[,Truth], family="binomial",
                                alpha=elasticnet_cv_res$bestTune$alpha, 
                                lambda=elasticnet_cv_res$bestTune$lambda)
      
      model_predict<- predict(elasticnet_model, as.matrix(as.data.frame(test_dt[,-c("Truth")])))
      prediction_labels<- factor(model_predict>=prob_cutoff, levels=c("FALSE", "TRUE"))
      
    } else if(model=="knn"){
      
      knn_cv_res<- hyperparam_tune(train_dt, model="knn")
      print(knn_cv_res[[1]]$bestTune$k)
      
      # After choosing nvars, select the best features
      best_subset_res<- regsubsets(Truth~., train_dt, nvmax=knn_cv_res[[3]])
      best_subset_res_summary<- summary(best_subset_res)
      select_features_counts<- apply(best_subset_res_summary$which[,-1], 2, sum)
      select_features<- names(select_features_counts)[select_features_counts>0]
      select_features<- c(select_features, "Truth")
      train_select_dt<- train_dt[,..select_features]
      
      knn_model<- knn3(Truth~., train_select_dt, k=knn_cv_res[[1]]$bestTune$k)
      model_predict<- predict(knn_model, test_dt, type="prob")
      prediction_labels<- factor(model_predict[,2]>=prob_cutoff, levels=c("FALSE", "TRUE"))
      
    } else if(model=="rf"){
      
      rf_cv_res<- hyperparam_tune(train_dt, model="rf")
      print(rf_cv_res$bestTune)
      
      # Consider tuning mtrees parameter
      rf_model<- randomForest(
        Truth~., train_dt, mtry=rf_cv_res$bestTune$mtry, ntree=rf_cv_res$bestTune$ntree)
      
      model_predict<- predict(rf_model, test_dt, type="prob")
      prediction_labels<- factor(model_predict[,2]>=prob_cutoff, levels=c("FALSE", "TRUE"))
      
    } else if(model=="gbm"){
      
      gbm_cv_res<- hyperparam_tune(train_dt, model="gbm")
      print(gbm_cv_res$bestTune)
      
      # The outcome needs to be 0 or 1
      train_dt[,Truth_int:=as.integer(as.logical(Truth))]
      test_dt[,Truth_int:=as.integer(as.logical(Truth))]
      
      gbm_model<- gbm(Truth_int~., train_dt, distribution="bernoulli", 
                      n.trees=gbm_cv_res$bestTune$n.trees, 
                      interaction.depth=gbm_cv_res$bestTune$interaction.depth, 
                      shrinkage=gbm_cv_res$bestTune$shrinkage,
                      n.minobsinnode=gbm_cv_res$bestTune$n.minobsinnode)
      
      model_predict<- predict(gbm_model, test_dt, type="response", 
                              n.trees=gbm_cv_res$bestTune$n.trees)

      prediction_labels<- factor(model_predict>=prob_cutoff, levels=c("FALSE", "TRUE"))
      
    } else if(model=="svm"){
      
      svm_cv_res<- hyperparam_tune(train_dt, model="svm")
      print(svm_cv_res$bestTune)
      
      # This is a "support vector classifier" with linear kernel, but use SVM naming
      svm_model<- svm(Truth~., train_dt, kernel="linear", cost=svm_cv_res$bestTune$C)
      
      # For some reason couldn't get probability=TRUE to give me probs, even with type="prob"
      model_predict<- predict(svm_model, test_dt)
      prediction_labels<- factor(model_predict, levels=c("FALSE", "TRUE"))
      
    } else {
      stop("Unsupported model.")
    }
    
    contigency_table<- table(test_dt[,Truth], prediction_labels, dnn=c("Truth", "Predicted"))
    cv_perf[i,1]<- get_accuracy(contigency_table)
    cv_perf[i,2]<- get_precision(contigency_table)
    cv_perf[i,3]<- get_recall(contigency_table)
  }
  colnames(cv_perf)<- c("Accuracy", "Precision", "Recall")
  return(cv_perf)
}

