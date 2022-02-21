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


#' Merges a simulated truth dataset with a satmut_utils call dataset
#' 
#' @param train_dt data.table training data
#' @param truth_dt data.table containing truth variants
#' @return data.table modeling input
merge_truth_dt<- function(train_dt, truth_dt){
  
  copy_train_dt<- copy(train_dt)
  copy_train_dt[,VAR_ID:=paste(`#CHROM`, POS, REF, ALT, sep=":")]
    
  copy_truth_dt<- copy(truth_dt)
  copy_truth_dt[,VAR_ID:=paste(`#CHROM`, POS, REF, ALT, sep=":")]
  
  copy_truth_dt[,CAO:=as.integer(sapply(strsplit(
    sapply(strsplit(INFO, ";"), "[", 1), "=", fixed=T), "[", 2))]
  
  copy_truth_dt[,CAF:=as.numeric(sapply(strsplit(
    sapply(strsplit(INFO, ";"), "[", 2), "=", fixed=T), "[", 2))]
  
  # Some variants may have not been edited; remove them from final truth dataset
  copy_truth_filt_dt<- copy_truth_dt[CAO>0]
  copy_truth_filt_dt[,log10_CAF:=log10(CAF)]
  
  copy_train_postproc_dt<- postprocess_vcf_summary_dt(
    in_dt=copy_train_dt, mapping_dt=NULL, tile_positions_dt=NULL, 
    key_features=c("VAR_ID"))
  
  copy_train_postproc_merge_dt<- merge(
    copy_truth_filt_dt[,.(VAR_ID, CAO, CAF, log10_CAF)], 
    copy_train_postproc_dt, by=c("VAR_ID"), all=TRUE)
  
  copy_train_postproc_merge_dt[,Truth:="FALSE"]
  copy_train_postproc_merge_dt[!is.na(CAO.x), Truth:="TRUE"]
  copy_train_postproc_merge_dt[,Truth:=factor(Truth, levels=cont_table_factor_levels)]
  
  return(copy_train_postproc_merge_dt)
  
}


#' Prepares a training dataset for modeling by coercing types
#' 
#' @param train_dt data.table training data
#' @param character modeling_features features to select
#' @param truth_dt data.table containing truth variants, or NULL if a dataset to filter
#' @param caf_thresh numeric maximum frequency threshold
#' @return data.table modeling input
prep_train_dt<- function(train_dt, modeling_features, truth_dt=NULL, caf_thresh=0.3){
  
  copy_train_dt<- copy(train_dt)
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
  
  # Keep track of VAR_ID and TYPE after one-hot encoding the factors
  full_features<- modeling_features
  
  if(!"Sample"%in%modeling_features){
    full_features<- c(full_features, "Sample")
  }
  
  if(!"VAR_ID"%in%modeling_features){
    full_features<- c(full_features, "VAR_ID")
  }
  
  if(!"TYPE"%in%modeling_features){
    full_features<- c(full_features, "TYPE")
  }
  
  # Filter observed variants based on region and frequency
  # Remove any outlier variants with higher frequency (background or systematic errors)
  # Ensure all modeled variants are in the coding region so they have valid MATCHES_MUT_SIG
  copy_train_filt_dt<- copy_train_dt[LOCATION=="CDS" & CAF < 0.3, ..full_features]
  
  if(!is.null(truth_dt)){
    
    copy_truth_dt<- copy(truth_dt)
    
    copy_truth_dt[,CAO:=as.integer(sapply(strsplit(
      sapply(strsplit(INFO, ";"), "[", 1), "=", fixed=T), "[", 2))]
    
    copy_truth_dt[,VAR_ID:=paste(`#CHROM`, POS, REF, ALT, sep=":")]
    
    # Some variants may have not been edited; remove them from final truth dataset
    copy_truth_filt_dt<- copy_truth_dt[CAO>0]
    
    copy_train_filt_dt[,Truth:=FALSE]
    copy_train_filt_dt[VAR_ID%in%copy_truth_filt_dt[,VAR_ID], Truth:=TRUE]
    copy_train_filt_dt[,Truth:=factor(Truth, levels=cont_table_factor_levels)]
    
    # Use model matrix to one-hot encode factors
    copy_train_filt_mm<- model.matrix(Truth~., data=copy_train_filt_dt[,-c("VAR_ID", "Sample")])
    
    # Convert back to a data.table
    copy_train_filt_input_dt<- as.data.table(copy_train_filt_mm[,-1])
    copy_train_filt_input_dt[,Truth:=copy_train_filt_dt[,Truth]]
    
    return(copy_train_filt_input_dt)
  }
  
  # Otherwise we have a dataset to filter and there will be no Truth column
  # Use model matrix to one-hot encode factors
  copy_train_filt_mm<- model.matrix(~., data=copy_train_filt_dt[,-c("VAR_ID", "Sample")])
  
  # Convert back to a data.table and remove the intercept
  copy_train_filt_input_dt<- cbind(as.data.table(
    copy_train_filt_mm[,-1]), copy_train_filt_dt[,.(VAR_ID, TYPE, Sample)])
  
  return(copy_train_filt_input_dt)
  
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
      k=c(3,5,7,9)))
    
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
    
    # Would need to tune ntrees parameter manually
    # for now just use default
    rf_cv_res<- train(Truth~., in_dt, method="rf", tuneGrid=expand.grid(mtry=seq(nvmin, nvmax)), 
                      na.action=na.roughfix)
    
    return(rf_cv_res)
    
  } else if(model=="gbm"){
    
    # ntrees should not be too high for GBM as this can lead to over-fitting
    # also low values for interaction.depth (stumps) typically work well
   
    gbm_cv_res<- train(Truth~., in_dt, method="gbm", tuneGrid=expand.grid(
      interaction.depth=seq(1,5), n.trees=c(50,100), 
      shrinkage=c(0.01, 0.1, 0.2), n.minobsinnode=10), verbose=FALSE)
    
    return(gbm_cv_res)
    
  } else if(model=="svm"){
    
    # Tune different costs
    svm_cv_res<- train(Truth~., in_dt, method="svmLinear", 
                       tuneGrid=expand.grid(C=c(0.001, 0.01, 0.1, 1, 10)))
    return(svm_cv_res)
    
  }
}


#' Hyperparameter tuning of nvars/nfeatures for several models
#' 
#' @param in_dt data.table input data (should contain both train, test) with Truth column
#' @param model character one of c("elasticnet", "knn", "rf", "gbm", "svm"). Default rf.
#' @param prob_cutoff numeric probability cutoff to use for binary classification. Default 0.5.
#' @param kfolds integer number of folds for CV. Default 10. Uses one fold training set for tuning.
#' @param hyperparam_prop numeric proportion of each fold's training dataset to use for hyperparam tuning
#' @return character vector of features to use, selected with best subset for optimal nvars/nfeatures
run_nested_cv<- function(in_dt,  model="rf", prob_cutoff=0.5, kfolds=10, hyperparam_prop=0.2){
  
  set.seed(9)
  cv_dt<- copy(in_dt)
  kfolds_list<- create_kfolds(cv_dt, kfolds)
  
  # Test models and keep track of performance and hyperparams
  cv_perf<- matrix(nrow=kfolds, ncol=3)
  cv_hyperparams<- NULL
  
  for(i in 1:kfolds){
    
    train_dt<- cv_dt[!kfolds_list[[i]],]
    test_dt<- cv_dt[kfolds_list[[i]],]
    
    train_nrows<- nrow(train_dt)
    
    # Use a smaller proportion of each fold's training dataset for hyperparam tuning
    hyperparam_dt<- train_dt[sample(1:train_nrows, round(train_nrows*hyperparam_prop))]
    
    if(model=="elasticnet"){
      
      elasticnet_cv_res<- hyperparam_tune(hyperparam_dt, model="elasticnet")
      best_tune_dt<- as.data.table(elasticnet_cv_res$bestTune)
      best_tune_dt[,Fold:=i]
      cv_hyperparams<- rbind(cv_hyperparams, best_tune_dt)
      
      elasticnet_model<- glmnet(x=as.matrix(as.data.frame(train_dt[,-c("Truth")])), 
                                y=train_dt[,Truth], family="binomial",
                                alpha=elasticnet_cv_res$bestTune$alpha, 
                                lambda=elasticnet_cv_res$bestTune$lambda)
      
      model_predict<- predict(elasticnet_model, as.matrix(as.data.frame(test_dt[,-c("Truth")])))
      prediction_labels<- factor(model_predict>=prob_cutoff, levels=c("FALSE", "TRUE"))
      
    } else if(model=="knn"){
      
      knn_cv_res<- hyperparam_tune(hyperparam_dt, model="knn")
      best_tune_dt<- as.data.table(knn_cv_res[[1]]$bestTune)
      # For kNN also include the best nvars selected by manual CV
      best_tune_dt[,nvars:=knn_cv_res[[3]]]
      best_tune_dt[,Fold:=i]
      cv_hyperparams<- rbind(cv_hyperparams, best_tune_dt)
      
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
      
      rf_cv_res<- hyperparam_tune(hyperparam_dt, model="rf")
      best_tune_dt<- as.data.table(rf_cv_res$bestTune)
      best_tune_dt[,Fold:=i]
      cv_hyperparams<- rbind(cv_hyperparams, best_tune_dt)
      
      # Consider tuning mtrees parameter
      rf_model<- randomForest(
        Truth~., train_dt, mtry=rf_cv_res$bestTune$mtry, na.action=na.roughfix)
      
      model_predict<- predict(rf_model, test_dt, type="prob")
      prediction_labels<- factor(model_predict[,2]>=prob_cutoff, levels=c("FALSE", "TRUE"))
      
    } else if(model=="gbm"){
      
      gbm_cv_res<- hyperparam_tune(hyperparam_dt, model="gbm")
      best_tune_dt<- as.data.table(gbm_cv_res$bestTune)
      best_tune_dt[,Fold:=i]
      cv_hyperparams<- rbind(cv_hyperparams, best_tune_dt)
      
      # The outcome needs to be 0 or 1
      train_dt[,Truth_int:=as.integer(as.logical(Truth))]
      test_dt[,Truth_int:=as.integer(as.logical(Truth))]
      
      gbm_model<- gbm(Truth_int~., train_dt[,-c("Truth")], distribution="bernoulli", 
                      n.trees=gbm_cv_res$bestTune$n.trees, 
                      interaction.depth=gbm_cv_res$bestTune$interaction.depth, 
                      shrinkage=gbm_cv_res$bestTune$shrinkage,
                      n.minobsinnode=gbm_cv_res$bestTune$n.minobsinnode)
      
      model_predict<- predict(gbm_model, test_dt, type="response", 
                              n.trees=gbm_cv_res$bestTune$n.trees)

      prediction_labels<- factor(model_predict>=prob_cutoff, levels=c("FALSE", "TRUE"))
      
    } else if(model=="svm"){
      
      svm_cv_res<- hyperparam_tune(hyperparam_dt, model="svm")
      best_tune_dt<- as.data.table(svm_cv_res$bestTune)
      best_tune_dt[,Fold:=i]
      cv_hyperparams<- rbind(cv_hyperparams, best_tune_dt)
      
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
  
  final_res<- list(cv_perf, cv_hyperparams)
  return(final_res)
  
}

