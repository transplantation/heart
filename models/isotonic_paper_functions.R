### In this file, we included all functions we wrote and used in the
### article: A Two-Stage Machine Learning Approach to Predict Heart 
### Transplantation Survival Probabilities over Time with a Monotonic 
### Probability Constraint

   
## cat_changer() function re-groups values in a categorical variable
## data_set: a data object that contains the varaible that we want to re-group
## var: the variable that we want to re-group
## val_old: some old levels in the variable
## val_new: some new levels in the variable
## val_old and val_new should be ONE-TO-ONE corresponds to each other
## Other levels that are not specified in val_old are put into "OTHER"
## Return value: a data object after re-groupping levels in the given variable
cat_changer <- function(data_set,var,val_old,val_new){
  temp_var <- dplyr::pull(data_set,var)
  cum_index <-c()
  for (i in 1:length(val_old)){
    index <- which(data_set[,var]==val_old[i])
    temp_var[index] <- val_new[i]
    cum_index <- c(cum_index,index)
  }
  na.index <- which(is.na(data_set[,var]))
  temp_var[-sort(c(cum_index, na.index))] <- "OTHER"
  data_set[,var] <- temp_var
  return(data_set)
}


## detect_terms() detects if the given term is in a vector
## x: a vector where we want to detect a given term
## term: a given term we want to detect
## Return Value: "YES" or "NO", "YES" means the given term appears at least once in the vector; "NO" means it doesn't appear in the vector
detect_terms <- function(x,term){
  result <- sapply(x, function(y) any(gsub(" ", "", strsplit(y, ",")[[1]])==term))
  if (all(is.na(result))) return(NA)
  result <- ifelse(is.na(result), FALSE, result)
  if (any(result==TRUE)) return("YES")
  return("NO")
}


## col_missing_function() counts number of nas in each column in a given data object
## input_data: a array, including a matrix, dataframe
## Return Values: a data frame object that reconds percentage of missing values in each column in the input data
col_missing_function <- function(input_data){
  # first, we count nas in each column
  na_count_col <- apply(input_data, 2, function(y) length(which(is.na(y)))/nrow(input_data)) 
  # we saved the object into a data_frame
  na_count_col <- data.frame(na_count_col) 
  return(na_count_col)
}


## dummy_maker() performs One-hot encoding (creates dummy variables) for categorical variables in a given data object
## input_data: a data object
## char_var: a vector of all independent categorical variables in the data
## Return Values: a data object after One-hot encoding is applied
dummy_maker <- function(input_data,char_var){
  for (i in 1:ncol(input_data)){
    if(names(input_data[i]) %in% char_var){
      # Use createDummyvars() function to create dummy variables for each categorical variable
      # The definition of createDummyvars() can be found in this file
      temp <- createDummyvars(input_data[i])
      names(temp) <- paste(names(input_data[i]),levels(as.factor(input_data[,i])),sep="_")
      input_data <- cbind(input_data,temp)
      input_data[,ncol(input_data)]<-NULL}
  }
  # We removed the dependent dummy variable in each categorical variable
  input_data <- input_data[-which(names(input_data) %in% char_var)]
  return(input_data)
}


## createDummyvars() function creates dummy variables for a given categorical variable
## data0: a column (vector, categoric) in the data
## Return Values: a data object after one-hot encoding is applied to the given variable
createDummyvars <- function(data0){
  dd<-as.data.frame(table(data0))
  dum_data<-as.data.frame(matrix(0, ncol = nrow(dd), nrow = nrow(data0)))
  names(dum_data) <- dd[,1]
  for(i in 1:ncol(dum_data)){
    dum_data[which(data0==names(dum_data)[i]),i]<-1
    dum_data[i]<-as.factor(dum_data[,i])
  }
  return(dum_data)
}


## class_generator_bino() function creates binary response variable
## gstatus: GRAFT FAILED (1=YES)
## gtime: GRAFT LIFESPAN-Days From Transplant to Failure/Death/Last Follow-Up
## p_unit: time period, the unit is year
## predict_length: the length of days, 365 here if the unit is 1 year 
## Return Value: a vector (our response/dependent variable), 0: death, 1:survival, NA: unknow status
class_generator_bino <- function(gstatus,gtime,p_unit,predict_length){
  p_unit <- as.numeric(p_unit)
  predict_length <- as.numeric(predict_length)
  if(gtime < p_unit*predict_length){
    if(is.na(gstatus)){return(NA)}else{
      if(gstatus==0){return(NA)}
      if(gstatus==1){return(0)}  # death
    }
  }else{
    return(1)  # survival
  }
}


## FFS_bin() function performs Fast Feature Selection (FFS) for a given data object
## ii: a number that indicates which response variable (corresponds to the time point) should be used
## In our study, ii=0: first month, ii=1: 1st year, ii=2: 2nd year, ..., up to ii=10: 10th year 
## df: a data object (this should be the finalized data after one-hot encoding is applied)
## ids: a list object that we save corresponding ids for all traing and test data objects
## exclud: a vector of variables excluded when finding important variables
## seed: random seed used in the function 
## Return Value: a vector of all variables selected, saved it as a list
FFS_bin <- function(ii,df,ids,exclud,seed=110){
  library(Biocomb)
  set.seed(seed)
  df1 <- df[df$ID %in% ids[[paste("ID_train",ii, sep="")]],c(names(df)[
    !names(df) %in% exclud])]
  disc <- "MDL"
  threshold <- 0.001
  attrs.nominal <- numeric()
  FF_vars <- select.fast.filter(df1, disc.method=disc, threshold=threshold, attrs.nominal=attrs.nominal)
  FF_vars$Information.Gain<-NULL
  FF_vars$NumberFeature<-NULL
  names(FF_vars) <- "variables"
  newlist <- list(FF_vars)
  names(newlist) <- paste("Year", ii, sep="")
  # colnames(FFSV)<-c("variables")
  #return(FF_vars)
  return(newlist)
}


## Lasso_bin() function performs LASSO Feature Selection for a given data object
## ii: a number that indicates which response variable (corresponds to the time point) should be used
## In our study, ii=0: first month, ii=1: 1st year, ii=2: 2nd year, ..., up to ii=10: 10th year 
## df: a data object (this should be the finalized data after one-hot encoding is applied)
## ids: a list object that we save corresponding ids for all traing and test data objects
## exclud: a vector of variables excluded when finding important variables
## seed: random seed used in the function 
## folds, trace, alpha: parameters used in cv.glmnet() function, details can be found in the package: glmnet
## Return Value: a vector of all variables selected, saved it as a list
Lasso_bin <- function(ii,df,ids,exclud,folds=5,trace=F,alpha=1,seed=110){
  set.seed(seed)
  yvar <- paste("year", ii, sep="")
  df1 <- df[df$ID %in% ids[[paste("ID_train",ii,sep="")]],c(names(df)[
    !names(df) %in% exclud],yvar)]
  
  dff <- df1[!names(df1) %in% yvar]
  
  for(i in 1:ncol(dff)){
    dff[i] <- as.numeric(dff[,i])
  }
  x <- data.matrix(dff)
  glmnet1 <- glmnet::cv.glmnet(x=x,y=as.factor(df1[,yvar]),type.measure='auc',nfolds=folds,alpha=alpha, family="binomial")
  co <- coef(glmnet1,s = "lambda.1se")
  inds <- which(co[,1]!=0)
  variables <- row.names(co)[inds]
  variables <- as.data.frame(variables[!(variables %in% '(Intercept)')])
  colnames(variables) <- c("variables")
  newlist <- list(variables)
  names(newlist) <- paste("Year", ii, sep="")
  return(newlist)
}


## RF_bin() function performs Feature Selection using Randon Forest Algorithm in the package: Boruta for a given data object
## ii: a number that indicates which response variable (corresponds to the time point) should be used
## In our study, ii=0: first month, ii=1: 1st year, ii=2: 2nd year, ..., up to ii=10: 10th year 
## df: a data object (this should be the finalized data after one-hot encoding is applied)
## ids: a list object that we save corresponding ids for all traing and test data objects
## exclud: a vector of variables excluded when finding important variables
## seed: random seed used in the function 
## Return Value: a vector of all variables selected, saved it as a list
RF_bin <- function(ii, df, ids, exclud, seed=110){
  set.seed(seed)
  TARGET <- paste("year", ii, sep="")
  df1 <- df[df$ID %in% ids[[paste("ID_train",ii,sep="")]],c(names(df)[
    !names(df) %in% exclud],TARGET)]
  
  dataset <- df1[!names(df1) %in% TARGET]
  dataset$TARGET <- as.factor(df1[,TARGET])
  
  Random_Forrest.train <- Boruta::Boruta(TARGET~., data = dataset, doTrace = 2)
  Random_Forrest_fs <- as.data.frame(Random_Forrest.train$finalDecision)
  names(Random_Forrest_fs)[1] <- paste("criteria")
  Random_Forrest_imp <- as.data.frame.table(subset(Random_Forrest_fs, Random_Forrest_fs[,1] == "Confirmed"))
  names(Random_Forrest_imp)[1] <- paste("variables")
  names(Random_Forrest_imp)[2] <- c("version")
  R_Forrest <- as.data.frame(Random_Forrest_imp$variables)
  colnames(R_Forrest) <- c("variables")
  newlist <- list(R_Forrest)
  names(newlist) <- paste("Year", ii, sep="")
  return(newlist)
}


## all_bin() combines all varaibles selected from FFS, LASSO, and Random Forest
## x: time, x=0: 1st month, x=1: 1st year, ...x=10: 10th year
## df: a list object where all variables from each selection algorithm are saved
## Return Value: a vector of all variables combined, saved it as a list
all_bin <- function(x, df){
  all_bind <- data.frame(variables = union(union(df$FFS[[x]][[1]], df$LASSO[[x]][[1]]), df$RF[[x]][[1]]))
  newlist <- list(all_bind)
  names(newlist) <- paste("Year", (x-1), sep="")
  return(newlist)
}


## pred_func() performs down-sampling with replacement on the training data
## and trains models given a machine learning algorithm and the down-sampling samples
## it also evaluates the performance measures: AUC, Sensitivity, Specificity, Accuracy on the holdout set
## This function is used with parSapply() function in the package: snow
## All functions used in this function should be imported locally since multi-cores are used
## ii: the ii-th down-sampling sample is used
## traindata: the training data object
## hold_out: the hold_out / test data object
## TARGET: the response / dependent variable in the data
## formul: formula used in the machine learning algorithm, refer to the function train() in the package caret
## var_numeric: a vector of names for numerical variables
## assigned_seed: set a random seed
## methods_input="log": machine learning algorithm used
## methods_input can be one of the following: "log", "rf_bag", "gbm_boost", "cart_bag"
## fold_no,repeat_no: parameters used in trainControl() function in the package: caret
## The default fold_no is 5 and repeat_no is 3
## Return Values: a list contains the following:
## 1. Performance measures from the holdout object
## 2. Predicted Survival Probabilities for the patients in the holdout object
## 3. AUC value for the down-sampling sample
pred_func <- function(ii,traindata,hold_out,TARGET,formul,var_numeric,assigned_seed,methods_input="log",fold_no,repeat_no){
  library(caret)
  library(AUC)
  library(MASS)
  
  ## RUS_func() performs down-sampling with replacement algorithm
  ## input_data: a data frame object (should be the training data)
  ## TARGET: the response / dependent variable in the data
  ## Return Value: a down-sampling sample
  RUS_func <- function(input_data,TARGET){
    Train_Two <- input_data[ which(input_data[TARGET]=="Two"), ]
    Train_One <- input_data[ which(input_data[TARGET]=="One"), ]
    if(nrow(Train_Two)<=nrow(Train_One)){
      sample_size<-nrow(Train_Two)
      Train_One <- Train_One[sample(nrow(Train_One), sample_size, replace=T), ]
      Train_Two <- Train_Two[sample(nrow(Train_Two), sample_size, replace=T), ]
    }else{
      sample_size<-nrow(Train_One)
      Train_One <- Train_One[sample(nrow(Train_One), sample_size, replace=T), ]
      Train_Two <- Train_Two[sample(nrow(Train_Two), sample_size, replace=T), ]
    }
    input_data<-rbind(Train_One,Train_Two)
    
    return(input_data)
  }
  
  set.seed((assigned_seed+ii))
  
  traindata[TARGET] <- as.factor(ifelse(traindata[TARGET]==0, "One", "Two"))
  hold_out[TARGET] <- as.factor(ifelse(hold_out[TARGET]==0, "One", "Two"))
  
  ## we obtain the bootstrap down-sampling data using RUS_func() function and use it to train the model 
  traindata <- RUS_func(traindata,TARGET)
  
  ## we scale numerical variables in the training data and holdout data
  for (i in 1:ncol(traindata)){
    if (colnames(traindata)[i]%in%var_numeric){
      temp_normalized <- scale(traindata[,i])
      traindata[,i] <- as.numeric(temp_normalized)
      hold_out[,i] <- (hold_out[,i]-attributes(temp_normalized)$`scaled:center`)/attributes(temp_normalized)$`scaled:scale`
    }
  }
  
  control_ <- trainControl(method = "repeatedcv", number=fold_no,  
                           repeats = repeat_no, classProbs = TRUE)
  
  if(methods_input=="log"){
    result_model <- train(formul, data=traindata, method="glm", family="binomial",
                          trControl = control_, metric="Accuracy")
  }
  
  
  if(methods_input=="rf_bag"){
    #mtry <- sqrt(ncol(traindata))
    #tunegrid <- expand.grid(.mtry=mtry)
    result_model <- train(formul,  data=traindata, method="rf",
                          trControl = control_, tuneLength=10, metric="Accuracy")
  }
  
  if(methods_input=="gbm_boost"){
    result_model <- train(formul, data=traindata, method="gbm",
                          trControl = control_, tuneLength=10, metric="Accuracy")
  }
  
  if(methods_input=="cart_bag"){
    result_model <- train(formul, data=traindata, method="treebag", family="binomial",
                          trControl = control_, tuneLength=10, metric="Accuracy")
  }
  
  resul_raw <- as.data.frame(matrix(NA, ncol = 3, nrow = nrow(hold_out)))
  colnames(resul_raw) <- c("TARGET", methods_input, "Probability")
  resul_raw$TARGET <- hold_out[TARGET]
  
  train_raw <- as.data.frame(matrix(NA, ncol = 3, nrow = nrow(traindata)))
  colnames(train_raw) <- c("TARGET", methods_input, "Probability")
  train_raw$TARGET <- traindata[TARGET]
  
  resul_pred_perf<-as.data.frame(matrix(NA, ncol = 1, nrow = 4))
  colnames(resul_pred_perf)<-c(methods_input)
  rownames(resul_pred_perf)<-c("auc","sen","spec","accu")
  train_auc <- NA
  
  train_raw$Probability <- predict(result_model, newdata=traindata, type="prob")[,2]
  train_raw[methods_input] <- predict(result_model, newdata=traindata, type="raw")
  train_auc <- AUC::auc(roc(train_raw$Probability, traindata[,TARGET]))
  resul_raw[methods_input] <- predict(result_model, newdata=hold_out, type="raw")
  resul_raw$Probability <- predict(result_model, newdata=hold_out, type="prob")[,2]
  resul_pred_perf[1,1] <- AUC::auc(roc(resul_raw$Probability,hold_out[,TARGET]))
  resul_pred_perf[2,1] <- caret::sensitivity(resul_raw[,methods_input],hold_out[,TARGET])
  resul_pred_perf[3,1] <- caret::specificity(resul_raw[,methods_input],hold_out[,TARGET])
  resul_pred_perf[4,1] <- (as.data.frame(confusionMatrix(resul_raw[,methods_input], hold_out[,TARGET])$overall))[1,]
  
  return(list(Performance=resul_pred_perf, Predicted=resul_raw, AUC=train_auc))
}


#The following two ________________________________________________________________________

# check_range: check if numerical values are in the corresponding range
check_range <- function(df, vars, Num_Range){
  out_range <- rep(NA, 13)
  for (i in 1:length(vars)){
    m <- min(df[,vars[i]])
    M <- max(df[,vars[i]])
    out_range[i] <- ifelse((m<Num_Range[i,1]|M>Num_Range[i,2]), 1, 0)
  }
  return(out_range)
}
#________________________________________________________________________

# check_levels: check if categorical values are in our data
check_levels <- function(df, vars, Cat_Levels){
  CAT_match <- rep(NA, length(vars))
  for (i in 1:length(vars)){
    temp_level <- levels(as.factor(df[,vars[i]]))
    CAT_match[i] <- ifelse(all(temp_level%in%Cat_Levels[[i]]), 0, 1)
  }
  return(CAT_match)
}


#________________________________________________________________________

# creating dummy variables for input data (this is for shinny app)
dummy_maker_app<-function(input_data, char_var){
  for (i in 1:ncol(input_data)){
    if(names(input_data[i]) %in% char_var){
      temp<-createDummyvars(input_data[i])
      names(temp)<-paste(names(input_data[i]),levels(as.factor(input_data[,i])),sep="_")
      
      input_data<-cbind(input_data,temp)
      if (length(levels(as.factor(input_data[,i])))!=1){
        input_data[,ncol(input_data)]<-NULL
      }else{
        input_data[,ncol(input_data)]<-factor(input_data[,ncol(input_data)], levels = c(0,1))
      }
    }
  }
  input_data<-input_data[-which(names(input_data) %in% char_var)]
  return(input_data)
}

#________________________________________________________________________

# change all "U" to Unknown (this is for shinny app)
find_U <- function(x){
  col_index <- which(x=="U")
  if (length(col_index)!=0) x[col_index] <- "UNKOWN"
  return(x)
}