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
cat_changer<- function(data_set,var,val_old,val_new){
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


######################################################
############monotonic survival prediction#############
############monotonic survival prediction#############
######################################################
survivals_cal<-function(patients){
  
  #cat("\014") # Clear the console
  
  source("https://raw.githubusercontent.com/transplantation/heart/master/models/isotonic_paper_functions.R")
  library(shiny)
  # read csv file
  #file1 <- read.csv("https://raw.githubusercontent.com/transplantation/heart/master/models/example_data.csv", stringsAsFactors = FALSE)
  file1 <- patients
  
  # load important variables (UNOS) and corresponding information
  variables <- read.csv("https://raw.githubusercontent.com/transplantation/heart/master/models/important_variables.csv")
  # Load the center (mean) and spread (standard devation) for numerical variables
  Num_Scales <- readRDS(gzcon(url("https://github.com/transplantation/heart/raw/master/models/Num_Scales.rds")))
  # Load the levels for categroical variables
  Cat_Levels <- readRDS(gzcon(url("https://github.com/transplantation/heart/raw/master/models/Cat_Levels.rds")))
  # Load the features from LASSO selection 
  features <- readRDS(gzcon(url("https://github.com/transplantation/heart/raw/master/models/LASSO_features.rds")))
  
  log_models<-read.csv("https://raw.githubusercontent.com/transplantation/heart/master/models/log_models.csv",row.names = 1)
  
  is_error<-"NO"
  #error for no of rows equal zero
  error_no<-""
  # check if the variables are provided
  error_vars<-""
  # error if the numerical variables are out of range
  error_range<-""
  # error if the value of categorical variables are valid
  error_cat<-""
  
  if (nrow(file1)==0){error_no<- c("  Your file is empty  ")
  is_error<-"YES"
  }
  
  # the rest will run only if the is_error=="NO"
  if(is_error=="NO"){
    
    # we just work with first 10 patients' data
    if (nrow(file1)>10) file1 <- file1[1:10,]
    patients_no<-nrow(file1)
    
    num_vars <- as.character(variables$Name[variables$Type=="NUM"])
    cat_vars <- as.character(variables$Name[variables$Type=="CAT"])
    cat_vars[which(cat_vars=="ETH_MAT")] <- "ETHCAT_DON" 
    index <- which(!(c(num_vars, cat_vars)%in%colnames(file1)))
    
    
    
    for (i in 1:length(num_vars)){
      file1[,num_vars[i]] <- as.numeric(file1[,num_vars[i]])
    }
    
    
    if (length(index)>0){
      # return("One or more variables are not found. Please build file based on the Example CSV.")
      #__
      error_vars<- c("  One or more variables are not found. Please build file based on the Example CSV.  ")
      is_error<-"YES"
    }
    
    # the rest will run only if the is_error=="NO"
    if(is_error=="NO"){
      Num_Range <- matrix(c(5, 10, 0, 0, 0, 10, 0, 120, 120, 10, 120, 10, 0, 100, 45, 100, 3000, 3000, 45, 200, 230, 230, 45, 230, 720, 30) ,ncol=2, nrow=13, byrow=F)
      
      # ISCHTIME
      file1$ISCHTIME <- file1$ISCHTIME*60
      
      # a function designed for checking if the numerical variables selected by user are within a range
      out_range <- check_range(file1, num_vars, Num_Range)
      
      
      if (sum(out_range)>0){
        error_range<-paste("One or more values for", num_vars[which(out_range==1)],
                           "are outside of the following range: (", Num_Range[which(out_range==1),1],
                           ",", Num_Range[which(out_range==1),2],").")
        is_error<-"YES"
      }else{
        error_range<- c("  All the numerical values are within the expected range. ")
      }
      
      # the rest will run only if the is_error=="NO"
      if(is_error=="NO"){
        
        # Now, change categorical variables to our features
        # ANCEF
        find_UNKNOWN <- which(file1$ANCEF=="UNKNOWN")
        if (length(find_UNKNOWN)!=0) file1$ANCEF[find_UNKNOWN] <- "UNKOWN"
        
        # CARD_SURG
        find_UNKNOWN <- which(file1$CARD_SURG=="UNKNOWN")
        if (length(find_UNKNOWN)!=0) file1$CARD_SURG[find_UNKNOWN] <- "UNKOWN"
        
        # DIAG
        val_old <- c(1000,1001,1002,1003,1004,1005,1006,1049,1007,1200)
        val_new <- c("DILATED_MYOPATHY_IDI","DILATED_MYOPATHY_OTH","DILATED_MYOPATHY_OTH","DILATED_MYOPATHY_OTH","DILATED_MYOPATHY_OTH","DILATED_MYOPATHY_OTH","DILATED_MYOPATHY_OTH","DILATED_MYOPATHY_OTH","DILATED_MYOPATHY_ISC","CORONARY")
        file1 <- cat_changer(file1,var="DIAG",val_old,val_new)
        
        # CHEST_XRAY_DON
        val_old <- c(0,1,2,3,4,5,998,999)
        val_new <- c(NA,NA,"NOR","AB","AB","ABboth",NA,NA)
        file1 <- cat_changer(file1,var="CHEST_XRAY_DON",val_old,val_new)
        
        # CMV_DON & EBV_SEROSTATUS
        val_old <- c("C","I","N","ND","P","U")
        val_new <- c(NA,NA,"Neg",NA,"POS",NA)
        file1 <- cat_changer(file1,var="CMV_DON",val_old,val_new)
        file1 <- cat_changer(file1,var="EBV_SEROSTATUS",val_old,val_new)
        
        # COD_CAD_DON
        val_old <- c(1,2,3,4,999,"Unknown")
        val_new <- c("ANOXIA","CEREBROVASCULAR_STROKE","HEAD_TRAUMA","OTHER","OTHER",NA)
        file1 <- cat_changer(file1,var="COD_CAD_DON",val_old,val_new)
        
        # CORONARY_ANGIO
        val_old <- c(1,2,3)
        val_new <- c("NO", "YES", "YES")
        file1 <- cat_changer(file1,var="CORONARY_ANGIO",val_old,val_new)
        
        # DIAB
        val_old <- c(1,2,3,4,5,998)
        val_new<-c("no","one","two","OTHER","OTHER",NA)
        file1 <- cat_changer(file1,var="DIAB",val_old,val_new)
        
        # DIAG 
        val_old <- c(1000,1001,1002,1003,1004,1005,1006,1049,1007,1200)
        val_new <- c("DILATED_MYOPATHY_IDI","DILATED_MYOPATHY_OTH","DILATED_MYOPATHY_OTH","DILATED_MYOPATHY_OTH","DILATED_MYOPATHY_OTH","DILATED_MYOPATHY_OTH","DILATED_MYOPATHY_OTH","DILATED_MYOPATHY_OTH","DILATED_MYOPATHY_ISC","CORONARY")
        file1 <- cat_changer(file1,var="DIAG",val_old,val_new)
        file1 <- cat_changer(file1,var="TCR_DGN",val_old,val_new)
        file1 <- cat_changer(file1,var="THORACIC_DGN",val_old,val_new)
        
        # DOPAMINE
        find_UNKNOWN <- which(file1$DOPAMINE=="UNKNOWN")
        if (length(find_UNKNOWN)!=0) file1$DOPAMINE[find_UNKNOWN] <- "UNKOWN"
        
        # EDUCATION
        val_old <- c(1,2,3,4,5,6,996,998)
        val_new <- c("a","a","b","c","d","d",NA,NA)
        file1 <- cat_changer(file1,var="EDUCATION",val_old,val_new)
        
        # ETH_MAT
        file1$ETH_MAT <- NA
        for(i in 1:nrow(file1)){
          if(!is.na(file1$ETHCAT[i])){
            if(!is.na(file1$ETHCAT_DON[i])){
              if(file1$ETHCAT_DON[i]==file1$ETHCAT[i]){
                file1$ETH_MAT[i] <- "Y"
              }else{
                file1$ETH_MAT[i]<-"N"
              }
            }
          }
        }
        
        # ETHCAT, ETHCAT_DON
        val_old <- c(1,2,4,5,6,7,9,998)
        val_new <- c("w","b","h","o","o","o","o",NA)
        file1 <- cat_changer(file1,var="ETHCAT",val_old,val_new)
        file1 <- cat_changer(file1,var="ETHCAT_DON",val_old,val_new)
        
        # FUNC_STAT_TCR, FUNC_STAT_TRR
        val_old <- c(1,2,3,996,998,2010,2020,2030,2040,2050,2060,2070,2080,2090,2100)
        val_new <- c("A","B","B",NA,NA,"C","C","C","C","D","D","D","E","E","E")
        file1 <- cat_changer(file1,var="FUNC_STAT_TRR",val_old,val_new)
        file1 <- cat_changer(file1,var="FUNC_STAT_TCR",val_old,val_new)
        
        # HEPARIN
        find_UNKNOWN <- which(file1$HEPARIN=="UNKNOWN")
        if (length(find_UNKNOWN)!=0) file1$HEPARIN[find_UNKNOWN] <- "UNKOWN"
        
        # HLAMIS
        val_old <- 0:6
        val_new <- c("a","a","a","b","c","f","e")
        file1 <- cat_changer(file1,var="HLAMIS",val_old,val_new)
        
        # INIT_STAT
        val_old <- c(2010,2020,2030,2090,2999)
        val_new <- c("ONE","ONE","TWO","ONE","OTHER")
        file1 <- cat_changer(file1,var="INIT_STAT",val_old,val_new)
        
        # LAST_INACT_REASON
        val_old <- seq(1,14)
        val_new <- c("ONE", "ONE","ONE","ONE","ONE", "ONE","TWO","ONE","TWO", "ONE","ONE","ONE","ONE","ONE")
        file1 <- cat_changer(file1,var="LAST_INACT_REASON",val_old,val_new)
        
        # PRI_PAYMENT_TCR, PRI_PAYMENT_TRR
        val_old <- seq(1,14)
        val_new<-c("pv","pbma","pbmcffs","pbmoth","pbmoth","pbmoth","pbmoth","OTHER","OTHER","OTHER","OTHER","OTHER","OTHER","OTHER")
        file1 <- cat_changer(file1,var="PRI_PAYMENT_TCR",val_old,val_new)
        file1 <- cat_changer(file1,var="PRI_PAYMENT_TRR",val_old,val_new)
        
        # PRIOR_CARD_SURG_TYPE_TCR
        val_old <- seq(1,31)
        val_new <- c("CABG","VALV", "CABG","OTHER", "CABG","VALV", "CABG","OTHER", "CABG","VALV",
                     "CABG","OTHER", "CABG","VALV", "CABG","OTHER", "CABG","VALV", "CABG","OTHER",
                     "CABG","VALV", "CABG","OTHER", "CABG","VALV", "CABG","OTHER", "CABG","VALV", "CABG")
        file1 <- cat_changer(file1,var="PRIOR_CARD_SURG_TYPE_TCR",val_old,val_new)
        
        # PROC_TY_HR
        val_old <- c(1,2)
        val_new <- c("Bicaval","Traditional")
        file1 <- cat_changer(file1,var="PROC_TY_HR",val_old,val_new)
        
        # REGION
        val_old <- seq(1,11)
        val_new <- c("NE","NE","SE","SE","W","W","MW","MW","NE","MW","SE")
        file1 <- cat_changer(file1,var="REGION",val_old,val_new)
        
        # SHARE_TY"
        val_old <- c(3,4)
        val_new <- c("LOCAL","REGIONAL")
        file1 <- cat_changer(file1,var="SHARE_TY",val_old,val_new)
        
        # STERNOTOMY_TRR
        val_old <- c(1,2,3,998)
        val_new <- c("ONE", "MORE","MORE", NA)
        file1 <- cat_changer(file1,var="STERNOTOMY_TRR",val_old,val_new)
        ###
        if (nrow(file1)==1){
          file1 <- as.data.frame(t(apply(file1, 2, function(x) gsub("^$| ^", NA, x))), stringsAsFactors=FALSE)
        }else{
          file1 <- as.data.frame(apply(file1, 2, function(x) gsub("^$| ^", NA, x)), stringsAsFactors=FALSE)
        }
        
        file1[is.na(file1)] <- "UNKOWN"
        
        if (nrow(file1)==1){
          file1 <- as.data.frame(t(apply(file1, 2, find_U)))
        }else{
          file1 <- as.data.frame(apply(file1, 2, find_U))
        }
        
        file1[,num_vars] <- apply(file1[,num_vars],2,as.numeric)
        
        
        ####
        # One-hot encoding for the data
        cat_vars <- c(cat_vars, "ETH_MAT")
        for (i in 1:length(cat_vars)){
          file1[,cat_vars[i]] <- gsub(" ", "", file1[,cat_vars[i]],fixed=TRUE)
          file1[,cat_vars[i]] <- factor(file1[,cat_vars[i]], levels=Cat_Levels[[i]])
        }
        
        #new_variables <- c("ETH_MAT", "ANCEF", "DOPAMINE", "HEPARIN", "CARD_SURG")
        CAT_match <- check_levels(file1, cat_vars, Cat_Levels)
        
        if (sum(CAT_match)>0){
          return(paste("One or more values for", cat_vars[which(CAT_match==1)],
                       "are invalid."))
        }else{
          renderText("All the categorical values are valid. Now. ")
        }
        
        
        if (sum(CAT_match)>0){
          error_cat<-paste("One or more values for", cat_vars[which(CAT_match==1)],
                           "are invalid.")
          is_error<-"YES"
          
        }else{
          error_cat<-paste(" All categorical values are valid. ")
        }
        if(is_error=="NO"){
          
          temp.dum <- dummy_maker_app(file1, cat_vars)
          
          
          all_variables <- c()
          for (i in 1:11){
            all_variables <- c(all_variables, as.character(features[[i]][[1]]))
          }
          all_variables <- unique(all_variables)
          all_num_index <- which(all_variables%in%num_vars)
          all_num <- all_variables[all_num_index]  # num variables from LASSO
          all_cat <- all_variables[-all_num_index] # cat variables from LASSO
          
          file1.dum <- matrix(NA, ncol=length(all_variables), nrow=patients_no) 
          colnames(file1.dum) <- all_variables
          file1.dum <- as.data.frame(file1.dum)
          
          file1.dum[,all_num] <- file1[,all_num]
          
          for (i in 1:length(all_cat)){
            find_vars <- which(all_cat[i]%in%colnames(temp.dum))
            if (length(find_vars)==0){
              file1.dum[,all_cat[i]] <- 0 # zero vector
            }else{
              file1.dum[,all_cat[i]] <- temp.dum[,all_cat[i]]
            }
          }
          
          all_year_data <- rep(list(NA),11)
          for (i in 1:11){
            temp_vars <- as.character(features[[i]][[1]])
            temp_data <- file1.dum[,temp_vars]
            temp_index <- which(temp_vars%in%all_num)
            temp_cat_index <- setdiff(1:ncol(temp_data), temp_index)
            for (j in 1:length(temp_cat_index)){
              temp_data[,temp_cat_index[j]] <- factor(temp_data[,temp_cat_index[j]])
            }
            temp_scale <- Num_Scales[[i]]
            for (j in 1:length(temp_index)){
              var_name <- temp_vars[temp_index[j]]
              find_index <- which(temp_scale$Var_names==var_name)
              temp_data[,var_name] <- (temp_data[,var_name]-temp_scale[find_index,2])/temp_scale[find_index,3]
            }
            all_year_data[[i]] <- temp_data
          }
          
          # the following is for prediction but I use the matrix
          
          names(all_year_data) <- paste0("Year",0:10)
          survivals <- matrix(NA, ncol=11, nrow=patients_no)
          survival_Probability <- matrix(NA, ncol=11, nrow=patients_no)
          
          # predictions<-matrix(NA, ncol = 11, nrow = patients_no)
          
          for (i in 1:11){
            # Insert one column with all 1s for the intercept
            X <- cbind(INTERCEPT=rep(1,nrow(file1.dum)), all_year_data[[i]])
            # Make sure the data type is numeric
            if (nrow(X)==1){
              X <- as.matrix(t(apply(X, 2, as.numeric)))}else{
                X <- as.matrix(apply(X, 2, as.numeric))
              }
            # Get the coefficient matrix for the corresponding time point
            BETA_X <- log_models[i,colnames(X)]
            # Find the fitted value for each patient 
            Y <- as.matrix(BETA_X)%*%t(as.matrix(X))
            # Find the survival probability for each patient
            survivals[,i] <- exp(Y)/(1+exp(Y))
          }
          
          survivals_reverse <- matrix(NA, ncol=11, nrow=nrow(file1.dum))
          for (i in 1:11){
            survivals_reverse[,i] <- survivals[,(12-i)]
          }
          
          # apply isotonic regression to the probability matrix
          survivals_isotonic <- t(apply(survivals_reverse, 1, function(x) isoreg(x)$yf))
          
          for (i in 1:11){
            survival_Probability[,i] <- survivals_isotonic[,(12-i)]
          } 
          
          
          # the survival probability matrix after isotonic regression is applied
          survival_Probability <- cbind(1:patients_no, survival_Probability)
          colnames(survival_Probability) <- c("Patient_No"," Month_1 Survival Prob.", paste(" Year_", 1:10, " Survival Prob.",sep=""))
          
          survival_Probability <- as.data.frame(survival_Probability) 
          outcome<-list()
          outcome$survival_Probability<-survival_Probability
          outcome$patients_no<-patients_no
          outcome$error_no<-error_no
          outcome$error_vars<-error_vars
          outcome$error_range<-error_range
          outcome$error_cat<-error_cat
          
          return(outcome)
        }
      }
    }
  }
  if(is_error=="YES"){
    outcome<-list()
    outcome$survival_Probability<-print(message)
    outcome$patients_no<-0
    outcome$error_no<-error_no
    outcome$error_vars<-error_vars
    outcome$error_range<-error_range
    outcome$error_cat<-error_cat
    return(outcome)
  }
  
}

















