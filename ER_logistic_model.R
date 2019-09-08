# use the created negative and positive set of sequences
# perform random forest classification

########################################################################################################################
library(aod)
library(ROCR)
library(caret)
library(DMwR)
library(ROSE)
library(VGAM)
library(LiblineaR)
library(MLmetrics)
library(PRROC)
library(purrr)
library(dplyr)
########################################################################################################################
# Functions
dataset_to_numericFeature_factor_label <- function(my_dataset){
  dataset <- my_dataset
  aaclass <- sapply(dataset, class)
  stopifnot(all(aaclass %in% "character"))
  aacolsnum <- colnames(dataset)[1:(ncol(dataset) - 1)]
  dataset[aacolsnum] <- sapply(dataset[aacolsnum], as.numeric)
  colnames(dataset)[ncol(dataset)] <- "label_set"
  dataset <- transform(dataset,  label_set = as.factor(label_set))
  return(dataset)
}
########################################################################################################################
calc_auprc <- function(model, data){
  
  index_Pos <- data$label_set == "Pos"
  index_Neg<- data$label_set == "Neg"
  
  predictions <- predict(model, data, type = "prob")
  
  pr.curve(predictions$Pos[index_Pos], predictions$Pos[index_Neg], curve = TRUE)
  
}
########################################################################################################################
test_roc <- function(model, data) {
  library(pROC)
  roc_obj <- roc(data$label, 
                 predict(model, data, type = "prob")[, "pos"],
                 levels = c("neg", "pos"))
  ci(roc_obj)
}
########################################################################################################################


#########################################################################################################
# feature sets:

feature_set_sum_LLR_1000bp_4
# feature sets:
#sum LR
feature_set_sum_LLR_1000bp_4 <- rbind(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[4]], 
                                      Positive_set_seq_list_chopped_scoreMat_list_1000bp[[4]])
colnames(feature_set_sum_LLR_1000bp_4)[13] <- "NKX3_1"
# adjacency
feature_set_Adjmat_1000bp_4 <- rbind(Negative_set_seq_list_chopped_adjanScore_list_1000bp[[4]], 
                                     Positive_set_seq_list_chopped_adjanScore_list_1000bp[[4]])
label_set <- c(rep("Neg", nrow(Negative_set_seq_list_chopped_scoreMat_list_1000bp[[4]])), 
               rep("Pos", nrow(Positive_set_seq_list_chopped_scoreMat_list_1000bp[[4]])))

# changing the name of NKX3-1 to NKX3_1 so it works in the classification models
aa <- strsplit(colnames(feature_set_Adjmat_1000bp_4), split = "_vs_")
for(i in 1:length(aa)){
  aa[[i]][aa[[i]] == "NKX3-1"] <- "NKX3_1"
}
for(i in 1:length(aa)){
  aa[[i]] <- paste(aa[[i]][1], aa[[i]][2], sep = "_vs_")
}
colnames(feature_set_Adjmat_1000bp_4) <- unlist(aa)
# overlap
feature_set_Overlapmat_1000bp_4 <- rbind(Negative_set_seq_list_chopped_OverlapScore_list_1000bp[[4]], 
                                         Positive_set_seq_list_chopped_OverlapScore_list_1000bp[[4]])
colnames(feature_set_Overlapmat_1000bp_4) <- unlist(aa)


my_learning_datasets_1000bp_4 <- list()
# sum LR only
my_learning_datasets_1000bp_4[[1]] <- as.data.frame(cbind(feature_set_sum_LLR_1000bp_4,
                                                          label_set), 
                                                    stringsAsFactors =F)
my_learning_datasets_1000bp_4[[1]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_4[[1]])
# adjacancy only
my_learning_datasets_1000bp_4[[2]] <- as.data.frame(cbind(feature_set_Adjmat_1000bp_4,
                                                          label_set), 
                                                    stringsAsFactors =F)
my_learning_datasets_1000bp_4[[2]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_4[[2]])

# Overlap only
my_learning_datasets_1000bp_4[[3]] <- as.data.frame(cbind(feature_set_Overlapmat_1000bp_4,
                                                          label_set), 
                                                    stringsAsFactors =F)
my_learning_datasets_1000bp_4[[3]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_4[[3]])
# LLR and adj
my_learning_datasets_1000bp_4[[4]] <- as.data.frame(cbind(cbind(feature_set_sum_LLR_1000bp_4,
                                                                feature_set_Adjmat_1000bp_4),
                                                          label_set), 
                                                    stringsAsFactors =F)
my_learning_datasets_1000bp_4[[4]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_4[[4]])
# LLR and overlap
my_learning_datasets_1000bp_4[[5]] <- as.data.frame(cbind(cbind(feature_set_sum_LLR_1000bp_4,
                                                                feature_set_Overlapmat_1000bp_4),
                                                          label_set), 
                                                    stringsAsFactors =F)
my_learning_datasets_1000bp_4[[5]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_4[[5]])
# LLR and adj and overlap
my_learning_datasets_1000bp_4[[6]] <- as.data.frame(cbind(cbind(feature_set_sum_LLR_1000bp_4,
                                                                cbind(feature_set_Adjmat_1000bp_4, 
                                                                      feature_set_Overlapmat_1000bp_4)),
                                                          label_set), 
                                                    stringsAsFactors =F)
colnames(my_learning_datasets_1000bp_4[[6]])[28:405] <- paste(colnames(my_learning_datasets_1000bp_4[[6]])[28:405], "adj", sep = "_at_")
colnames(my_learning_datasets_1000bp_4[[6]])[406:783] <- paste(colnames(my_learning_datasets_1000bp_4[[6]])[406:783], "ovl", sep = "_at_")
my_learning_datasets_1000bp_4[[6]]  <- dataset_to_numericFeature_factor_label(my_learning_datasets_1000bp_4[[6]])

names(my_learning_datasets_1000bp_4) <- c("Sum_LR","Adjacency", "Overlap", 
                                          "Sum_LR_plus_Adjacency", "Sum_LR_plus_Overlap", 
                                          "Sum_LR_plus_Adjacency_plus_Overlap")

my_data_partitions[[7]] <- createDataPartition(y = my_learning_datasets_1000bp_4[[1]][, ncol(my_learning_datasets_1000bp_4[[1]])], times = 1, p = 0.75,list = F)
##############
save(list = c("my_learning_datasets_1000bp_4",
              "my_data_partitions"), file = "Classification_input_set4.RData")

my_learning_datasets_results_RF_kappa_down_1000bp_4 <- list()

length(my_learning_datasets_1000bp_4)
for(i in 1:1){
  aanzv <- nearZeroVar(my_learning_datasets_1000bp_4[[i]], saveMetrics= F)
  if(length(aanzv) > 0){
    aa_all_data <- my_learning_datasets_1000bp_4[[i]][,-aanzv]
  }else{
    aa_all_data <- my_learning_datasets_1000bp_4[[i]]
  }
  
  
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
    aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
    aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
    if(length(aahighlyCorDescr) > 0){
      aa_all_data <- aa_all_data[,-aahighlyCorDescr]
    }
  }
  print("#########################################################################################")
  print(paste0("number of features in original feature set number ", i," named ", 
               names(my_learning_datasets_1000bp_4)[i] ,":"))
  print(ncol(my_learning_datasets_1000bp_4[[i]]))
  print("number of features after removing near zero variance and highly correlated features ( > 0.9)")
  print(ncol(aa_all_data))
  print("#########################################################################################")
  aa_train_data <- aa_all_data[my_data_partitions[[7]], ]
  aa_test_data <-  aa_all_data[-my_data_partitions[[7]], ]
  
  
  aactrl <- trainControl(method = "repeatedcv",
                         number = 5,
                         repeats = 3,
                         classProbs = TRUE,
                         summaryFunction = mnLogLoss,
                         savePredictions = TRUE,
                         ## new option here:
                         sampling = "down", 
                         search = 'random')
  print("aaRF_down")
  aaRF_down <- train(label_set ~ ., data = aa_train_data,
                     method = "svmRadial",
                     preProc = c("center", "scale"),
                     #nbagg = 50,
                     metric = "logLoss",
                     trControl = aactrl,
                     tuneLength  = 10)
  my_learning_datasets_results_RF_kappa_down_1000bp_4[[i]] <- aaRF_down  
}

save(list = c("my_learning_datasets_results_RF_kappa_down_1000bp_4"),
     file = "classification_results_set_4.RData")
load("classification_results_set_5.RData")

aa_gglist_2 <- list()

for(i in 1:length(my_learning_datasets_1000bp_4)){
  aanzv1 <- nearZeroVar(my_learning_datasets_1000bp_4[[i]], saveMetrics= F)
  if(length(aanzv1) > 0){
    aa_all_data1 <- my_learning_datasets_1000bp_4[[i]][,-aanzv1]
  }else{
    aa_all_data1 <- my_learning_datasets_1000bp_4[[i]]
  }
  
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data1[, 1:(ncol(aa_all_data1) - 1)]) > 1){
    aadescrCor1 <- cor(aa_all_data1[, 1:(ncol(aa_all_data1) - 1)])
    aahighlyCorDescr1 <- findCorrelation(aadescrCor1, cutoff = .9)
    if(length(aahighlyCorDescr1) > 0){
      aa_all_data1 <- aa_all_data1[,-aahighlyCorDescr1]
    }
  }
  
  aa_train_data1 <- aa_all_data1[my_data_partitions[[7]], ]
  aa_test_data1 <-  aa_all_data1[-my_data_partitions[[7]], ]
  
  aaindex_pos1 <- aa_test_data1$label_set == "Pos"
  aaindex_neg1 <- aa_test_data1$label_set == "Neg"
  
  aarfdown1 <- predict(my_learning_datasets_results_RF_kappa_down_1000bp_4[[i]], 
                       aa_test_data1, type = "prob")
  
  aapreds_list <- list(aarfdown1$Pos)
  # List of actual values (same for all)
  m <- length(aapreds_list)
  aactuals_list <- rep(list(aa_test_data1$label_set), m)
  # Plot the ROC curves
  aapred <- prediction(aapreds_list, aactuals_list)
  aauc <- performance(aapred, measure = "auc")
  aarocs <- performance(aapred, "tpr", "fpr")
  print(aauc@y.values)
  plot(aarocs, col = as.list(1:m), main = paste0("Test Set ROC Curves for feature set: ",
                                                 names(my_learning_datasets_1000bp_4)[i]))
  abline(a = 0, b = 1, col = 6, lty = 3, lwd = 0.7)
  legend(x = "bottomright", cex = 0.6,
         legend = c(paste0("RF_down -Set4 ", "AUC: " , format(aauc@y.values[[1]], digits = 2))),
         fill = 1:m)
  aa_model_list <- list(RF_down1 = my_learning_datasets_results_RF_kappa_down_1000bp_4[[i]])
  
  
  aamodel_list_pr <- aa_model_list %>%
    map(.f = calc_auprc, data = (aa_test_data1))
  aamodel_list_pr_auc <- aamodel_list_pr %>%
    map(function(the_mod) the_mod$auc.integral)
  aamodel_list_pr_auc <- unlist(aamodel_list_pr_auc)
  
  aaresults_list_pr <- list(NA)
  num_mod <- 1
  
  for(the_pr in aamodel_list_pr){
    
    aaresults_list_pr[[num_mod]] <- data_frame(recall = the_pr$curve[, 1],
                                               precision = the_pr$curve[, 2],
                                               model = names(aamodel_list_pr)[num_mod])
    
    num_mod <- num_mod + 1
    
  }
  
  aaresults_df_pr <- bind_rows(aaresults_list_pr)
  
  aacustom_col <- colorspace::rainbow_hcl(length(aa_model_list))
  
  aa_gglist_2[[i]] <- ggplot(aes(x = recall, y = precision, group = model), data = aaresults_df_pr) +
    geom_line(aes(color = model), size = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = aacustom_col, labels = paste(names(aa_model_list) ,"AUC: " ,
                                                             format(aamodel_list_pr_auc, digits = 2))) +
    geom_abline(intercept = sum(aa_test_data1$label_set == "Pos")/nrow(aa_test_data1),
                slope = 0, color = "gray", size = 1) +
    ggtitle(paste0("PR curve ", names(my_learning_datasets_1000bp_4)[i]))+
    theme_bw()
  
}
names(aa_gglist_2) <-  names(my_learning_datasets_1000bp_4)
aa_gglist_2[[7]]

varImp(my_learning_datasets_results_RF_kappa_down_1000bp_5[[4]])
################################################
####### get performance on training
################################################
aa_gglist_3 <- list()

for(i in 1:length(my_learning_datasets_1000bp_4)){
  aanzv1 <- nearZeroVar(my_learning_datasets_1000bp_4[[i]], saveMetrics= F)
  if(length(aanzv1) > 0){
    aa_all_data1 <- my_learning_datasets_1000bp_4[[i]][,-aanzv1]
  }else{
    aa_all_data1 <- my_learning_datasets_1000bp_4[[i]]
  }
  
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data1[, 1:(ncol(aa_all_data1) - 1)]) > 1){
    aadescrCor1 <- cor(aa_all_data1[, 1:(ncol(aa_all_data1) - 1)])
    aahighlyCorDescr1 <- findCorrelation(aadescrCor1, cutoff = .9)
    if(length(aahighlyCorDescr1) > 0){
      aa_all_data1 <- aa_all_data1[,-aahighlyCorDescr1]
    }
  }
  
  aa_train_data1 <- aa_all_data1[my_data_partitions[[7]], ]
  aa_test_data1 <-  aa_all_data1[-my_data_partitions[[7]], ]
  
  aaindex_pos1 <- aa_test_data1$label_set == "Pos"
  aaindex_neg1 <- aa_test_data1$label_set == "Neg"
  
  aarfdown1 <- predict(my_learning_datasets_results_RF_kappa_down_1000bp_4[[i]], 
                       aa_train_data1, type = "prob")
  
  aapreds_list <- list(aarfdown1$Pos)
  # List of actual values (same for all)
  m <- length(aapreds_list)
  aactuals_list <- rep(list(aa_train_data1$label_set), m)
  # Plot the ROC curves
  aapred <- prediction(aapreds_list, aactuals_list)
  aauc <- performance(aapred, measure = "auc")
  aarocs <- performance(aapred, "tpr", "fpr")
  print(aauc@y.values)
  plot(aarocs, col = as.list(1:m), main = paste0("Train Set ROC Curves for feature set: ", names(my_learning_datasets_1000bp_4)[i]))
  abline(a = 0, b = 1, col = 6, lty = 3, lwd = 0.7)
  legend(x = "bottomright", cex = 0.6,
         legend = c(paste0("RF_down-Set4 ", "AUC: " , format(aauc@y.values[[1]], digits = 2))),
         fill = 1:m)
  aa_model_list <- list(RF_down1 = my_learning_datasets_results_RF_kappa_down_1000bp_4[[i]])
  
  
  aamodel_list_pr <- aa_model_list %>%
    map(.f = calc_auprc, data = (aa_train_data1))
  aamodel_list_pr_auc <- aamodel_list_pr %>%
    map(function(the_mod) the_mod$auc.integral)
  aamodel_list_pr_auc <- unlist(aamodel_list_pr_auc)
  
  aaresults_list_pr <- list(NA)
  num_mod <- 1
  
  for(the_pr in aamodel_list_pr){
    
    aaresults_list_pr[[num_mod]] <- data_frame(recall = the_pr$curve[, 1],
                                               precision = the_pr$curve[, 2],
                                               model = names(aamodel_list_pr)[num_mod])
    
    num_mod <- num_mod + 1
    
  }
  
  aaresults_df_pr <- bind_rows(aaresults_list_pr)
  
  aacustom_col <- colorspace::rainbow_hcl(length(aa_model_list))
  
  aa_gglist_3[[i]] <- ggplot(aes(x = recall, y = precision, group = model), data = aaresults_df_pr) +
    geom_line(aes(color = model), size = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = aacustom_col, labels = paste(names(aa_model_list) ,"AUC: " 
                                                             ,format(aamodel_list_pr_auc, digits = 2))) +
    geom_abline(intercept = sum(aa_train_data1$label_set == "Pos")/nrow(aa_train_data1),
                slope = 0, color = "gray", size = 1) +
    ggtitle(paste0("Training PR curve ", names(my_learning_datasets_1000bp_4)[i]))+
    theme_bw()
  
}
names(aa_gglist_3) <-  names(my_learning_datasets_1000bp_4)





data_partition_mut_GEMSTAT <- list()
aa_train_valid <- createDataPartition(y = my_learning_datasets_1000bp_4[[1]][, ncol(my_learning_datasets_1000bp_4[[1]])],
                          times = 1, p = 0.80,list = F)
aa_test <- rownames(my_learning_datasets_1000bp_4[[1]])[setdiff(c(1:nrow(my_learning_datasets_1000bp_4[[1]])), aa_train_valid)]
table(unlist(lapply(strsplit(aa_test, split="_"), "[[", 1)))
aa_train_validSet <- my_learning_datasets_1000bp_4[[1]][aa_train_valid,]
aa_trainind <- createDataPartition(y = aa_train_validSet[, ncol(aa_train_validSet)],
                                      times = 1, p = 0.75,list = F)
aa_train <- rownames(aa_train_validSet)[aa_trainind]
table(unlist(lapply(strsplit(aa_train, split="_"), "[[", 1)))

aa_validind <- setdiff(c(1:nrow(aa_train_validSet)), aa_trainind)
aa_valid <- rownames(aa_train_validSet)[aa_validind]
table(unlist(lapply(strsplit(aa_valid, split="_"), "[[", 1)))
data_partition_mut_GEMSTAT[[1]] <- aa_train
data_partition_mut_GEMSTAT[[2]] <- aa_valid
data_partition_mut_GEMSTAT[[3]] <- aa_test
names(data_partition_mut_GEMSTAT) <- c("train", "validation", "test")

table(unlist(lapply(strsplit(data_partition_mut_GEMSTAT[[1]], split="_"), "[[", 1)))
table(unlist(lapply(strsplit(data_partition_mut_GEMSTAT[[2]], split="_"), "[[", 1)))
table(unlist(lapply(strsplit(data_partition_mut_GEMSTAT[[3]], split="_"), "[[", 1)))

# train a RF model on training and validation sets together, test the results on the test set
# create rf train and test data
my_learning_datasets_1000bp_4_ForGEMSTAT <- my_learning_datasets_1000bp_4[c(1, 4, 7)]

save(list = c("my_learning_datasets_1000bp_4_ForGEMSTAT", "data_partition_mut_GEMSTAT"), 
     file = "Classification_set4_GEMSTAT.RData")

my_learning_datasets_1000bp_4_ForGEMSTAT_results <- list()
for(i in 1:length(my_learning_datasets_1000bp_4_ForGEMSTAT)){
  aanzv <- nearZeroVar(my_learning_datasets_1000bp_4_ForGEMSTAT[[i]], saveMetrics= F)
  if(length(aanzv) > 0){
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]][,-aanzv]
  }else{
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]]
  }
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
    aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
    aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
    if(length(aahighlyCorDescr) > 0){
      aa_all_data <- aa_all_data[,-aahighlyCorDescr]
    }
  }
  print("#########################################################################################")
  print(paste0("number of features in original feature set number ", i," named ", 
               names(my_learning_datasets_1000bp_4_ForGEMSTAT)[i] ,":"))
  print(ncol(my_learning_datasets_1000bp_4_ForGEMSTAT[[i]]))
  print("number of features after removing near zero variance and highly correlated features ( > 0.9)")
  print(ncol(aa_all_data))
  print("#########################################################################################")
  aatrainind <- match(c(data_partition_mut_GEMSTAT[[1]], data_partition_mut_GEMSTAT[[2]]),
                      rownames(my_learning_datasets_1000bp_4_ForGEMSTAT[[1]]))
  aagridLen <- c(5, 30, 50)
  aa_train_data <- aa_all_data[aatrainind, ]
  aa_test_data <-  aa_all_data[-aatrainind, ]
  aactrl <- trainControl(method = "repeatedcv",
                         number = 6,
                         repeats = 3,
                         classProbs = TRUE,
                         summaryFunction = mnLogLoss,
                         savePredictions = TRUE,
                         ## new option here:
                         sampling = "down", 
                         search = 'random')
  aaRF_down <- train(label_set ~ ., data = aa_train_data,
                     method = "rf",
                     #preProc = c("center", "scale"),
                     ntree = 1000,
                     metric = "logLoss",
                     trControl = aactrl,
                     tuneLength  = aagridLen[i])
  my_learning_datasets_1000bp_4_ForGEMSTAT_results[[i]] <- aaRF_down
}

save(list = c("my_learning_datasets_1000bp_4_ForGEMSTAT_results"),
     file = "classification_results_set_4_GEMSTAT_logloss.RData")



my_learning_datasets_1000bp_4_ForGEMSTAT_results_2 <- list()
for(i in 1:length(my_learning_datasets_1000bp_4_ForGEMSTAT)){
  aanzv <- nearZeroVar(my_learning_datasets_1000bp_4_ForGEMSTAT[[i]], saveMetrics= F)
  if(length(aanzv) > 0){
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]][,-aanzv]
  }else{
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]]
  }
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
    aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
    aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
    if(length(aahighlyCorDescr) > 0){
      aa_all_data <- aa_all_data[,-aahighlyCorDescr]
    }
  }
  print("#########################################################################################")
  print(paste0("number of features in original feature set number ", i," named ", 
               names(my_learning_datasets_1000bp_4_ForGEMSTAT)[i] ,":"))
  print(ncol(my_learning_datasets_1000bp_4_ForGEMSTAT[[i]]))
  print("number of features after removing near zero variance and highly correlated features ( > 0.9)")
  print(ncol(aa_all_data))
  print("#########################################################################################")
  aatrainind <- match(c(data_partition_mut_GEMSTAT[[1]], data_partition_mut_GEMSTAT[[2]]),
                      rownames(my_learning_datasets_1000bp_4_ForGEMSTAT[[1]]))
  aagridLen <- c(5, 30, 50)
  aa_train_data <- aa_all_data[aatrainind, ]
  aa_test_data <-  aa_all_data[-aatrainind, ]
  aactrl <- trainControl(method = "repeatedcv",
                         number = 6,
                         repeats = 3,
                         classProbs = TRUE,
                         summaryFunction = defaultSummary,
                         savePredictions = TRUE,
                         ## new option here:
                         sampling = "down", 
                         search = 'random')
  aaRF_down <- train(label_set ~ ., data = aa_train_data,
                     method = "rf",
                     #preProc = c("center", "scale"),
                     ntree = 1000,
                     metric = "Kappa",
                     trControl = aactrl,
                     tuneLength  = aagridLen[i])
  my_learning_datasets_1000bp_4_ForGEMSTAT_results_2[[i]] <- aaRF_down
}
save(list = c("my_learning_datasets_1000bp_4_ForGEMSTAT_results_2"),
     file = "classification_results_set_4_GEMSTAT_Kappa.RData")



names(my_learning_datasets_1000bp_4_ForGEMSTAT_results_2) <- c("Sum_LLR", "Sum_LR_plus_Adjacency", "kmer5")

aa_gglist_4 <- list()

for(i in 1:length(my_learning_datasets_1000bp_4_ForGEMSTAT)){
  aanzv <- nearZeroVar(my_learning_datasets_1000bp_4_ForGEMSTAT[[i]], saveMetrics= F)
  if(length(aanzv) > 0){
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]][,-aanzv]
  }else{
    aa_all_data <- my_learning_datasets_1000bp_4_ForGEMSTAT[[i]]
  }
  
  #  Identifying Correlated Predictors
  if(ncol(aa_all_data[, 1:(ncol(aa_all_data) - 1)]) > 1){
    aadescrCor <- cor(aa_all_data[, 1:(ncol(aa_all_data) - 1)])
    aahighlyCorDescr <- findCorrelation(aadescrCor, cutoff = .9)
    if(length(aahighlyCorDescr) > 0){
      aa_all_data <- aa_all_data[,-aahighlyCorDescr]
    }
  }
  
  aatrainind <- match(c(data_partition_mut_GEMSTAT[[1]], data_partition_mut_GEMSTAT[[2]]),
                      rownames(my_learning_datasets_1000bp_4_ForGEMSTAT[[1]]))
  aagridLen <- c(5, 30, 50)
  aa_train_data <- aa_all_data[aatrainind, ]
  aa_test_data <-  aa_all_data[-aatrainind, ]
  
  aaindex_pos <- aa_test_data$label_set == "Pos"
  aaindex_neg <- aa_test_data$label_set == "Neg"
  
  aarfdown <- predict(my_learning_datasets_1000bp_4_ForGEMSTAT_results_2[[i]], 
                      aa_test_data, type = "prob")
  
  aapreds_list <- list(aarfdown$Pos)
  # List of actual values (same for all)
  m <- length(aapreds_list)
  aactuals_list <- rep(list(aa_test_data$label_set), m)
  # Plot the ROC curves
  aapred <- prediction(aapreds_list, aactuals_list)
  aauc <- performance(aapred, measure = "auc")
  aarocs <- performance(aapred, "tpr", "fpr")
  print(aauc@y.values)
  plot(aarocs, col = as.list(1:m), main = paste0("Test Set ROC Curves for feature set: ",
                                                 names(my_learning_datasets_1000bp_4_ForGEMSTAT_results_2)[i]))
  abline(a = 0, b = 1, col = 6, lty = 3, lwd = 0.7)
  legend(x = "bottomright", cex = 0.6,
         legend = c(paste0("RF_down_Kappa ", "AUC: " , format(aauc@y.values[[1]], digits = 2))),
         fill = 1:m)
  aa_model_list <- list(RF_down = my_learning_datasets_1000bp_4_ForGEMSTAT_results_2[[i]])
  
  
  aamodel_list_pr <- aa_model_list %>%
    map(.f = calc_auprc, data = (aa_test_data))
  aamodel_list_pr_auc <- aamodel_list_pr %>%
    map(function(the_mod) the_mod$auc.integral)
  aamodel_list_pr_auc <- unlist(aamodel_list_pr_auc)
  
  aaresults_list_pr <- list(NA)
  num_mod <- 1
  
  for(the_pr in aamodel_list_pr){
    
    aaresults_list_pr[[num_mod]] <- data_frame(recall = the_pr$curve[, 1],
                                               precision = the_pr$curve[, 2],
                                               model = names(aamodel_list_pr)[num_mod])
    
    num_mod <- num_mod + 1
    
  }
  
  aaresults_df_pr <- bind_rows(aaresults_list_pr)
  
  aacustom_col <- colorspace::rainbow_hcl(length(aa_model_list))
  
  aa_gglist_4[[i]] <- ggplot(aes(x = recall, y = precision, group = model), data = aaresults_df_pr) +
    geom_line(aes(color = model), size = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = aacustom_col, labels = paste(names(aa_model_list) ,"AUC: " ,
                                                             format(aamodel_list_pr_auc, digits = 2))) +
    geom_abline(intercept = sum(aa_test_data$label_set == "Pos")/nrow(aa_test_data),
                slope = 0, color = "gray", size = 1) +
    ggtitle(paste0("PR curve ", names(my_learning_datasets_1000bp_4_ForGEMSTAT_results_2)[i]))+
    theme_bw()
  
}
names(aa_gglist_4) <-  names(my_learning_datasets_1000bp_4_ForGEMSTAT_results_2)
aa_gglist_4[[2]]
print(my_learning_datasets_1000bp_4_ForGEMSTAT_results_2[[2]])
