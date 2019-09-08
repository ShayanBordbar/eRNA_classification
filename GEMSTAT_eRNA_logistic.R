library("PRROC")
library(rjson)
#######################################################################################################################
#######################################################################################################################
#####################################        FUNCTIONS         ########################################################
#######################################################################################################################
#######################################################################################################################
read_parameters_GEMSTAT_indiv <- function(par_file_name, 
                                          .convert_to_mat=T){
  # par_file_name : is the name of the pararmeter file
  # .convert_to_mat : is whether to convert from json to a matrix like structure
  my_par <- fromJSON(file = par_file_name)
  if(typeof(my_par$qbtm) != "list"){
    my_par$qbtm <- list(my_par$qbtm)
  }
  
  if(.convert_to_mat){
    # TFs
    nu_TFs <- length(my_par$tfs)
    TF_annot <- numeric(nu_TFs)
    TF_binding <- numeric(nu_TFs)
    TF_alpha <- numeric(nu_TFs)
    names(TF_annot) <- names(my_par$tfs)
    names(TF_binding) <- names(my_par$tfs)
    names(TF_alpha) <- names(my_par$tfs)
    for(cur_tf in 1:nu_TFs){
      TF_annot[cur_tf] <- my_par$tfs[[cur_tf]]$annot_thresh
      TF_binding[cur_tf] <- my_par$tfs[[cur_tf]]$maxbind
      TF_alpha[cur_tf] <- my_par$tfs[[cur_tf]]$alpha_a
    } 
    # intercations
    Coop_par <- unlist(my_par$inter)
    # qBTM
    qbtm_par <- unlist(my_par$qbtm)
    # beta
    beta_par <- my_par$enh[[1]]$beta
    # logistic
    logistic_par <- unlist(my_par$log_Reg[[1]])
    
    return(list(annot = TF_annot, binding = TF_binding,
                alpha = TF_alpha, coop = Coop_par, qbtm=qbtm_par,
                beta=beta_par, logistic = logistic_par))
  }else{
    return(my_par)
  }
}
#######################################################################################################################
#######################################################################################################################
read_parameters_GEMSTAT_ensemble <- function(par_directory, convert_to_mat){
  library(rjson)
  # function to read parameters of an ensemble
  # par_directory : is the name of the ensemble parameter directory
  # convert_to_mat : is whether to convert from json to a matrix like structure
  
  aa_tpar_f <- list.files(path = par_directory, 
                          full.names = T)
  parameter_holder_list <- list()
  for(cur_par in 1:length(aa_tpar_f)){
    parameter_holder_list[[cur_par]] <- read_parameters_GEMSTAT_indiv(par_file_name = aa_tpar_f[cur_par],
                                                                      .convert_to_mat=convert_to_mat)
  }
  aa_tpar_fn <- list.files(path = par_directory, 
                           full.names = F)
  aa_tpar_fn_sp <- unlist(lapply(strsplit(aa_tpar_fn, split = "\\."), "[[", 1))
  names(parameter_holder_list) <- aa_tpar_fn_sp
  if(convert_to_mat){
    annot_mat <- do.call(rbind,  lapply(parameter_holder_list, "[[", 1))
    binding_mat <- do.call(rbind,  lapply(parameter_holder_list, "[[", 2))
    alpha_mat <- do.call(rbind,  lapply(parameter_holder_list, "[[", 3))
    coop_mat <- do.call(rbind,  lapply(parameter_holder_list, "[[", 4))
    qbtms <- do.call(c,  lapply(parameter_holder_list, "[[", 5))
    betas <- do.call(c,  lapply(parameter_holder_list, "[[", 6))
    logistic_mat <- do.call(rbind,  lapply(parameter_holder_list, "[[", 7))
    return(list(annot = annot_mat, binding=binding_mat, alpha=alpha_mat, coop=coop_mat,
           qbtm = qbtms, beta = betas, logistic=logistic_mat))
  }else{
    return(parameter_holder_list)
  }
}
#######################################################################################################################
#######################################################################################################################
# example
aaa <- read_parameters_GEMSTAT_ensemble(par_directory = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_1/Trained_par",
                                        convert_to_mat = T)
#######################################################################################################################
#######################################################################################################################
# read output of ensemble, compute AUROC, AUPRC and have plots ready
read_output_train_test_GEMSTAT_indiv <- function(output_file, .plot=F){
  
  output_all <- read.table(file = output_file ,
                           header = T, stringsAsFactors = F)
  output_gt <- output_all[seq(1, nrow(output_all), 2),]
  output_model <- output_all[seq(2, nrow(output_all), 2),]
  if(.plot){
    boxplot(output_model$X1 ~ output_gt$X1,
            xlab = "label", ylab = "predicted exp", 
            main = output_file)
  }
  my_ROC <- roc.curve(scores.class0 = output_model$X1[output_gt$X1 == 1],
                      scores.class1 = output_model$X1[output_gt$X1 == 0], 
                      curve = TRUE)
  my_PRC <- pr.curve(scores.class0 = output_model$X1[output_gt$X1 == 1],
                     scores.class1 = output_model$X1[output_gt$X1 == 0],
                     curve = TRUE)
  return(list(GT=output_gt,
              pred=output_model,
              ROC_curve=my_ROC,
              PRC_curve = my_PRC))
  
}
#######################################################################################################################
#######################################################################################################################

read_output_train_test_GEMSTAT_ensemble <- function(output_dir, should_plot=F, validation=F){
  # validation: if True looks for and analyzes validation results as well
  all_files <- list.files(output_dir, full.names = T)
  test_files <- list.files(output_dir, pattern = "*_test.txt", full.names = T)
  valid_files <- list.files(output_dir, pattern = "*_valid.txt", full.names = T)
  train_files <- setdiff(all_files, c(test_files, valid_files))
  
  all_files_n <- list.files(output_dir, full.names = F)
  test_files_n <- list.files(output_dir, pattern = "*_test.txt", full.names = F)
  test_files_n_sp <- unlist(lapply(strsplit(test_files_n, split = "\\."), "[[", 1))
  
  valid_files_n <- list.files(output_dir, pattern = "*_valid.txt", full.names = F)
  valid_files_n_sp <- unlist(lapply(strsplit(valid_files_n, split = "\\."), "[[", 1))
  
  train_files_n <- setdiff(all_files_n, c(test_files_n, valid_files_n))
  train_files_n_sp <- unlist(lapply(strsplit(train_files_n, split = "\\."), "[[", 1))
  train_Results <- list()
  test_Results <- list()
  valid_Results <- list()
  stopifnot(length(test_files_n) == length(train_files_n))
  for(cur_tr in 1:length(train_files)){
    train_Results[[cur_tr]] <- read_output_train_test_GEMSTAT_indiv(output_file = train_files[cur_tr],
                                                                    .plot = should_plot)
  }
  names(train_Results) <- train_files_n_sp
  
  for(cur_ts in 1:length(test_files)){
    test_Results[[cur_ts]] <- read_output_train_test_GEMSTAT_indiv(output_file = test_files[cur_ts],
                                                                    .plot = should_plot)
  }
  names(test_Results) <- test_files_n_sp
  
  if(validation){
    for(cur_va in 1:length(valid_files)){
      valid_Results[[cur_va]] <- read_output_train_test_GEMSTAT_indiv(output_file = valid_files[cur_va],
                                                                     .plot = should_plot)
    }
    names(valid_Results) <- valid_files_n_sp
    valid_ROC_list <- lapply(valid_Results, "[[", 3)
    valid_PRC_list <- lapply(valid_Results, "[[", 4)
    Valid_ROC <-  unlist(lapply(valid_ROC_list, "[[", 2))
    Valid_PRC <-  unlist(lapply(valid_PRC_list, "[[", 2))
    names(Valid_ROC) <- valid_files_n_sp
    names(Valid_PRC) <- valid_files_n_sp
  }
  
  Train_ROC_list <- lapply(train_Results, "[[", 3)
  Train_PRC_list <- lapply(train_Results, "[[", 4)
  
  Test_ROC_list <- lapply(test_Results, "[[", 3)
  Test_PRC_list <- lapply(test_Results, "[[", 4)
  
  stopifnot(length(Train_ROC_list) == length(Test_ROC_list))
  
  Train_ROC <-  unlist(lapply(Train_ROC_list, "[[", 2))
  names(Train_ROC) <- train_files_n_sp
  Test_ROC <-  unlist(lapply(Test_ROC_list, "[[", 2))
  names(Test_ROC) <- test_files_n_sp
  stopifnot(length(Train_ROC) == length(Test_ROC))
  Train_PRC <- unlist(lapply(Train_PRC_list, "[[", 2))
  names(Train_PRC) <- train_files_n_sp
  Test_PRC <- unlist(lapply(Test_PRC_list, "[[", 2))
  names(Test_PRC) <- test_files_n_sp
  
  par(mfrow = c(1, 2), mar= c(4,4,4,4))
  plot(Test_ROC, Train_ROC, main = "AUROC", xlab = "test", ylab = "training", pch = 16, cex = 0.8)
  plot(Test_PRC, Train_PRC, main = "AUPRC", xlab = "test", ylab = "training", pch = 16, cex = 0.8)
  if(validation){
    return(list(train_results =train_Results,
                test_results=test_Results,
                train_ROC = Train_ROC,
                test_ROC=Test_ROC,
                train_PRC=Train_PRC, 
                test_PRC=Test_PRC, 
                valid_results = valid_Results,
                valid_ROC = Valid_ROC,
                valid_PRC=Valid_PRC))
  }
  return(list(train_results =train_Results,
              test_results=test_Results,
              train_ROC = Train_ROC,
              test_ROC=Test_ROC,
              train_PRC=Train_PRC, 
              test_PRC=Test_PRC))
}
#######################################################################################################################
#######################################################################################################################
# example
aaa <- read_output_train_test_GEMSTAT_ensemble(output_dir = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_1/Out", 
                                               should_plot = F)
#######################################################################################################################
#######################################################################################################################

GEMSTAT_create_KD_paramFile <- function(param_file_address, .TF_names ,TF_to_KD, Coop_to_KD, ouput_address){
  # this function creates a parameter file for KD experiemnt given a trained param file as well as the names of the TFs to KD
  # param_file_address : character indicating the WT par file address
  # .TF_names : char vector contatining the name of all TFs
  # TF_to_KD : name of the TFs that are going to be KD
  # Coop_to_KD : a matrix where each row indicates an interaction to KD, it has two columns. entries are the index of the TF in TF_names vector
  # ouput_address : is the address where the output parameter file should be written to
  if(length(Coop_to_KD) == 2){
    Coop_to_KD <- matrix(Coop_to_KD, nrow = 1)
  }
  stopifnot(is.character(param_file_address),
            length(.TF_names) > 0,
            (length(TF_to_KD) > 0 | length(Coop_to_KD) > 0),
            is.character(ouput_address), 
            (is.matrix(Coop_to_KD) | length(Coop_to_KD) == 0),
            (max(Coop_to_KD) <= length(.TF_names) | length(Coop_to_KD) == 0), 
            (all(TF_to_KD %in% .TF_names) | length(TF_to_KD) == 0))
  
  cur_par_file <- fromJSON(file = param_file_address)
  if(typeof(cur_par_file$qbtm) != "list"){
    cur_par_file$qbtm <- list(cur_par_file$qbtm)
  }
  if(length(TF_to_KD) > 0){
    for(cur_tf in 1:length(TF_to_KD)){
      cur_par_file$tfs[[TF_to_KD[cur_tf]]]$annot_thresh <- 0
    }
  }
  if(length(Coop_to_KD) >= 2){
    all_coops <- names(cur_par_file$inter)
    for(cur_row in 1:nrow(Coop_to_KD)){
      ref_name <- paste(.TF_names[Coop_to_KD[cur_row, ]], collapse = ":")
      alt_name <- paste(.TF_names[rev(Coop_to_KD[cur_row, ])], collapse = ":")
      cur_existance <- c(ref_name, alt_name) %in% all_coops
      stopifnot(sum(cur_existance) == 1 | (sum(cur_existance) == 2 & ref_name==alt_name))
      cur_par_file$inter[[which(all_coops %in% c(ref_name, alt_name))]] <- 1
    }
  }
  cat(toJSON(cur_par_file), append = F, 
      file = ouput_address)
}
#######################################################################################################################
#######################################################################################################################
#example
aa_names <- c("ESR1_2", "ESR1_3","FOXA1","JUN_1","LEF1","NFIB","NKX3_1","NR2F2","NR3C1","NR5A2",
              "PBX1","PGR","PPARD","PPARG","RARA","RUNX1","SP1","YBX1")
aacoop <- rbind(c(1, 1), c(1, 2), c(1, 3), c(1, 8),
                c(1, 12), c(1,13), c(1, 14), c(1, 18),
                c(8, 8), c(8, 16))
GEMSTAT_create_KD_paramFile(param_file_address = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Trained_par/log_par_1_1.txt.Filter",
                            .TF_names = aa_names,
                            TF_to_KD = c("ESR1_2", "RARA"),
                            Coop_to_KD = numeric(0),
                            ouput_address = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Test_kd.par"
)
#######################################################################################################################
#######################################################################################################################
GEMSTAT_create_KD_paramFile_ensemble <- function(param_addresses, 
                                                 KD_exp_list, 
                                                 TF_names,
                                                 dest_directory){
  # to create parameters for KD experiments for multiple models, it also creates the necessary directories
  # param_addresses: is a char vector, each entry is the address for a parameter file
  # KD_exp_list: is a named list. each entry is a list with two entries: TF, COOP. TF is a 
  #  character vector containing the name of the TFs to be Knocked down. COOP is a matrix where each
  #  row indicates an interaction to KD, it has two columns. entries are the index of the TF in TF_names vector
  # TF_names: character vector containing the names of the TFs in the model
  # dest_directory: the directory to work in
  stopifnot(is.character(param_addresses), 
            length(param_addresses) > 0,
            is.list(KD_exp_list),
            length(names(KD_exp_list)) > 0,
            all(unlist(lapply(KD_exp_list, length)) %in% 2),
            length(TF_names) > 0, 
            is.character(TF_names), 
            is.character(dest_directory), 
            length(dest_directory) == 1)
  
  library(rjson)
  outside_directory <- getwd()
  on.exit(setwd(outside_directory))
  setwd(dest_directory)
  
  if (! dir.exists("KnockDown")){
    system("mkdir KnockDown")
  }
  KD_directory <- paste0(dest_directory,"/KnockDown")
  setwd(KD_directory)
  for(cur_kd in 1:length(KD_exp_list)){
    if (! dir.exists(names(KD_exp_list)[cur_kd])){
      system(paste0("mkdir ", names(KD_exp_list)[cur_kd]))
      system(paste0("mkdir ", names(KD_exp_list)[cur_kd], "/Par"))
      system(paste0("mkdir ", names(KD_exp_list)[cur_kd], "/Out"))
    }
    cur_dir <- paste0(KD_directory, "/", names(KD_exp_list)[cur_kd])
    for(cur_model in 1:length(param_addresses)){
      aaname <- strsplit(param_addresses[cur_model], split = "\\/")[[1]]
      aaname <- aaname[length(aaname)]
      cur_address <- paste0(cur_dir, "/Par/",aaname,"_KD")
      
      GEMSTAT_create_KD_paramFile(param_file_address = param_addresses[cur_model],
                                  .TF_names = TF_names,
                                  TF_to_KD = KD_exp_list[[cur_kd]]$TF,
                                  Coop_to_KD = KD_exp_list[[cur_kd]]$COOP,
                                  ouput_address = cur_address)
    } # end of loop over models
  } # end of loop over KD experiments
  
}
#######################################################################################################################
#######################################################################################################################
#example
aa_names <- c("ESR1_2", "ESR1_3","FOXA1","JUN_1","LEF1","NFIB","NKX3_1","NR2F2","NR3C1","NR5A2",
              "PBX1","PGR","PPARD","PPARG","RARA","RUNX1","SP1","YBX1")
aacoop <- rbind(c(1, 1), c(1, 2), c(1, 3), c(1, 8),
                c(1, 12), c(1,13), c(1, 14), c(1, 18),
                c(8, 8), c(8, 16))
aaKDexp <- list()
for(i in 1:5){
  aaKDexp[[i]] <- list()
  aaKDexp[[i]]$TF <- aa_names[i]
  aaKDexp[[i]]$COOP <- aacoop[i,]
}
names(aaKDexp) <- c(1:length(aaKDexp))
aa_address <- c("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Trained_par/log_par_1_1.txt.Filter",
                "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Trained_par/log_par_10_1.txt.Filter",
                "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Trained_par/log_par_9_1.txt.Filter")
aadest_directory <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble"
GEMSTAT_create_KD_paramFile_ensemble(param_addresses = aa_address, 
                                     KD_exp_list = aaKDexp, 
                                     TF_names=aa_names,
                                     dest_directory = aadest_directory)
#######################################################################################################################
#######################################################################################################################
GEMSTAT_read_KD <- function(KD_directory, kd_exp_names, WT_output){
  # WT_output : is the output of read_output_train_test_GEMSTAT_ensemble function for the full ensemble from which these KDs where constructed
  # KD_directory: is the KnockDown directory
  # kd_exp_names: is a character vector containing the names of KD experiments we want
  outside_directory <- getwd()
  on.exit(setwd(outside_directory))
  setwd(KD_directory)
  output_list <- list()
  model_expmat_list <- list()
  for(cur_exp in 1:length(kd_exp_names)){
    output_list[[cur_exp]] <- list()
    cur_outs <- list.files(path = paste0(kd_exp_names[cur_exp], "/Out/"),
                           full.names = T)
    cur_outs_n <- list.files(path = paste0(kd_exp_names[cur_exp], "/Out/"),
                           full.names = F)
    cur_outs_sp <- unlist(lapply(lapply(strsplit(cur_outs_n, split= "_"), "[", c(1,3,4,5)), paste, collapse= "_"))
    for(cur_model in 1:length(cur_outs_sp)){
      output_list[[cur_exp]][[cur_model]] <- read_output_train_test_GEMSTAT_indiv(output_file = cur_outs[cur_model],
                                                                      .plot = F)
      if(cur_exp == 1){ #initialize the matrix for each model
        model_expmat_list[[cur_model]] <- matrix(nrow = nrow(output_list[[cur_exp]][[cur_model]]$pred),
                                                 ncol = (length(kd_exp_names) + 2))
        colnames(model_expmat_list[[cur_model]]) <- c("GT", "WT" ,kd_exp_names)
        rownames(model_expmat_list[[cur_model]]) <- rownames(output_list[[cur_exp]][[cur_model]]$pred)
        model_expmat_list[[cur_model]][, 1] <- output_list[[cur_exp]][[cur_model]]$GT$X1
        model_expmat_list[[cur_model]][, 2] <- WT_output$train_results[[cur_outs_sp[cur_model]]]$pred$X1
      }
      model_expmat_list[[cur_model]][, cur_exp + 2] <- output_list[[cur_exp]][[cur_model]]$pred$X1
    } # end of loop over models
    if(cur_exp == 1){
      names(model_expmat_list) <- cur_outs_sp
    }
    names(output_list[[cur_exp]]) <- cur_outs_sp
  } # end of loop over KD experiments
  names(output_list) <- kd_exp_names
  return(list(Output_raw_list =output_list, 
         output_mat_list =model_expmat_list))
}
list(GT=output_gt,
     pred=output_model,
     ROC_curve=my_ROC,
     PRC_curve = my_PRC)
"output_log_par_81_1_KD.txt"
E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$train_results$output_par_1_1$pred$X1
#######################################################################################################################
#######################################################################################################################

aa <- GEMSTAT_read_KD(KD_directory = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/KnockDown",
                      kd_exp_names = as.character(c(1:4)),
                      WT_output = E_RNA_GEMSTAT_Ensemble_Outlist[[7]])
heatmap.2(x = aa$output_mat_list$output_par_1033_1[,2:6], Rowv = T, 
          Colv = T, dendrogram = "both", trace = "none")
#######################################################################################################################
#######################################################################################################################
# next function starts here

#######################################################################################################################
#######################################################################################################################

# GEMSTAT logistic eRNA
########
# aalab1 <- unlist(lapply((strsplit(aa_seqnames1, split = "_")), "[[", 2))
# aalab2 <- unlist(lapply((strsplit(aa_seqnames2, split = "_")), "[[", 2))
# 
# aalab1[aalab1 == "pos"] <- 1
# aalab1[aalab1 == "neg"] <- 0
# aalab1 <- as.numeric(aalab1)
# aalab1 <- matrix(aalab1, nrow = length(aalab1))
# aalab2[aalab2 == "pos"] <- 1
# aalab2[aalab2 == "neg"] <- 0
# aalab2 <- as.numeric(aalab2)
# aalab2 <- matrix(aalab2, nrow = length(aalab2))
# rownames(aalab1) <- aa_seqnames1
# rownames(aalab2) <- aa_seqnames2
# 
# aalab1 <- aalab1[1:(which(rownames(aalab1)=="test_neg_7_1") - 1),]
# aalab1 <- matrix(aalab1, nrow = length(aalab1))
# aalab2 <- aalab2[1:(which(rownames(aalab2)=="test_neg_5_1") - 1),]
# aalab2 <- matrix(aalab2, nrow = length(aalab2))
# rownames(aalab1) <- aa_seqnames1[1:length(aalab1)]
# rownames(aalab2) <- aa_seqnames2[1:length(aalab2)]
# 
# aalab1[aalab1 == 0] <- 0.258
# aalab1[aalab1 == 1] <- 0.741
# cat(c("Rows", paste0("1", "\n")), sep = "\t",
#     append = F, file = "weight_1_train.txt")
# for(i in 1:length(aalab1)){
#   cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
#                                   rep("\n", as.integer(i != length(aalab1)))), 
#       append = T, file = "weight_1_train.txt", sep = "\t" )
# }
# 
# 
# 
# aalab2[aalab2 == 0] <- 0.349
# aalab2[aalab2 == 1] <- 0.650
# cat(c("Rows", paste0("1", "\n")), sep = "\t",
#     append = F, file = "weight_2_train.txt")
# for(i in 1:length(aalab2)){
#   cat(rownames(aalab2)[i], paste0(aalab2[i, 1],
#                                   rep("\n", as.integer(i != length(aalab2)))), 
#       append = T, file = "weight_2_train.txt", sep = "\t" )
# }
################# # create inputs for GEMSTAT -> test case
MotifWriter(TF.motifs.Expanded_new_pseudo_padded[1:5],
            pseudo = 0.001, output.File.Name = "motif_list_NN_5")

# TF expression
aaTFexp <- rep(1, length(TF.motifs.Expanded_new_pseudo_padded[1:5]))
names(aaTFexp) <- names(TF.motifs.Expanded_new_pseudo_padded)[1:5]
cat(c("Rows", paste0("1", "\n")), sep = "\t", append = F, file = "NN_TF_exp_5.tab")
for(i in 1:length(aaTFexp)){
  cat(names(aaTFexp)[i], paste0(aaTFexp[i], rep("\n", as.integer(i != length(aaTFexp)))), 
      append = T, file = "NN_TF_exp_5.tab", sep = "\t" )
}
par_ff_lb_ub_coop_Writer_multi_enh(.TFnames=names(TF.motifs.Expanded_new_pseudo_padded),
                                   .annotation_thresh=numeric(0),
                                   .annotation_range=numeric(0),
                                   .initial_bind_w=numeric(0),
                                   .bind_w_range=numeric(0),
                                   .initial_alpha=numeric(0),
                                   .initial_log_reg_bias = numeric(0),
                                   .initial_log_reg_coeff = numeric(0),
                                   .initial_log_reg_bias_range = numeric(0),
                                   .initial_log_reg_coeff_range = numeric(0),
                                   .alpha_range=numeric(0),
                                   .coop_tf_mat=numeric(0),
                                   .initial_coop_weight=numeric(0),
                                   .coop_weight_range=numeric(0),
                                   .coop_type=character(0),
                                   .coop_dist=integer(0),
                                   .coop_orientation=integer(0),
                                   one_qbtm_per_enh=F,
                                   .initial_qBTM=numeric(0),
                                   .qBTMrange=numeric(0),
                                   one_beta_per_enh=F,
                                   .initial_pi_beta=numeric(0), 
                                   .pi_beta_upper=numeric(0),
                                   .pi_beta_lower=numeric(0),
                                   .annotation_thresh_ff=integer(0),
                                   .initial_bind_w_ff=integer(0),
                                   .initial_alpha_ff=integer(0),
                                   .initial_coop_weight_ff=integer(0),
                                   .initial_qBTM_ff=integer(0),
                                   .initial_pi_beta_ff=integer(0),
                                   .filename_start=character(0),
                                   .filename_ff=character(0),
                                   .filename_upper=character(0),
                                   .filename_lower=character(0),
                                   .file_name_coop=character(0),
                                   .initial_log_reg_bias_ff = 1,
                                   .initial_log_reg_coeff_ff = 1,
                                   .nu_enhacners=5,
                                   ensemble_mode = F,
                                   .nu_samples = 2,
                                   .annotation_thresh_ff_ens = c(0,0,0),
                                   .initial_coop_weight_ff_ens = c(0,0),
                                   .initial_qBTM_ff_ens = 0,
                                   .initial_log_reg_bias_ff_ens = 0,
                                   .initial_log_reg_coeff_ff_ens = 0,
                                   Write_to_file = T,
                                   create_folders = T,
                                   TF_role_decided = F,
                                   factor_role_vector=integer(0),  
                                   .create_bounds=T,
                                   .create_par=T,
                                   .create_ff=T,
                                   .create_coop=T,
                                   logistic_params = T)
#################
# # read GEMSTAT results
# 
# aa_GEM_1 <- read.table(file = "~/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/e_RNA_test/test_1_train.out",
#                        header = T, stringsAsFactors = F)
# aa_GEM_2 <- read.table(file = "~/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/e_RNA_test/test_1_train_weighted.out",
#                        header = T, stringsAsFactors = F)
# aa_GEM_3 <- read.table(file = "~/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/e_RNA_test/logistic/test_1_train_weighted.out",
#                        header = T, stringsAsFactors = F)
# aa_GEM_4 <- read.table(file = "~/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/e_RNA_test/logistic/train_logisitic_weighted.out",
#                        header = T, stringsAsFactors = F)
# aa_GEM_5 <- read.table(file = "~/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/e_RNA_test/logistic/Test2/train_logisitic_weighted_2.out",
#                        header = T, stringsAsFactors = F)
# 
# aa_GEM_1$Rows[1:10]
# aa1_1 <- aa_GEM_1[seq(2, nrow(aa_GEM_1), 2),]
# aa1_gt <- aa_GEM_1[seq(1, nrow(aa_GEM_1), 2),]
# boxplot(aa1_1$X1 ~ aa1_gt$X1)
# aa1_2 <- aa_GEM_2[seq(2, nrow(aa_GEM_2), 2),]
# aa1_gt2 <- aa_GEM_2[seq(1, nrow(aa_GEM_2), 2),]
# boxplot(aa1_2$X1 ~ aa1_gt2$X1)
# aa1_3 <- aa_GEM_3[seq(2, nrow(aa_GEM_3), 2),]
# aa1_gt3 <- aa_GEM_3[seq(1, nrow(aa_GEM_3), 2),]
# boxplot(aa1_3$X1 ~ aa1_gt3$X1)
# aa1_4 <- aa_GEM_4[seq(2, nrow(aa_GEM_4), 2),]
# aa1_gt4 <- aa_GEM_4[seq(1, nrow(aa_GEM_4), 2),]
# boxplot(aa1_4$X1 ~ aa1_gt4$X1)
# aa1_5 <- aa_GEM_5[seq(2, nrow(aa_GEM_5), 2),]
# aa1_gt5 <- aa_GEM_5[seq(1, nrow(aa_GEM_5), 2),]
# boxplot(aa1_5$X1 ~ aa1_gt5$X1)
#######################################################################################################################
#######################################################################################################################
################# ################# Trying smaller sets of sequences and TFs
################# CASE 1
# write GEMSTAT input with 40 training and 20 test sequences. each containing equal number of positives and negatives

aatpos <- sample(Positive_set_seq_list_char_1000bp[[4]], size = 30, replace = F)
aatneg <- sample(Negative_set_seq_list_char_1000bp[[4]], size = 30, replace = F)
names(aatpos) <- paste("pos", c(1:length(aatpos)), sep = "_")
names(aatneg) <- paste("neg", c(1:length(aatneg)), sep = "_")
GEMSTAT_ToyTrainSeq <- sample(x = c(aatpos[1:20], aatneg[1:20]), size = 40, replace = F)
GEMSTAT_ToyTestSeq <- c(aatpos[21:30], aatneg[21:30])
DNAStringSet
# Sequence
library(ShortRead)
writeFasta(DNAStringSet(GEMSTAT_ToyTrainSeq), width = 1000,
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_1/Sequence/GEMSTAT_ToyTrainSeq.fa")
writeFasta(DNAStringSet(GEMSTAT_ToyTestSeq), width = 1000, 
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_1/Sequence/GEMSTAT_ToyTestSeq.fa")

# Motif
aa_names <- sort(c("ESR1_2", "SP1", "PGR", "FOXA1", "GATA3", "JUN_2",  "MYC",   "PBX1", "PPARG", "RARA"))


aa_coop_name[[1]] <- rbind(c("ESR1_2", "PGR"), 
                           c("ESR1_2", "PPARG"),
                           c("ESR1_2", "RARA"), 
                           c("ESR1_2", "SP1"),
                           c("ESR1_2", "FOXA1"))


MotifWriter(TF.motifs.Expanded_pseudo_count[aa_names],
            pseudo = 0.001, 
            output.File.Name = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_1/Motif/motifs_T1")


# TF expression
aaTFexp <- rep(1, length(aa_names))
names(aaTFexp) <- aa_names
cat(c("Rows", paste0("1", "\n")), sep = "\t", append = F,
    file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_1/TF_exp/TF_exp_1.tab")
for(i in 1:length(aaTFexp)){
  cat(names(aaTFexp)[i], paste0(aaTFexp[i], rep("\n", as.integer(i != length(aaTFexp)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_1/TF_exp/TF_exp_1.tab",
      sep = "\t" )
}
# labels
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_ToyTrainSeq), split = "_")), "[[", 1))

aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_ToyTrainSeq)
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_1/Labels/Label_train_1.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_1/Labels/Label_train_1.txt", sep = "\t" )
}
# label for test
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_ToyTestSeq), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_ToyTestSeq)
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_1/Labels/Label_test_1.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_1/Labels/Label_test_1.txt", sep = "\t" )
}

#weights (in this case all equal)
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_ToyTrainSeq), split = "_")), "[[", 1))

aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_ToyTrainSeq)

aalab1[aalab1 == 0] <- 0.5
aalab1[aalab1 == 1] <- 0.5
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_1/weight_1_train.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_1/weight_1_train.txt", sep = "\t" )
}
# weights for test just in case of error
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_ToyTestSeq), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_ToyTestSeq)
aalab1[aalab1 == 0] <- 0.5
aalab1[aalab1 == 1] <- 0.5
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_1/weight_1_test.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_1/weight_1_test.txt", sep = "\t" )
}


# coop
aacoop <- rbind(c(1, 2), c(1, 3), c(1, 4), c(1,9), c(1, 10))
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_1")
par_ff_lb_ub_coop_Writer_multi_enh(.TFnames=aa_names,
                                   .annotation_thresh=numeric(0),
                                   .annotation_range=numeric(0),
                                   .initial_bind_w=numeric(0),
                                   .bind_w_range=numeric(0),
                                   .initial_alpha=numeric(0),
                                   .initial_log_reg_bias = numeric(0),
                                   .initial_log_reg_coeff = numeric(0),
                                   .initial_log_reg_bias_range = numeric(0),
                                   .initial_log_reg_coeff_range = numeric(0),
                                   .alpha_range=numeric(0),
                                   .coop_tf_mat=aacoop,
                                   .initial_coop_weight=numeric(0),
                                   .coop_weight_range=numeric(0),
                                   .coop_type=character(0),
                                   .coop_dist=integer(0),
                                   .coop_orientation=integer(0),
                                   one_qbtm_per_enh=F,
                                   .initial_qBTM=numeric(0),
                                   .qBTMrange=numeric(0),
                                   one_beta_per_enh=F,
                                   .initial_pi_beta=numeric(0), 
                                   .pi_beta_upper=numeric(0),
                                   .pi_beta_lower=numeric(0),
                                   .annotation_thresh_ff=integer(0),
                                   .initial_bind_w_ff=integer(0),
                                   .initial_alpha_ff=integer(0),
                                   .initial_coop_weight_ff=integer(0),
                                   .initial_qBTM_ff=integer(0),
                                   .initial_pi_beta_ff=integer(0),
                                   .filename_start=character(0),
                                   .filename_ff=character(0),
                                   .filename_upper=character(0),
                                   .filename_lower=character(0),
                                   .file_name_coop=character(0),
                                   .initial_log_reg_bias_ff = 1,
                                   .initial_log_reg_coeff_ff = 1,
                                   .nu_enhacners=40,
                                   ensemble_mode = F,
                                   .nu_samples = 2,
                                   .annotation_thresh_ff_ens = rep(0, 10),
                                   .initial_coop_weight_ff_ens = rep(0, 5),
                                   .initial_qBTM_ff_ens = 0,
                                   .initial_log_reg_bias_ff_ens = 0,
                                   .initial_log_reg_coeff_ff_ens = 0,
                                   Write_to_file = T,
                                   create_folders = T,
                                   TF_role_decided = F,
                                   factor_role_vector=integer(0),  
                                   .create_bounds=T,
                                   .create_par=T,
                                   .create_ff=T,
                                   .create_coop=T,
                                   logistic_params = T)

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")

# read training res
aa_GEM_10 <- read.table(file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_1/Out/output_case_1_train.out",
                        header = T, stringsAsFactors = F)

aa_GEM_10$Rows[1:10]
aa1_1 <- aa_GEM_10[seq(2, nrow(aa_GEM_10), 2),]
aa1_gt <- aa_GEM_10[seq(1, nrow(aa_GEM_10), 2),]
boxplot(aa1_1$X1 ~ aa1_gt$X1, xlab = "label", ylab = "predicted exp", main = "case 1: train")

aaroc<-roc.curve(scores.class0 = aa1_1$X1[aa1_gt$X1 == 1], scores.class1 = aa1_1$X1[aa1_gt$X1 == 0], curve = TRUE)
aapr<-pr.curve(scores.class0 = aa1_1$X1[aa1_gt$X1 == 1], scores.class1 = aa1_1$X1[aa1_gt$X1 == 0], curve = TRUE)
plot(aaroc, main = "case 1 train ROC")
plot(aapr, main = "case 1 train PR")

# read test res
aa_GEM_11 <- read.table(file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_1/Out/output_case_1_test.out",
                        header = T, stringsAsFactors = F)

aa_GEM_11$Rows[1:10]
aa1_1 <- aa_GEM_11[seq(2, nrow(aa_GEM_11), 2),]
aa1_gt <- aa_GEM_11[seq(1, nrow(aa_GEM_11), 2),]
boxplot(aa1_1$X1 ~ aa1_gt$X1, xlab = "label", 
        ylab = "predicted exp", main = "case 1: test")
# get PR and ROC curves

# aarfdown1 <- predict(my_learning_datasets_results_RF_kappa_down_1000bp_4[[i]], 
#                      aa_test_data1, type = "prob")

# aapreds_list <- list(aa1_1$X1)
# # List of actual values (same for all)
# m <- length(aapreds_list)
# aactuals_list <- rep(list(factor(aa1_gt$X1)), m)
# # Plot the ROC curves
# aapred <- prediction(aapreds_list, aactuals_list)
# aauc <- performance(aapred, measure = "auc")
# aarocs <- performance(aapred, "tpr", "fpr")
# print(aauc@y.values)
# plot(aarocs, col = as.list(1:m),
#      main = paste0("case 1 Test Set ROC"))
# abline(a = 0, b = 1, col = 6, lty = 3, lwd = 0.7)
# legend(x = "bottomright", cex = 0.8,
#        legend = c(paste0("GEMSTAT ", "AUC: " , format(aauc@y.values[[1]], digits = 2))),
#        fill = 1:m)

aaroc<-roc.curve(scores.class0 = aa1_1$X1[aa1_gt$X1 == 1], scores.class1 = aa1_1$X1[aa1_gt$X1 == 0], curve = TRUE)
aapr<-pr.curve(scores.class0 = aa1_1$X1[aa1_gt$X1 == 1], scores.class1 = aa1_1$X1[aa1_gt$X1 == 0], curve = TRUE)
plot(aaroc, main = "case 1 test ROC")
plot(aapr, main = "case 1 test PR")
################################################################################################################
################################################################################################################
############### CASE 2
# Choosing a larger test and training set
# test case 2: only increasing the number of training and test points: 100 training and 100 test. see what happens

# write GEMSTAT input with 40 training and 20 test sequences. each containing equal number of positives and negatives

aatpos <- sample(Positive_set_seq_list_char_1000bp[[4]], size = 100, replace = F)
aatneg <- sample(Negative_set_seq_list_char_1000bp[[4]], size = 100, replace = F)
names(aatpos) <- paste("pos", c(1:length(aatpos)), sep = "_")
names(aatneg) <- paste("neg", c(1:length(aatneg)), sep = "_")
GEMSTAT_ToyTrainSeq_2 <- sample(x = c(aatpos[1:50], aatneg[1:50]), size = 100, replace = F)
GEMSTAT_ToyTestSeq_2 <- c(aatpos[51:100], aatneg[51:100])

# Sequence
library(ShortRead)
writeFasta(DNAStringSet(GEMSTAT_ToyTrainSeq_2), width = 1000,
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_2/Sequence/GEMSTAT_ToyTrainSeq.fa")
writeFasta(DNAStringSet(GEMSTAT_ToyTestSeq_2), width = 1000, 
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_2/Sequence/GEMSTAT_ToyTestSeq.fa")

# Motif
aa_names <- sort(c("ESR1_2", "SP1", "PGR", "FOXA1", "GATA3", "JUN_2",  "MYC",   "PBX1", "PPARG", "RARA"))


# aa_coop_name[[1]] <- rbind(c("ESR1_2", "PGR"), 
#                            c("ESR1_2", "PPARG"),
#                            c("ESR1_2", "RARA"), 
#                            c("ESR1_2", "SP1"),
#                            c("ESR1_2", "FOXA1"))


MotifWriter(TF.motifs.Expanded_pseudo_count[aa_names],
            pseudo = 0.001, 
            output.File.Name = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_2/Motif/motifs")


# TF expression
aaTFexp <- rep(1, length(aa_names))
names(aaTFexp) <- aa_names
cat(c("Rows", paste0("1", "\n")), sep = "\t", append = F,
    file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_2/TF_exp/TF_exp.tab")
for(i in 1:length(aaTFexp)){
  cat(names(aaTFexp)[i], paste0(aaTFexp[i], rep("\n", as.integer(i != length(aaTFexp)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_2/TF_exp/TF_exp.tab",
      sep = "\t" )
}
# labels
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_ToyTrainSeq_2), split = "_")), "[[", 1))

aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_ToyTrainSeq_2)
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_2/Labels/Label_train.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_2/Labels/Label_train.txt", sep = "\t" )
}
# label for test
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_ToyTestSeq_2), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_ToyTestSeq_2)
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_2/Labels/Label_test_1.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_2/Labels/Label_test_1.txt", sep = "\t" )
}

#weights (in this case all equal)
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_ToyTrainSeq_2), split = "_")), "[[", 1))

aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_ToyTrainSeq_2)

aalab1[aalab1 == 0] <- 0.5
aalab1[aalab1 == 1] <- 0.5
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_2/weight_1_train.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_2/weight_1_train.txt", sep = "\t" )
}
# weights for test just in case of error
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_ToyTestSeq_2), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_ToyTestSeq_2)
aalab1[aalab1 == 0] <- 0.5
aalab1[aalab1 == 1] <- 0.5
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_2/weight_1_test.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_2/weight_1_test.txt", sep = "\t" )
}


# coop
aacoop <- rbind(c(1, 2), c(1, 3), c(1, 4), c(1,9), c(1, 10))
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_2")
par_ff_lb_ub_coop_Writer_multi_enh(.TFnames=aa_names,
                                   .annotation_thresh=numeric(0),
                                   .annotation_range=numeric(0),
                                   .initial_bind_w=numeric(0),
                                   .bind_w_range=numeric(0),
                                   .initial_alpha=numeric(0),
                                   .initial_log_reg_bias = numeric(0),
                                   .initial_log_reg_coeff = numeric(0),
                                   .initial_log_reg_bias_range = numeric(0),
                                   .initial_log_reg_coeff_range = numeric(0),
                                   .alpha_range=numeric(0),
                                   .coop_tf_mat=aacoop,
                                   .initial_coop_weight=numeric(0),
                                   .coop_weight_range=numeric(0),
                                   .coop_type=character(0),
                                   .coop_dist=integer(0),
                                   .coop_orientation=integer(0),
                                   one_qbtm_per_enh=F,
                                   .initial_qBTM=numeric(0),
                                   .qBTMrange=numeric(0),
                                   one_beta_per_enh=F,
                                   .initial_pi_beta=numeric(0), 
                                   .pi_beta_upper=numeric(0),
                                   .pi_beta_lower=numeric(0),
                                   .annotation_thresh_ff=integer(0),
                                   .initial_bind_w_ff=integer(0),
                                   .initial_alpha_ff=integer(0),
                                   .initial_coop_weight_ff=integer(0),
                                   .initial_qBTM_ff=integer(0),
                                   .initial_pi_beta_ff=integer(0),
                                   .filename_start=character(0),
                                   .filename_ff=character(0),
                                   .filename_upper=character(0),
                                   .filename_lower=character(0),
                                   .file_name_coop=character(0),
                                   .initial_log_reg_bias_ff = 1,
                                   .initial_log_reg_coeff_ff = 1,
                                   .nu_enhacners=100,
                                   ensemble_mode = F,
                                   .nu_samples = 2,
                                   .annotation_thresh_ff_ens = rep(0, 10),
                                   .initial_coop_weight_ff_ens = rep(0, 5),
                                   .initial_qBTM_ff_ens = 0,
                                   .initial_log_reg_bias_ff_ens = 0,
                                   .initial_log_reg_coeff_ff_ens = 0,
                                   Write_to_file = T,
                                   create_folders = T,
                                   TF_role_decided = F,
                                   factor_role_vector=integer(0),  
                                   .create_bounds=T,
                                   .create_par=T,
                                   .create_ff=T,
                                   .create_coop=T,
                                   logistic_params = T)

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")

# read training res
aa_GEM_12 <- read.table(file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_2/Out/output_case_2_train.out",
                        header = T, stringsAsFactors = F)

aa_GEM_12$Rows[1:10]
aa1_1 <- aa_GEM_12[seq(2, nrow(aa_GEM_12), 2),]
aa1_gt <- aa_GEM_12[seq(1, nrow(aa_GEM_12), 2),]
boxplot(aa1_1$X1 ~ aa1_gt$X1, xlab = "label", ylab = "predicted exp", main = "case 2: train")

aaroc<-roc.curve(scores.class0 = aa1_1$X1[aa1_gt$X1 == 1], scores.class1 = aa1_1$X1[aa1_gt$X1 == 0], curve = TRUE)
aapr<-pr.curve(scores.class0 = aa1_1$X1[aa1_gt$X1 == 1], scores.class1 = aa1_1$X1[aa1_gt$X1 == 0], curve = TRUE)
plot(aaroc, main = "case 2 train ROC")
plot(aapr, main = "case 2 train PR")

# read test results
aa_GEM_13 <- read.table(file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_Case_2/Out/output_case_2_test.out",
                        header = T, stringsAsFactors = F)

aa_GEM_13$Rows[1:10]
aa1_1 <- aa_GEM_13[seq(2, nrow(aa_GEM_13), 2),]
aa1_gt <- aa_GEM_13[seq(1, nrow(aa_GEM_13), 2),]
boxplot(aa1_1$X1 ~ aa1_gt$X1, xlab = "label", ylab = "predicted exp", main = "case 2: test")
# get PR and ROC curves
aaroc<-roc.curve(scores.class0 = aa1_1$X1[aa1_gt$X1 == 1], scores.class1 = aa1_1$X1[aa1_gt$X1 == 0], curve = TRUE)
aapr<-pr.curve(scores.class0 = aa1_1$X1[aa1_gt$X1 == 1], scores.class1 = aa1_1$X1[aa1_gt$X1 == 0], curve = TRUE)
plot(aaroc, main = "case 2 test ROC")
plot(aapr, main = "case 2 test PR")

###################################################################################################
# create a small ensemble just for test
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Test_ensemble")
aa_names <- sort(c("ESR1_2", "SP1", "PGR", "FOXA1", "GATA3", "JUN_2",  "MYC",   "PBX1", "PPARG", "RARA"))
aacoop <- rbind(c(1, 2), c(1, 3), c(1, 4), c(1,9), c(1, 10))

aamin_LLR <- numeric(length(TF.motifs.Expanded_LLR_to_pval_pseudo))
aamaxL_LLR <- numeric(length(TF.motifs.Expanded_LLR_to_pval_pseudo))
for(i in 1:length(TF.motifs.Expanded_LLR_to_pval_pseudo)){
  aamaxL_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][1,1] / 1000
  aa <- which(TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][,2] > 0.0001)[1]
  if(aa == 1){
    aamin_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][aa, 1] / 1000
  }else{
    aamin_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][(aa-1), 1] / 1000
  }
}
names(aamin_LLR) <- names(TF.motifs.Expanded_LLR_to_pval_pseudo)
names(aamaxL_LLR) <- names(TF.motifs.Expanded_LLR_to_pval_pseudo)
aathr <- 1- (aamin_LLR/aamaxL_LLR)

par_ff_lb_ub_coop_Writer_multi_enh(.TFnames=aa_names,
                                   .annotation_thresh=aathr[aa_names],
                                   .annotation_range=numeric(0),
                                   .initial_bind_w=numeric(0),
                                   .bind_w_range=numeric(0),
                                   .initial_alpha=numeric(0),
                                   .initial_log_reg_bias = numeric(0),
                                   .initial_log_reg_coeff = numeric(0),
                                   .initial_log_reg_bias_range = numeric(0),
                                   .initial_log_reg_coeff_range = numeric(0),
                                   .alpha_range=numeric(0),
                                   .coop_tf_mat=aacoop,
                                   .initial_coop_weight=numeric(0),
                                   .coop_weight_range=numeric(0),
                                   .coop_type=character(0),
                                   .coop_dist=integer(0),
                                   .coop_orientation=integer(0),
                                   one_qbtm_per_enh=F,
                                   .initial_qBTM=numeric(0),
                                   .qBTMrange=numeric(0),
                                   one_beta_per_enh=F,
                                   .initial_pi_beta=numeric(0), 
                                   .pi_beta_upper=numeric(0),
                                   .pi_beta_lower=numeric(0),
                                   .annotation_thresh_ff=integer(0),
                                   .initial_bind_w_ff=integer(0),
                                   .initial_alpha_ff=integer(0),
                                   .initial_coop_weight_ff=integer(0),
                                   .initial_qBTM_ff=1,
                                   .initial_pi_beta_ff=integer(0),
                                   .filename_start=character(0),
                                   .filename_ff=character(0),
                                   .filename_upper=character(0),
                                   .filename_lower=character(0),
                                   .file_name_coop=character(0),
                                   .initial_log_reg_bias_ff = 1,
                                   .initial_log_reg_coeff_ff = 1,
                                   .nu_enhacners=100,
                                   ensemble_mode = T,
                                   .nu_samples = 1,
                                   .annotation_thresh_ff_ens = rep(0, 10),
                                   .initial_coop_weight_ff_ens = c(0,1,0,0,0),
                                   .initial_qBTM_ff_ens = 1,
                                   .initial_log_reg_bias_ff_ens = 1,
                                   .initial_log_reg_coeff_ff_ens = 0,
                                   Write_to_file = T,
                                   create_folders = T,
                                   TF_role_decided = F,
                                   factor_role_vector=integer(0),  
                                   .create_bounds=T,
                                   .create_par=T,
                                   .create_ff=T,
                                   .create_coop=T,
                                   logistic_params = T,
                                   .initial_bind_w_ff_ens = c(1,0,0,0,0,0,0,0,0,0),
                                   .initial_alpha_ff_ens = c(0,0,0,0,0,1,0,0,0,0))
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")

###################################################################################################
# create first actual ensemble
# data_partition_mut_GEMSTAT[[1]] 
#data_partition_mut_GEMSTAT[[2]]
# data_partition_mut_GEMSTAT[[3]]
GEMSTAT_Ensemble_train_SeqList <- list()
GEMSTAT_Ensemble_test_SeqList <- list()
data_partition_mut_GEMSTAT_trimmed <- list()
data_partition_mut_GEMSTAT_trimmed[[1]] <- unlist(lapply(lapply(strsplit(data_partition_mut_GEMSTAT[[1]], split = "_"), "[", c(1, 2)), paste, collapse = "_"))
data_partition_mut_GEMSTAT_trimmed[[2]] <- unlist(lapply(lapply(strsplit(data_partition_mut_GEMSTAT[[2]], split = "_"), "[", c(1, 2)), paste, collapse = "_"))
data_partition_mut_GEMSTAT_trimmed[[3]] <- unlist(lapply(lapply(strsplit(data_partition_mut_GEMSTAT[[3]], split = "_"), "[", c(1, 2)), paste, collapse = "_"))
aalab <- unlist(lapply(strsplit(data_partition_mut_GEMSTAT_trimmed[[1]], split="_"), "[[", 1))
aa_train_ind_pos <- data_partition_mut_GEMSTAT_trimmed[[1]][sample(x = which(aalab %in% "pos"), size = 120, replace = F)]
aa_train_ind_neg <- data_partition_mut_GEMSTAT_trimmed[[1]][sample(x = which(aalab %in% "neg"), size = 120, replace = F)]

aamatchp <- match(aa_train_ind_pos, names(Positive_set_seq_list_char_1000bp[[4]]))
aamatchn <- match(aa_train_ind_neg, names(Negative_set_seq_list_char_1000bp[[4]]))

aatpos <- Positive_set_seq_list_char_1000bp[[4]][aamatchp]
aatneg <- Negative_set_seq_list_char_1000bp[[4]][aamatchn]
GEMSTAT_Ensemble_train_SeqList[[1]] <- sample(x = c(aatpos[1:120], aatneg[1:120]), size = 240, replace = F)

aalab <- unlist(lapply(strsplit(data_partition_mut_GEMSTAT_trimmed[[3]], split="_"), "[[", 1))
aa_test_ind_pos <- data_partition_mut_GEMSTAT_trimmed[[3]][sample(x = which(aalab %in% "pos"), size = 89, replace = F)]
aa_test_ind_neg <- data_partition_mut_GEMSTAT_trimmed[[3]][sample(x = which(aalab %in% "neg"), size = 333, replace = F)]

aamatchp <- match(aa_test_ind_pos, names(Positive_set_seq_list_char_1000bp[[4]]))
aamatchn <- match(aa_test_ind_neg, names(Negative_set_seq_list_char_1000bp[[4]]))

aatpos <- Positive_set_seq_list_char_1000bp[[4]][aamatchp]
aatneg <- Negative_set_seq_list_char_1000bp[[4]][aamatchn]
GEMSTAT_Ensemble_test_SeqList[[1]] <- sample(x = c(aatpos, aatneg), size = 422, replace = F)

#DNAStringSet
# Sequence
library(ShortRead)
writeFasta(DNAStringSet(GEMSTAT_Ensemble_train_SeqList[[1]]), width = 1000,
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_1/Sequence/train_seq.fa")
writeFasta(DNAStringSet(GEMSTAT_Ensemble_test_SeqList[[1]]), width = 1000, 
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_1/Sequence/test_seq.fa")

# Motif
aa_names <- c("ESR1_2","FOXA1","JUN_1","NFIB","NR2F2","NR3C1","NR5A2",
              "PBX1","POU5F1","PPARD","RARA","RUNX1","SP1","YBX1")

# coop set1
# "ESR1_2", "FOXA1"
# "ESR1_2", "NR2F2"
# "ESR1_2", "SP1"
# "NR2F2", "PPARD"
# "NR2F2", "YBX1"
# "POU5F1", "SP1"
# "PPARD", "YBX1"
# "RUNX1", "YBX1"

MotifWriter(TF.motifs.Expanded_pseudo_count[aa_names],
            pseudo = 0.001, 
            output.File.Name = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_1/motifs")


# TF expression
aaTFexp <- rep(1, length(aa_names))
names(aaTFexp) <- aa_names
cat(c("Rows", paste0("1", "\n")), sep = "\t", append = F,
    file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_1/TF_exp.tab")
for(i in 1:length(aaTFexp)){
  cat(names(aaTFexp)[i], paste0(aaTFexp[i], rep("\n", as.integer(i != length(aaTFexp)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_1/TF_exp.tab",
      sep = "\t" )
}
# labels
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_train_SeqList[[1]]), split = "_")), "[[", 1))

aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_train_SeqList[[1]])
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_1/Labels/Label_train.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_1/Labels/Label_train.txt", sep = "\t" )
}
# label for test
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_test_SeqList[[1]]), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_test_SeqList[[1]])
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_1/Labels/Label_test.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_1/Labels/Label_test.txt", sep = "\t" )
}

#weights (in this case all equal)
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_train_SeqList[[1]]), split = "_")), "[[", 1))

aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_train_SeqList[[1]])

aalab1[aalab1 == 0] <- 0.5
aalab1[aalab1 == 1] <- 0.5
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_1/weight_train.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_1/weight_train.txt", sep = "\t" )
}
# weights for test just in case of error
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_test_SeqList[[1]]), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_test_SeqList[[1]])
aalab1[aalab1 == 0] <- 0.21
aalab1[aalab1 == 1] <- 0.79
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_1/weight_test.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_1/weight_test.txt", sep = "\t" )
}


# coop
aa_names <- c("ESR1_2","FOXA1","JUN_1","NFIB","NR2F2","NR3C1","NR5A2",
              "PBX1","POU5F1","PPARD","RARA","RUNX1","SP1","YBX1")

# coop set1
# "ESR1_2", "FOXA1"
# "ESR1_2", "NR2F2"
# "ESR1_2", "SP1"
# "NR2F2", "PPARD"
# "NR2F2", "YBX1"
# "POU5F1", "SP1"
# "PPARD", "YBX1"
# "RUNX1", "YBX1"

aacoop <- rbind(c(1, 2), c(1, 5), c(1, 13), c(5,10),
                c(5, 14), c(9, 13), c(10, 14), c(12, 14))

aamin_LLR <- numeric(length(TF.motifs.Expanded_LLR_to_pval_pseudo))
aamaxL_LLR <- numeric(length(TF.motifs.Expanded_LLR_to_pval_pseudo))
for(i in 1:length(TF.motifs.Expanded_LLR_to_pval_pseudo)){
  aamaxL_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][1,1] / 1000
  aa <- which(TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][,2] > 0.0001)[1]
  if(aa == 1){
    aamin_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][aa, 1] / 1000
  }else{
    aamin_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][(aa-1), 1] / 1000
  }
}
names(aamin_LLR) <- names(TF.motifs.Expanded_LLR_to_pval_pseudo)
names(aamaxL_LLR) <- names(TF.motifs.Expanded_LLR_to_pval_pseudo)
aathr <- 1- (aamin_LLR/aamaxL_LLR)

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_1")
par_ff_lb_ub_coop_Writer_multi_enh(.TFnames=aa_names,
                                   .annotation_thresh=aathr[aa_names],
                                   .annotation_range=numeric(0),
                                   .initial_bind_w=numeric(0),
                                   .bind_w_range=numeric(0),
                                   .initial_alpha=numeric(0),
                                   .initial_log_reg_bias = numeric(0),
                                   .initial_log_reg_coeff = numeric(0),
                                   .initial_log_reg_bias_range = numeric(0),
                                   .initial_log_reg_coeff_range = numeric(0),
                                   .alpha_range=numeric(0),
                                   .coop_tf_mat=aacoop,
                                   .initial_coop_weight=numeric(0),
                                   .coop_weight_range=numeric(0),
                                   .coop_type=character(0),
                                   .coop_dist=integer(0),
                                   .coop_orientation=integer(0),
                                   one_qbtm_per_enh=F,
                                   .initial_qBTM=numeric(0),
                                   .qBTMrange=numeric(0),
                                   one_beta_per_enh=F,
                                   .initial_pi_beta=numeric(0), 
                                   .pi_beta_upper=numeric(0),
                                   .pi_beta_lower=numeric(0),
                                   .annotation_thresh_ff=integer(0),
                                   .initial_bind_w_ff=integer(0),
                                   .initial_alpha_ff=integer(0),
                                   .initial_coop_weight_ff=integer(0),
                                   .initial_qBTM_ff=1,
                                   .initial_pi_beta_ff=integer(0),
                                   .filename_start=character(0),
                                   .filename_ff=character(0),
                                   .filename_upper=character(0),
                                   .filename_lower=character(0),
                                   .file_name_coop=character(0),
                                   .initial_log_reg_bias_ff = 1,
                                   .initial_log_reg_coeff_ff = 1,
                                   .nu_enhacners=240,
                                   ensemble_mode = T,
                                   .nu_samples = 1,
                                   .annotation_thresh_ff_ens = rep(0, 14),
                                   .initial_coop_weight_ff_ens = c(0,1,0,0,0,1,0,0),
                                   .initial_qBTM_ff_ens = 1,
                                   .initial_log_reg_bias_ff_ens = 1,
                                   .initial_log_reg_coeff_ff_ens = 0,
                                   Write_to_file = T,
                                   create_folders = T,
                                   TF_role_decided = F,
                                   factor_role_vector=integer(0),  
                                   .create_bounds=T,
                                   .create_par=T,
                                   .create_ff=T,
                                   .create_coop=T,
                                   logistic_params = T,
                                   .initial_bind_w_ff_ens = c(0,1,0,0,0,0,0,0,0,1),
                                   .initial_alpha_ff_ens = c(0,0,1,0,0,0,0,0,1,0))


# write GEMSTAT jobs 
cat("#!/bin/bash\n", file = "Ensemble_Exp1_job.sh", append = F)
aa_par_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_1/Parameters")
aa_par_files_trimmed <- unlist(lapply(strsplit(aa_par_files, "\\."), "[[", 1))
for(i in 1:length(aa_par_files)){
  cat(c("~/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/src/seq2expr -s Sequence/train_seq.fa",
        "-e Labels/Label_train.txt",
        "-m motifs.wtmx -f TF_exp.tab", 
        paste0("-fo Out/output_", aa_par_files_trimmed[i], ".txt"), 
        "-o DIRECT -c Coop/coop.par ",
        paste0("-p Parameters/", aa_par_files[i]),
        "-lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -na 2 -onebeta 1",
        "-train_weights weight_train.txt > ", 
        paste0("Log/log_", aa_par_files_trimmed[i], ".txt\n")),
      sep = " ", append = T, file = "Ensemble_Exp1_job.sh")
}

#write test script job for hal
#cat("#!/bin/bash\n", file = "Ensemble_Exp1_job_test.sh", append = F)
aa_par_trained_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_1/Trained_par")
aa_par_trained_files_trimmed <- unlist(lapply(strsplit(aa_par_trained_files, "\\."), "[[", 1))
for(i in 1:length(aa_par_trained_files)){
  cat(c("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr -s Sequence/test_seq.fa",
        "-e Labels/Label_test.txt",
        "-m motifs.wtmx -f TF_exp.tab", 
        paste0("-fo Out/output_", aa_par_trained_files_trimmed[i], "_test.txt"), 
        "-o DIRECT -c Coop/coop.par ",
        paste0("-p Trained_par/", aa_par_trained_files[i]),
        "-lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -na 0 -onebeta 1",
        "-train_weights weight_test.txt > ", 
        paste0("Log/log_", aa_par_trained_files_trimmed[i], "_test.txt\n")),
      sep = " ", append = T, file = "Ensemble_Exp1_job_test.sh")
}

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")
########################################################################################################################

E_RNA_GEMSTAT_Ensemble_Outlist <- list()
E_RNA_GEMSTAT_Ensemble_Parlist <- list()
# read train-test results
E_RNA_GEMSTAT_Ensemble_Outlist[[1]] <- read_output_train_test_GEMSTAT_ensemble(output_dir = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_1/Out",
                                                                               should_plot = F)
###### read parameters
E_RNA_GEMSTAT_Ensemble_Parlist[[1]] <- read_parameters_GEMSTAT_ensemble(par_directory = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_1/Trained_par/",
                                                                        convert_to_mat = T)

library(gplots)
aa_perf_srt <-names(E_RNA_GEMSTAT_Ensemble_Outlist[[1]]$test_ROC)[sort(E_RNA_GEMSTAT_Ensemble_Outlist[[1]]$test_ROC, decreasing = T, index.return=T)$ix]
aa_perf_srt2 <- unlist(lapply(lapply(strsplit(aa_perf_srt, split="_"), "[", c(3,4,5)), paste, collapse="_"))
aa_perf_srt22 <- paste0("log_", aa_perf_srt2)
heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[1]]$binding[aa_perf_srt22, ]), Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none')
heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[1]]$alpha[aa_perf_srt22, ] + 1e-6), Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none')
heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[1]]$coop[aa_perf_srt22, ]), Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none', margins = c(10,10))
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[1]]$qbtm[aa_perf_srt22],ylab = "" ,main = "qBTM sorted by perf exp1")
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[1]]$beta[aa_perf_srt22],ylab = "" ,main = "beta sorted by perf exp1")

heatmap.2(E_RNA_GEMSTAT_Ensemble_Parlist[[1]]$logistic[aa_perf_srt22,], Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none', margins = c(10,10))
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[1]]$logistic[aa_perf_srt22, 1], ylab = "", main="logistic bias sorted by perf")
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[1]]$logistic[aa_perf_srt22, 2], ylab = "", main="logistic coeff sorted by perf")

########################################################################################################################
########################################################################################################################
# Experiemnt 2
# set2



#DNAStringSet
# Sequence
library(ShortRead)
writeFasta(DNAStringSet(GEMSTAT_Ensemble_train_SeqList[[1]]), width = 1000,
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_2/Sequence/train_seq.fa")
writeFasta(DNAStringSet(GEMSTAT_Ensemble_test_SeqList[[1]]), width = 1000, 
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_2/Sequence/test_seq.fa")

# Motif
aa_names <- c("AR","ESR1_1","FOXA1","JUN_1","LEF1","NFIB","NKX3_1",
              "NR2C1","NR2F2","PGR", "PPARD","PPARG","RARA","SP1","YBX1")

"ESR1_1_vs_ESR1_1"
"ESR_1_vs_FOXA1"
"ESR1_1_vs_NR2F2"
"ESR_1_vs_PGR"
"ESR_1_vs_PPARG"
"ESR_1_vs_SP1"
"ESR_1_vs_YBX1"
"NR2F2_vs_NR2F2"
names(TF.motifs.Expanded_pseudo_count)[names(TF.motifs.Expanded_pseudo_count) == "NKX3-1"] <- "NKX3_1"

MotifWriter(TF.motifs.Expanded_pseudo_count[aa_names],
            pseudo = 0.001, 
            output.File.Name = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_2/motifs")


# TF expression
aaTFexp <- rep(1, length(aa_names))
names(aaTFexp) <- aa_names
cat(c("Rows", paste0("1", "\n")), sep = "\t", append = F,
    file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_2/TF_exp.tab")
for(i in 1:length(aaTFexp)){
  cat(names(aaTFexp)[i], paste0(aaTFexp[i], rep("\n", as.integer(i != length(aaTFexp)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_2/TF_exp.tab",
      sep = "\t" )
}
# labels
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_train_SeqList[[1]]), split = "_")), "[[", 1))

aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_train_SeqList[[1]])
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_2/Labels/Label_train.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_2/Labels/Label_train.txt", sep = "\t" )
}
# label for test
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_test_SeqList[[1]]), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_test_SeqList[[1]])
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_2/Labels/Label_test.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_2/Labels/Label_test.txt", sep = "\t" )
}

#weights (in this case all equal)
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_train_SeqList[[1]]), split = "_")), "[[", 1))

aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_train_SeqList[[1]])

aalab1[aalab1 == 0] <- 0.5
aalab1[aalab1 == 1] <- 0.5
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_2/weight_train.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_2/weight_train.txt", sep = "\t" )
}
# weights for test just in case of error
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_test_SeqList[[1]]), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_test_SeqList[[1]])
aalab1[aalab1 == 0] <- 0.21
aalab1[aalab1 == 1] <- 0.79
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_2/weight_test.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_2/weight_test.txt", sep = "\t" )
}


# coop
aa_names <- c("AR","ESR1_1","FOXA1","JUN_1","LEF1","NFIB","NKX3_1",
              "NR2C1","NR2F2","PGR", "PPARD","PPARG","RARA","SP1","YBX1")

"ESR1_1_vs_ESR1_1"
"ESR_1_vs_FOXA1"
"ESR1_1_vs_NR2F2"
"ESR_1_vs_PGR"
"ESR_1_vs_PPARG"
"ESR_1_vs_SP1"
"ESR_1_vs_YBX1"
"NR2F2_vs_NR2F2"
aa_coop <- rbind(c(2,2), c(2,3) ,c(2, 9), c(2, 10),
                 c(2, 12),c(2, 14),c(2, 15), c(9, 9))

aamin_LLR <- numeric(length(TF.motifs.Expanded_LLR_to_pval_pseudo))
aamaxL_LLR <- numeric(length(TF.motifs.Expanded_LLR_to_pval_pseudo))
for(i in 1:length(TF.motifs.Expanded_LLR_to_pval_pseudo)){
  aamaxL_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][1,1] / 1000
  aa <- which(TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][,2] > 0.0001)[1]
  if(aa == 1){
    aamin_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][aa, 1] / 1000
  }else{
    aamin_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][(aa-1), 1] / 1000
  }
}
names(aamin_LLR) <- names(TF.motifs.Expanded_LLR_to_pval_pseudo)
names(aamaxL_LLR) <- names(TF.motifs.Expanded_LLR_to_pval_pseudo)
aathr <- 1- (aamin_LLR/aamaxL_LLR)
names(aathr)[names(aathr) == "NKX3-1"] <- "NKX3_1"

aa_bind_ff_ens <- c(0,0,1,0,0,0,0,0,1,0,0,1,0,0,1)
aa_alpha_ff_ens<- c(0,0,0,1,0,0,0,0,0,1,0,0,0,0,0)
names(aa_bind_ff_ens) <- aa_names
names(aa_alpha_ff_ens) <- aa_names

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_2")
par_ff_lb_ub_coop_Writer_multi_enh(.TFnames=aa_names,
                                   .annotation_thresh=aathr[aa_names],
                                   .annotation_range=numeric(0),
                                   .initial_bind_w=numeric(0),
                                   .bind_w_range=numeric(0),
                                   .initial_alpha=numeric(0),
                                   .initial_log_reg_bias = numeric(0),
                                   .initial_log_reg_coeff = numeric(0),
                                   .initial_log_reg_bias_range = numeric(0),
                                   .initial_log_reg_coeff_range = numeric(0),
                                   .alpha_range=numeric(0),
                                   .coop_tf_mat=aa_coop,
                                   .initial_coop_weight=numeric(0),
                                   .coop_weight_range=numeric(0),
                                   .coop_type=character(0),
                                   .coop_dist=integer(0),
                                   .coop_orientation=integer(0),
                                   one_qbtm_per_enh=F,
                                   .initial_qBTM=numeric(0),
                                   .qBTMrange=numeric(0),
                                   one_beta_per_enh=F,
                                   .initial_pi_beta=numeric(0), 
                                   .pi_beta_upper=numeric(0),
                                   .pi_beta_lower=numeric(0),
                                   .annotation_thresh_ff=integer(0),
                                   .initial_bind_w_ff=integer(0),
                                   .initial_alpha_ff=integer(0),
                                   .initial_coop_weight_ff=integer(0),
                                   .initial_qBTM_ff=1,
                                   .initial_pi_beta_ff=integer(0),
                                   .filename_start=character(0),
                                   .filename_ff=character(0),
                                   .filename_upper=character(0),
                                   .filename_lower=character(0),
                                   .file_name_coop=character(0),
                                   .initial_log_reg_bias_ff = 1,
                                   .initial_log_reg_coeff_ff = 1,
                                   .nu_enhacners=240,
                                   ensemble_mode = T,
                                   .nu_samples = 1,
                                   .annotation_thresh_ff_ens = rep(0, 14),
                                   .initial_coop_weight_ff_ens = c(0,0,1,0,0,1,0,0),
                                   .initial_qBTM_ff_ens = 1,
                                   .initial_log_reg_bias_ff_ens = 1,
                                   .initial_log_reg_coeff_ff_ens = 0,
                                   Write_to_file = T,
                                   create_folders = T,
                                   TF_role_decided = F,
                                   factor_role_vector=integer(0),  
                                   .create_bounds=T,
                                   .create_par=T,
                                   .create_ff=T,
                                   .create_coop=T,
                                   logistic_params = T,
                                   .initial_bind_w_ff_ens = aa_bind_ff_ens,
                                   .initial_alpha_ff_ens = aa_alpha_ff_ens)


# write GEMSTAT jobs
cat("#!/bin/bash\n", file = "Ensemble_Exp2_job.sh", append = F)
aa_par_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_2/Parameters")
aa_par_files_trimmed <- unlist(lapply(strsplit(aa_par_files, "\\."), "[[", 1))
for(i in 1:length(aa_par_files)){
  cat(c("~/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/src/seq2expr -s Sequence/train_seq.fa",
        "-e Labels/Label_train.txt",
        "-m motifs.wtmx -f TF_exp.tab", 
        paste0("-fo Out/output_", aa_par_files_trimmed[i], ".txt"), 
        "-o DIRECT -c Coop/coop.par ",
        paste0("-p Parameters/", aa_par_files[i]),
        "-lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -na 2 -onebeta 1",
        "-train_weights weight_train.txt > ", 
        paste0("Log/log_", aa_par_files_trimmed[i], ".txt\n")),
      sep = " ", append = T, file = "Ensemble_Exp2_job.sh")
}

#write test script job for hal
#cat("#!/bin/bash\n", file = "Ensemble_Exp2_job_test.sh", append = F)
aa_par_trained_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_2/Trained_par")
aa_par_trained_files_trimmed <- unlist(lapply(strsplit(aa_par_trained_files, "\\."), "[[", 1))
for(i in 1:length(aa_par_trained_files)){
  cat(c("~/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/src/seq2expr -s Sequence/test_seq.fa",
        "-e Labels/Label_test.txt",
        "-m motifs.wtmx -f TF_exp.tab", 
        paste0("-fo Out/output_", aa_par_trained_files_trimmed[i], "_test.txt"), 
        "-o DIRECT -c Coop/coop.par ",
        paste0("-p Trained_par/", aa_par_trained_files[i]),
        "-lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -na 0 -onebeta 1",
        "-train_weights weight_test.txt > ", 
        paste0("Log/log_", aa_par_trained_files_trimmed[i], "_test.txt\n")),
      sep = " ", append = T, file = "Ensemble_Exp2_job_test.sh")
}

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")
# read train-test results

# read train-test results
E_RNA_GEMSTAT_Ensemble_Outlist[[2]] <- read_output_train_test_GEMSTAT_ensemble(output_dir = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_2/Out",
                                                                               should_plot = F)
###### read parameters
E_RNA_GEMSTAT_Ensemble_Parlist[[2]] <- read_parameters_GEMSTAT_ensemble(par_directory = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_2/Trained_par/",
                                                                        convert_to_mat = T)

library(gplots)
sort(E_RNA_GEMSTAT_Ensemble_Outlist[[2]]$test_ROC, decreasing = T)[1:24]
aa_perf_srt <-names(E_RNA_GEMSTAT_Ensemble_Outlist[[2]]$test_ROC)[sort(E_RNA_GEMSTAT_Ensemble_Outlist[[2]]$test_ROC, decreasing = T, index.return=T)$ix]
aa_perf_srt2 <- unlist(lapply(lapply(strsplit(aa_perf_srt, split="_"), "[", c(3,4,5)), paste, collapse="_"))
aa_perf_srt22 <- paste0("log_", aa_perf_srt2)
heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[2]]$binding[aa_perf_srt22[c(c(1:24), c(101:124))], ]), Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none')
heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[2]]$alpha[aa_perf_srt22, ] + 1e-6), Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none')
heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[2]]$coop[aa_perf_srt22, ]), Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none', margins = c(10,10))
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[2]]$qbtm[aa_perf_srt22], main = "qbtm_sorted_by_perf_exp2", pch = 16)
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[2]]$beta[aa_perf_srt22], main = "beta_sorted_by_perf_exp2", pch = 16)

heatmap.2(E_RNA_GEMSTAT_Ensemble_Parlist[[2]]$logistic[aa_perf_srt22,], Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none', margins = c(10,10))
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[2]]$logistic[aa_perf_srt22, 1], ylab="", main = "logistic_bias_sorted_by_perf")
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[2]]$logistic[aa_perf_srt22, 2], ylab="", main = "logistic_coeff_sorted_by_perf")

########################################################################################################################
########################################################################################################################
# Experiemnt 3
# set3
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/")
system("mkdir Experiment_3")
system("mkdir Experiment_3/Sequence")
system("mkdir Experiment_3/Labels")

# Sequence
library(ShortRead)
writeFasta(DNAStringSet(GEMSTAT_Ensemble_train_SeqList[[1]]), width = 1000,
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_3/Sequence/train_seq.fa")
writeFasta(DNAStringSet(GEMSTAT_Ensemble_test_SeqList[[1]]), width = 1000, 
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_3/Sequence/test_seq.fa")

# Motif
# coop
aa_names <- c("ESR1_2","FOXA1","JUN_1","NFIB","NKX3_1","NR2F2","NR3C1","NR5A2",
              "PBX1","PGR","PPARD","PPARG","RARA","RUNX1","SP1","YBX1")
# coop set1
# "ESR1_2", "ESR1_2"
# "ESR1_2", "FOXA1"
# "ESR1_2", "NR2F2"
# "ESR1_2", "PGR"
# "ESR1_2", "SP1"
# "ESR1_2", "YBX1"
# "NR2F2", "YBX1"
aacoop <- rbind(c(1, 1), c(1, 2), c(1, 6),
                c(1, 10), c(1,15),
                c(1, 16), c(6, 16))
MotifWriter(TF.motifs.Expanded_pseudo_count[aa_names],
            pseudo = 0.001, 
            output.File.Name = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_3/motifs")


# TF expression
aaTFexp <- rep(1, length(aa_names))
names(aaTFexp) <- aa_names
cat(c("Rows", paste0("1", "\n")), sep = "\t", append = F,
    file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_3/TF_exp.tab")
for(i in 1:length(aaTFexp)){
  cat(names(aaTFexp)[i], paste0(aaTFexp[i], rep("\n", as.integer(i != length(aaTFexp)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_3/TF_exp.tab",
      sep = "\t" )
}
# labels
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_train_SeqList[[1]]), split = "_")), "[[", 1))

aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_train_SeqList[[1]])
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_3/Labels/Label_train.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_3/Labels/Label_train.txt", sep = "\t" )
}
# label for test
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_test_SeqList[[1]]), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_test_SeqList[[1]])
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_3/Labels/Label_test.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_3/Labels/Label_test.txt", sep = "\t" )
}

#weights (in this case all equal)
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_train_SeqList[[1]]), split = "_")), "[[", 1))

aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_train_SeqList[[1]])

aalab1[aalab1 == 0] <- 0.5
aalab1[aalab1 == 1] <- 0.5
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_3/weight_train.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_3/weight_train.txt", sep = "\t" )
}
# weights for test just in case of error
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_test_SeqList[[1]]), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_test_SeqList[[1]])
aalab1[aalab1 == 0] <- 0.21
aalab1[aalab1 == 1] <- 0.79
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_3/weight_test.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_3/weight_test.txt", sep = "\t" )
}


# coop
# coop set1
# "ESR1_2", "ESR1_2"
# "ESR1_2", "FOXA1"
# "ESR1_2", "NR2F2"
# "ESR1_2", "PGR"
# "ESR1_2", "SP1"
# "ESR1_2", "YBX1"
# "NR2F2", "YBX1"
aacoop <- rbind(c(1, 1), c(1, 2), c(1, 6),
                c(1, 10), c(1,15),
                c(1, 16), c(6, 16))

aamin_LLR <- numeric(length(TF.motifs.Expanded_LLR_to_pval_pseudo))
aamaxL_LLR <- numeric(length(TF.motifs.Expanded_LLR_to_pval_pseudo))
for(i in 1:length(TF.motifs.Expanded_LLR_to_pval_pseudo)){
  aamaxL_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][1,1] / 1000
  aa <- which(TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][,2] > 0.0001)[1]
  if(aa == 1){
    aamin_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][aa, 1] / 1000
  }else{
    aamin_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][(aa-1), 1] / 1000
  }
}
names(aamin_LLR) <- names(TF.motifs.Expanded_LLR_to_pval_pseudo)
names(aamaxL_LLR) <- names(TF.motifs.Expanded_LLR_to_pval_pseudo)
aathr <- 1- (aamin_LLR/aamaxL_LLR)
names(aathr)[names(aathr) == "NKX3-1"] <- "NKX3_1"

aa_bind_ff_ens <- c(0,0,0,0,0,1,1,0,1,1,0,1,0,0,0,1)
aa_alpha_ff_ens<- c(0,0,0,0,0,1,1,0,1,1,0,1,0,0,0,1)
names(aa_bind_ff_ens) <- aa_names
names(aa_alpha_ff_ens) <- aa_names

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_3")
par_ff_lb_ub_coop_Writer_multi_enh(.TFnames=aa_names,
                                   .annotation_thresh=aathr[aa_names],
                                   .annotation_range=numeric(0),
                                   .initial_bind_w=numeric(0),
                                   .bind_w_range=numeric(0),
                                   .initial_alpha=numeric(0),
                                   .initial_log_reg_bias = numeric(0),
                                   .initial_log_reg_coeff = 5,
                                   .initial_log_reg_bias_range = numeric(0),
                                   .initial_log_reg_coeff_range = c(1, 20),
                                   .alpha_range=numeric(0),
                                   .coop_tf_mat=aacoop,
                                   .initial_coop_weight=numeric(0),
                                   .coop_weight_range=numeric(0),
                                   .coop_type=character(0),
                                   .coop_dist=integer(0),
                                   .coop_orientation=integer(0),
                                   one_qbtm_per_enh=F,
                                   .initial_qBTM=numeric(0),
                                   .qBTMrange=numeric(0),
                                   one_beta_per_enh=F,
                                   .initial_pi_beta=numeric(0), 
                                   .pi_beta_upper=numeric(0),
                                   .pi_beta_lower=numeric(0),
                                   .annotation_thresh_ff=integer(0),
                                   .initial_bind_w_ff=integer(0),
                                   .initial_alpha_ff=integer(0),
                                   .initial_coop_weight_ff=integer(0),
                                   .initial_qBTM_ff=1,
                                   .initial_pi_beta_ff=integer(0),
                                   .filename_start=character(0),
                                   .filename_ff=character(0),
                                   .filename_upper=character(0),
                                   .filename_lower=character(0),
                                   .file_name_coop=character(0),
                                   .initial_log_reg_bias_ff = 1,
                                   .initial_log_reg_coeff_ff = 1,
                                   .nu_enhacners=240,
                                   ensemble_mode = T,
                                   .nu_samples = 1,
                                   .annotation_thresh_ff_ens = rep(0, 14),
                                   .initial_coop_weight_ff_ens = c(0,0,0,1,0,0,0),
                                   .initial_qBTM_ff_ens = 0,
                                   .initial_log_reg_bias_ff_ens = 0,
                                   .initial_log_reg_coeff_ff_ens = 0,
                                   Write_to_file = T,
                                   create_folders = T,
                                   TF_role_decided = F,
                                   factor_role_vector=integer(0),  
                                   .create_bounds=T,
                                   .create_par=T,
                                   .create_ff=T,
                                   .create_coop=T,
                                   logistic_params = T,
                                   .initial_bind_w_ff_ens = aa_bind_ff_ens,
                                   .initial_alpha_ff_ens = aa_alpha_ff_ens)


# write GEMSTAT jobs
#cat("#!/bin/bash\n", file = "Ensemble_Exp3_job_train_hal.job", append = F)
aa_par_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_3/Parameters")
aa_par_files_trimmed <- unlist(lapply(strsplit(aa_par_files, "\\."), "[[", 1))
for(i in 1:length(aa_par_files)){
  cat(c("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr -s Sequence/train_seq.fa",
        "-e Labels/Label_train.txt",
        "-m motifs.wtmx -f TF_exp.tab", 
        paste0("-fo Out/output_", aa_par_files_trimmed[i], ".txt"), 
        "-o DIRECT -c Coop/coop.par ",
        paste0("-p Parameters/", aa_par_files[i]),
        "-lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -na 2 -onebeta 1",
        "-train_weights weight_train.txt > ", 
        paste0("Log/log_", aa_par_files_trimmed[i], ".txt\n")),
      sep = " ", append = T, file = "Ensemble_Exp3_job_train_hal.job")
}

#write test script job for hal
#cat("#!/bin/bash\n", file = "Ensemble_Exp2_job_test.sh", append = F)
aa_par_trained_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_3/Trained_par")
aa_par_trained_files_trimmed <- unlist(lapply(strsplit(aa_par_trained_files, "\\."), "[[", 1))
for(i in 1:length(aa_par_trained_files)){
  cat(c("~/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/src/seq2expr -s Sequence/test_seq.fa",
        "-e Labels/Label_test.txt",
        "-m motifs.wtmx -f TF_exp.tab", 
        paste0("-fo Out/output_", aa_par_trained_files_trimmed[i], "_test.txt"), 
        "-o DIRECT -c Coop/coop.par ",
        paste0("-p Trained_par/", aa_par_trained_files[i]),
        "-lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -na 0 -onebeta 1",
        "-train_weights weight_test.txt > ", 
        paste0("Log/log_", aa_par_trained_files_trimmed[i], "_test.txt\n")),
      sep = " ", append = T, file = "Ensemble_Exp2_job_test.sh")
}

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")
# read train-test results

# read train-test results
E_RNA_GEMSTAT_Ensemble_Outlist[[3]] <- read_output_train_test_GEMSTAT_ensemble(output_dir = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_3/Out",
                                                                               should_plot = F)
###### read parameters
E_RNA_GEMSTAT_Ensemble_Parlist[[3]] <- read_parameters_GEMSTAT_ensemble(par_directory = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_3/Trained_par/",
                                                                        convert_to_mat = T)
c(c(1:400), c(7001:7400))
library(gplots)
sort(E_RNA_GEMSTAT_Ensemble_Outlist[[3]]$test_ROC, decreasing = T)
aa_perf_srt <-names(E_RNA_GEMSTAT_Ensemble_Outlist[[3]]$test_ROC)[sort(E_RNA_GEMSTAT_Ensemble_Outlist[[3]]$test_ROC, decreasing = T, index.return=T)$ix]
aa_perf_srt2 <- unlist(lapply(lapply(strsplit(aa_perf_srt, split="_"), "[", c(3,4,5)), paste, collapse="_"))
aa_perf_srt22 <- paste0("log_", aa_perf_srt2)
heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[3]]$binding[aa_perf_srt22[c(c(1:400), c(7001:7400))], ]), Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none')
heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[3]]$alpha[aa_perf_srt22[c(c(1:400), c(7001:7400))], ] + 1e-6), Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none')
heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[3]]$coop[aa_perf_srt22[c(c(1:400), c(7001:7400))], ]), Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none', margins = c(10,10))

plot(E_RNA_GEMSTAT_Ensemble_Parlist[[3]]$qbtm[aa_perf_srt22], main = "qbtm_sorted_by_perf_exp3", pch = 16)
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[3]]$beta[aa_perf_srt22[-4133]], main = "beta_sorted_by_perf_exp3", pch = 16)

heatmap.2(E_RNA_GEMSTAT_Ensemble_Parlist[[3]]$logistic[aa_perf_srt22,], Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none', margins = c(10,10))
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[3]]$logistic[aa_perf_srt22, 1], ylab="", main = "logistic_bias_sorted_by_perf_exp3")
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[3]]$logistic[aa_perf_srt22, 2], ylab="", main = "logistic_coeff_sorted_by_perf_exp3")
# plot the best models
which.max(E_RNA_GEMSTAT_Ensemble_Outlist[[3]]$test_ROC)
which.max(E_RNA_GEMSTAT_Ensemble_Outlist[[3]]$test_PRC)
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[3]]$train_results[[4818]]$ROC_curve, main = "ROC train bestModel exp3")
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[3]]$train_results[[4818]]$PRC_curve, main = "PRC train bestModel exp3")
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[3]]$test_results[[4818]]$ROC_curve, main = "ROC test bestModel exp3")
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[3]]$test_results[[4818]]$PRC_curve, main = "PRC test bestModel exp3")

########################################################################################################################
########################################################################################################################
# Experiemnt 4
# set3
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/")
system("mkdir Experiment_4")
system("mkdir Experiment_4/Sequence")
system("mkdir Experiment_4/Labels")

data_partition_mut_GEMSTAT_trimmed <- list()
data_partition_mut_GEMSTAT_trimmed[[1]] <- unlist(lapply(lapply(strsplit(data_partition_mut_GEMSTAT[[1]], split = "_"), "[", c(1, 2)), paste, collapse = "_"))
data_partition_mut_GEMSTAT_trimmed[[2]] <- unlist(lapply(lapply(strsplit(data_partition_mut_GEMSTAT[[2]], split = "_"), "[", c(1, 2)), paste, collapse = "_"))
data_partition_mut_GEMSTAT_trimmed[[3]] <- unlist(lapply(lapply(strsplit(data_partition_mut_GEMSTAT[[3]], split = "_"), "[", c(1, 2)), paste, collapse = "_"))
aalab <- unlist(lapply(strsplit(data_partition_mut_GEMSTAT_trimmed[[1]], split="_"), "[[", 1))
aa_train_ind_pos <- data_partition_mut_GEMSTAT_trimmed[[1]][sample(x = which(aalab %in% "pos"), size = 200, replace = F)]
aa_train_ind_neg <- data_partition_mut_GEMSTAT_trimmed[[1]][sample(x = which(aalab %in% "neg"), size = 600, replace = F)]

aamatchp <- match(aa_train_ind_pos, names(Positive_set_seq_list_char_1000bp[[4]]))
aamatchn <- match(aa_train_ind_neg, names(Negative_set_seq_list_char_1000bp[[4]]))

aatpos <- Positive_set_seq_list_char_1000bp[[4]][aamatchp]
aatneg <- Negative_set_seq_list_char_1000bp[[4]][aamatchn]
GEMSTAT_Ensemble_train_SeqList[[2]] <- sample(x = c(aatpos, aatneg), size = 800, replace = F)

aalab <- unlist(lapply(strsplit(data_partition_mut_GEMSTAT_trimmed[[3]], split="_"), "[[", 1))
aa_test_ind_pos <- data_partition_mut_GEMSTAT_trimmed[[3]][sample(x = which(aalab %in% "pos"), size = 89, replace = F)]
aa_test_ind_neg <- data_partition_mut_GEMSTAT_trimmed[[3]][sample(x = which(aalab %in% "neg"), size = 333, replace = F)]

aamatchp <- match(aa_test_ind_pos, names(Positive_set_seq_list_char_1000bp[[4]]))
aamatchn <- match(aa_test_ind_neg, names(Negative_set_seq_list_char_1000bp[[4]]))

aatpos <- Positive_set_seq_list_char_1000bp[[4]][aamatchp]
aatneg <- Negative_set_seq_list_char_1000bp[[4]][aamatchn]
GEMSTAT_Ensemble_test_SeqList[[2]] <- sample(x = c(aatpos, aatneg), size = 422, replace = F)



# Sequence
library(ShortRead)

writeFasta(DNAStringSet(GEMSTAT_Ensemble_train_SeqList[[2]]), width = 1000,
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_4/Sequence/train_seq.fa")
writeFasta(DNAStringSet(GEMSTAT_Ensemble_test_SeqList[[2]]), width = 1000, 
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_4/Sequence/test_seq.fa")
writeFasta(DNAStringSet(GEMSTAT_Ensemble_validation_SeqList[[1]]), width = 1000, 
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_4/Sequence/valid_seq.fa")

# Motif
# coop
aa_names <- c("ESR1_2","FOXA1","JUN_1","NFIB","NKX3_1","NR2F2","NR3C1","NR5A2",
              "PBX1","PGR","PPARD","PPARG","RARA","RUNX1","SP1","YBX1")
# coop set1
# "ESR1_2", "ESR1_2"
# "ESR1_2", "FOXA1"
# "ESR1_2", "NR2F2"
# "ESR1_2", "PGR"
# "ESR1_2", "SP1"
# "ESR1_2", "YBX1"
# "NR2F2", "YBX1"
aacoop <- rbind(c(1, 1), c(1, 2), c(1, 6),
                c(1, 10), c(1,15),
                c(1, 16), c(6, 16))
MotifWriter(TF.motifs.Expanded_pseudo_count[aa_names],
            pseudo = 0.001, 
            output.File.Name = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_4/motifs")


# TF expression
aaTFexp <- rep(1, length(aa_names))
names(aaTFexp) <- aa_names
cat(c("Rows", paste0("1", "\n")), sep = "\t", append = F,
    file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_4/TF_exp.tab")
for(i in 1:length(aaTFexp)){
  cat(names(aaTFexp)[i], paste0(aaTFexp[i], rep("\n", as.integer(i != length(aaTFexp)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_4/TF_exp.tab",
      sep = "\t" )
}
# labels
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_train_SeqList[[2]]), split = "_")), "[[", 1))

aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_train_SeqList[[2]])
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_4/Labels/Label_train.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_4/Labels/Label_train.txt", sep = "\t" )
}
# label for test
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_test_SeqList[[2]]), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_test_SeqList[[2]])
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_4/Labels/Label_test.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_4/Labels/Label_test.txt", sep = "\t" )
}

#weights (in this case not equal)
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_train_SeqList[[2]]), split = "_")), "[[", 1))

aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_train_SeqList[[2]])

aalab1[aalab1 == 0] <- 0.25
aalab1[aalab1 == 1] <- 0.75
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_4/weight_train.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_4/weight_train.txt", sep = "\t" )
}
# weights for test just in case of error
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_test_SeqList[[2]]), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_test_SeqList[[2]])
aalab1[aalab1 == 0] <- 0.21
aalab1[aalab1 == 1] <- 0.79
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_4/weight_test.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_4/weight_test.txt", sep = "\t" )
}


# coop
# coop set1
# "ESR1_2", "ESR1_2"
# "ESR1_2", "FOXA1"
# "ESR1_2", "NR2F2"
# "ESR1_2", "PGR"
# "ESR1_2", "SP1"
# "ESR1_2", "YBX1"
# "NR2F2", "YBX1"
aacoop <- rbind(c(1, 1), c(1, 2), c(1, 6),
                c(1, 10), c(1,15),
                c(1, 16), c(6, 16))

aamin_LLR <- numeric(length(TF.motifs.Expanded_LLR_to_pval_pseudo))
aamaxL_LLR <- numeric(length(TF.motifs.Expanded_LLR_to_pval_pseudo))
for(i in 1:length(TF.motifs.Expanded_LLR_to_pval_pseudo)){
  aamaxL_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][1,1] / 1000
  aa <- which(TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][,2] > 0.0001)[1]
  if(aa == 1){
    aamin_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][aa, 1] / 1000
  }else{
    aamin_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][(aa-1), 1] / 1000
  }
}
names(aamin_LLR) <- names(TF.motifs.Expanded_LLR_to_pval_pseudo)
names(aamaxL_LLR) <- names(TF.motifs.Expanded_LLR_to_pval_pseudo)
aathr <- 1- (aamin_LLR/aamaxL_LLR)
names(aathr)[names(aathr) == "NKX3-1"] <- "NKX3_1"

aa_bind_ff_ens <- c(0,0,0,0,0,1,1,0,1,1,0,1,0,0,0,1)
aa_alpha_ff_ens<- c(0,0,0,0,0,1,1,0,1,1,0,1,0,0,0,1)
names(aa_bind_ff_ens) <- aa_names
names(aa_alpha_ff_ens) <- aa_names

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_4")
par_ff_lb_ub_coop_Writer_multi_enh(.TFnames=aa_names,
                                   .annotation_thresh=aathr[aa_names],
                                   .annotation_range=numeric(0),
                                   .initial_bind_w=numeric(0),
                                   .bind_w_range=numeric(0),
                                   .initial_alpha=numeric(0),
                                   .initial_log_reg_bias = numeric(0),
                                   .initial_log_reg_coeff = 5,
                                   .initial_log_reg_bias_range = numeric(0),
                                   .initial_log_reg_coeff_range = c(1, 20),
                                   .alpha_range=numeric(0),
                                   .coop_tf_mat=aacoop,
                                   .initial_coop_weight=numeric(0),
                                   .coop_weight_range=numeric(0),
                                   .coop_type=character(0),
                                   .coop_dist=integer(0),
                                   .coop_orientation=integer(0),
                                   one_qbtm_per_enh=F,
                                   .initial_qBTM=numeric(0),
                                   .qBTMrange=numeric(0),
                                   one_beta_per_enh=F,
                                   .initial_pi_beta=numeric(0), 
                                   .pi_beta_upper=numeric(0),
                                   .pi_beta_lower=numeric(0),
                                   .annotation_thresh_ff=integer(0),
                                   .initial_bind_w_ff=integer(0),
                                   .initial_alpha_ff=integer(0),
                                   .initial_coop_weight_ff=integer(0),
                                   .initial_qBTM_ff=1,
                                   .initial_pi_beta_ff=integer(0),
                                   .filename_start=character(0),
                                   .filename_ff=character(0),
                                   .filename_upper=character(0),
                                   .filename_lower=character(0),
                                   .file_name_coop=character(0),
                                   .initial_log_reg_bias_ff = 1,
                                   .initial_log_reg_coeff_ff = 1,
                                   .nu_enhacners=240,
                                   ensemble_mode = T,
                                   .nu_samples = 1,
                                   .annotation_thresh_ff_ens = rep(0, 14),
                                   .initial_coop_weight_ff_ens = c(0,0,0,1,0,0,0),
                                   .initial_qBTM_ff_ens = 0,
                                   .initial_log_reg_bias_ff_ens = 0,
                                   .initial_log_reg_coeff_ff_ens = 0,
                                   Write_to_file = T,
                                   create_folders = T,
                                   TF_role_decided = F,
                                   factor_role_vector=integer(0),  
                                   .create_bounds=T,
                                   .create_par=T,
                                   .create_ff=T,
                                   .create_coop=T,
                                   logistic_params = T,
                                   .initial_bind_w_ff_ens = aa_bind_ff_ens,
                                   .initial_alpha_ff_ens = aa_alpha_ff_ens)

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")
# write GEMSTAT jobs
#cat("#!/bin/bash\n", file = "Ensemble_Exp3_job_train_hal.job", append = F)
aa_par_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_4/Parameters")
aa_par_files_trimmed <- unlist(lapply(strsplit(aa_par_files, "\\."), "[[", 1))
for(i in 1:length(aa_par_files)){
  cat(c("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr -s Sequence/train_seq.fa",
        "-e Labels/Label_train.txt",
        "-m motifs.wtmx -f TF_exp.tab", 
        paste0("-fo Out/output_", aa_par_files_trimmed[i], ".txt"), 
        "-o DIRECT -c Coop/coop.par ",
        paste0("-p Parameters/", aa_par_files[i]),
        "-lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -na 2 -onebeta 1",
        "-train_weights weight_train.txt > ", 
        paste0("Log/log_", aa_par_files_trimmed[i], ".txt\n")),
      sep = " ", append = T, file = "Ensemble_Exp4_job_train_hal.job")
}

#write test script job for hal
#cat("#!/bin/bash\n", file = "Ensemble_Exp2_job_test.sh", append = F)
aa_par_trained_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_3/Trained_par")
aa_par_trained_files_trimmed <- unlist(lapply(strsplit(aa_par_trained_files, "\\."), "[[", 1))
for(i in 1:length(aa_par_trained_files)){
  cat(c("~/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/src/seq2expr -s Sequence/test_seq.fa",
        "-e Labels/Label_test.txt",
        "-m motifs.wtmx -f TF_exp.tab", 
        paste0("-fo Out/output_", aa_par_trained_files_trimmed[i], "_test.txt"), 
        "-o DIRECT -c Coop/coop.par ",
        paste0("-p Trained_par/", aa_par_trained_files[i]),
        "-lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -na 0 -onebeta 1",
        "-train_weights weight_test.txt > ", 
        paste0("Log/log_", aa_par_trained_files_trimmed[i], "_test.txt\n")),
      sep = " ", append = T, file = "Ensemble_Exp2_job_test.sh")
}

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")
# read train-test results

# read train-test results
E_RNA_GEMSTAT_Ensemble_Outlist[[4]] <- read_output_train_test_GEMSTAT_ensemble(output_dir = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_4/Out",
                                                                               should_plot = F)
###### read parameters
E_RNA_GEMSTAT_Ensemble_Parlist[[4]] <- read_parameters_GEMSTAT_ensemble(par_directory = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_4/Trained_par/",
                                                                        convert_to_mat = T)

library(gplots)
sort(E_RNA_GEMSTAT_Ensemble_Outlist[[4]]$test_ROC, decreasing = T)[1:20]
aa_perf_srt <-names(E_RNA_GEMSTAT_Ensemble_Outlist[[4]]$test_ROC)[sort(E_RNA_GEMSTAT_Ensemble_Outlist[[4]]$test_ROC, decreasing = T, index.return=T)$ix]
aa_perf_srt2 <- unlist(lapply(lapply(strsplit(aa_perf_srt, split="_"), "[", c(3,4,5)), paste, collapse="_"))
aa_perf_srt22 <- paste0("log_", aa_perf_srt2)
heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[4]]$binding[aa_perf_srt22[c(c(1:400), c(7001:7400))], ]), Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none')
heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[4]]$alpha[aa_perf_srt22[c(c(1:400), c(7001:7400))], ] + 1e-6), Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none')
heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[4]]$coop[aa_perf_srt22[c(c(1:400), c(7001:7400))], ]), Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none', margins = c(10,10))
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[4]]$qbtm[aa_perf_srt22], main = "qbtm_sorted_by_perf_exp4", pch = 16, ylab = "")
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[4]]$beta[aa_perf_srt22], main = "beta_sorted_by_perf_exp4", pch = 16)

heatmap.2(E_RNA_GEMSTAT_Ensemble_Parlist[[4]]$logistic[aa_perf_srt22,], Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none', margins = c(10,10))
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[4]]$logistic[aa_perf_srt22, 1], ylab="", main = "logistic_bias_sorted_by_perf_exp4")
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[4]]$logistic[aa_perf_srt22, 2], ylab="", main = "logistic_coeff_sorted_by_perf_exp4")

# plot the best models
which.max(E_RNA_GEMSTAT_Ensemble_Outlist[[4]]$test_ROC)
which.max(E_RNA_GEMSTAT_Ensemble_Outlist[[4]]$test_PRC)
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[4]]$train_results[[3835]]$ROC_curve, main = "ROC train bestModel exp4")
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[4]]$train_results[[3835]]$PRC_curve, main = "PRC train bestModel exp4")
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[4]]$test_results[[3835]]$ROC_curve, main = "ROC test bestModel exp4")
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[4]]$test_results[[3835]]$PRC_curve, main = "PRC test bestModel exp4")

######################################################################################################
##########################################################################################
# experiment 5: using all features that came out to be important from RF model

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/")
system("mkdir Experiment_5")
system("mkdir Experiment_5/Sequence")
system("mkdir Experiment_5/Labels")
system("mkdir Experiment_5/Log")
system("mkdir Experiment_5/Out")
# Sequence
GEMSTAT_Ensemble_validation_SeqList <- list()
aalab <- unlist(lapply(strsplit(data_partition_mut_GEMSTAT_trimmed[[2]], split="_"), "[[", 1))
aa_val_ind_pos <- data_partition_mut_GEMSTAT_trimmed[[2]][sample(x = which(aalab %in% "pos"), size = 90, replace = F)]
aa_val_ind_neg <- data_partition_mut_GEMSTAT_trimmed[[2]][sample(x = which(aalab %in% "neg"), size = 334, replace = F)]

aamatchp <- match(aa_val_ind_pos, names(Positive_set_seq_list_char_1000bp[[4]]))
aamatchn <- match(aa_val_ind_neg, names(Negative_set_seq_list_char_1000bp[[4]]))

aatpos <- Positive_set_seq_list_char_1000bp[[4]][aamatchp]
aatneg <- Negative_set_seq_list_char_1000bp[[4]][aamatchn]
GEMSTAT_Ensemble_validation_SeqList[[1]] <- sample(x = c(aatpos, aatneg), size = length(c(aatpos, aatneg)), replace = F)

library(ShortRead)

writeFasta(DNAStringSet(GEMSTAT_Ensemble_train_SeqList[[2]]), width = 1000,
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Sequence/train_seq.fa")
writeFasta(DNAStringSet(GEMSTAT_Ensemble_test_SeqList[[2]]), width = 1000, 
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Sequence/test_seq.fa")
writeFasta(DNAStringSet(GEMSTAT_Ensemble_validation_SeqList[[1]]), width = 1000, 
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Sequence/valid_seq.fa")

# Motif
aa_names <- c("ESR1_2", "ESR1_3","FOXA1","JUN_1","LEF1","NFIB","NKX3_1","NR2F2","NR3C1","NR5A2",
              "PBX1","PGR","PPARD","PPARG","RARA","RUNX1","SP1","YBX1")
# coop 
# "ESR1_2", "ESR1_2"
# "ESR1_2", "ESR1_3"
# "ESR1_2", "FOXA1"
# "ESR1_2", "NR2F2"
# "ESR1_2", "PGR"
# "ESR1_2", "PPARD"
# "ESR1_2", "PPARG"
# "ESR1_2", "YBX1"
# "NR2F2", "NR2F2"
# "NR2F2", "YBX1"
aacoop <- rbind(c(1, 1), c(1, 2), c(1, 3), c(1, 8),
                c(1, 12), c(1,13), c(1, 14), c(1, 18),
                c(8, 8), c(8, 16))

MotifWriter(TF.motifs.Expanded_pseudo_count[aa_names],
            pseudo = 0.001, 
            output.File.Name = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/motifs")


# TF expression
aaTFexp <- rep(1, length(aa_names))
names(aaTFexp) <- aa_names
cat(c("Rows", paste0("1", "\n")), sep = "\t", append = F,
    file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/TF_exp.tab")
for(i in 1:length(aaTFexp)){
  cat(names(aaTFexp)[i], paste0(aaTFexp[i], rep("\n", as.integer(i != length(aaTFexp)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/TF_exp.tab",
      sep = "\t" )
}
# labels
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_train_SeqList[[2]]), split = "_")), "[[", 1))

aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_train_SeqList[[2]])
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Labels/Label_train.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Labels/Label_train.txt", sep = "\t" )
}

# label for validation
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_validation_SeqList[[1]]), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_validation_SeqList[[1]])
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Labels/Label_validation.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Labels/Label_validation.txt", sep = "\t" )
}

# label for test
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_test_SeqList[[2]]), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_test_SeqList[[2]])
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Labels/Label_test.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Labels/Label_test.txt", sep = "\t" )
}


#weights (in this case not equal)
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_train_SeqList[[2]]), split = "_")), "[[", 1))

aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_train_SeqList[[2]])

aalab1[aalab1 == 0] <- 0.25
aalab1[aalab1 == 1] <- 0.75
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/weight_train.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/weight_train.txt", sep = "\t" )
}

# weights for test just in case of error
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_test_SeqList[[2]]), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_test_SeqList[[2]])
aalab1[aalab1 == 0] <- 0.21
aalab1[aalab1 == 1] <- 0.79
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/weight_test.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/weight_test.txt", sep = "\t" )
}

# weights for validation just in case of error
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_validation_SeqList[[1]]), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_validation_SeqList[[1]])
aalab1[aalab1 == 0] <- 90/(334 + 90)
aalab1[aalab1 == 1] <- 334/(334 + 90)
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/weight_validation.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/weight_validation.txt", sep = "\t" )
}



aamin_LLR <- numeric(length(TF.motifs.Expanded_LLR_to_pval_pseudo))
aamaxL_LLR <- numeric(length(TF.motifs.Expanded_LLR_to_pval_pseudo))
for(i in 1:length(TF.motifs.Expanded_LLR_to_pval_pseudo)){
  aamaxL_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][1,1] / 1000
  aa <- which(TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][,2] > 0.0001)[1]
  if(aa == 1){
    aamin_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][aa, 1] / 1000
  }else{
    aamin_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][(aa-1), 1] / 1000
  }
}
names(aamin_LLR) <- names(TF.motifs.Expanded_LLR_to_pval_pseudo)
names(aamaxL_LLR) <- names(TF.motifs.Expanded_LLR_to_pval_pseudo)
aathr <- 1- (aamin_LLR/aamaxL_LLR)
names(aathr)[names(aathr) == "NKX3-1"] <- "NKX3_1"

aa_bind_ff_ens <- c(0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,1,0,0)
aa_alpha_ff_ens<- c(0,0,0,1,0,0,1,0,1,0,0,1,0,1,0,1,0,0)
names(aa_bind_ff_ens) <- aa_names
names(aa_alpha_ff_ens) <- aa_names

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5")
par_ff_lb_ub_coop_Writer_multi_enh(.TFnames=aa_names,
                                   .annotation_thresh=aathr[aa_names],
                                   .annotation_range=numeric(0),
                                   .initial_bind_w=numeric(0),
                                   .bind_w_range=numeric(0),
                                   .initial_alpha=numeric(0),
                                   .initial_log_reg_bias = numeric(0),
                                   .initial_log_reg_coeff = 5,
                                   .initial_log_reg_bias_range = numeric(0),
                                   .initial_log_reg_coeff_range = c(1, 20),
                                   .alpha_range=numeric(0),
                                   .coop_tf_mat=aacoop,
                                   .initial_coop_weight=numeric(0),
                                   .coop_weight_range=numeric(0),
                                   .coop_type=character(0),
                                   .coop_dist=integer(0),
                                   .coop_orientation=integer(0),
                                   one_qbtm_per_enh=F,
                                   .initial_qBTM=numeric(0),
                                   .qBTMrange=numeric(0),
                                   one_beta_per_enh=F,
                                   .initial_pi_beta=numeric(0), 
                                   .pi_beta_upper=numeric(0),
                                   .pi_beta_lower=numeric(0),
                                   .annotation_thresh_ff=integer(0),
                                   .initial_bind_w_ff=integer(0),
                                   .initial_alpha_ff=integer(0),
                                   .initial_coop_weight_ff=integer(0),
                                   .initial_qBTM_ff=1,
                                   .initial_pi_beta_ff=integer(0),
                                   .filename_start=character(0),
                                   .filename_ff=character(0),
                                   .filename_upper=character(0),
                                   .filename_lower=character(0),
                                   .file_name_coop=character(0),
                                   .initial_log_reg_bias_ff = 1,
                                   .initial_log_reg_coeff_ff = 1,
                                   .nu_enhacners=600,
                                   ensemble_mode = T,
                                   .nu_samples = 1,
                                   .annotation_thresh_ff_ens = rep(0, 14),
                                   .initial_coop_weight_ff_ens = c(0,0,0,0,1,0,0,0,0,0),
                                   .initial_qBTM_ff_ens = 0,
                                   .initial_log_reg_bias_ff_ens = 0,
                                   .initial_log_reg_coeff_ff_ens = 0,
                                   Write_to_file = T,
                                   create_folders = T,
                                   TF_role_decided = F,
                                   factor_role_vector=integer(0),  
                                   .create_bounds=T,
                                   .create_par=T,
                                   .create_ff=T,
                                   .create_coop=T,
                                   logistic_params = T,
                                   .initial_bind_w_ff_ens = aa_bind_ff_ens,
                                   .initial_alpha_ff_ens = aa_alpha_ff_ens)

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")
# write GEMSTAT jobs
#cat("#!/bin/bash\n", file = "Ensemble_Exp3_job_train_hal.job", append = F)
aa_par_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Parameters")
aa_par_files_trimmed <- unlist(lapply(strsplit(aa_par_files, "\\."), "[[", 1))
for(i in 1:length(aa_par_files)){
  cat(c("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr -s Sequence/train_seq.fa",
        "-e Labels/Label_train.txt",
        "-m motifs.wtmx -f TF_exp.tab", 
        paste0("-fo Out/output_", aa_par_files_trimmed[i], ".txt"), 
        "-o DIRECT -c Coop/coop.par ",
        paste0("-p Parameters/", aa_par_files[i]),
        "-lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -na 2 -onebeta 1",
        "-train_weights weight_train.txt > ", 
        paste0("Log/log_", aa_par_files_trimmed[i], ".txt\n")),
      sep = " ", append = T, file = "Ensemble_Exp5_job_train_hal.job")
}

#write test script job for hal
#cat("#!/bin/bash\n", file = "Ensemble_Exp2_job_test.sh", append = F)
aa_par_trained_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_3/Trained_par")
aa_par_trained_files_trimmed <- unlist(lapply(strsplit(aa_par_trained_files, "\\."), "[[", 1))
for(i in 1:length(aa_par_trained_files)){
  cat(c("~/Documents/GEMSTAT_git/my_GEM_fork/GEMSTAT/src/seq2expr -s Sequence/test_seq.fa",
        "-e Labels/Label_test.txt",
        "-m motifs.wtmx -f TF_exp.tab", 
        paste0("-fo Out/output_", aa_par_trained_files_trimmed[i], "_test.txt"), 
        "-o DIRECT -c Coop/coop.par ",
        paste0("-p Trained_par/", aa_par_trained_files[i]),
        "-lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -na 0 -onebeta 1",
        "-train_weights weight_test.txt > ", 
        paste0("Log/log_", aa_par_trained_files_trimmed[i], "_test.txt\n")),
      sep = " ", append = T, file = "Ensemble_Exp2_job_test.sh")
}

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")
# read train-test results

# read train-test results
E_RNA_GEMSTAT_Ensemble_Outlist[[5]] <- read_output_train_test_GEMSTAT_ensemble(output_dir = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Out",
                                                                               should_plot = F)
###### read parameters
E_RNA_GEMSTAT_Ensemble_Parlist[[5]] <- read_parameters_GEMSTAT_ensemble(par_directory = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Trained_par/",
                                                                        convert_to_mat = T)

library(gplots)
sort(E_RNA_GEMSTAT_Ensemble_Outlist[[5]]$test_ROC, decreasing = T)[1:20]
aa_perf_srt <-names(E_RNA_GEMSTAT_Ensemble_Outlist[[5]]$test_ROC)[sort(E_RNA_GEMSTAT_Ensemble_Outlist[[5]]$test_ROC, decreasing = T, index.return=T)$ix]
aa_perf_srt2 <- unlist(lapply(lapply(strsplit(aa_perf_srt, split="_"), "[", c(3,4,5)), paste, collapse="_"))
aa_perf_srt22 <- paste0("log_", aa_perf_srt2)
heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[5]]$binding[aa_perf_srt22, ]), Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none')
#[c(c(1:400), c(7001:7400))]
heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[5]]$alpha[aa_perf_srt22, ] + 1e-6), Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none')
heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[5]]$coop[aa_perf_srt22, ]), Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none', margins = c(10,10))
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[5]]$qbtm[aa_perf_srt22], 
     main = "qbtm_sorted_by_perf_exp5", pch = 16, ylab = "")
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[5]]$beta[aa_perf_srt22], 
     main = "beta_sorted_by_perf_exp5", pch = 16)

heatmap.2(E_RNA_GEMSTAT_Ensemble_Parlist[[5]]$logistic[aa_perf_srt22,], Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none', margins = c(10,10))
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[5]]$logistic[aa_perf_srt22, 1], ylab="", main = "logistic_bias_sorted_by_perf_exp5")
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[5]]$logistic[aa_perf_srt22, 2], ylab="", main = "logistic_coeff_sorted_by_perf_exp5")

# plot the best models
which.max(E_RNA_GEMSTAT_Ensemble_Outlist[[5]]$train_ROC)
which.max(E_RNA_GEMSTAT_Ensemble_Outlist[[5]]$train_PRC)
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[5]]$train_results[[1051]]$ROC_curve, main = "ROC train bestModel exp5")
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[5]]$train_results[[1051]]$PRC_curve, main = "PRC train bestModel exp5")
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[5]]$test_results[[1051]]$ROC_curve, main = "ROC test bestModel exp5")
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[5]]$test_results[[1051]]$PRC_curve, main = "PRC test bestModel exp5")

######################################################################################################
# experiment 6: CCV mode
CCV_training_partition_list
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/")
system("mkdir Experiment_6")
system("mkdir Experiment_6/Sequence")
system("mkdir Experiment_6/Labels")
system("mkdir Experiment_6/Log")
system("mkdir Experiment_6/Out")

data_partition_mut_GEMSTAT_trimmed_ccv <- list()
data_partition_mut_GEMSTAT_trimmed_ccv[[1]] <- unlist(lapply(lapply(strsplit(CCV_training_partition_list[[1]], split = "_"), "[", c(1, 2)), paste, collapse = "_"))
data_partition_mut_GEMSTAT_trimmed_ccv[[2]] <- unlist(lapply(lapply(strsplit(setdiff(rownames(my_learning_datasets_1000bp_4_ForGEMSTAT[[1]]),
                                                                                     CCV_training_partition_list[[1]]),
                                                                             split = "_"), "[", c(1, 2)), paste, collapse = "_"))
aalab <- unlist(lapply(strsplit(data_partition_mut_GEMSTAT_trimmed_ccv[[1]], split="_"), "[[", 1))
aa_train_ind_pos <- data_partition_mut_GEMSTAT_trimmed_ccv[[1]][sample(x = which(aalab %in% "pos"), size = sum(aalab %in% "pos"), replace = F)]
aa_train_ind_neg <- data_partition_mut_GEMSTAT_trimmed_ccv[[1]][sample(x = which(aalab %in% "neg"), size = sum(aalab %in% "neg"), replace = F)]

aamatchp <- match(aa_train_ind_pos, names(Positive_set_seq_list_char_1000bp[[4]]))
aamatchn <- match(aa_train_ind_neg, names(Negative_set_seq_list_char_1000bp[[4]]))

aatpos <- Positive_set_seq_list_char_1000bp[[4]][aamatchp]
aatneg <- Negative_set_seq_list_char_1000bp[[4]][aamatchn]
GEMSTAT_Ensemble_train_SeqList[[3]] <- sample(x = c(aatpos, aatneg), size = length(c(aatpos, aatneg)), replace = F)

aalab <- unlist(lapply(strsplit(data_partition_mut_GEMSTAT_trimmed_ccv[[2]], split="_"), "[[", 1))
aa_test_ind_pos <- data_partition_mut_GEMSTAT_trimmed_ccv[[2]][sample(x = which(aalab %in% "pos"), size = sum(aalab %in% "pos"), replace = F)]
aa_test_ind_neg <- data_partition_mut_GEMSTAT_trimmed_ccv[[2]][sample(x = which(aalab %in% "neg"), size = sum(aalab %in% "neg"), replace = F)]

aamatchp <- match(aa_test_ind_pos, names(Positive_set_seq_list_char_1000bp[[4]]))
aamatchn <- match(aa_test_ind_neg, names(Negative_set_seq_list_char_1000bp[[4]]))

aatpos <- Positive_set_seq_list_char_1000bp[[4]][aamatchp]
aatneg <- Negative_set_seq_list_char_1000bp[[4]][aamatchn]
GEMSTAT_Ensemble_test_SeqList[[3]] <- sample(x = c(aatpos, aatneg),
                                             size = length(c(aatpos, aatneg)), replace = F)

# Sequence
library(ShortRead)

writeFasta(DNAStringSet(GEMSTAT_Ensemble_train_SeqList[[3]]), width = 1000,
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_6/Sequence/train_seq.fa")
writeFasta(DNAStringSet(GEMSTAT_Ensemble_test_SeqList[[3]]), width = 1000, 
           file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_6/Sequence/test_seq.fa")

# Motif
aa_names <- c("ESR1_2", "ESR1_3","FOXA1","JUN_1","LEF1","NFIB","NKX3_1","NR2F2","NR3C1","NR5A2",
              "PBX1","PGR","PPARD","PPARG","RARA","RUNX1","SP1","YBX1")
# coop 
# "ESR1_2", "ESR1_2"
# "ESR1_2", "ESR1_3"
# "ESR1_2", "FOXA1"
# "ESR1_2", "NR2F2"
# "ESR1_2", "PGR"
# "ESR1_2", "PPARD"
# "ESR1_2", "PPARG"
# "ESR1_2", "YBX1"
# "NR2F2", "NR2F2"
# "NR2F2", "YBX1"
aacoop <- rbind(c(1, 1), c(1, 2), c(1, 3), c(1, 8),
                c(1, 12), c(1,13), c(1, 14), c(1, 18),
                c(8, 8), c(8, 16))

MotifWriter(TF.motifs.Expanded_pseudo_count[aa_names],
            pseudo = 0.001, 
            output.File.Name = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_6/motifs")


# TF expression
aaTFexp <- rep(1, length(aa_names))
names(aaTFexp) <- aa_names
cat(c("Rows", paste0("1", "\n")), sep = "\t", append = F,
    file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_6/TF_exp.tab")
for(i in 1:length(aaTFexp)){
  cat(names(aaTFexp)[i], paste0(aaTFexp[i], rep("\n", as.integer(i != length(aaTFexp)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_6/TF_exp.tab",
      sep = "\t" )
}
# labels
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_train_SeqList[[3]]), split = "_")), "[[", 1))

aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_train_SeqList[[3]])
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_6/Labels/Label_train.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_6/Labels/Label_train.txt", sep = "\t" )
}
# label for test
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_test_SeqList[[3]]), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_test_SeqList[[3]])
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_6/Labels/Label_test.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_6/Labels/Label_test.txt", sep = "\t" )
}

#weights (in this case not equal)
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_train_SeqList[[3]]), split = "_")), "[[", 1))

aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_train_SeqList[[3]])

aalab1[aalab1 == 0] <- 336/(336+1393)
aalab1[aalab1 == 1] <- 1393/(336+1393)
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_6/weight_train.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_6/weight_train.txt", sep = "\t" )
}
# weights for test just in case of error
aalab1 <- unlist(lapply((strsplit(names(GEMSTAT_Ensemble_test_SeqList[[3]]), split = "_")), "[[", 1))
aalab1[aalab1 == "pos"] <- 1
aalab1[aalab1 == "neg"] <- 0
aalab1 <- as.numeric(aalab1)
aalab1 <- matrix(aalab1, nrow = length(aalab1))
rownames(aalab1) <- names(GEMSTAT_Ensemble_test_SeqList[[3]])
aalab1[aalab1 == 0] <- 113/(113 + 276)
aalab1[aalab1 == 1] <- 276/(113 + 276)
cat(c("Rows", paste0("1", "\n")), sep = "\t",
    append = F, file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_6/weight_test.txt")
for(i in 1:length(aalab1)){
  cat(rownames(aalab1)[i], paste0(aalab1[i, 1],
                                  rep("\n", as.integer(i != length(aalab1)))), 
      append = T,
      file = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_6/weight_test.txt", sep = "\t" )
}



aamin_LLR <- numeric(length(TF.motifs.Expanded_LLR_to_pval_pseudo))
aamaxL_LLR <- numeric(length(TF.motifs.Expanded_LLR_to_pval_pseudo))
for(i in 1:length(TF.motifs.Expanded_LLR_to_pval_pseudo)){
  aamaxL_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][1,1] / 1000
  aa <- which(TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][,2] > 0.0001)[1]
  if(aa == 1){
    aamin_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][aa, 1] / 1000
  }else{
    aamin_LLR[i] <- TF.motifs.Expanded_LLR_to_pval_pseudo[[i]][(aa-1), 1] / 1000
  }
}
names(aamin_LLR) <- names(TF.motifs.Expanded_LLR_to_pval_pseudo)
names(aamaxL_LLR) <- names(TF.motifs.Expanded_LLR_to_pval_pseudo)
aathr <- 1- (aamin_LLR/aamaxL_LLR)
names(aathr)[names(aathr) == "NKX3-1"] <- "NKX3_1"

aa_bind_ff_ens <- c(0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,1,0,0)
aa_alpha_ff_ens<- c(0,0,0,1,0,0,1,0,1,0,0,1,0,1,0,1,0,0)
names(aa_bind_ff_ens) <- aa_names
names(aa_alpha_ff_ens) <- aa_names

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_6")
par_ff_lb_ub_coop_Writer_multi_enh(.TFnames=aa_names,
                                   .annotation_thresh=aathr[aa_names],
                                   .annotation_range=numeric(0),
                                   .initial_bind_w=numeric(0),
                                   .bind_w_range=numeric(0),
                                   .initial_alpha=numeric(0),
                                   .initial_log_reg_bias = numeric(0),
                                   .initial_log_reg_coeff = 5,
                                   .initial_log_reg_bias_range = numeric(0),
                                   .initial_log_reg_coeff_range = c(1, 20),
                                   .alpha_range=numeric(0),
                                   .coop_tf_mat=aacoop,
                                   .initial_coop_weight=numeric(0),
                                   .coop_weight_range=numeric(0),
                                   .coop_type=character(0),
                                   .coop_dist=integer(0),
                                   .coop_orientation=integer(0),
                                   one_qbtm_per_enh=F,
                                   .initial_qBTM=numeric(0),
                                   .qBTMrange=numeric(0),
                                   one_beta_per_enh=F,
                                   .initial_pi_beta=numeric(0), 
                                   .pi_beta_upper=numeric(0),
                                   .pi_beta_lower=numeric(0),
                                   .annotation_thresh_ff=integer(0),
                                   .initial_bind_w_ff=integer(0),
                                   .initial_alpha_ff=integer(0),
                                   .initial_coop_weight_ff=integer(0),
                                   .initial_qBTM_ff=1,
                                   .initial_pi_beta_ff=integer(0),
                                   .filename_start=character(0),
                                   .filename_ff=character(0),
                                   .filename_upper=character(0),
                                   .filename_lower=character(0),
                                   .file_name_coop=character(0),
                                   .initial_log_reg_bias_ff = 1,
                                   .initial_log_reg_coeff_ff = 1,
                                   .nu_enhacners=600,
                                   ensemble_mode = T,
                                   .nu_samples = 1,
                                   .annotation_thresh_ff_ens = rep(0, 14),
                                   .initial_coop_weight_ff_ens = c(0,0,0,0,1,0,0,0,0,0),
                                   .initial_qBTM_ff_ens = 0,
                                   .initial_log_reg_bias_ff_ens = 0,
                                   .initial_log_reg_coeff_ff_ens = 0,
                                   Write_to_file = T,
                                   create_folders = T,
                                   TF_role_decided = F,
                                   factor_role_vector=integer(0),  
                                   .create_bounds=T,
                                   .create_par=T,
                                   .create_ff=T,
                                   .create_coop=T,
                                   logistic_params = T,
                                   .initial_bind_w_ff_ens = aa_bind_ff_ens,
                                   .initial_alpha_ff_ens = aa_alpha_ff_ens)

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")
# write GEMSTAT jobs
aa_par_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_6/Parameters")
aa_par_files_trimmed <- unlist(lapply(strsplit(aa_par_files, "\\."), "[[", 1))
for(i in 1:length(aa_par_files)){
  cat(c("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr -s Sequence/train_seq.fa",
        "-e Labels/Label_train.txt",
        "-m motifs.wtmx -f TF_exp.tab", 
        paste0("-fo Out/output_", aa_par_files_trimmed[i], ".txt"), 
        "-o DIRECT -c Coop/coop.par ",
        paste0("-p Parameters/", aa_par_files[i]),
        "-lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -na 2 -onebeta 1",
        "-train_weights weight_train.txt > ", 
        paste0("Log/log_", aa_par_files_trimmed[i], ".txt\n")),
      sep = " ", append = T, file = "Ensemble_Exp6_job_train_hal.job")
}

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor")
# read train-test results

# read train-test results
E_RNA_GEMSTAT_Ensemble_Outlist[[6]] <- read_output_train_test_GEMSTAT_ensemble(output_dir = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_6/Out",
                                                                               should_plot = F)
###### read parameters
E_RNA_GEMSTAT_Ensemble_Parlist[[6]] <- read_parameters_GEMSTAT_ensemble(par_directory = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_6/Trained_par/",
                                                                        convert_to_mat = T)

library(gplots)
sort(E_RNA_GEMSTAT_Ensemble_Outlist[[6]]$test_ROC, decreasing = T)[1:20]
aa_perf_srt <-names(E_RNA_GEMSTAT_Ensemble_Outlist[[6]]$test_ROC)[sort(E_RNA_GEMSTAT_Ensemble_Outlist[[6]]$test_ROC, decreasing = T, index.return=T)$ix]
aa_perf_srt2 <- unlist(lapply(lapply(strsplit(aa_perf_srt, split="_"), "[", c(3,4,5)), paste, collapse="_"))
aa_perf_srt22 <- paste0("log_", aa_perf_srt2)
heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[6]]$binding[aa_perf_srt22[c(c(1:400), c(7001:7400))], ]), Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none')
heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[6]]$alpha[aa_perf_srt22[c(c(1:400), c(7001:7400))], ] + 1e-6), Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none')
heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[6]]$coop[aa_perf_srt22[c(c(1:400), c(7001:7400))], ]), Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none', margins = c(10,10))
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[6]]$qbtm[aa_perf_srt22], main = "qbtm_sorted_by_perf_exp6", pch = 16, ylab = "")
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[6]]$beta[aa_perf_srt22], main = "beta_sorted_by_perf_exp6", pch = 16)

heatmap.2(E_RNA_GEMSTAT_Ensemble_Parlist[[6]]$logistic[aa_perf_srt22,], Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none', margins = c(10,10))
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[6]]$logistic[aa_perf_srt22, 1], 
     ylab="", main = "logistic_bias_sorted_by_perf_exp6")
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[6]]$logistic[aa_perf_srt22, 2],
     ylab="", main = "logistic_coeff_sorted_by_perf_exp6")

# plot the best models
which.max(E_RNA_GEMSTAT_Ensemble_Outlist[[6]]$train_ROC)
which.max(E_RNA_GEMSTAT_Ensemble_Outlist[[6]]$train_PRC)
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[6]]$train_results[[4587]]$ROC_curve, main = "ROC train bestModel exp6")
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[6]]$train_results[[4587]]$PRC_curve, main = "PRC train bestModel exp6")
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[6]]$test_results[[4587]]$ROC_curve, main = "ROC test bestModel exp6")
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[6]]$test_results[[4587]]$PRC_curve, main = "PRC test bestModel exp6")

########## Testing something exactly similar to exp 5, except fixing the alpha param of half-ER to be 1
aa_par_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Parameters")
aa_par_files_trimmed <- unlist(lapply(strsplit(aa_par_files, "\\."), "[[", 1))
for(i in 1:length(aa_par_files)){
  cat(c("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr -s Sequence/train_seq.fa",
        "-e Labels/Label_train.txt",
        "-m motifs.wtmx -f TF_exp.tab", 
        paste0("-fo Out_hER/output_", aa_par_files_trimmed[i], ".txt"), 
        "-o DIRECT -c Coop/coop.par ",
        paste0("-p Parameters/", aa_par_files[i]),
        "-lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix_hER/ff.par -classifier -na 2 -onebeta 1",
        "-train_weights weight_train.txt > ", 
        paste0("Log_hER/log_", aa_par_files_trimmed[i], ".txt\n")),
      sep = " ", append = T, file = "Ensemble_Exp5_job_train_hal_hER.job")
}

# read train-test results
E_RNA_GEMSTAT_Ensemble_Outlist[[7]] <- read_output_train_test_GEMSTAT_ensemble(output_dir = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Out_hER/",
                                                                               should_plot = F, validation = T)
###### read parameters
E_RNA_GEMSTAT_Ensemble_Parlist[[7]] <- read_parameters_GEMSTAT_ensemble(par_directory = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Trained_par_hER/",
                                                                        convert_to_mat = T)

library(gplots)
sort(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$test_PRC, decreasing = T)[1:20]
aa_perf_srt <-names(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$test_PRC)[sort(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$test_PRC,
                                                                       decreasing = T, index.return=T)$ix]
aa_perf_srt2 <- unlist(lapply(lapply(strsplit(aa_perf_srt, split="_"), "[", c(3,4,5)), paste, collapse="_"))
aa_perf_srt22 <- paste0("log_", aa_perf_srt2)
# [c(c(1:400), c(7001:7400))]
heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[7]]$binding[aa_perf_srt22, ]), Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none')
heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[7]]$alpha[aa_perf_srt22, ] + 1e-6), Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none')
heatmap.2(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[7]]$coop[aa_perf_srt22, ]), Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none', margins = c(10,10))

aa_all_param <- cbind(log10(E_RNA_GEMSTAT_Ensemble_Parlist[[7]]$binding[aa_perf_srt22, ]), 
                      log10(E_RNA_GEMSTAT_Ensemble_Parlist[[7]]$alpha[aa_perf_srt22, ] + 1e-2), 
                      log10(E_RNA_GEMSTAT_Ensemble_Parlist[[7]]$coop[aa_perf_srt22, ]), 
                      E_RNA_GEMSTAT_Ensemble_Parlist[[7]]$qbtm[aa_perf_srt22], 
                      log10(E_RNA_GEMSTAT_Ensemble_Parlist[[7]]$beta[aa_perf_srt22]))
colnames(aa_all_param)[1:18] <- paste(colnames(aa_all_param)[1:18], "_binding")
colnames(aa_all_param)[19:36] <- paste(colnames(aa_all_param)[19:36], "_activation")
colnames(aa_all_param)[47] <- "qbtm"
colnames(aa_all_param)[48] <- "beta"

heatmap.2(aa_all_param, Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none', margins = c(9,5))


plot(E_RNA_GEMSTAT_Ensemble_Parlist[[7]]$qbtm[aa_perf_srt22],
     main = "qbtm_sorted_by_perf_exp5_hER", pch = 16, ylab = "")
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[7]]$beta[aa_perf_srt22],
     main = "beta_sorted_by_perf_exp5_hER", pch = 16)

heatmap.2(E_RNA_GEMSTAT_Ensemble_Parlist[[7]]$logistic[aa_perf_srt22,],
          Rowv = F, Colv = T,
          dendrogram = 'column', trace = 'none', margins = c(10,10))
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[7]]$logistic[aa_perf_srt22, 1], 
     ylab="", main = "logistic_bias_sorted_by_perf_exp5_hER")
plot(E_RNA_GEMSTAT_Ensemble_Parlist[[7]]$logistic[aa_perf_srt22, 2],
     ylab="", main = "logistic_coeff_sorted_by_perf_exp5_hER")

# plot the best models
which.max(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$test_ROC)
which.max(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$test_PRC)
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$train_results[[3835]]$ROC_curve, main = "ROC train bestModel exp5_hER")
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$train_results[[3835]]$PRC_curve, main = "PRC train bestModel exp5_hER")
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$test_results[[3835]]$ROC_curve, main = "ROC test bestModel exp5_hER")
plot(E_RNA_GEMSTAT_Ensemble_Outlist[[7]]$test_results[[3835]]$PRC_curve, main = "PRC test bestModel exp5_hER")

# aava <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Out_hER_valid/", full.names = F)
# setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Out_hER_valid/")
# aava_sp <- strsplit(aava, split="_")
# for(i in 1:length(aava_sp)){
#   aava_sp[[i]][6] <- "valid.txt"
# }
# aava_sp2 <- unlist(lapply(aava_sp, paste, collapse = "_"))
# for(i in 1:length(aava)){
#   system(paste0("mv ", aava[i], " ", aava_sp2[i]))
# }
# perform coop for experiment 5_hER (7)

aa_names <- c("ESR1_2", "ESR1_3","FOXA1","JUN_1","LEF1","NFIB","NKX3_1","NR2F2","NR3C1","NR5A2",
              "PBX1","PGR","PPARD","PPARG","RARA","RUNX1","SP1","YBX1")
aacoop <- rbind(c(1, 1), c(1, 2), c(1, 3), c(1, 8),
                c(1, 12), c(1,13), c(1, 14), c(1, 18),
                c(8, 8), c(8, 16))
aaKDexp <- list()
for(i in 1:length(aa_names)){
  aaKDexp[[i]] <- list()
  aaKDexp[[i]]$TF <- aa_names[i]
  aaKDexp[[i]]$COOP <- numeric(0)
}
for(i in 1:nrow(aacoop)){
  aaKDexp[[length(aa_names) + i]] <- list()
  aaKDexp[[length(aa_names) + i]]$TF <- character(0)
  aaKDexp[[length(aa_names) + i]]$COOP <- aacoop[i, ]
}
names(aaKDexp) <- c(1:length(aaKDexp))

aas1 <- sort(aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$valid_PRC, decreasing = T, index.return=T)$ix
aaadd1 <- rownames(aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$annot)[aas1[1:20]]
aa_address1 <- paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Trained_par_hER/", aaadd1, ".txt.Filter")

aasum <- aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$valid_ROC + aaE_RNA_GEMSTAT_Ensemble_Outlist_7_filtered$valid_PRC
aas2 <- sort(aasum, decreasing = T, index.return=T)$ix
aaadd2 <- rownames(aaE_RNA_GEMSTAT_Ensemble_Parlist_7_filtered$annot)[aas2[1:20]]
aa_address2 <- paste0("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/Trained_par_hER/", aaadd2, ".txt.Filter")
aa_address <- unique(c(aa_address1, aa_address2))

aadest_directory <- "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5"
GEMSTAT_create_KD_paramFile_ensemble(param_addresses = aa_address, 
                                     KD_exp_list = aaKDexp, 
                                     TF_names=aa_names,
                                     dest_directory = aadest_directory)
# write hal job for KDs
aa_par_trained_files <- list.files(path = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/KnockDown",
                                   recursive = T)
aa_par_trained_files_trimmed_1 <- unlist(lapply(strsplit(aa_par_trained_files, "\\/"), "[[", 1))
aa_par_trained_files_trimmed_2 <- unlist(lapply(strsplit(unlist(lapply(strsplit(aa_par_trained_files, "\\/"), "[[", 3)), split = "\\."), "[[", 1))

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5")
for(i in 1:length(aa_par_trained_files)){
  cat(c("/shared-mounts/sinhas-storage1/tabebor2/ER_Project/GEMSTAT/src/seq2expr -s Sequence/train_seq.fa",
        "-e Labels/Label_train.txt",
        "-m motifs.wtmx -f TF_exp.tab", 
        paste0("-fo KnockDown/", aa_par_trained_files_trimmed_1[i], "/Out/output_", aa_par_trained_files_trimmed_2[i], "_KD.txt"), 
        "-o DIRECT -c Coop/coop.par ",
        paste0("-p KnockDown/", aa_par_trained_files[i]),
        "-lower_bound bounds/lower.par -upper_bound bounds/upper.par -ff free_fix/ff.par -classifier -na 0 -onebeta 1",
        "-train_weights weight_train.txt\n"),
      sep = " ", append = T, file = "Ensemble_Exp5_hER_job_KD.sh")
}

# read KD results
aa_names <- c("ESR1_2", "ESR1_3","FOXA1","JUN_1","LEF1","NFIB","NKX3_1","NR2F2","NR3C1","NR5A2",
              "PBX1","PGR","PPARD","PPARG","RARA","RUNX1","SP1","YBX1")
aacoop <- rbind(c(1, 1), c(1, 2), c(1, 3), c(1, 8),
                c(1, 12), c(1,13), c(1, 14), c(1, 18),
                c(8, 8), c(8, 16))
aacoop_names <- character(length = nrow(aacoop))
for(i in 1:nrow(aacoop)){
  aacoop_names[i] <- paste0(aa_names[aacoop[i, 1]],"__" ,aa_names[aacoop[i, 2]])
}

E_RNA_GEMSTAT_Ensemble_KD_list <- list()
E_RNA_GEMSTAT_Ensemble_KD_list[[1]] <- GEMSTAT_read_KD(KD_directory = "~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/KnockDown",
                                                       kd_exp_names = as.character(c(1:28)),
                                                       WT_output = E_RNA_GEMSTAT_Ensemble_Outlist[[7]])
for(i in 1:length( E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)){
  colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]])[1:2] <- c("GT", "WT")
  colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]])[3:20] <- aa_names
  colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]])[21:30] <- aacoop_names
}
aarowside <- E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list$output_par_420_1[, 1]
aarowside[aarowside==1] <- "yellow"
aarowside[aarowside==0] <- "red"
heatmap.2(x = E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list$output_par_420_1[,1:30], Rowv = T, 
          Colv = T, dendrogram = "both", trace = "none", margins = c(8,5)
          , breaks = c(0, 0.05, 0.1, 0.2, 0.30, 0.50, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5)
          #,RowSideColors = aarowside
          )

setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/E_RNA/GEMSTAT_eRNA_logistic/Ensemble/Experiment_5/KD_plots/")
for(j in 1:length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)){
  jpeg(filename = paste0(names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)[j], "_KD.jpeg"), 
       width = 20, height = 17, units = "cm", res = 300)
  par(mfrow = c(5, 6), mar = c(1,1,1,1))
  for(i in 1:ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[j]])){
    plot( E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[j]][, 2], 
          E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[j]][, i], 
          ylab="KD", xlab = "WT", xaxt = "n", yaxt = "n",
          main = colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[j]])[i])
    
  }
  dev.off()
}
setwd("~/Documents/Shayan/BioInf/EstrogenReceptor/")

aa_wilcox <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list), 
                    ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]]) - 2))
rownames(aa_wilcox) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)
colnames(aa_wilcox) <- colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]])[3:30]
aa_wilcox_list <- list()
for(i in 1:length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)){
  aa_wilcox_list[[i]] <- list()
  for(j in 3:ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]])){
    aa_wilcox_list[[i]][[j-2]] <-  wilcox.test(x = E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, 2],
                                             y = E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, j],
                                             paired = T)
    aa_wilcox[i, j-2] <- aa_wilcox_list[[i]][[j-2]]$p.value
  }
  names(aa_wilcox_list[[i]]) <- colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]])[3:30]
}
names(aa_wilcox_list) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)

aa_wilcoxmlog10 <- -log10(aa_wilcox)
aa <- order(aa_wilcoxmlog10[1,], decreasing = T)
colnames(aa_wilcoxmlog10)[aa]

aa_wilcox_ex <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list), 
                    ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]]) - 2))
rownames(aa_wilcox_ex) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)
colnames(aa_wilcox_ex) <- colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]])[3:30]
aa_wilcox_list_ex <- list()
for(i in 1:length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)){
  aa_wilcox_list_ex[[i]] <- list()
  for(j in 3:ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]])){
    aa_wilcox_list_ex[[i]][[j-2]] <-  wilcox.test(x = E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, 2],
                                               y = E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, j],
                                               paired = T, alternative="two.sided", conf.int = T)
    aa_wilcox_ex[i, j-2] <- aa_wilcox_list_ex[[i]][[j-2]]$p.value
  }
  names(aa_wilcox_list_ex[[i]]) <- colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]])[3:30]
}
names(aa_wilcox_list_ex) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)



aa_wilcox_list_ex$output_par_1033_1$LEF1$estimate
for(i in 1:length(aa_wilcox_list_ex$output_par_1033_1)){
  print(names(aa_wilcox_list_ex$output_par_1033_1)[i])
  print(aa_wilcox_list_ex$output_par_1033_1[[i]]$statistic)
  print(aa_wilcox_list_ex$output_par_1033_1[[i]]$p.value)
  print(aa_wilcox_list_ex$output_par_1033_1[[i]]$conf.int)
  print(aa_wilcox_list_ex$output_par_1033_1[[i]]$estimate)
}

aa_diff_meanabs <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list), 
                       ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]]) - 2))
rownames(aa_diff_meanabs) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)
colnames(aa_diff_meanabs) <- colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]])[3:30]

aa_diff_maxabs <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list), 
                     ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]]) - 2))
rownames(aa_diff_maxabs) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)
colnames(aa_diff_maxabs) <- colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]])[3:30]

aa_diff_meanabsgt1 <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list), 
                         ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]]) - 2))
rownames(aa_diff_meanabsgt1) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)
colnames(aa_diff_meanabsgt1) <- colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]])[3:30]

aa_diff_nuabsgt1 <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list), 
                             ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]]) - 2))
rownames(aa_diff_nuabsgt1) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)
colnames(aa_diff_nuabsgt1) <- colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]])[3:30]

aa_diff_list_ex <- list()
for(i in 1:length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)){
  #print(sum(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, 2] == 0))
  aa_diff_list_ex[[i]] <- (E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]] - E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, 2])/E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, 2]
  
  for(j in 3:ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]])){
    
    aa_diff_meanabs[i, j-2] <- mean(abs(aa_diff_list_ex[[i]][,j]))
    aa_diff_maxabs[i, j-2] <- max(abs(aa_diff_list_ex[[i]][,j]))
    aa_diff_meanabsgt1[i, j-2] <- mean(abs(aa_diff_list_ex[[i]][abs(aa_diff_list_ex[[i]][, j]) > 0.05 ,j]))
    aa_diff_nuabsgt1[i, j-2] <- sum(abs(aa_diff_list_ex[[i]][, j]) > 0.05)/nrow(aa_diff_list_ex[[i]])
  }
  #names(aa_diff_list_ex[[i]]) <- colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]])[3:30]
}
names(aa_diff_list_ex) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)

par(mfrow = c(1,1), mar = c(8,4,2,2))
boxplot.matrix(aa_diff_meanabs, las = 2, outline= F, xaxt = 's', main = "mean absolute fold change")
boxplot.matrix(aa_diff_nuabsgt1, las = 2, outline= F, xaxt = 's', main = "fraction of enhancers affected > 5%")
boxplot.matrix(aa_diff_meanabsgt1, las = 2, outline= F, xaxt = 's', main = "mean absolute fold change among altered > 5%")
boxplot.matrix(aa_diff_maxabs, las = 2, outline= F, main= "max absolute fold change")



aa_diff_meanabs_order <- aa_diff_meanabs
aa_diff_meanabsgt1_order <- aa_diff_meanabsgt1
aa_diff_maxabs_order <- aa_diff_maxabs
aa_diff_nuabsgt1_order <- aa_diff_nuabsgt1
for(i in 1:nrow(aa_diff_meanabs)){
  aa_diff_meanabs_order[i, ] <- rev(order(aa_diff_meanabs[i,]))
  aa_diff_meanabsgt1_order[i, ] <- rev(order(aa_diff_meanabsgt1[i, ]))
  aa_diff_maxabs_order[i, ] <- rev(order(aa_diff_maxabs[i, ]))
  aa_diff_nuabsgt1_order[i, ] <- rev(order(aa_diff_nuabsgt1[i, ]))
}

par(mfrow = c(1,1), mar = c(8,4,2,2))
boxplot.matrix(aa_diff_meanabs_order, las = 2, outline= T, xaxt = 's', main="rank of mean absolute fold change")
boxplot.matrix(aa_diff_nuabsgt1_order, las = 2, outline= T, xaxt = 's',main = "rank of #enhancers affected > 5%")
boxplot.matrix(aa_diff_meanabsgt1_order, las = 2, outline= T, xaxt = 's', main = "rank of mean absolute fold change among altered > 5%")
boxplot.matrix(aa_diff_maxabs_order, las = 2, outline= T, main= "rank of max absolute fold change")



aa_diff_meanabs_sep <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list), 
                          ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]]) - 2) * 2)
rownames(aa_diff_meanabs_sep) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)
aa <- (cbind(paste(colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]])[3:30], "pos", sep = "_"),
             paste(colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]])[3:30], "neg", sep = "_")))
colnames(aa_diff_meanabs_sep) <- matrix(t(aa), ncol = 1, byrow = T)[,1]

aa_diff_maxabs_sep <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list), 
                         ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]]) - 2)* 2)
rownames(aa_diff_maxabs_sep) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)
colnames(aa_diff_maxabs_sep) <- matrix(t(aa), ncol = 1, byrow = T)[,1]

aa_diff_meanabsgt1_sep <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list), 
                             ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]]) - 2)*2)
rownames(aa_diff_meanabsgt1_sep) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)
colnames(aa_diff_meanabsgt1_sep) <- matrix(t(aa), ncol = 1, byrow = T)[,1]

aa_diff_nuabsgt1_sep <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list), 
                           ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]]) - 2)*2)
rownames(aa_diff_nuabsgt1_sep) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)
colnames(aa_diff_nuabsgt1_sep) <- matrix(t(aa), ncol = 1, byrow = T)[,1]

aa_diff_list_ex <- list()
for(i in 1:length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)){
  #print(sum(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, 2] == 0))
  aa_diff_list_ex[[i]] <- (E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]] - E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, 2])/E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, 2]
  
  for(j in 3:ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]])){
    aaone <- E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, 1] == 1
    aazero <- E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, 1] == 0
    aa_diff_meanabs_sep[i, (2*(j-2) - 1)] <- mean(abs(aa_diff_list_ex[[i]][aaone,j]))
    aa_diff_maxabs_sep[i, (2*(j-2) - 1)] <- max(abs(aa_diff_list_ex[[i]][aaone,j]))
    aa_diff_meanabsgt1_sep[i, (2*(j-2) - 1)] <- mean(abs(aa_diff_list_ex[[i]][(abs(aa_diff_list_ex[[i]][, j]) > 0.05) & aaone ,j]))
    aa_diff_nuabsgt1_sep[i, (2*(j-2) - 1)] <- sum((abs(aa_diff_list_ex[[i]][, j]) > 0.05) & aaone)/sum(aaone)
    aa_diff_meanabs_sep[i, (2*(j-2) )] <- mean(abs(aa_diff_list_ex[[i]][aazero,j]))
    aa_diff_maxabs_sep[i, (2*(j-2))] <- max(abs(aa_diff_list_ex[[i]][aazero,j]))
    aa_diff_meanabsgt1_sep[i, (2*(j-2))] <- mean(abs(aa_diff_list_ex[[i]][(abs(aa_diff_list_ex[[i]][, j]) > 0.05) & aazero ,j]))
    aa_diff_nuabsgt1_sep[i, (2*(j-2))] <- sum((abs(aa_diff_list_ex[[i]][, j]) > 0.05) & aazero)/sum(aazero)
  }
  #names(aa_diff_list_ex[[i]]) <- colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]])[3:30]
}
names(aa_diff_list_ex) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)

par(mfrow = c(1,1), mar = c(8,4,2,2))
boxplot.matrix(aa_diff_meanabs_sep, las = 2, outline= F, xaxt = 's', main = "mean absolute fold change", col = rep(c("red", "blue"), 28))
boxplot.matrix(aa_diff_nuabsgt1_sep, las = 2, outline= F, xaxt = 's', main = "fraction of enhancers affected > 5%", col = rep(c("red", "blue"), 28))
boxplot.matrix(aa_diff_meanabsgt1_sep, las = 2, outline= F, xaxt = 's', main = "mean absolute fold change among altered > 5%", col = rep(c("red", "blue"), 28))
boxplot.matrix(aa_diff_maxabs_sep, las = 2, outline= F, main= "max absolute fold change", col = rep(c("red", "blue"), 28))


aa_diff_mean_sep <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list), 
                              ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]]) - 2) * 2)
rownames(aa_diff_mean_sep) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)
aa <- (cbind(paste(colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]])[3:30], "pos", sep = "_"),
             paste(colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]])[3:30], "neg", sep = "_")))
colnames(aa_diff_mean_sep) <- matrix(t(aa), ncol = 1, byrow = T)[,1]

aa_diff_max_sep <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list), 
                             ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]]) - 2)* 2)
rownames(aa_diff_max_sep) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)
colnames(aa_diff_max_sep) <- matrix(t(aa), ncol = 1, byrow = T)[,1]

aa_diff_meangt1_sep <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list), 
                                 ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]]) - 2)*2)
rownames(aa_diff_meangt1_sep) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)
colnames(aa_diff_meangt1_sep) <- matrix(t(aa), ncol = 1, byrow = T)[,1]

aa_diff_nugt1_sep <- matrix(nrow = length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list), 
                               ncol = (ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[1]]) - 2)*2)
rownames(aa_diff_nugt1_sep) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)
colnames(aa_diff_nugt1_sep) <- matrix(t(aa), ncol = 1, byrow = T)[,1]

aa_diff_list_ex <- list()
for(i in 1:length(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)){
  #print(sum(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, 2] == 0))
  aa_diff_list_ex[[i]] <- (E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]] - E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, 2])/E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, 2]
  
  for(j in 3:ncol(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]])){
    aaone <- E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, 1] == 1
    aazero <- E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]][, 1] == 0
    aa_diff_mean_sep[i, (2*(j-2) - 1)] <- mean((aa_diff_list_ex[[i]][aaone,j]))
    aavv <- aa_diff_list_ex[[i]][aaone,j]
    aa_diff_max_sep[i, (2*(j-2) - 1)] <- max(abs(aavv)) * sign(aavv[which.max(abs(aavv))])
    aa_diff_meangt1_sep[i, (2*(j-2) - 1)] <- mean(aa_diff_list_ex[[i]][(abs(aa_diff_list_ex[[i]][, j]) > 0.05) & aaone ,j])
    aa_diff_nugt1_sep[i, (2*(j-2) - 1)] <- sum((abs(aa_diff_list_ex[[i]][, j]) > 0.05) & aaone)/sum(aaone)
    aa_diff_mean_sep[i, (2*(j-2) )] <- mean((aa_diff_list_ex[[i]][aazero,j]))
    aavv <- aa_diff_list_ex[[i]][aazero,j]
    aa_diff_max_sep[i, (2*(j-2))] <- max(abs(aavv)) * sign(aavv[which.max(abs(aavv))])
    aa_diff_meangt1_sep[i, (2*(j-2))] <- mean(aa_diff_list_ex[[i]][(abs(aa_diff_list_ex[[i]][, j]) > 0.05) & aazero ,j])
    aa_diff_nugt1_sep[i, (2*(j-2))] <- sum((abs(aa_diff_list_ex[[i]][, j]) > 0.05) & aazero)/sum(aazero)
  }
  #names(aa_diff_list_ex[[i]]) <- colnames(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list[[i]])[3:30]
}
names(aa_diff_list_ex) <- names(E_RNA_GEMSTAT_Ensemble_KD_list[[1]]$output_mat_list)

par(mfrow = c(1,1), mar = c(8,4,2,2))
boxplot.matrix(aa_diff_mean_sep, las = 2, outline= F, xaxt = 's',
               main = "mean fold change", col = rep(c("red", "blue"), 28))
abline(h = seq(-1, 1, 0.05), col = 3, lty = 3, lwd = 0.8)
abline(v = seq(0.5, 60.5, 2), col = 3,lty = 3, lwd = 0.8 )
boxplot.matrix(aa_diff_nugt1_sep, las = 2, outline= F, xaxt = 's', 
               main = "fraction of enhancers affected > 5%", col = rep(c("red", "blue"), 28))
abline(h = seq(-1, 1, 0.05), col = 3, lty = 3, lwd = 0.8)
abline(v = seq(0.5, 60.5, 2), col = 3,lty = 3, lwd = 0.8 )

boxplot.matrix(aa_diff_meangt1_sep, las = 2, outline= F,
               xaxt = 's', main = "mean  fold change among altered > 5%",
               col = rep(c("red", "blue"), 28)
               ,ylim = c(-0.5, 0.5)
               )
abline(h = seq(-5, 5, 0.05), col = 3, lty = 3, lwd = 0.8)
abline(v = seq(0.5, 60.5, 2), col = 3,lty = 3, lwd = 0.8 )

boxplot.matrix(aa_diff_max_sep, las = 2, outline= F,
               main= "max absolute fold change", col = rep(c("red", "blue"), 28)
               ,ylim = c(-1, 2)
               )
abline(h = seq(-10, 10, 0.05), col = 3, lty = 3, lwd = 0.8)
abline(v = seq(0.5, 60.5, 2), col = 3,lty = 3, lwd = 0.8 )














