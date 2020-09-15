# This program builds a variant of the XRF.R code from causalToolbox package (https://github.com/soerenkuenzel/causalToolbox/blob/master/R/XRF.R)
# using functions or adapted functions from the grf package.

library('dplyr')
library('readr')
library('grf')
library('tidyr')
library('magrittr')
library('ggplot2')




default_sample_fraction <- 0.5
default_n_trees <- 5000
to_cluster <- F
drop_duplicates <- T
set.seed(392809)

filename <- "grf/final_forests_missingness"

if (to_cluster==F){
  filename <- filename %>% paste0("_no_cluster")
}

if (drop_duplicates==T){
  filename <- filename %>% paste0("_no_dupes")
}



# load data ----
setwd("/export/projects/migrated-BAM/Match_AnalysisFiles_May2016/Results/2020")
source('grf/grf_globals.R')
master_pool <- load_master_dataset()



# helper functions ----

impute_block_baselines <- function(input_dataset,
                                   list_of_controls,
                                   type = 'median'){
  
  input_dataset <- input_dataset %>% mutate_at(controls, as.numeric)
  
  for (block in unique(input_dataset$blocknum)){
    # for each variable...
    for (baseline_var in list_of_controls){
      # get block statistic for this variable
      if (type %in% 'median'){
        temp_block_stat <- master_pool %>% filter(blocknum %in% block) %>% pull(baseline_var) %>% median(na.rm=T)
        
      }
      if (type %in% 'mean'){
        temp_block_stat <- master_pool %>% filter(blocknum %in% block) %>% pull(baseline_var) %>% mean(na.rm=T)
      }  
      
      # find all ppl in block missing this variable, replace value
      input_dataset[which(input_dataset$blocknum %in% block &
                            is.na(input_dataset[,baseline_var])), baseline_var] <- temp_block_stat
    }
  }
  
  return(input_dataset)
  
}


makeDummies <- function(dataframe,
                        flipped=F) {
  dataframe$d_positive <- ifelse(dataframe$grf_tau_hat > 0, 1, 0)
  dataframe$d_negative <- ifelse(dataframe$grf_tau_hat <= 0, 1, 0)
  
  dataframe$d_positive <- ifelse(dataframe$grf_tau_hat > 0, 1, 0)
  dataframe$d_negative <- ifelse(dataframe$grf_tau_hat <= 0, 1, 0)
  
  #qcut function taken from Hadley Wickham answer on stackoverflow
  qcut <- function(x, n) {
    cut(x, quantile(x, seq(0, 1, length = n + 1)), labels = seq_len(n),
        include.lowest = TRUE)
  }
  
  dataframe$tau_quartile <- qcut(dataframe$grf_tau_hat, 4) %>% as.numeric()
  
  if (flipped ==T){
    dataframe$tau_quartile <- abs(4 - dataframe$tau_quartile) + 1
  }
  
  dummy_df <- data.frame(dummies::dummy(dataframe$tau_quartile, sep="_"))
  
  #note: dummy function has issues with column names; overwrite them manually here
  colnames(dummy_df) <- c("d_quart1", "d_quart2", "d_quart3", "d_quart4")
  dataframe <- cbind(dataframe, dummy_df)
  return(dataframe)
  
}



scalar1 <- function(x) {x / sqrt(sum(x^2, na.rm=T))}

add_pairwise_interactions <- function(data, list_of_controls){
  
  temp_df <- data
  
  for (control_1 in list_of_controls){
    for (control_2 in list_of_controls){
      # if control1==control2, pass
      if (control_1 == control_2) next
      temp_df[,paste0(control_1, "_x_", control_2)] <- scalar1(temp_df[,control_1]) * scalar1(temp_df[,control_2])
      
    }
  }
  
  return(temp_df)
  
}



avg_effect <- function(Y,Y.hat,W,W.hat, weights,predictions,subset) {
  
  if (is.null(subset)) {
    subset <- 1:length(forest$Y.hat)
  }
  if (class(subset) == "logical" & length(subset) == length(forest$Y.hat)) {
    subset <- which(subset)
  }
  
  subset.W.orig <- W[subset]
  subset.W.hat <- W.hat[subset]
  subset.Y.orig <- Y[subset]
  subset.Y.hat <- Y.hat[subset]
  tau.hat.pointwise <- predictions[subset]
  subset.weights <- weights[subset]
  
  control.idx <- which(subset.W.orig == 0)
  treated.idx <- which(subset.W.orig == 1)
  
  tau.hat.pointwise <- predictions[subset]
  tau.avg.raw <- weighted.mean(tau.hat.pointwise, subset.weights)
  
  Y.hat.0 <- subset.Y.hat - subset.W.hat * tau.hat.pointwise
  Y.hat.1 <- subset.Y.hat + (1 - subset.W.hat) * tau.hat.pointwise
  
  gamma.control.raw <- 1/(1 - subset.W.hat[control.idx])
  gamma.treated.raw <- 1/subset.W.hat[treated.idx]
  gamma <- rep(0, length(subset.W.orig))
  gamma[control.idx] <- gamma.control.raw/sum(subset.weights[control.idx] * 
                                                gamma.control.raw) * sum(subset.weights)
  gamma[treated.idx] <- gamma.treated.raw/sum(subset.weights[treated.idx] * 
                                                gamma.treated.raw) * sum(subset.weights)  
  dr.correction.all <- subset.W.orig * gamma * (subset.Y.orig - 
                                                  Y.hat.1) - (1 - subset.W.orig) * gamma * (subset.Y.orig - 
                                                                                              Y.hat.0)
  dr.correction <- weighted.mean(dr.correction.all, subset.weights)
  
  tau.avg <- tau.avg.raw + dr.correction
  
  return(c(estimate = tau.avg))
}


make_single_X_RF <- function(dataset,
                             controls,
                             outcome,
                             outcome_label,
                             flip_outcome=FALSE,
                             cluster=to_cluster,
                             add_all_interactions = F,
                             remove_duplicates=drop_duplicates,
                             tune_num_reps = 1000,
                             tune_num_draws = 2000,
                             tune_num_trees = 200){
  
  
  # Data Cleaning from build_grf_models
  
  if (remove_duplicates==T){
    # get duplicate sids
    duplicate_sids <- dataset %>% 
      group_by(sid) %>% filter(n()>1) %>% 
      distinct(sid) %>% pull(sid)
    
    dataset <- dataset %>% filter(!sid %in% duplicate_sids)
  }
  
  
  # keep observations complete on outcome and that has an inv prob weight
  working_df <- dataset %>% drop_na(outcome, inv_prob_weight)
  
  n.obs <- nrow(working_df)
  
  X <- working_df[,controls]
  W <- working_df %>% pull(dmatch)
  Y <- working_df %>% pull(outcome)
  sample_weights <- working_df %>% pull(inv_prob_weight)
  cluster_ids <- working_df %>% pull(sid) # we have pooled the data, and with duplicates, we should cluster on sid
  
  
  # X-Learner (Based on causalToolbox package + weights)
  
  
  
  # First stage --------------------------------------------------------------
  Y_0 <- Y[W == 0]
  Y_1 <- Y[W == 1]
  
  X_0 <- X[W == 0, ]
  X_1 <- X[W == 1, ]
  
  sample_weights_0 <- sample_weights[W == 0]
  sample_weights_1 <- sample_weights[W == 1]
  
  m_0 <- regression_forest(X_0, Y_0, 
                           tune.parameters = c("sample.fraction", "mtry",
                                               "honesty.fraction", "honesty.prune.leaves",
                                               "alpha", "imbalance.penalty"),
                           min.node.size = 5,
                           sample.weights = sample_weights_0#,
                           #clusters = cluster_ids,
                           # sample.fraction = default_sample_fraction,
                           # num.trees = default_n_trees,
                           # tune.num.trees = tune_num_trees,
                           # tune.num.draws = tune_num_draws,
                           # tune.num.reps = tune_num_reps
                           )
  
  
  m_1 <- regression_forest(X_1, Y_1, 
                           tune.parameters = c("sample.fraction", "mtry",
                                               "honesty.fraction", "honesty.prune.leaves",
                                               "alpha", "imbalance.penalty"),
                           min.node.size = 5,
                           sample.weights = sample_weights_1#,
                           #clusters = cluster_ids,
                           # sample.fraction = default_sample_fraction,
                           # num.trees = default_n_trees,
                           # tune.num.trees = tune_num_trees,
                           # tune.num.draws = tune_num_draws,
                           # tune.num.reps = tune_num_reps
                           )
  
  
  Y.hat <- Y*NA
  Y.hat[W == 0] <- predict(m_0) %>% pull(predictions)
  Y.hat[W == 1] <- predict(m_1) %>% pull(predictions)
  
  # Second Stage -------------------------------------------------------------
  r_0 <-  predict(m_1, X_0) %>% pull(predictions) - Y_0
  r_1 <- Y_1 - (predict(m_0, X_1) %>% pull(predictions))
  
  m_tau_0 <- regression_forest(X_0, r_0, 
                               tune.parameters = c("sample.fraction", "mtry",
                                                   "honesty.fraction", "honesty.prune.leaves",
                                                   "alpha", "imbalance.penalty"),
                               min.node.size = 5,
                               sample.weights = sample_weights_0#,
                               #clusters = cluster_ids,
                               # sample.fraction = default_sample_fraction,
                               # num.trees = default_n_trees,
                               # tune.num.trees = tune_num_trees,
                               # tune.num.draws = tune_num_draws,
                               # tune.num.reps = tune_num_reps
                               )
  
  
  m_tau_1 <- regression_forest(X_1, r_1, 
                               tune.parameters = c("sample.fraction", "mtry",
                                                   "honesty.fraction", "honesty.prune.leaves",
                                                   "alpha", "imbalance.penalty"),
                               min.node.size = 5,
                               sample.weights = sample_weights_1#,
                               #clusters = cluster_ids,
                               # sample.fraction = default_sample_fraction,
                               # num.trees = default_n_trees,
                               # tune.num.trees = tune_num_trees,
                               # tune.num.draws = tune_num_draws,
                               # tune.num.reps = tune_num_reps
                               )
  
  
  
  # Prop score estimation ----------------------------------------------------
  m_prop <-
    regression_forest(X, W, 
                      tune.parameters = c("sample.fraction", "mtry",
                                          "honesty.fraction", "honesty.prune.leaves",
                                          "alpha", "imbalance.penalty"),
                      min.node.size = 5,
                      sample.weights = sample_weights#,
                      #clusters = cluster_ids,
                      # sample.fraction = default_sample_fraction,
                      # num.trees = default_n_trees,
                      # tune.num.trees = tune_num_trees,
                      # tune.num.draws = tune_num_draws,
                      # tune.num.reps = tune_num_reps
                      ) 
  
  prop_scores <- predict(m_prop) %>% pull(predictions)
  
  predictions_oob <- prop_scores * (predict(m_tau_0,X) %>% pull(predictions)) +
    (1 - prop_scores)  * (predict(m_tau_1,X) %>% pull(predictions))
  
  
  
  # Calibration Plot --------------------------------------------
  mean.pred <- weighted.mean(predictions_oob, sample_weights)
  DF <- data.frame(target = unname(Y - Y.hat), 
                   mean.forest.prediction = unname(W - prop_scores) * 
                     mean.pred,
                   differential.forest.prediction = unname(W - prop_scores) * (predictions_oob - mean.pred))
  
  # changed "observation.weight" to "sample_weights"
  best.linear.predictor <- lm(target ~ mean.forest.prediction + 
                                differential.forest.prediction + 0, weights = sample_weights, 
                              data = DF)
  # blp.summary <- lmtest::coeftest(best.linear.predictor, vcov = sandwich::vcovCL, 
  #                                 type = "HC3", cluster = clusters)
  # changed this because we are not clustering (can change type of HC errors )
  blp.summary <- lmtest::coeftest(best.linear.predictor,
                                  vcov=sandwich::vcovHC(best.linear.predictor, type="HC1"))
  
  attr(blp.summary, "method") <- paste("Best linear fit using forest predictions (on held-out data)", 
                                       "as well as the mean forest prediction as regressors, along", 
                                       "with one-sided heteroskedasticity-robust (HC3) SEs", 
                                       sep = "\n")
  dimnames(blp.summary)[[2]][4] <- gsub("[|]", "", dimnames(blp.summary)[[2]][4])
  blp.summary[, 4] <- ifelse(blp.summary[, 3] < 0, 1 - blp.summary[, 
                                                                   4]/2, blp.summary[, 4]/2)
  forest_calibration <- blp.summary %>% broom::tidy()
  
  # output should be a list with
  # - avg causal effects
  group_tau_avgs <- c('overall'= mean.pred#avg_tx_effect_overall,
                      # 'treated'=avg_tx_effect_treated,
                      # 'control'=avg_tx_effect_control
                      )
  
  # - the predictions + variance in a df with the sids
  if (flip_outcome == FALSE) {
    tau_df <- tibble(sid = cluster_ids,
                     grf_tau_hat = predictions_oob) %>% 
      makeDummies()
  }
  
  if (flip_outcome == TRUE) {
    tau_df <- tibble(sid = cluster_ids,
                     grf_tau_hat = predictions_oob) %>% 
      makeDummies(flipped = T)
  }  
  
  # now that we have the quartile dummies, we can get quartile level effects
  ate.highest_25 <- avg_effect(Y = Y,Y.hat = Y.hat ,W = W, W.hat = prop_scores, weights = sample_weights, predictions = predictions_oob,
                               subset = tau_df$tau_quartile == 4)
  
  ate.bottom_25 <- avg_effect(Y = Y,Y.hat = Y.hat ,W = W, W.hat = prop_scores, weights = sample_weights, predictions = predictions_oob,
                              subset = tau_df$tau_quartile == 1)
  ate.quartile_2 <- avg_effect(Y = Y,Y.hat = Y.hat ,W = W, W.hat = prop_scores, weights = sample_weights, predictions = predictions_oob,
                               subset = tau_df$tau_quartile == 2)
  ate.quartile_3 <- avg_effect(Y = Y,Y.hat = Y.hat ,W = W, W.hat = prop_scores, weights = sample_weights, predictions = predictions_oob,
                               subset = tau_df$tau_quartile == 3)
  
  ate.bottom_75 <- avg_effect(Y = Y,Y.hat = Y.hat ,W = W, W.hat = prop_scores, weights = sample_weights, predictions = predictions_oob,
                              subset = tau_df$tau_quartile != 4)
  ate.high <- avg_effect(Y = Y,Y.hat = Y.hat ,W = W, W.hat = prop_scores, weights = sample_weights, predictions = predictions_oob,
                         subset = tau_df$tau_quartile %in% c(3,4))
  
  
  ate.low <- avg_effect(Y = Y,Y.hat = Y.hat ,W = W, W.hat = prop_scores, weights = sample_weights, predictions = predictions_oob,
                        subset = tau_df$tau_quartile %in% c(1,2))
  
  subsample_tau_avgs <- list('highest_quartile' = ate.highest_25,
                             'bottom_three_quartiles' = ate.bottom_75,
                             'above_median'=ate.high,
                             'below_median'=ate.low,
                             'quartile_1' = ate.bottom_25,
                             'quartile_2' = ate.quartile_2,
                             'quartile_3' = ate.quartile_3,
                             'quartile_4' = ate.highest_25)
  
  # - 'augmented' master df (the above df merged with the full master dataset)
  augmented_df <- working_df %>% left_join(tau_df, by='sid')
  
  # - quartile heterogeneity table
  # - plot: prediction vs rank plot
  
  
  output_list <- list('tau_df' = tau_df,
                      'augmented_df' = augmented_df,
                      #'quartile_heterogeneity_table' = quartile_heterogeneity_table,
                      #'tau_rank_plot'=tau_rank_plot,
                      'group_tau_avgs'=group_tau_avgs,
                      'subsample_tau_avgs'=subsample_tau_avgs,
                      #'forest_object'=tau.forest,
                      'calibration_test'=forest_calibration,
                      'n_observations' = n.obs,
                      'tuning_output' = tau.forest$tuning.output)
  
  return(output_list)  
}


block_treatment_probability_df <- master_pool %>% group_by(blocknum, study) %>% summarize(p_block = mean(dmatch)) %>% ungroup()

singleton_blocks <- master_pool %>% group_by(blocknum, study) %>% filter(n()==1) %>% ungroup() %>% pull(blocknum)

master_pool <- master_pool %>%
  left_join(block_treatment_probability_df, by=c('study', 'blocknum')) %>% 
  mutate(inv_prob_weight = dmatch/p_block + (1-dmatch)/(1-p_block)) %>% 
  filter(! blocknum %in% singleton_blocks)

# W_i = D/p_block + (1-D)/(1-p_block)



make_single_X_RF(master_pool, controls_sans_missingness_dummies, outcomes_of_interest[1], outcomes_of_interest_labels[1])




make_xrf_for_each_outcome <- function(input_dataset,
                                      input_controls,
                                      outcomes,
                                      outcome_labels,
                                      flipped_outcomes=c()){
  final_forest_list <- list()
  
  for (i in 1:length(outcomes_of_interest)){
    
    current_outcome <- outcomes[i]
    current_outcome_label <- outcome_labels[i]
    cat("Running: ", current_outcome)
    if (current_outcome %in% flipped_outcomes){
      single_forest_output <- make_single_X_RF(dataset = input_dataset,
                                               controls = input_controls,
                                               outcome = current_outcome,
                                               outcome_label = current_outcome_label,
                                               flip_outcome = T)
    } else {
      single_forest_output <- make_single_X_RF(dataset = input_dataset,
                                               controls = input_controls,
                                               outcome = current_outcome,
                                               outcome_label = current_outcome_label)
    }
    
    
    
    final_forest_list[[current_outcome_label]] <- single_forest_output
    
  }
  
  return(final_forest_list)
}

# # function arguments
# input_dataset <- master_pool_imputed
# input_controls <- controls
# outcomes <- outcomes_of_interest
# outcome_labels <- outcomes_of_interest

# 
# final_forests_imputed <- make_forest_for_each_outcome(input_dataset = master_pool_imputed,
#                              controls,
#                              outcomes_of_interest,
#                              outcomes_of_interest_labels)


start_time <- Sys.time()
final_forests_xrf <- make_forest_for_each_outcome(master_pool,
                                                  controls_sans_missingness_dummies,
                                                  outcomes_of_interest,
                                                  outcomes_of_interest_labels,
                                                  flipped_outcomes)
end_time <- Sys.time()

# saveRDS(final_forests_imputed, file='grf/final_forests_imputed.rds')
saveRDS(final_forests_xrf, file=paste0(filename, ".rds"))





