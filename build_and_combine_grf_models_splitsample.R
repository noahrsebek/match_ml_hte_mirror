library('dplyr')
library('readr')
library('grf')
library('tidyr')
library('magrittr')
library('ggplot2')
library('furrr')
library('purrr')


n_splitsample_iterations <- 2
splitsample_fraction <- 0.5
default_n_trees <- 50
to_cluster <- F
drop_duplicates <- T



filename <- "grf/final_forests_missingness"

if (to_cluster==F){
  filename <- filename %>% paste0("_no_cluster")
}

if (drop_duplicates==T){
  filename <- filename %>% paste0("_no_dupes")
}



# helper functions ----
source(file = 'helper_functions.R')


# load data ----
setwd("/export/projects/migrated-BAM/Match_AnalysisFiles_May2016/Results/2020")
master_pool <- load_master_dataset()
setwd('grf')
source('grf_globals.R')

# main ML functions
make_single_causal_forest <- function(dataset,
                                      controls,
                                      outcome,
                                      outcome_label,
                                      flip_outcome=FALSE,
                                      cluster=to_cluster,
                                      add_all_interactions = F,
                                      remove_duplicates=drop_duplicates,
                                      tune_num_reps = 50,
                                      tune_num_draws = 1000,
                                      tune_num_trees = 200,
                                      #splitsample_frac = splitsample_fraction
                                      splitsample_df,
                                      tuned_parameters = c("sample.fraction", "mtry",
                                                           "honesty.fraction", "honesty.prune.leaves",
                                                           "alpha", "imbalance.penalty") 
                                      ){
  
  
  if (remove_duplicates==T){
    # get duplicate sids
    duplicate_sids <- dataset %>% 
      group_by(sid) %>% filter(n()>1) %>% 
      distinct(sid) %>% pull(sid)
    
    dataset <- dataset %>% filter(!sid %in% duplicate_sids)
  }
  
  
  # keep observations complete on outcome and that has an inv prob weight
  working_df <- dataset %>% drop_na(all_of(outcome)) #%>%  mutate(row_id = row_number())
  splitsample_df <- splitsample_df %>% drop_na(all_of(outcome))
  
  n.obs <- nrow(working_df)
  
  X <- working_df[,controls]
  W <- working_df %>% pull(dmatch)
  Y <- working_df %>% pull(tidyselect::all_of(outcome))
  sample_weights <- working_df %>% pull(inv_prob_weight)
  
  X_splitsample <- splitsample_df[,controls]
  Y_splitsample <- splitsample_df %>% pull(tidyselect::all_of(outcome))
  W_splitsample <- splitsample_df %>% pull(dmatch)
  sample_weights_splitsample <- splitsample_df %>% pull(inv_prob_weight)
  
  cluster_ids <- working_df %>% pull(sid) # we'll use these if we have duplicates n the data
  cluster_ids_splitsample <- splitsample_df %>% pull(sid)
  
  if (add_all_interactions==T){
    X <- X %>% add_pairwise_interactions(controls)
    X_splitsample <- X_splitsample %>% add_pairwise_interactions(controls)
  }
  
  if (cluster==T){
    # fit Y forest
    
    
    # fit W forest
    
    # fit tau forest
    tau.forest <- causal_forest(X, Y, W,
                                tune.parameters = c("sample.fraction", "mtry",
                                                    "honesty.fraction", "honesty.prune.leaves",
                                                    "alpha", "imbalance.penalty"),
                                min.node.size = 5,
                                sample.weights = sample_weights,
                                clusters = cluster_ids,
                                #sample.fraction = default_sample_fraction,
                                num.trees = default_n_trees,
                                tune.num.trees = tune_num_trees,
                                tune.num.draws = tune_num_draws,
                                tune.num.reps = tune_num_reps)
  }
  if (cluster==F){
    
    # fit W forest
    forest.W <- regression_forest(X, W, tune.parameters = tuned_parameters, min.node.size = 5)
    W.hat_splitsample <- predict(forest.W, newdata=X_splitsample)$predictions
    W.hat <- predict(forest.W)$predictions
    
    
    # fit Y forest
    forest.Y <- regression_forest(X, Y, tune.parameters = tuned_parameters, min.node.size = 5)
    Y.hat_splitsample <- predict(forest.Y, newdata=X_splitsample)$predictions
    Y.hat <- predict(forest.Y)$predictions
    
    
    tau.forest <- causal_forest(X, Y, W,
                                tune.parameters = tuned_parameters,
                                min.node.size = 5,
                                sample.weights = sample_weights,
                                W.hat = W.hat,
                                Y.hat = Y.hat,
                                #clusters = cluster_ids,
                                #sample.fraction = default_sample_fraction,
                                num.trees = default_n_trees,
                                tune.num.trees = tune_num_trees,
                                tune.num.draws = tune_num_draws,
                                tune.num.reps = tune_num_reps
                                )
  }

  
  tau.hat.oob <- predict(tau.forest,
                         estimate.variance = T)
  
  tau.hat.oob.splitsample <- predict(tau.forest,
                         newdata = X_splitsample,
                         estimate.variance = T)
  
  predictions_oob <- tau.hat.oob$predictions
  variance <- tau.hat.oob$variance.estimates
  standard_dev <- sqrt(variance)
  
  predictions_oob_splitsample <- tau.hat.oob.splitsample$predictions
  variance_splitsample <- tau.hat.oob.splitsample$variance.estimates
  standard_dev_splitsample <- sqrt(variance_splitsample)
    
  
  
  # get forest-wide avg tau's/tx effects (and the standard errpr)
  # NOTE:this is *not* with the split sample
  #avg_tx_effect_overall <- average_treatment_effect(tau.forest, target.sample = 'all')
  avg_tx_effect_overall <- avg_effect(Y = Y_splitsample, Y.hat = Y.hat_splitsample, W = W_splitsample, W.hat = W.hat_splitsample,
             weights = sample_weights_splitsample, predictions = predictions_oob_splitsample,
             subset = NULL)
  
  # test calibration of forest
  forest_calibration <- test_calibration(tau.forest) %>% broom::tidy()
  
  # output should be a list with
  # - avg causal effects
  group_tau_avgs <- c('overall'=avg_tx_effect_overall)
  
  # - the predictions + variance in a df with the sids
  if (flip_outcome == F) {
    tau_df_splitsample <- tibble(sid = cluster_ids_splitsample,
                     grf_tau_hat = predictions_oob_splitsample,
                     grf_variance = variance_splitsample,
                     grf_sd = standard_dev_splitsample) %>% 
      makeDummies()
    
    tau_df <- tibble(sid = cluster_ids,
                                 grf_tau_hat = predictions_oob,
                                 grf_variance = variance,
                                 grf_sd = standard_dev) %>% 
      makeDummies()
  }

  if (flip_outcome == T) {
    tau_df <- tibble(sid = cluster_ids,
                     grf_tau_hat = predictions_oob,
                     grf_variance = variance,
                     grf_sd = standard_dev) %>% 
      makeDummies(flipped = T)
    
    tau_df_splitsample <- tibble(sid = cluster_ids_splitsample,
                                 grf_tau_hat = predictions_oob_splitsample,
                                 grf_variance = variance_splitsample,
                                 grf_sd = standard_dev_splitsample) %>% 
      makeDummies(flipped = T)
  }  
  # now that we have the quartile dummies, we can get quartile level effects
  # ate.highest_25 <- average_treatment_effect(tau.forest,
  #                                            subset = tau_df$tau_quartile == 4)
  # ate.bottom_25 <- average_treatment_effect(tau.forest,
  #                                           subset = tau_df$tau_quartile == 1)
  # ate.quartile_2 <- average_treatment_effect(tau.forest,
  #                                            subset = tau_df$tau_quartile == 2)
  # ate.quartile_3 <- average_treatment_effect(tau.forest,
  #                                            subset = tau_df$tau_quartile == 3)
  # 
  # ate.bottom_75 <- average_treatment_effect(tau.forest,
  #                                           subset = tau_df$tau_quartile != 4)
  # ate.high <- average_treatment_effect(tau.forest,
  #                                      subset = tau_df$tau_quartile %in% c(3,4))
  # ate.low <- average_treatment_effect(tau.forest,
  #                                     subset = tau_df$tau_quartile %in% c(1,2))
  
  
  ate.highest_25 <- avg_effect(Y = Y_splitsample, Y.hat = Y.hat_splitsample, W = W_splitsample, W.hat = W.hat_splitsample,
                               weights = sample_weights_splitsample, predictions = predictions_oob_splitsample,
                               subset = tau_df_splitsample$tau_quartile == 4)
  
  ate.bottom_25 <- avg_effect(Y = Y_splitsample, Y.hat = Y.hat_splitsample, W = W_splitsample, W.hat = W.hat_splitsample,
                              weights = sample_weights_splitsample, predictions = predictions_oob_splitsample,
                              subset = tau_df_splitsample$tau_quartile == 1)
  
  ate.quartile_2 <- avg_effect(Y = Y_splitsample, Y.hat = Y.hat_splitsample, W = W_splitsample, W.hat = W.hat_splitsample,
                               weights = sample_weights_splitsample, predictions = predictions_oob_splitsample,
                               subset = tau_df_splitsample$tau_quartile == 2)
  
  ate.quartile_3 <- avg_effect(Y = Y_splitsample, Y.hat = Y.hat_splitsample, W = W_splitsample, W.hat = W.hat_splitsample,
                               weights = sample_weights_splitsample, predictions = predictions_oob_splitsample,
                               subset = tau_df_splitsample$tau_quartile == 3)
  
  ate.bottom_75 <- avg_effect(Y = Y_splitsample, Y.hat = Y.hat_splitsample, W = W_splitsample, W.hat = W.hat_splitsample,
                              weights = sample_weights_splitsample, predictions = predictions_oob_splitsample,
                              subset = tau_df_splitsample$tau_quartile != 4)
  
  ate.high <- avg_effect(Y = Y_splitsample, Y.hat = Y.hat_splitsample, W = W_splitsample, W.hat = W.hat_splitsample,
                         weights = sample_weights_splitsample, predictions = predictions_oob_splitsample,
                         subset = tau_df_splitsample$tau_quartile %in% c(3,4))
  
  
  ate.low <- avg_effect(Y = Y_splitsample, Y.hat = Y.hat_splitsample, W = W_splitsample, W.hat = W.hat_splitsample,
                        weights = sample_weights_splitsample, predictions = predictions_oob_splitsample,
                        subset = tau_df_splitsample$tau_quartile %in% c(1,2))
  
  
  
  subsample_tau_avgs <- list('highest_quartile' = ate.highest_25,
                             'bottom_three_quartiles' = ate.bottom_75,
                             'above_median'=ate.high,
                             'below_median'=ate.low,
                             'quartile_1' = ate.bottom_25,
                             'quartile_2' = ate.quartile_2,
                             'quartile_3' = ate.quartile_3,
                             'quartile_4' = ate.highest_25)
  
  # - 'augmented' master df (the above df merged with the full master dataset)
  augmented_df <- working_df %>% left_join(tau_df,#_splitsample,
                                           by='sid')
  
  
  
  # quartile baseline table
  quartile_heterogeneity_table <- augmented_df %>% make_quartile_baseline_table(subsample_tau_avgs)
  
  # subsample ate table
  subsample_ate_table <- subsample_tau_avgs %>% make_ate_summary_table(group_tau_avgs)
  
  # # subsample difference table
  # subsample_difference_table <- subsample_tau_avgs %>% make_subsample_difference_table()
  
  # calibration test (Jon D) output (no plot)
  
  
  output_list <- list(
    # 'tau_df' = tau_df,
    # 'augmented_df' = augmented_df,
    #'quartile_heterogeneity_table' = quartile_heterogeneity_table,
    #'tau_rank_plot'=tau_rank_plot,
    #'group_tau_avgs'=group_tau_avgs,
    #'subsample_tau_avgs'=subsample_tau_avgs,
    'subsample_ate_table' = subsample_ate_table,
    #'subsample_difference_table' = subsample_difference_table,
    #'forest_object'=tau.forest,
    'calibration_test'=forest_calibration,
    #'n_observations' = n.obs,
    #'tuning_output' = tau.forest$tuning.output,
    'quartile_heterogeneity_table' = quartile_heterogeneity_table)
  
  return(output_list)  
}



make_causal_forest_for_each_outcome <- function(input_dataset,
                                                splitsample_df,
                                                input_controls,
                                                outcomes,
                                                outcome_labels,
                                                flipped_outcomes=c()){
  final_forest_list <- list()
  
  for (i in 1:length(outcomes)){
    
    current_outcome <- outcomes[i]
    current_outcome_label <- outcome_labels[i]
    cat("Running: ", current_outcome)
    if (current_outcome %in% flipped_outcomes){
      single_forest_output <- make_single_causal_forest(dataset = input_dataset,
                                                        controls = input_controls,
                                                        outcome = current_outcome,
                                                        outcome_label = current_outcome_label,
                                                        flip_outcome = T,
                                                        splitsample_df = splitsample_df)
    } else {
      single_forest_output <- make_single_causal_forest(dataset = input_dataset,
                                                        controls = input_controls,
                                                        outcome = current_outcome,
                                                        outcome_label = current_outcome_label,
                                                        splitsample_df = splitsample_df)
    }
    
    
    
    final_forest_list[[current_outcome_label]] <- single_forest_output
    
  }
  
  return(final_forest_list)
}




# running models ----
start_time <- Sys.time()
# setting up furr plan
plan(multisession, workers = 10)


# setting u + running w sequential seeds

input_dataset <- master_pool
input_controls <- controls_sans_missingness_dummies


all_outcome_splitsample_lists <- list()

for (i in 1:length(outcomes_of_interest[1:2])){
  # make run splitsample forest function with this outcome
  
  current_outcome_of_interest <- outcomes_of_interest[i]
  current_outcome_label <- outcomes_of_interest_labels[i]
  
  run_splitsample_forest <- function(seed,
                                     current_outcome = current_outcome_of_interest,
                                     outcome_label = current_outcome_label,
                                     flipped_outcome_list = flipped_outcomes){
    # set seed
    set.seed(seed)
    
    # run single causal forest function
    if (current_outcome %in% flipped_outcome_list){
      single_forest_output <- make_single_causal_forest(dataset = input_dataset,
                                                        controls = input_controls,
                                                        outcome = current_outcome,
                                                        outcome_label = outcome_label,
                                                        flip_outcome = T)
    } else {
      single_forest_output <- make_single_causal_forest(dataset = input_dataset,
                                                        controls = input_controls,
                                                        outcome = current_outcome,
                                                        outcome_label = outcome_label)
    }
    
    return(single_forest_output) 
  }
  
  # run it as many times as we want using 
  splitsample_iterations <- future_map(1:n_splitsample_iterations, run_splitsample_forest)
  
  all_outcome_splitsample_lists[[i]] <- splitsample_iterations
}

end_time <- Sys.time()





# saveRDS(final_forests_missingness, file=paste0(filename, ".rds"))


