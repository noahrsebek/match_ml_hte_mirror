library('dplyr')
library('readr')
library('grf')
library('tidyr')
library('magrittr')
library('ggplot2')



n_splitsample_iterations <- 10
splitsample_fraction <- 0.5








default_sample_fraction <- 0.5
default_n_trees <- 5000
to_cluster <- F
drop_duplicates <- T
set.seed(20200529)

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


makeDummies <- function(dataframe,
                        flipped=F){
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



make_single_causal_forest <- function(dataset,
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
  
  
  if (add_all_interactions==T){
    X <- X %>% add_pairwise_interactions(controls)
  }
  
  if (cluster==T){
    tau.forest <- causal_forest(X, Y, W,
                                tune.parameters = c("sample.fraction", "mtry",
                                                    "honesty.fraction", "honesty.prune.leaves",
                                                    "alpha", "imbalance.penalty"),
                                min.node.size = 5,
                                sample.weights = sample_weights,
                                clusters = cluster_ids,
                                sample.fraction = default_sample_fraction,
                                num.trees = default_n_trees,
                                tune.num.trees = tune_num_trees,
                                tune.num.draws = tune_num_draws,
                                tune.num.reps = tune_num_reps)
  }
  if (cluster==F){
    tau.forest <- causal_forest(X, Y, W,
                                tune.parameters = c("sample.fraction", "mtry",
                                                    "honesty.fraction", "honesty.prune.leaves",
                                                    "alpha", "imbalance.penalty"),
                                min.node.size = 5,
                                sample.weights = sample_weights,
                                #clusters = cluster_ids,
                                sample.fraction = default_sample_fraction,
                                num.trees = default_n_trees,
                                tune.num.trees = tune_num_trees,
                                tune.num.draws = tune_num_draws,
                                tune.num.reps = tune_num_reps)
  }

  
  tau.hat.oob <- predict(tau.forest, estimate.variance = T)
  predictions_oob <- tau.hat.oob$predictions
  variance <- tau.hat.oob$variance.estimates
  standard_dev <- sqrt(variance)
  
  # get forest-wide avg tau's/tx effects (and the standard errpr)
  avg_tx_effect_overall <- average_treatment_effect(tau.forest, target.sample = 'all')
  avg_tx_effect_treated <- average_treatment_effect(tau.forest, target.sample = 'treated')
  avg_tx_effect_control <- average_treatment_effect(tau.forest, target.sample = 'control')
  
  # test calibration of forest
  forest_calibration <- test_calibration(tau.forest) %>% broom::tidy()
  
  # output should be a list with
  # - avg causal effects
  group_tau_avgs <- c('overall'=avg_tx_effect_overall,
                      'treated'=avg_tx_effect_treated,
                      'control'=avg_tx_effect_control)
  
  # - the predictions + variance in a df with the sids
  if (flip_outcome == F) {
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
  }  
  # now that we have the quartile dummies, we can get quartile level effects
  ate.highest_25 <- average_treatment_effect(tau.forest,
                                             subset = tau_df$tau_quartile == 4)
  ate.bottom_25 <- average_treatment_effect(tau.forest,
                                            subset = tau_df$tau_quartile == 1)
  ate.quartile_2 <- average_treatment_effect(tau.forest,
                                             subset = tau_df$tau_quartile == 2)
  ate.quartile_3 <- average_treatment_effect(tau.forest,
                                             subset = tau_df$tau_quartile == 3)
  
  ate.bottom_75 <- average_treatment_effect(tau.forest,
                                            subset = tau_df$tau_quartile != 4)
  ate.high <- average_treatment_effect(tau.forest,
                                       subset = tau_df$tau_quartile %in% c(3,4))
  ate.low <- average_treatment_effect(tau.forest,
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

# OVERALL DATA WORK ----
# - adding inv probability weights
# - rand block mean imputation
# - 
block_treatment_probability_df <- master_pool %>% group_by(blocknum, study) %>% summarize(p_block = mean(dmatch)) %>% ungroup()

singleton_blocks <- master_pool %>% group_by(blocknum, study) %>% filter(n()==1) %>% ungroup() %>% pull(blocknum)

master_pool <- master_pool %>%
  left_join(block_treatment_probability_df, by=c('study', 'blocknum')) %>% 
  mutate(inv_prob_weight = dmatch/p_block + (1-dmatch)/(1-p_block)) %>% 
  filter(! blocknum %in% singleton_blocks)

# W_i = D/p_block + (1-D)/(1-p_block)




# question: % of missing values for these variables?
master_pool[,controls] %>% is.na() %>% colSums()


# let's impute these variables
#master_pool_imputed <- master_pool %>% impute_block_baselines(controls, type='mean')








make_forest_for_each_outcome <- function(input_dataset,
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
      single_forest_output <- make_single_causal_forest(dataset = input_dataset,
                                                        controls = input_controls,
                                                        outcome = current_outcome,
                                                        outcome_label = current_outcome_label,
                                                        flip_outcome = T)
    } else {
      single_forest_output <- make_single_causal_forest(dataset = input_dataset,
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


start_time <- Sys.time()
final_forests_missingness <- make_forest_for_each_outcome(master_pool,
                                                          controls_sans_missingness_dummies,
                                                          outcomes_of_interest,
                                                          outcomes_of_interest_labels,
                                                          flipped_outcomes)
end_time <- Sys.time()


# saveRDS(final_forests_missingness, file=paste0(filename, ".rds"))


