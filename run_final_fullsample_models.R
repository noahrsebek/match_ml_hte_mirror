set.seed(20201111)

# loading packages and setting global vars ----
library('dplyr')
library('readr')
library('grf')
library('tidyr')
library('magrittr')
library('ggplot2')
library('furrr')
library('purrr')
library(forcats)

default_n_trees <- 10000
to_cluster <- F



# load data ----
setwd("/export/projects/migrated-BAM/Match_AnalysisFiles_May2016/Results/2020/grf")
source(file = 'helper_functions.R')
source(file = 'mltx_globals.R')



master_pool <- load_master_dataset()

# main ML functions ----
make_single_causal_forest <- function(dataset,
                                      controls,
                                      outcome,
                                      outcome_label,
                                      tune_num_reps = 50,
                                      tune_num_draws = 1000,
                                      tune_num_trees = 200,
                                      #splitsample_frac = splitsample_fraction
                                      splitsample_df,
                                      tuned_parameters = c("sample.fraction", "mtry",
                                                           "honesty.fraction", "honesty.prune.leaves",
                                                           "alpha", "imbalance.penalty")){
  
  
  
  # keep observations complete on outcome and that has an inv prob weight
  working_df <- dataset %>% drop_na(all_of(outcome))
  
  n.obs <- nrow(working_df)
  
  X <- working_df[,controls]
  W <- working_df %>% pull(dmatch)
  Y <- working_df %>% pull(tidyselect::all_of(outcome))
  sample_weights <- working_df %>% pull(inv_prob_weight)
  
  cluster_ids <- working_df %>% pull(sid) # we'll use these if we have duplicates n the data
  
  
  
  
  # fit W forest
  forest.W <- regression_forest(X, W, tune.parameters = tuned_parameters, min.node.size = 5)
  W.hat <- predict(forest.W)$predictions
  
  
  # fit Y forest
  forest.Y <- regression_forest(X, Y, tune.parameters = tuned_parameters, min.node.size = 5)
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

  
  
  tau.hat.oob <- predict(tau.forest,
                         estimate.variance = T)
  
  
  predictions_oob <- tau.hat.oob$predictions
  variance <- tau.hat.oob$variance.estimates
  standard_dev <- sqrt(variance)
  

  
  
  
  # get forest-wide avg tau's/tx effects (and the standard errpr)
  # NOTE:this is *not* with the split sample
  #avg_tx_effect_overall <- average_treatment_effect(tau.forest, target.sample = 'all')
  avg_tx_effect_overall <- c('estimate'=mean(predictions_oob))
  
  
  # test calibration of forest
  #W.hat_splitsample <- W.hat_splitsample
  
  mean.pred <- weighted.mean(predictions_oob, sample_weights)
  # DF <- data.frame(target = unname(Y - Y.hat), 
  #                  mean.forest.prediction = unname(W - W.hat) * mean.pred,
  #                  differential.forest.prediction = unname(W - W.hat) * (predictions_oob - mean.pred))
  
  
  DF <- data.frame(target = unname(Y),
                   Yhat = Y.hat,
                   mean.forest.prediction = unname(W - W.hat) * mean.pred,
                   differential.forest.prediction = unname(W - W.hat) * (predictions_oob - mean.pred))
  
  
  
  # changed "observation.weight" to "sample_weights"
  best.linear.predictor <- lm(target ~ Yhat + mean.forest.prediction + 
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
  forest_calibration <- blp.summary %>%
    broom::tidy() %>% 
    mutate(lower_CI_full = estimate - 1.96*std.error,
           upper_CI_full = estimate + 1.96*std.error) %>% 
    select(-statistic, "In-sample P-val"=p.value) %>% filter(! term %in% 'Yhat')
  
  #forest_calibration <- test_calibration(tau.forest) %>% broom::tidy()
  
  # output should be a list with
  # - avg causal effects
  group_tau_avgs <- c('overall'=avg_tx_effect_overall)
  
  # - the predictions + variance in a df with the sids
  tau_df <- tibble(sid = cluster_ids,
                   grf_tau_hat = predictions_oob,
                   grf_variance = variance,
                   grf_sd = standard_dev) %>% 
    makeDummies()
  
  # now that we have the quartile dummies, we can get quartile level effects
  ate.highest_25 <- avg_effect_simple(predictions = predictions_oob,
                                      subset = tau_df$tau_quartile == 4)
  
  ate.bottom_25 <- avg_effect_simple(predictions = predictions_oob,
                                     subset = tau_df$tau_quartile == 1)
  
  ate.quartile_2 <- avg_effect_simple(predictions = predictions_oob,
                                      subset = tau_df$tau_quartile == 2)
  
  ate.quartile_3 <- avg_effect_simple(predictions = predictions_oob,
                                      subset = tau_df$tau_quartile == 3)
  
  ate.bottom_75 <- avg_effect_simple(predictions = predictions_oob,
                                     subset = tau_df$tau_quartile != 4)
  
  ate.high <- avg_effect_simple(predictions = predictions_oob,
                                subset = tau_df$tau_quartile %in% c(3,4))
  
  
  ate.low <- avg_effect_simple(predictions = predictions_oob,
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
  
  
  # quartile baseline table
  # quartile_heterogeneity_table <- augmented_df %>% make_quartile_baseline_table(subsample_tau_avgs)
  # median_heterogeneity_table <- augmented_df %>% make_above_below_median_baseline_table(subsample_tau_avgs)
  quartile_heterogeneity_table <- tau_df %>% make_quartile_summary_table()
  median_heterogeneity_table <- tau_df %>% make_above_below_median_summary_table()
  
  # subsample ate table
  subsample_ate_table <- subsample_tau_avgs %>% make_ate_summary_table(group_tau_avgs)
  
  # calibration test (Jon D) output (no plot)
  calibration_tests_naive <- make_naive_calibration_tests(tau_df, outcome)
  calibration_tests_naive[[1]] <- calibration_tests_naive[[1]] %>% 
    select(-p.value)
  
  calibration_tests_naive[[2]] <- calibration_tests_naive[[2]] %>% 
    select(-p.value)
  
  
  # isr itt table by quartile
  isr_itt_tables_quartile <- make_isr_itt_quantile_tables(tau_df, n_calibration_quantiles = 4)

  # isr itt table by above/below median
  isr_itt_tables_median <- make_isr_itt_quantile_tables(tau_df, n_calibration_quantiles = 2)
  
    
  # tau rank plot
  tau_rank_plot <- tau_df %>% make_tau_rank_order_plot(outcome_label)
  
  # calibration plot
  calibration_plot <- tau_df %>% make_calibration_plot(outcome)
  
  output_list <- list(
    'tau_df' = tau_df,
    # 'augmented_df' = augmented_df,
    #'quartile_heterogeneity_table' = quartile_heterogeneity_table,
    'tau_rank_plot'=tau_rank_plot,
    'calibration_plot'=calibration_plot,
    #'group_tau_avgs'=group_tau_avgs,
    #'subsample_tau_avgs'=subsample_tau_avgs,
    'subsample_ate_table' = subsample_ate_table,
    #'subsample_difference_table' = subsample_difference_table,
    #'forest_object'=tau.forest,
    'calibration_test'=forest_calibration,
    'calibration_tests_naive_quantile_dummies' = calibration_tests_naive[[1]],
    'calibration_tests_naive_linear' = calibration_tests_naive[[2]],
    'quartile_heterogeneity_table' = quartile_heterogeneity_table,
    'median_heterogeneity_table' = median_heterogeneity_table,
    'ISR_ITT_quartile_table_wave1' = isr_itt_tables_quartile[[1]],
    'ISR_ITT_quartile_table_wave2' = isr_itt_tables_quartile[[2]],
    'ISR_ITT_median_table_wave1' = isr_itt_tables_median[[1]],
    'ISR_ITT_median_table_wave2' = isr_itt_tables_median[[2]]
  )
  
  return(output_list)  
}



make_causal_forest_for_each_outcome <- function(input_dataset,
                                                input_controls,
                                                outcomes,
                                                outcome_labels,
                                                flipped_outcomes=c()){
  final_forest_list <- list()
  
  for (i in 1:length(outcomes)){
    
    current_outcome <- outcomes[i]
    current_outcome_label <- outcome_labels[i]
    cat(":: Running: ", current_outcome, "~ \n")
    
    single_forest_output <- make_single_causal_forest(dataset = input_dataset,
                                                      controls = input_controls,
                                                      outcome = current_outcome,
                                                      outcome_label = current_outcome_label)
    
    
    
    
    final_forest_list[[current_outcome_label]] <- single_forest_output
    
  }
  
  return(final_forest_list)
}





# OVERALL DATA WORK ----
# - adding inv probability weights

block_treatment_probability_df <- master_pool %>% group_by(blocknum, study) %>% summarize(p_block = mean(dmatch)) %>% ungroup()

singleton_blocks <- master_pool %>% group_by(blocknum, study) %>% filter(n()==1) %>% ungroup() %>% pull(blocknum)

master_pool <- master_pool %>%
  left_join(block_treatment_probability_df, by=c('study', 'blocknum')) %>% 
  mutate(inv_prob_weight = dmatch/p_block + (1-dmatch)/(1-p_block)) %>% 
  filter(! blocknum %in% singleton_blocks)




# deduping and prepping data for split
duplicate_sids <- master_pool %>% 
  group_by(sid) %>% filter(n()>1) %>% 
  distinct(sid) %>% pull(sid)

master_pool <- master_pool %>% filter(!sid %in% duplicate_sids,
                                      !is.na(inv_prob_weight))



# running models ----
all_fullsample_cf_models <- make_causal_forest_for_each_outcome(master_pool,
                                                               controls_with_missingness_dummies,
                                                               outcomes_of_interest,
                                                               outcomes_of_interest_labels)



all_fullsample_cf_models %>% write_rds("all_fullsample_cf_models.Rds")

