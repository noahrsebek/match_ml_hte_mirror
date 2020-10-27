set.seed(10132020)

# loading packages and data ----

library('dplyr')
library('readr')
library('grf')
library('tidyr')
library('magrittr')
library('ggplot2')
library('furrr')
library('purrr')


#default_n_trees <- 5000
drop_duplicates <- T


bivariate_pte_comparisons <- tribble(
  ~outcome1, ~outcome2,
  "mathxil_z_post1_np", "mathgpa_post1",
  "mathxil_z_post1_np", "mathfail_post1",
  "mathxil_z_post1_np", "mathfailpercent_post1",
  "mathxil_z_post1_np", "nonmathgpa_post1",
  "mathgpa_post1", "nonmathgpa_post1",
  "mathfail_post1", "nonmathgpa_post1",
  "mathfailpercent_post1", "nonmathgpa_post1"
)


# load data ----
setwd("/export/projects/migrated-BAM/Match_AnalysisFiles_May2016/Results/2020/grf")
source(file = 'helper_functions.R')
source(file = 'grf_globals.R')
source(file = 'bart_helper_functions.R')

# overriding flipped outcomes (for the time being)
flipped_outcomes <- c()



master_pool <- load_master_dataset()



block_treatment_probability_df <- master_pool %>% group_by(blocknum, study) %>% summarize(p_block = mean(dmatch)) %>% ungroup()

singleton_blocks <- master_pool %>% group_by(blocknum, study) %>% filter(n()==1) %>% ungroup() %>% pull(blocknum)

master_pool <- master_pool %>%
  left_join(block_treatment_probability_df, by=c('study', 'blocknum')) %>% 
  mutate(inv_prob_weight = dmatch/p_block + (1-dmatch)/(1-p_block)) %>% 
  filter(! blocknum %in% singleton_blocks)


# # we want each core to run ALL outcomes
# for each iteration:
#   randomly split sample (block/treatment stratfication)
#   for each method:
#     for each outcome:
#       run method(outcome)  




# deduping and prepping data for split
duplicate_sids <- master_pool %>% 
  group_by(sid) %>% filter(n()>1) %>% 
  distinct(sid) %>% pull(sid)

master_pool <- master_pool %>% filter(!sid %in% duplicate_sids)

master_pool <- master_pool %>%  drop_na(inv_prob_weight) %>%  mutate(row_id = row_number())

# # imputing baselines
# master_pool <- master_pool %>% impute_block_baselines(controls_sans_missingness_dummies)

# editing our code to it runs with BART ----
library(bcf)

# dataset <- master_pool
# controls <- controls_sans_missingness_dummies
# outcome = outcomes_of_interest[1]



make_single_bcf <- function(dataset, controls, outcome){
  n_burn = 2000
  n_sim = 2000
  
  
  # keep observations complete on outcome
  working_df <- dataset %>% drop_na(all_of(outcome)) #%>%  mutate(row_id = row_number())
  
  working_df <- working_df[c(controls, outcome, 'dmatch', 'inv_prob_weight', 'sid')] %>% drop_na()
  
  n.obs <- nrow(working_df)
  
  X <- working_df[,controls] %>% as.matrix()
  W <- working_df %>% pull(dmatch)
  Y <- working_df %>% pull(tidyselect::all_of(outcome))
  #sample_weights <- working_df %>% pull(inv_prob_weight)
  sample_weights <- rep(1, length(Y)) # we don't need to do weighting; make some dummy weights that are all
  sids <- working_df %>% pull(sid)
  
  # estimate pi hat with grf regression tree
  forest.W <- regression_forest(X, W,
                                tune.parameters = c("sample.fraction", "mtry",
                                                    "honesty.fraction", "honesty.prune.leaves",
                                                    "alpha", "imbalance.penalty"),
                                min.node.size = 5)
  W.hat <- predict(forest.W)$predictions
  
  # fit bayesian causal forest
  bcf_fit <- bcf(Y, W, X, X, W.hat, nburn=n_burn, nsim=n_sim,
                 include_pi = 'both', verbose=T)
  
  # pull posterior estimates, take mean of each
  tau.hat_df <- bcf_fit$tau
  tau.hat <- tau.hat_df %>% colMeans()
  
  Y.hat <- Y - W*tau.hat
  
  # get overall stuff done  
  # - plots
  # - calibration tests
  # - naive calibration tests
  
  
  
  # make tau df
  tau_df <- tibble(sid = sids,
                   grf_tau_hat = tau.hat) %>% makeDummies()
  
  # make tau rank order plot
  tau_rank_order_plot <- make_tau_rank_order_plot(tau_df, n.obs,
                                                  names(outcome_and_label_list)[which(outcome_and_label_list == outcome)])
  
  
  # calibration tests and plot
  naive_calibration_tests <- make_naive_calibration_tests(tau_df, outcome)
  naive_linear <- naive_calibration_tests[[2]]
  naive_quantile <- naive_calibration_tests[[1]]
  calibration_plot <- make_calibration_plot(tau_df, outcome)
  calibration_table <- make_calibration_table(tau.hat, Y.hat)
  

  # average effect table
  avg_tx_effect_overall <- avg_effect(Y = Y, Y.hat = Y.hat, W = W, W.hat = pi_hat,
                                      weights = sample_weights,
                                      predictions = tau.hat,
                                      subset = NULL)
  group_tau_avgs <- c('overall'=avg_tx_effect_overall)
  subsample_tau_avgs <- make_average_effect_list(tau.hat, Y.hat)
  
  subsample_ate_table <- subsample_tau_avgs %>% make_ate_summary_table(group_tau_avgs)
  
  
  # # augmented df
  # augmented_df <- working_df %>% left_join(tau_df, by='sid')
  
    
  # now, we loop over the other realizations, calculating:
  # - calibration tests
  # - average treatment effects
  # and turning those ensembles of estimates into confidence intervals
  naive_linear_mashup <- NULL
  naive_quantile_mashup <- NULL
  calibration_mashup <- NULL
  subsample_avgs_mashup <- NULL
  
  for (i in 1:n_sim){
    
    temp_tau.hat <- tau.hat_df[i,]
    temp_tau_df <- tibble(sid = sids,
                          grf_tau_hat = tau.hat) %>% makeDummies()
    temp_Y.hat <- Y - W*temp_tau.hat
    
    # calibration tests and plot
    temp_naive_calibration_tests <- suppressMessages(make_naive_calibration_tests(temp_tau_df, outcome))
    
    temp_naive_linear <- temp_naive_calibration_tests[[2]] %>%
      select(term, estimate)
    temp_naive_quantile <- temp_naive_calibration_tests[[1]] %>%
      select(calibration_quantiles, estimate)
    
    temp_calibration_table <- make_calibration_table(temp_tau.hat, temp_Y.hat) %>% select(term, estimate)
    
    
    # average effect table
    avg_tx_effect_overall <- avg_effect(Y = Y,
                                        Y.hat = temp_Y.hat,
                                        W = W, W.hat = pi_hat,
                                        weights = sample_weights,
                                        predictions = temp_tau.hat,
                                        subset = NULL)
    temp_group_tau_avgs <- c('overall'=avg_tx_effect_overall)
    temp_subsample_tau_avgs <- make_average_effect_list(temp_tau.hat, temp_Y.hat)
    
    temp_subsample_ate_table <- temp_subsample_tau_avgs %>% make_ate_summary_table(temp_group_tau_avgs)
    
    
    # now bind rows to the 4 tables
    naive_linear_mashup <- bind_rows(naive_linear_mashup, temp_naive_linear)
    naive_quantile_mashup <- bind_rows(naive_quantile_mashup, temp_naive_quantile)
    calibration_mashup <- bind_rows(calibration_mashup, temp_calibration_table)
    subsample_avgs_mashup <- bind_rows(subsample_avgs_mashup, temp_subsample_ate_table)
    
    
    
  }
  

  
  # now turn the mashup tables into confidence intervals, and add into the 'main' table we estimated before the loop
  naive_linear <- naive_linear %>% left_join(add_CIs_to_mashup_tables(naive_linear_mashup))
  naive_quantile <- naive_quantile %>% left_join(add_CIs_to_mashup_tables(naive_quantile_mashup))
  calibration_table <- calibration_table %>% left_join(add_CIs_to_mashup_tables(calibration_mashup))
  subsample_ate_table <- subsample_ate_table %>% left_join(add_CIs_to_mashup_tables(subsample_avgs_mashup))
  
  
  bart_output = list(
    'tau_rank_order_plot' = tau_rank_order_plot,
    'calibration_test' = calibration_table,
    'ate_table' = subsample_ate_table,
    'calibration_tests_naive_quantile_dummies' = naive_quantile,
    'calibration_tests_naive_linear' = naive_linear,
    'calibration_plot' = calibration_plot
  )
  
  return(bart_output)
}

make_single_bcf(master_pool, controls_sans_missingness_dummies, outcomes_of_interest[1])



make_bcf_for_each_outcome <- function(input_dataset,
                                      input_controls,
                                      outcomes,
                                      outcome_labels){
  final_forest_list <- list()
  
  for (i in 1:length(outcomes)){
    
    current_outcome <- outcomes[i]
    current_outcome_label <- outcome_labels[i]
    cat("Running: ", current_outcome)
    
    single_forest_output <- make_single_bcf(dataset = input_dataset,
                                             controls = input_controls,
                                             outcome = current_outcome)
    
    
    
    
    final_forest_list[[current_outcome_label]] <- single_forest_output
    
  }
  
  return(final_forest_list)
}


all_bcf_output <- make_bcf_for_each_outcome(master_pool,
                          controls_sans_missingness_dummies,
                          outcomes_of_interest,
                          outcomes_of_interest_labels)

saveRDS(all_bcf_output, 'all_bcf_models.Rds')
