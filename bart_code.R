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
  sample_weights <- working_df %>% pull(inv_prob_weight)
  sids <- working_df %>% pull(sid)
  
  pi_hat <- lm(dmatch ~ ., data=working_df[c(controls, 'dmatch')]) %>% predict()
  
  bcf_fit <- bcf(Y, W, X, X, pi_hat, nburn=n_burn, nsim=n_sim,
                 include_pi = 'both', verbose=T)
  tau.hat <- bcf_fit$tau %>% colMeans()
  
  y_hat <- bcf_fit$yhat
  
  
  # calibration test
  mean.pred <- weighted.mean(tau.hat, sample_weights)
  
  
  
  
  # tau df
  tau_df <- tibble(sid = sids,
                   grf_tau_hat = tau.hat) %>% makeDummies()
  
  # tau rank order plot
  tau_rank_order_plot <- make_tau_rank_order_plot(tau_df, n.obs,
                                                  names(outcome_and_label_list)[which(outcome_and_label_list == outcome)])
  
  
  # calibration tests and plot
  naive_calibration_tests <- make_naive_calibration_tests(tau_df, outcome)
  
  calibration_plot <- make_calibration_plot(tau_df, outcome)
  
  
  # augmented df
  augmented_df <- working_df %>% left_join(tau_df, by='sid')
  
  # quartile heterogeneity table (not crucial rn)
  
  
  # bart output
  # - naive calibration tests
  # - calibraton plot
  # - rank order plot
  # - quartile baseline table
  bart_output = list(
    'tau_rank_order_plot' = tau_rank_order_plot,
    'calibration_tests_naive_quantile_dummies' = naive_calibration_tests[[1]],
    'calibration_tests_naive_linear' = naive_calibration_tests[[2]],
    'calibration_plot' = calibration_plot
  )
  
  return(bart_output)
}

#make_single_bcf(master_pool, controls_sans_missingness_dummies, outcomes_of_interest[1])



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
