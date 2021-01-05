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
    dplyr::select(-statistic, "Full-sample P-val"=p.value) %>% filter(! term %in% 'Yhat')

  
  output_list <- list('calibration_test'=forest_calibration)
  
  return(output_list)  
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
make_single_causal_forest(master_pool, "graduated_ontime", "Graduated on-time",
                          controls = controls_with_missingness_dummies)






run_splitsample_models <- function(input_seed){
  
  # set seed
  set.seed(input_seed)
  
  # run cf for each outcome
  causal_forest_models <- make_single_causal_forest(master_pool, "graduated_ontime", "Graduated on-time",
                                                    controls = controls_with_missingness_dummies)$calibration_test
  
  
  return(causal_forest_models)
}



# run in parallel
n_iterations = 1000
seeds <- round(1000000*runif(n_iterations*2)) %>%
  unique() %>% .[1:n_iterations]


start_time <- Sys.time()
# setting up furr plan
plan(multisession, workers = 5)

all_splitsample_models <- future_map(seeds, run_splitsample_models,
                                     .options = future_options(scheduling=Inf, seed=T))

end_time <- Sys.time()

end_time - start_time

# note: despite the names of the above (adapted) functions, this is NOT a spitsample method; it is fullsample (name it as such)
all_splitsample_models %>% readr::write_rds('fullsample_graduation_1k_runs.Rds')


