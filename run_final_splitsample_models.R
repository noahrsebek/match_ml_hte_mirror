overall_seeds = c(09212020, 09242020, 09252020, 09302020, 05101995, 12301990, 05111995, 02271993)

for (overall_seed in overall_seeds){
  
  set.seed(overall_seed)
  
  # loading packages and setting global vars ----
  library('dplyr')
  library('readr')
  library('grf')
  library('tidyr')
  library('magrittr')
  library('ggplot2')
  library('furrr')
  library('purrr')
  
  
  n_splitsample_iterations <- 125
  splitsample_fraction <- 0.5
  default_n_trees <- 5000
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
                                        flip_outcome=FALSE,
                                        cluster=to_cluster,
                                        add_all_interactions = F,
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
      # fit W forest
      forest.W <- regression_forest(X, W, tune.parameters = tuned_parameters, min.node.size = 5, clusters = cluster_ids_splitsample)
      W.hat_splitsample <- predict(forest.W, newdata=X_splitsample)$predictions
      W.hat <- predict(forest.W)$predictions
      
      
      # fit Y forest
      forest.Y <- regression_forest(X, Y, tune.parameters = tuned_parameters, min.node.size = 5, clusters = cluster_ids_splitsample)
      Y.hat_splitsample <- predict(forest.Y, newdata=X_splitsample)$predictions
      Y.hat <- predict(forest.Y)$predictions
      
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
    avg_tx_effect_overall <- c('estimate'=mean(predictions_oob_splitsample))
    
    
    # test calibration of forest
    #W.hat_splitsample <- W.hat_splitsample
    
    mean.pred <- weighted.mean(predictions_oob_splitsample, sample_weights_splitsample)
    # DF <- data.frame(target = unname(Y_splitsample - Y.hat_splitsample), 
    #                  mean.forest.prediction = unname(W_splitsample - W.hat_splitsample) * mean.pred,
    #                  differential.forest.prediction = unname(W_splitsample - W.hat_splitsample) * (predictions_oob_splitsample - mean.pred))
    
    
    DF <- data.frame(target = unname(Y_splitsample),
                     Yhat = Y.hat_splitsample,
                     mean.forest.prediction = unname(W_splitsample - W.hat_splitsample) * mean.pred,
                     differential.forest.prediction = unname(W_splitsample - W.hat_splitsample) * (predictions_oob_splitsample - mean.pred))
    
    
    
    # changed "observation.weight" to "sample_weights"
    best.linear.predictor <- lm(target ~ Yhat + mean.forest.prediction + 
                                  differential.forest.prediction + 0, weights = sample_weights_splitsample, 
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
      mutate(lower_CI = estimate - 2.24*std.error,
             upper_CI = estimate + 2.24*std.error) %>% 
      select(-std.error, -statistic) %>% filter(! term %in% 'Yhat')
    
    #forest_calibration <- test_calibration(tau.forest) %>% broom::tidy()
    
    # output should be a list with
    # - avg causal effects
    group_tau_avgs <- c('overall'=avg_tx_effect_overall)
    
    # - the predictions + variance in a df with the sids
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
    
      # now that we have the quartile dummies, we can get quartile level effects
      ate.highest_25 <- avg_effect_simple(predictions = predictions_oob_splitsample,
                                          subset = tau_df_splitsample$tau_quartile == 4)
      
      ate.bottom_25 <- avg_effect_simple(predictions = predictions_oob_splitsample,
                                         subset = tau_df_splitsample$tau_quartile == 1)
      
      ate.quartile_2 <- avg_effect_simple(predictions = predictions_oob_splitsample,
                                          subset = tau_df_splitsample$tau_quartile == 2)
      
      ate.quartile_3 <- avg_effect_simple(predictions = predictions_oob_splitsample,
                                          subset = tau_df_splitsample$tau_quartile == 3)
      
      ate.bottom_75 <- avg_effect_simple(predictions = predictions_oob_splitsample,
                                         subset = tau_df_splitsample$tau_quartile != 4)
      
      ate.high <- avg_effect_simple(predictions = predictions_oob_splitsample,
                                    subset = tau_df_splitsample$tau_quartile %in% c(3,4))
      
      
      ate.low <- avg_effect_simple(predictions = predictions_oob_splitsample,
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
    #augmented_df <- working_df %>% left_join(tau_df, by='sid')
    augmented_df <- splitsample_df %>% left_join(tau_df_splitsample, by='sid')
    
    
    
    # # quartile baseline table
    # quartile_heterogeneity_table <- augmented_df %>% make_quartile_baseline_table(subsample_tau_avgs)
    # 
    # # subsample ate table
    # subsample_ate_table <- subsample_tau_avgs %>% make_ate_summary_table(group_tau_avgs)
    
    # calibration test (Jon D) output (no plot)
    calibration_tests_naive <- make_naive_calibration_tests(tau_df_splitsample, outcome)
    calibration_tests_naive[[1]] <- calibration_tests_naive[[1]] %>%
      mutate(lower_CI = estimate - 2.24*std.error,
             upper_CI = estimate + 2.24*std.error) %>% 
      select(-std.error)
    
    calibration_tests_naive[[2]] <- calibration_tests_naive[[2]] %>%
      mutate(lower_CI = estimate - 2.24*std.error,
             upper_CI = estimate + 2.24*std.error) %>% 
      select(-std.error)
    
    
    output_list <- list(
      # 'tau_df' = tau_df,
      # 'augmented_df' = augmented_df,
      #'quartile_heterogeneity_table' = quartile_heterogeneity_table,
      #'tau_rank_plot'=tau_rank_plot,
      #'group_tau_avgs'=group_tau_avgs,
      #'subsample_tau_avgs'=subsample_tau_avgs,
      #'subsample_ate_table' = subsample_ate_table,
      #'subsample_difference_table' = subsample_difference_table,
      #'forest_object'=tau.forest,
      'calibration_test'=forest_calibration,
      'calibration_tests_naive_quantile_dummies' = calibration_tests_naive[[1]],
      'calibration_tests_naive_linear' = calibration_tests_naive[[2]]#,
      #'n_observations' = n.obs,
      #'tuning_output' = tau.forest$tuning.output,
      #'tau_df_splitsample' = tau_df_splitsample,
      # 'quartile_heterogeneity_table' = quartile_heterogeneity_table
      )
    
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
      cat(":: Running: ", current_outcome, "~ \n")
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
  
  
  
  
  

  make_single_bart <- function(dataset,
                               controls,
                               outcome,
                               outcome_label,
                               splitsample_df){
    
  
    # keep observations complete on outcome and that has an inv prob weight
    working_df <- dataset %>% drop_na(all_of(outcome)) 
    splitsample_df <- splitsample_df %>% drop_na(all_of(outcome))
    
    n.obs <- nrow(working_df)
    
    X <- working_df[,controls]
    W <- working_df %>% pull(dmatch)
    Y <- working_df %>% pull(tidyselect::all_of(outcome))
    #sample_weights <- working_df %>% pull(inv_prob_weight)
    # make equal weights; some functions require weights but we dont need inv prob weights for BART
    #  because we include prop score estimates in the estimation procedure
    sample_weights <- rep(1, length(Y)) 
    
    X_splitsample <- splitsample_df[,controls]
    Y_splitsample <- splitsample_df %>% pull(tidyselect::all_of(outcome))
    W_splitsample <- splitsample_df %>% pull(dmatch)
    sample_weights_splitsample <- rep(1, length(Y_splitsample))
    
    cluster_ids <- working_df %>% pull(sid) # we'll use these if we have duplicates n the data
    cluster_ids_splitsample <- splitsample_df %>% pull(sid)
    
    
    # fit W forest
    forest.W <- regression_forest(X, W,
                                  tune.parameters = c("sample.fraction", "mtry",
                                                      "honesty.fraction", "honesty.prune.leaves",
                                                      "alpha", "imbalance.penalty"),
                                  min.node.size = 5)
    W.hat_splitsample <- predict(forest.W, newdata=X_splitsample)$predictions
    W.hat <- predict(forest.W)$predictions
    
    
    # making BART models and getting predictions
    # - make sure to include both W and propensity score as features
    X_bart <- bind_cols(X, "W.hat"=W.hat,"W"=W)
    X_bart_splitsample <- bind_cols(X_splitsample, "W.hat"=W.hat_splitsample,"W"=W_splitsample)
    
    # see if outcome is continuous
    if (length(unique(Y))==2){
      # discrete: use logit bart (lbart; or probit bart pbart)
      bart_fit <- mc.gbart(as.matrix(X_bart), Y, mc.cores = 2, nice=5, type='pbart')
      
      Y.hat_splitsample <- predict(bart_fit, as.matrix(X_bart_splitsample), mc.cores = 2, nice=5)
      Y.hat_splitsample <- Y.hat_splitsample$yhat.test %>% colMeans()
      
      X_bart_Y0_splitsample <- X_bart_splitsample %>% mutate(W=0)
      X_bart_Y1_splitsample <- X_bart_splitsample %>% mutate(W=1)
      
      bart_preds_Y0_split <- predict(bart_fit, as.matrix(X_bart_Y0_splitsample), mc.cores=2, nice=5)
      bart_preds_Y1_split <- predict(bart_fit, as.matrix(X_bart_Y1_splitsample), mc.cores=2, nice=5)
      
      bart_preds_Y0_split <- bart_preds_Y0_split$yhat.test %>% colMeans()
      bart_preds_Y1_split <- bart_preds_Y1_split$yhat.test %>% colMeans()
      
      
    } else {
      # continuous
      bart_fit <- mc.gbart(as.matrix(X_bart), Y, mc.cores = 2, nice=2)
      Y.hat_splitsample <- predict(bart_fit, as.matrix(X_bart_splitsample), mc.cores = 2, nice=5) %>% colMeans()
      
      
      # X_bart_Y0 <- X_bart %>% mutate(W=0)
      # X_bart_Y1 <- X_bart %>% mutate(W=1)
      # 
      # bart_preds_Y0 <- predict(bart_fit, as.matrix(X_bart_Y0), mc.cores=2, nice=5) %>% colMeans()
      # bart_preds_Y1 <- predict(bart_fit, as.matrix(X_bart_Y1), mc.cores=2, nice=5) %>% colMeans()
      # predictions <- bart_preds_Y1 - bart_preds_Y0
      
      
      X_bart_Y0_splitsample <- X_bart_splitsample %>% mutate(W=0)
      X_bart_Y1_splitsample <- X_bart_splitsample %>% mutate(W=1)
      
      bart_preds_Y0_split <- predict(bart_fit, as.matrix(X_bart_Y0_splitsample), mc.cores=2, nice=5) %>% colMeans()
      bart_preds_Y1_split <- predict(bart_fit, as.matrix(X_bart_Y1_splitsample), mc.cores=2, nice=5) %>% colMeans()
      
    }
    
    predictions_oob_splitsample <- bart_preds_Y1_split - bart_preds_Y0_split
    
    
    
    
    
    # get forest-wide avg tau's/tx effects (and the standard errpr)
    # avg_tx_effect_overall <- avg_effect(Y = Y_splitsample, Y.hat = Y.hat_splitsample, W = W_splitsample, W.hat = W.hat_splitsample,
    #                                     weights = sample_weights_splitsample, predictions = predictions_oob_splitsample,
    #                                     subset = NULL)
    avg_tx_effect_overall <- c('estimate'=mean(predictions_oob_splitsample))
    
    
    # test calibration of forest
    #W.hat_splitsample <- W.hat_splitsample
    
    mean.pred <- weighted.mean(predictions_oob_splitsample, sample_weights_splitsample)
    # DF <- data.frame(target = unname(Y_splitsample - Y.hat_splitsample), 
    #                  mean.forest.prediction = unname(W_splitsample - W.hat_splitsample) * mean.pred,
    #                  differential.forest.prediction = unname(W_splitsample - W.hat_splitsample) * (predictions_oob_splitsample - mean.pred))
    
    
    DF <- data.frame(target = unname(Y_splitsample),
                     Yhat = Y.hat_splitsample,
                     mean.forest.prediction = unname(W_splitsample - W.hat_splitsample) * mean.pred,
                     differential.forest.prediction = unname(W_splitsample - W.hat_splitsample) * (predictions_oob_splitsample - mean.pred))
    
    
    
    # changed "observation.weight" to "sample_weights"
    best.linear.predictor <- lm(target ~ Yhat + mean.forest.prediction + 
                                  differential.forest.prediction + 0, weights = sample_weights_splitsample, 
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
      mutate(lower_CI = estimate - 2.24*std.error,
             upper_CI = estimate + 2.24*std.error) %>% 
      select(-std.error, -statistic)
    
    #forest_calibration <- test_calibration(tau.forest) %>% broom::tidy()
    
    # output should be a list with
    # - avg causal effects
    group_tau_avgs <- c('overall'=avg_tx_effect_overall)
    
    # - the predictions + variance in a df with the sids
    tau_df_splitsample <- tibble(sid = cluster_ids_splitsample,
                                 grf_tau_hat = predictions_oob_splitsample) %>% 
      makeDummies()
    
    
    
    
    
    ate.highest_25 <- avg_effect_simple(predictions = predictions_oob_splitsample,
                                        subset = tau_df_splitsample$tau_quartile == 4)
    
    ate.bottom_25 <- avg_effect_simple(predictions = predictions_oob_splitsample,
                                       subset = tau_df_splitsample$tau_quartile == 1)
    
    ate.quartile_2 <- avg_effect_simple(predictions = predictions_oob_splitsample,
                                        subset = tau_df_splitsample$tau_quartile == 2)
    
    ate.quartile_3 <- avg_effect_simple(predictions = predictions_oob_splitsample,
                                        subset = tau_df_splitsample$tau_quartile == 3)
    
    ate.bottom_75 <- avg_effect_simple(predictions = predictions_oob_splitsample,
                                       subset = tau_df_splitsample$tau_quartile != 4)
    
    ate.high <- avg_effect_simple(predictions = predictions_oob_splitsample,
                                  subset = tau_df_splitsample$tau_quartile %in% c(3,4))
    
    
    ate.low <- avg_effect_simple(predictions = predictions_oob_splitsample,
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
    #augmented_df <- working_df %>% left_join(tau_df, by='sid')
    augmented_df <- splitsample_df %>% left_join(tau_df_splitsample, by='sid')
    
    
    
    # quartile baseline table
    quartile_heterogeneity_table <- augmented_df %>% make_quartile_baseline_table(subsample_tau_avgs)
    names(quartile_heterogeneity_table)[2:5] <- paste0("Q", 1:4)
    
    # subsample ate table
    subsample_ate_table <- subsample_tau_avgs %>% make_ate_summary_table(group_tau_avgs)
    
    # # subsample difference table
    # subsample_difference_table <- subsample_tau_avgs %>% make_subsample_difference_table()
    
    # calibration test (Jon D) output (no plot)
    calibration_tests_naive <- make_naive_calibration_tests(tau_df_splitsample, outcome)
    calibration_tests_naive[[1]] <- calibration_tests_naive[[1]] %>%
      mutate(lower_CI = estimate - 2.24*std.error,
             upper_CI = estimate + 2.24*std.error) %>% 
      select(-std.error)
    
    calibration_tests_naive[[2]] <- calibration_tests_naive[[2]] %>%
      mutate(lower_CI = estimate - 2.24*std.error,
             upper_CI = estimate + 2.24*std.error) %>% 
      select(-std.error)
    
    
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
      'calibration_tests_naive_quantile_dummies' = calibration_tests_naive[[1]],
      'calibration_tests_naive_linear' = calibration_tests_naive[[2]],
      #'n_observations' = n.obs,
      #'tuning_output' = tau.forest$tuning.output,
      'quartile_heterogeneity_table' = quartile_heterogeneity_table#,
      #'tau_df_splitsample' = tau_df_splitsample
    )
    
    return(output_list)  
  }
  
  
  
  make_bart_for_each_outcome <- function(input_dataset,
                                         splitsample_df,
                                         input_controls,
                                         outcomes,
                                         outcome_labels,
                                         flipped_outcomes=c()){
    final_forest_list <- list()
    
    for (i in 1:length(outcomes)){
      
      current_outcome <- outcomes[i]
      current_outcome_label <- outcome_labels[i]
      cat(":: Running: ", current_outcome, "~ \n")
      
      
      single_forest_output <- try(suppressMessages(make_single_bart(dataset = input_dataset,
                                                                    controls = input_controls,
                                                                    outcome = current_outcome,
                                                                    outcome_label = current_outcome_label,
                                                                    splitsample_df = splitsample_df)))
      
      if(inherits(single_forest_output, 'try-error')) next
      
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
  
  
  # running models ----
  seeds <- round(1000000*runif(n_splitsample_iterations*2)) %>%
    unique() %>% .[1:n_splitsample_iterations]
  
  input_dataset <- master_pool
  input_controls <- controls_with_missingness_dummies
  
  
  
  run_splitsample_models <- function(input_seed){
    
    # set seed
    set.seed(input_seed)
    
    estimation_df <- input_dataset %>% group_by(blocknum, dmatch) %>% sample_frac(splitsample_fraction)
    holdout_df <- input_dataset %>% anti_join(estimation_df, by='row_id') 
    
    # run cf for each outcome
    causal_forest_models <- make_causal_forest_for_each_outcome(estimation_df,
                                                         splitsample_df = holdout_df,
                                                         input_controls,
                                                         outcomes = outcomes_of_interest,
                                                         outcomes_of_interest_labels)
    
    # # run BART for each outcome
    # bart_models <- make_bart_for_each_outcome(estimation_df,
    #                                            splitsample_df = holdout_df,
    #                                            input_controls,
    #                                            outcomes = outcomes_of_interest,
    #                                            outcomes_of_interest_labels,
    #                                            flipped_outcomes)
    
    
    
    models <- list('causal_forest_models' = causal_forest_models)#,
                   #'bart_models' = bart_models)
    
    return(models)
  }
  
  
  
  # run in parallel
  
  start_time <- Sys.time()
  # setting up furr plan
  plan(multisession, workers = 13)
  all_splitsample_models <- future_map(seeds, run_splitsample_models, .options = future_options(scheduling=Inf))
  
  end_time <- Sys.time()
  
  end_time - start_time
  
  all_splitsample_models %>% saveRDS(paste0('all_splitsample_models_final',overall_seed,'.Rds'))
  
}


grf_rds_files <- list.files()[stringr::str_detect(list.files(), 'all_splitsample_models_final')]


all_grf_splitsamples <- NULL

for (grffile in grf_rds_files){
  temp_grffile <- readRDS(grffile)
  all_grf_splitsamples <- all_grf_splitsamples %>% append(temp_grffile)
}


# initiate with the first iteration, and then loop through 2:all grf models
combined_models <- all_grf_splitsamples[[1]]$causal_forest_models

for (i in 2:length(all_grf_splitsamples)){
  current_splitsample_iteration <- all_grf_splitsamples[[i]]$causal_forest_models
  n_outcomes <- length(current_splitsample_iteration)
  
  # for each outcome...
  for (j in 1:n_outcomes){
    # get each dataframe
    n_tables <- length(current_splitsample_iteration[[j]])
    
    for (k in 1:n_tables){
      # rbind the kth table to the same table in the combined model
      combined_models[[j]][[k]] <- bind_rows(combined_models[[j]][[k]],
                                             current_splitsample_iteration[[j]][[k]])
    }
  }
  
}



# now combine all these tables and keep the medians + sd's
for (j in 1:n_outcomes){
  for (k in 1:n_tables){
    combined_models[[j]][[k]] <-
      combined_models[[j]][[k]] %>%
      group_by_at(1) %>%
      summarise_all(list(median = median#,
                         #sd = sd
      )) %>%
      ungroup()
  }
}

saveRDS(combined_models, "combined_final_grf_splitsample_models.Rds")

