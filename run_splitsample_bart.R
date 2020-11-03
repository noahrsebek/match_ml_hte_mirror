
overall_seeds = c(09212020, 09242020, 09252020, 09302020, 05101995, 12301990, 05111995, 02271993)

for (overall_seed in overall_seeds){
  
  set.seed(overall_seed)
  #set.seed(09242020)
  #set.seed(09252020)
  
  # loading packages and setting global vars ----
  library('dplyr')
  library('readr')
  library('grf')
  library('tidyr')
  library('magrittr')
  library('ggplot2')
  library('furrr')
  library('purrr')
  library(BART)
  
  
  n_splitsample_iterations <- 25#125
  splitsample_fraction <- 0.5
  to_cluster <- F
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
  
  # overriding macros from grf_globals file
  # --> outcomes of interest (and labels)
  outcomes_of_interest <- c('mathxil_z_post1_np',
                            "math_gpa_full",
                            "math_failures_full",
                            #"mathfailpercent_post1",
                            "nonmathgpa_post1",
                            "readxil_z_post1_np",
                            "graduated_ontime",
                            "graduated_ever",
                            #'treat_post2',
                            "gpa11_math",
                            "eleventh_grade_math_z"
  )
  
  flipped_outcomes <- c()
  
  outcomes_of_interest_labels <- c('Math test score (Z)',
                                   'Math GPA',
                                   "Math Course Failures",
                                   #"Percent of Math Courses Failed",
                                   "Non-math GPA",
                                   "Reading test score (Z)",
                                   "Graduated on-time",
                                   "Graduated ever",
                                   #'Participated in Study 1 Year 2',
                                   "11th Grade Math GPA",
                                   "11th Grade Math Test Score (Z)")
  
  
  master_pool <- load_master_dataset()
  
  # main ML functions ----
  make_single_bart <- function(dataset,
                               controls,
                               outcome,
                               outcome_label,
                               flip_outcome=FALSE,
                               splitsample_df){
    
    
    
    # keep observations complete on outcome and that has an inv prob weight
    working_df <- dataset %>% drop_na(all_of(outcome)) 
    splitsample_df <- splitsample_df %>% drop_na(all_of(outcome))
    
    n.obs <- nrow(working_df)
    
    X <- working_df[,controls]
    W <- working_df %>% pull(dmatch)
    Y <- working_df %>% pull(tidyselect::all_of(outcome))
    #sample_weights <- working_df %>% pull(inv_prob_weight)
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
      bart_fit <- mc.gbart(as.matrix(X_bart), Y, mc.cores = 2, nice=5, type='lbart')
      
      Y.hat_splitsample <- predict(bart_fit, as.matrix(X_bart_splitsample), mc.cores = 2, nice=5)
      Y.hat_splitsample <- Y.hat_splitsample$yhat.test %>% colMeans()
      
      X_bart_YO_splitsample <- X_bart_splitsample %>% mutate(W=0)
      X_bart_Y1_splitsample <- X_bart_splitsample %>% mutate(W=1)
      
      bart_preds_Y0 <- predict(bart_fit, as.matrix(X_bart_YO_splitsample), mc.cores=2, nice=5)
      bart_preds_Y1 <- predict(bart_fit, as.matrix(X_bart_Y1_splitsample), mc.cores=2, nice=5)
      
      bart_preds_Y0 <- bart_preds_Y0$yhat.test %>% colMeans()
      bart_preds_Y1 <- bart_preds_Y1$yhat.test %>% colMeans()
      
      
    } else {
      # contiuous
      bart_fit <- mc.gbart(as.matrix(X_bart), Y, mc.cores = 4, nice=2)
      Y.hat_splitsample <- predict(bart_fit, as.matrix(X_bart_splitsample), mc.cores = 2, nice=2) %>% colMeans()
      
      X_bart_YO_splitsample <- X_bart_splitsample %>% mutate(W=0)
      X_bart_Y1_splitsample <- X_bart_splitsample %>% mutate(W=1)
      
      bart_preds_Y0 <- predict(bart_fit, as.matrix(X_bart_YO_splitsample), mc.cores=2, nice=2) %>% colMeans()
      bart_preds_Y1 <- predict(bart_fit, as.matrix(X_bart_Y1_splitsample), mc.cores=4, nice=2) %>% colMeans()
      
    }
    
    predictions_oob_splitsample <- bart_preds_Y1 - bart_preds_Y0
    
    
    
    
    
    # get forest-wide avg tau's/tx effects (and the standard errpr)
    avg_tx_effect_overall <- avg_effect(Y = Y_splitsample, Y.hat = Y.hat_splitsample, W = W_splitsample, W.hat = W.hat_splitsample,
                                        weights = sample_weights_splitsample, predictions = predictions_oob_splitsample,
                                        subset = NULL)
    
    
    # test calibration of forest
    W.hat_splitsample <- W.hat_splitsample
    
    mean.pred <- weighted.mean(predictions_oob_splitsample, sample_weights_splitsample)
    DF <- data.frame(target = unname(Y_splitsample - Y.hat_splitsample), 
                     mean.forest.prediction = unname(W_splitsample - W.hat_splitsample) * mean.pred,
                     differential.forest.prediction = unname(W_splitsample - W.hat_splitsample) * (predictions_oob_splitsample - mean.pred))
    
    # changed "observation.weight" to "sample_weights"
    best.linear.predictor <- lm(target ~ mean.forest.prediction + 
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
      broom::tidy() %>% select(-statistic) %>% 
      mutate(lower_CI = estimate - 2.24*std.error,
             upper_CI = estimate + 2.24*std.error)
    
    #forest_calibration <- test_calibration(tau.forest) %>% broom::tidy()
    
    # output should be a list with
    # - avg causal effects
    group_tau_avgs <- c('overall'=avg_tx_effect_overall)
    
    # - the predictions + variance in a df with the sids
    tau_df_splitsample <- tibble(sid = cluster_ids_splitsample,
                                 grf_tau_hat = predictions_oob_splitsample) %>% 
      makeDummies()
    
    
    
    
    
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
    #augmented_df <- working_df %>% left_join(tau_df, by='sid')
    augmented_df <- splitsample_df %>% left_join(tau_df_splitsample, by='sid')
    
    
    
    # quartile baseline table
    quartile_heterogeneity_table <- augmented_df %>% make_quartile_baseline_table(subsample_tau_avgs)
    
    # subsample ate table
    subsample_ate_table <- subsample_tau_avgs %>% make_ate_summary_table(group_tau_avgs)
    
    # # subsample difference table
    # subsample_difference_table <- subsample_tau_avgs %>% make_subsample_difference_table()
    
    # calibration test (Jon D) output (no plot)
    calibration_tests_naive <- make_naive_calibration_tests(tau_df_splitsample, outcome)
    
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
      'quartile_heterogeneity_table' = quartile_heterogeneity_table,
      'tau_df_splitsample' = tau_df_splitsample)
    
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

        single_forest_output <- suppressMessages(make_single_bart(dataset = input_dataset,
                                                          controls = input_controls,
                                                          outcome = current_outcome,
                                                          outcome_label = current_outcome_label,
                                                          splitsample_df = splitsample_df))
      
      final_forest_list[[current_outcome_label]] <- single_forest_output
      
    }
    
    return(final_forest_list)
  }
  
  
  
  
  
  # OVERALL DATA WORK ----
 
  
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
  
  master_pool <- master_pool %>% mutate(row_id = row_number())
  
  # add in new grade outcomes, add these to the list
  new_grade_vars <- read_csv("/export/projects/migrated-BAM/Match_AnalysisFiles_May2016/Results/2020/new_grade_variables.csv")
  new_grade_vars <- new_grade_vars %>% select(study, sid,
                                              math_failures_full, math_failures_no_selections,
                                              math_gpa_full,
                                              nonmathcore_gpa_all, nonmathcore_fails_all)
  
  master_pool <- master_pool %>% left_join(new_grade_vars, by=c('sid', 'study'))
  
  # running models ----
  #seeds <- round(1000000*runif(n_splitsample_iterations))
  # get n_splitsample_iterations worth of unique seeds
  seeds <- round(1000000*runif(n_splitsample_iterations*2)) %>%
    unique() %>% .[1:n_splitsample_iterations]
  
  input_dataset <- master_pool
  input_controls <- controls_sans_missingness_dummies
  
  
  
  run_splitsample_models <- function(input_seed){
    
    # set seed
    set.seed(input_seed)
    
    estimation_df <- input_dataset %>% group_by(blocknum, dmatch) %>% sample_frac(splitsample_fraction)
    holdout_df <- input_dataset %>% anti_join(estimation_df, by='row_id') 
    
    # run cf for each outcome
    bart_model <- make_bart_for_each_outcome(estimation_df,
                                                         splitsample_df = holdout_df,
                                                         input_controls,
                                                         outcomes = outcomes_of_interest,
                                                         outcomes_of_interest_labels,
                                                         flipped_outcomes)
    
    
    # do tau comparisons here
    cross_pte_comparisons <- list("bart_model" = bart_model %>% make_cross_pte_comparison_table())
    
    
    # drop tau df's from lists before saving
    for (outcome in outcomes_of_interest_labels){
      bart_model[[outcome]]$tau_df_splitsample <- NULL
    }
    
    
    models <- list('bart_model' = bart_model,
                   'extra_analyses' = list(cross_pte_comparisons))
    
    return(models)
  }
  
  
  
  # run in parallel
  
  start_time <- Sys.time()
  # setting up furr plan
  plan(multisession, workers = 13)
  all_splitsample_models <- future_map(seeds, run_splitsample_models, .options = future_options(scheduling=Inf))
  
  end_time <- Sys.time()
  
  end_time - start_time
  
  all_splitsample_models %>% saveRDS(paste0('all_splitsample_bart_',overall_seed,'.Rds'))
  
}




# combine models ----


models1 <- readRDS('all_splitsample_models1.Rds')
models2 <- readRDS('all_splitsample_models2.Rds')
models3 <- readRDS('all_splitsample_models3.Rds')


all_splitsample_models <- append(models1, models2) %>% append(models3)



combined_models <- all_splitsample_models[[1]]
n_models <- length(combined_models)


for (i in 1:length(all_splitsample_models)){
  # for each splitsample...
  current_splitsample_iteration <- all_splitsample_models[[i]]
  
  # go through each model and append each dataframe to all the other ones
  for (m in 1:n_models){
    # get mth model
    current_model <- current_splitsample_iteration[[m]]
    n_outcomes <- length(current_model)
    # for each outcome...
    for (j in 1:n_outcomes){
      # get each dataframe
      n_tables <- length(current_model[[j]])
      
      for (k in 1:n_tables){
        # rbind the kth table to the same table in the combined model
        combined_models[[m]][[j]][[k]] <- bind_rows(combined_models[[m]][[j]][[k]],
                                                    current_model[[j]][[k]])
      }
    }
  }
}


# now combine all these tables and keep the medians + sd's
for (m in 1:n_models){
  n_outcomes <- length(combined_models[[m]])
  for (j in 1:n_outcomes){
    n_tables <- length(combined_models[[m]][[j]])
    for (k in 1:n_tables){
      combined_models[[m]][[j]][[k]] <-
        combined_models[[m]][[j]][[k]] %>%
        group_by_at(1) %>%
        summarise_all(list(median = median,
                           sd = sd)) %>%
        ungroup()
    }
  }
}

saveRDS(combined_models, "combined_models.Rds")
