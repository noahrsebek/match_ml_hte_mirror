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


get_group_avg_tau_conf_interval <- function(group1, group2){
  group1[["estimate"]] - group2[["estimate"]] +
    c(-1, 1) * qnorm(0.975) * sqrt(group1[["std.err"]]^2 + group2[["std.err"]]^2)
}

clean_estimate_w_95_conf_intervals <- function(est, se){
  paste0(round( est , 3), " +/- ", round( qnorm(0.975) * se , 3))
}


make_quartile_baseline_table <- function(input_augmented_df,
                                         subsample_tau_avgs,
                                         baselines = table_baselines,
                                         baselines_labels = table_baselines_labels){
  quartile_baseline_table <- list()
  
  tau_avgs <- subsample_tau_avgs
  tau_avg_df <- tribble(~tau_quartile, ~avg,
                        1, tau_avgs$quartile_1[['estimate']],
                        2, tau_avgs$quartile_2[['estimate']],
                        3, tau_avgs$quartile_3[['estimate']],
                        4, tau_avgs$quartile_4[['estimate']])
  
  for (quartile in 1:4){
    
    # get N in quartile and get avg tau hat 
    n_quart <- input_augmented_df %>% filter(tau_quartile %in% quartile) %>% nrow()
    avg_tau_hat <- input_augmented_df %>% filter(tau_quartile %in% quartile) %>% pull(grf_tau_hat) %>% mean(na.rm=T)
    #avg_tau_hat <- tau_avg_df %>% filter(tau_quartile %in% quartile) %>% pull(avg)
    
    # make column name
    col_name <- paste0("$\\hat{\\tau}$ Quartile ", quartile)
    temp_quart_column <- c(avg_tau_hat, n_quart)
    
    for (bl_var in baselines){
      mean_val <- input_augmented_df %>% filter(tau_quartile %in% quartile) %>% pull(bl_var) %>% mean(na.rm=T)
      temp_quart_column <- c(temp_quart_column, mean_val)
    }
    quartile_baseline_table[[col_name]] <- temp_quart_column
  }
  
  quartile_baseline_table <- quartile_baseline_table %>% as_tibble()
  quartile_baseline_table <- bind_cols(
    tibble(Baseline = c("\\textit{Mean} $\\hat{\\tau}$",
                        "\\textit{N}",
                        baselines_labels)),
    quartile_baseline_table)
  
  return(quartile_baseline_table)
}



make_ate_summary_table <- function(subsample_tau_avgs,
                                   group_tau_avgs){
  
  subsample_ATEs_mean_se <- subsample_tau_avgs
  
  quartile_95_ci_col <- c()
  quartile_se_col <- c()
  for (i in 4:1){
    current_quartile <- paste0('quartile_', i)
    current_ATE <- subsample_ATEs_mean_se[[current_quartile]]
    quartile_95_ci_col <- c(quartile_95_ci_col,
                            current_ATE[1])
    quartile_se_col <- c(quartile_se_col, current_ATE[[2]])
  }
  
  
  
  
  # Quartile
  ate_95ci_table_pte_quartiles <- tibble('Sample'  = paste("Individual PTE Quartile", 4:1),
                                         "Avg. Treatment Effect" = quartile_95_ci_col,
                                         "Standard Error" = quartile_se_col)
  
  # Overall
  ate_95ci_table_overall <- tibble("Sample" = "Whole Sample",
                                   "Avg. Treatment Effect" = group_tau_avgs['overall.estimate'],
                                   "Standard Error" = group_tau_avgs['overall.std.err'])
  
  # Other subgroups
  other_subgroups <- c("above_median", "below_median", "bottom_three_quartiles")
  other_subgroup_labels <- c("Top 2 PTE Quartiles", "Bottom 2 PTE Quartiles", "Bottom 3 PTE Quartiles")
  
  other_subgroup_95_ci_col <- c()
  other_subgroup_se_col <- c()
  
  for (group in other_subgroups){
    current_ATE <- subsample_ATEs_mean_se[[group]]
    other_subgroup_95_ci_col <- c(other_subgroup_95_ci_col,
                                  current_ATE[1])
    other_subgroup_se_col <- c(other_subgroup_se_col, current_ATE[2])
  }
  
  ate_95ci_table_subgroup <- tibble('Sample'  = other_subgroup_labels,
                                    "Avg. Treatment Effect" = other_subgroup_95_ci_col,
                                    "Standard Error" = other_subgroup_se_col)
  
  
  
  bind_rows(ate_95ci_table_overall,
            ate_95ci_table_pte_quartiles,
            ate_95ci_table_subgroup)
}


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
                                      splitsample_frac = splitsample_fraction
                                      ){
  
  
  if (remove_duplicates==T){
    # get duplicate sids
    duplicate_sids <- dataset %>% 
      group_by(sid) %>% filter(n()>1) %>% 
      distinct(sid) %>% pull(sid)
    
    dataset <- dataset %>% filter(!sid %in% duplicate_sids)
  }
  
  
  # keep observations complete on outcome and that has an inv prob weight
  working_df <- dataset %>% drop_na(outcome, inv_prob_weight) %>% sample_frac(splitsample_frac)
  
  n.obs <- nrow(working_df)
  
  X <- working_df[,controls]
  W <- working_df %>% pull(dmatch)
  Y <- working_df %>% pull(tidyselect::all_of(outcome))
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
                                #sample.fraction = default_sample_fraction,
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
                                #sample.fraction = default_sample_fraction,
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
  
  
  # quartile baseline table
  quartile_heterogeneity_table <- augmented_df %>% make_quartile_baseline_table(subsample_tau_avgs)
  
  # subsample ate table
  subsample_ate_table <- subsample_tau_avgs %>% make_ate_summary_table(group_tau_avgs)
  
  # subsample difference table
  
  # 
  
  
  output_list <- list(
    # 'tau_df' = tau_df,
    # 'augmented_df' = augmented_df,
    #'quartile_heterogeneity_table' = quartile_heterogeneity_table,
    #'tau_rank_plot'=tau_rank_plot,
    'group_tau_avgs'=group_tau_avgs,
    'subsample_tau_avgs'=subsample_tau_avgs,
    'subsample_ate_table' = subsample_ate_table,
    #'forest_object'=tau.forest,
    'calibration_test'=forest_calibration,
    'n_observations' = n.obs,
    'tuning_output' = tau.forest$tuning.output,
    'quartile_heterogeneity_table' = quartile_heterogeneity_table)
  
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


