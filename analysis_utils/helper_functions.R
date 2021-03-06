

makeDummies <- function(dataframe,
                        flipped=F){
  dataframe$d_positive <- ifelse(dataframe$grf_tau_hat > 0, 1, 0)
  dataframe$d_negative <- ifelse(dataframe$grf_tau_hat <= 0, 1, 0)
  
  dataframe$d_positive <- ifelse(dataframe$grf_tau_hat > 0, 1, 0)
  dataframe$d_negative <- ifelse(dataframe$grf_tau_hat <= 0, 1, 0)
  
  #qcut function taken from Hadley Wickham answer on stackoverflow
  qcut <- function(x, n) {
    cut(x, quantile(x, seq(0, 1, length = n + 1), na.rm = T), labels = seq_len(n),
        include.lowest = TRUE)
  }
  
  dataframe$tau_quartile <- qcut(dataframe$grf_tau_hat, 4) %>% as.numeric()
  dataframe$tau_median_group <- qcut(dataframe$grf_tau_hat, 2) %>% as.numeric()
  
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
    #avg_tau_hat <- input_augmented_df %>% filter(tau_quartile %in% quartile) %>% pull(grf_tau_hat) %>% mean(na.rm=T)
    avg_tau_hat <- tau_avg_df %>% filter(tau_quartile %in% quartile) %>% pull(avg)
    
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


make_above_below_median_baseline_table <- function(input_augmented_df,
                                                   subsample_tau_avgs,
                                                   baselines = table_baselines,
                                                   baselines_labels = table_baselines_labels){
  quantile_baseline_table <- list()
  
  tau_avgs <- subsample_tau_avgs
  tau_avg_df <- tribble(~tau_median_group, ~avg,
                        2, tau_avgs$above_median[['estimate']],
                        1, tau_avgs$below_median[['estimate']])
  
  for (quantile in 1:2){
    
    # get N in quantile and get avg tau hat 
    n_quart <- input_augmented_df %>% filter(tau_median_group %in% quantile) %>% nrow()
    #avg_tau_hat <- input_augmented_df %>% filter(tau_median_group %in% quantile) %>% pull(grf_tau_hat) %>% mean(na.rm=T)
    avg_tau_hat <- tau_avg_df %>% filter(tau_median_group %in% quantile) %>% pull(avg)
    
    # make column name
    col_name <- paste0("$\\hat{\\tau}$ quantile ", quantile)
    temp_quart_column <- c(avg_tau_hat, n_quart)
    
    for (bl_var in baselines){
      mean_val <- input_augmented_df %>% filter(tau_median_group %in% quantile) %>% pull(bl_var) %>% mean(na.rm=T)
      temp_quart_column <- c(temp_quart_column, mean_val)
    }
    quantile_baseline_table[[col_name]] <- temp_quart_column
  }
  
  quantile_baseline_table <- quantile_baseline_table %>% as_tibble()
  quantile_baseline_table <- bind_cols(
    tibble(Baseline = c("\\textit{Mean} $\\hat{\\tau}$",
                        "\\textit{N}",
                        baselines_labels)),
    quantile_baseline_table)
  
  return(quantile_baseline_table)
}



make_above_below_median_summary_table <- function(tau_df){
  # get isr data
  isr_results_2014 <- read_csv("/export/projects/migrated-BAM/Match_AnalysisFiles_May2016/Results/isr_2014.csv")
  isr_results_2015 <- read_csv("/export/projects/migrated-BAM/Match_AnalysisFiles_May2016/Results/isr_2015.csv")
  source('isr_outcomes_and_labels_list.R')
  
  
  
  
  add_wave_1_label <- function(x){paste0(x, '_wave1')}
  add_wave_2_label <- function(x){paste0(x, '_wave2')}
  
  
  # wave 1
  isr_results_2014 <- isr_results_2014 %>%
    rename_with(add_wave_1_label, all_of(isr_outcomes_2014))
  
  isr_outcomes_2014 <- isr_outcomes_2014 %>% add_wave_1_label()
  
  # wave 2
  isr_results_2015 <- isr_results_2015 %>%
    rename_with(add_wave_2_label, all_of(isr_outcomes_2015))
  
  isr_outcomes_2015 <- isr_outcomes_2015 %>% add_wave_2_label()
  
  
  
  
  all_isr_outcomes <- c(isr_outcomes_2014, isr_outcomes_2015)
  all_isr_outcome_labels <- c(isr_outcome_labels_2014, isr_outcome_labels_2015)
  
  all_table_vars <- c(outcomes_of_interest, additional_table_outcomes, table_baselines, all_isr_outcomes)
  all_table_var_labels <- c(outcomes_of_interest_labels, additional_table_outcome_labels, table_baselines_labels, all_isr_outcome_labels)
  
  # merging master dataset with ISR
  table_master_df <- master_pool %>%
    left_join(isr_results_2014 %>% select(sid, all_of(isr_outcomes_2014)), by='sid') %>% 
    left_join(isr_results_2015 %>% select(sid, all_of(isr_outcomes_2015)), by='sid') %>% 
    select(sid, all_of(all_table_vars))
  
  
  # adding in arbitrary tau quantiles to this dataset
  table_master_df <- tau_df %>% as_tibble() %>%
    left_join(table_master_df, by='sid')
  
  
  quantile_baseline_table <- list()
  
  for (quantile in 1:2){
    
    # get N in quantile and get avg tau hat 
    n_quant <- table_master_df %>% filter(tau_median_group %in% quantile) %>% nrow()
    avg_tau_hat <- table_master_df %>% filter(tau_median_group %in% quantile) %>% pull(grf_tau_hat) %>% mean(na.rm=T)
    # avg_tau_hat <- tau_avg_df %>% filter(tau_median_group %in% quantile) %>% pull(avg)
    
    # make column name
    col_name <- paste0("$\\hat{\\tau}$ quantile ", quantile)
    temp_quant_column <- c(avg_tau_hat, n_quant)
    
    for (var in all_table_vars){
      mean_val <- table_master_df %>% filter(tau_median_group %in% quantile) %>% pull(var) %>% mean(na.rm=T)
      temp_quant_column <- c(temp_quant_column, mean_val)
    }
    quantile_baseline_table[[col_name]] <- temp_quant_column
  }
  
  
  quantile_baseline_table <- quantile_baseline_table %>% as_tibble()
  quantile_baseline_table <- bind_cols(
    tibble(Baseline = c("\\textit{Mean} $\\hat{\\tau}$",
                        "\\textit{N}",
                        all_table_var_labels)),
    quantile_baseline_table)
  
  
  return(quantile_baseline_table)
}

make_quartile_summary_table <- function(tau_df){
  # get isr data
  isr_results_2014 <- read_csv("/export/projects/migrated-BAM/Match_AnalysisFiles_May2016/Results/isr_2014.csv")
  isr_results_2015 <- read_csv("/export/projects/migrated-BAM/Match_AnalysisFiles_May2016/Results/isr_2015.csv")
  source('isr_outcomes_and_labels_list.R')
  
  
  
  
  add_wave_1_label <- function(x){paste0(x, '_wave1')}
  add_wave_2_label <- function(x){paste0(x, '_wave2')}
  
  
  # wave 1
  isr_results_2014 <- isr_results_2014 %>%
    rename_with(add_wave_1_label, all_of(isr_outcomes_2014))
  
  isr_outcomes_2014 <- isr_outcomes_2014 %>% add_wave_1_label()
  
  # wave 2
  isr_results_2015 <- isr_results_2015 %>%
    rename_with(add_wave_2_label, all_of(isr_outcomes_2015))
  
  isr_outcomes_2015 <- isr_outcomes_2015 %>% add_wave_2_label()
  
  
  
  
  all_isr_outcomes <- c(isr_outcomes_2014, isr_outcomes_2015)
  all_isr_outcome_labels <- c(isr_outcome_labels_2014, isr_outcome_labels_2015)
  
  all_table_vars <- c(outcomes_of_interest, additional_table_outcomes,table_baselines, all_isr_outcomes)
  all_table_var_labels <- c(outcomes_of_interest_labels, additional_table_outcome_labels, table_baselines_labels, all_isr_outcome_labels)
  
  # merging master dataset with ISR
  table_master_df <- master_pool %>%
    left_join(isr_results_2014 %>% select(sid, all_of(isr_outcomes_2014)), by='sid') %>% 
    left_join(isr_results_2015 %>% select(sid, all_of(isr_outcomes_2015)), by='sid') %>% 
    select(sid, all_of(all_table_vars))
  
  
  # adding in arbitrary tau quantiles to this dataset
  table_master_df <- tau_df %>% as_tibble() %>%
    left_join(table_master_df, by='sid')
  
  
  quantile_baseline_table <- list()
  
  for (quantile in 1:4){
    
    # get N in quantile and get avg tau hat 
    n_quant <- table_master_df %>% filter(tau_quartile %in% quantile) %>% nrow()
    avg_tau_hat <- table_master_df %>% filter(tau_quartile %in% quantile) %>% pull(grf_tau_hat) %>% mean(na.rm=T)
    # avg_tau_hat <- tau_avg_df %>% filter(tau_quartile %in% quantile) %>% pull(avg)
    
    # make column name
    col_name <- paste0("$\\hat{\\tau}$ quantile ", quantile)
    temp_quant_column <- c(avg_tau_hat, n_quant)
    
    for (var in all_table_vars){
      mean_val <- table_master_df %>% filter(tau_quartile %in% quantile) %>% pull(var) %>% mean(na.rm=T)
      temp_quant_column <- c(temp_quant_column, mean_val)
    }
    quantile_baseline_table[[col_name]] <- temp_quant_column
  }
  
  
  quantile_baseline_table <- quantile_baseline_table %>% as_tibble()
  quantile_baseline_table <- bind_cols(
    tibble(Baseline = c("\\textit{Mean} $\\hat{\\tau}$",
                        "\\textit{N}",
                        all_table_var_labels)),
    quantile_baseline_table)
  
  
  return(quantile_baseline_table)
}



avg_effect_simple <- function(predictions, subset){
  # taking simple avg of tau's of a given subsample
  est <- predictions[subset] %>% mean(na.rm=T)
  return(c('estimate'=est))
}


avg_effect <- function(Y,Y.hat,W,W.hat, weights,predictions,subset) {
  
  if (is.null(subset)) {
    subset <- 1:length(Y.hat)
  }
  if (class(subset) == "logical" & length(subset) == length(Y.hat)) {
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





make_ate_summary_table_with_se <- function(subsample_tau_avgs,
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



make_ate_summary_table <- function(subsample_tau_avgs,
                                   group_tau_avgs){
  
  subsample_ATEs_mean_se <- subsample_tau_avgs
  
  quartile_95_ci_col <- c()
  for (i in 4:1){
    current_quartile <- paste0('quartile_', i)
    current_ATE <- subsample_ATEs_mean_se[[current_quartile]]
    quartile_95_ci_col <- c(quartile_95_ci_col,
                            current_ATE[1])
    #quartile_se_col <- c(quartile_se_col, round(current_ATE[[2]], 3))
  }
  
  
  
  
  # Quartile
  ate_95ci_table_pte_quartiles <- tibble('Sample'  = paste("Individual PTE Quartile", 4:1),
                                         "Average Treatment Effect" = quartile_95_ci_col#,
                                         #"Standard Error" = quartile_se_col
  )
  
  # Overall
  ate_95ci_table_overall <- tibble("Sample" = "Whole Sample",
                                   "Average Treatment Effect" = group_tau_avgs['overall.estimate']#,
                                   #input_forest$group_tau_avgs['overall.std.err']),
                                   #"Standard Error" = round(input_forest$group_tau_avgs['overall.std.err'], 3)
  )
  
  # Other subgroups
  other_subgroups <- c("above_median", "below_median", "bottom_three_quartiles")
  other_subgroup_labels <- c("Top 2 PTE Quartiles", "Bottom 2 PTE Quartiles", "Bottom 3 PTE Quartiles")
  
  other_subgroup_95_ci_col <- c()
  other_subgroup_se_col <- c()
  
  for (group in other_subgroups){
    current_ATE <- subsample_ATEs_mean_se[[group]]
    other_subgroup_95_ci_col <- c(other_subgroup_95_ci_col,
                                  current_ATE[1])
    #other_subgroup_se_col <- c(other_subgroup_se_col, round(current_ATE[[2]], 3))
  }
  
  ate_95ci_table_subgroup <- tibble('Sample'  = other_subgroup_labels,
                                    "Average Treatment Effect" = other_subgroup_95_ci_col)
  
  
  
  bind_rows(ate_95ci_table_overall,
            ate_95ci_table_pte_quartiles,
            ate_95ci_table_subgroup)
}



make_subsample_difference_table <- function(subsample_tau_avgs){
  
  highest_quartile_vs_rest_ci <- get_group_avg_tau_conf_interval(
    subsample_tau_avgs$highest_quartile,
    subsample_tau_avgs$bottom_three_quartiles)
  
  above_median_vs_below_ci <- get_group_avg_tau_conf_interval(
    subsample_tau_avgs$above_median, 
    subsample_tau_avgs$below_median)
  
  top_vs_bottom_ci <- get_group_avg_tau_conf_interval(
    subsample_tau_avgs$quartile_4,
    subsample_tau_avgs$quartile_1)
  
  
  
  tribble(
    ~"Group 1", ~"Group 2", ~"95\\% Confidence Interval of Difference (low)", ~"95\\% CI (high)",
    "Quartile 4", "Bottom 3 Quartiles", highest_quartile_vs_rest_ci[1], highest_quartile_vs_rest_ci[2],
    "Top 2 Quartiles (3 \\& 4)", "Bottom 2 Quartiles (1 \\& 2)", above_median_vs_below_ci[1], above_median_vs_below_ci[2],
    "Quartile 4", "Quartile 1", top_vs_bottom_ci[1], top_vs_bottom_ci[2]
  )
}


make_naive_calibration_tests <- function(tau_df, outcome_of_interest,
                                         master=master_pool){
  
  
  # -> add quantiles
  qcut <- function(x, n) {
    cut(x, quantile(x, seq(0, 1, length = n + 1)), labels = seq_len(n),
        include.lowest = TRUE)
  }
  
  n_calibration_quantiles <- 4
  
  calibration_df <- tau_df %>% mutate(calibration_quantiles = grf_tau_hat %>% qcut(n_calibration_quantiles))
  
  
  
  master <- master %>% left_join(calibration_df, by='sid') %>% mutate(blocknum = factor(blocknum))
  
  
  itt_vars_calibration <- c("d13andunder","d14","d15","d16","d17andover","dlearningdisabled","dfreelunch",
                            "dblack","dhispanic","dother","dgrade9","dgrade10","gpa_pre_zeros",
                            "numAs_pre","numBs_pre","numCs_pre","numDs_pre","numFs_pre","missing_gpa_pre",
                            "days_absent_pre_zeros","missing_attend_pre","mathxil_z_pre_np_zeros","readxil_z_pre_np_zeros",
                            "mathxil_z_pre_missing","readxil_z_pre_missing","oss_dis_pre_zeros","incidents_pre_zeros",
                            "any_arrests_pre","violent_pre","property_pre","drug_pre", "blocknum")
  
  # add 
  quantile_itts <- lm(paste0(outcome_of_interest, " ~ ",
                             "dmatch:calibration_quantiles + calibration_quantiles +",
                             paste(itt_vars_calibration, collapse='+')),
                      data=master) %>% broom::tidy() %>% select(-statistic) %>% 
    filter(term %>% startsWith('dmatch:')) %>% 
    mutate( term = term %>% stringr::str_replace('dmatch:calibration_quantiles', "")) %>% 
    rename(calibration_quantiles = term)
  
  # get avg pte for each bucket
  avg_pte_df <- calibration_df %>% group_by(calibration_quantiles) %>% summarise(pte = mean(grf_tau_hat))
  
  
  
  
  # get the statistic from the linear test; this will be the plot subtitle
  calibration_test_statistic <- lm(paste0(outcome_of_interest, " ~ ",
                                          "dmatch + grf_tau_hat + grf_tau_hat*dmatch +",
                                          paste0(itt_vars_calibration, collapse="+")),
                                   data=master) %>% broom::tidy() %>% select(-statistic) %>% 
    filter(term %in% "dmatch:grf_tau_hat")
  
  
  calibration_test_note <- paste0("PTE x Treatment interaction estimate is ",
                                  calibration_test_statistic$estimate %>% round(3), 
                                  " (", calibration_test_statistic$std.error %>% round(3), "), with a p-value of ",
                                  calibration_test_statistic$p.value %>% round(3))
  
  
  
  
  
  return(list(quantile_itts,
              calibration_test_statistic))
  
}



make_calibration_plot <- function(tau_df, outcome_of_interest,
                                  master=master_pool){

  
  # -> add ventiles
  qcut <- function(x, n) {
    cut(x, quantile(x, seq(0, 1, length = n + 1)), labels = seq_len(n),
        include.lowest = TRUE)
  }
  
  n_calibration_quantiles <- 20
  
  calibration_df <- tau_df %>% mutate(calibration_quantiles = grf_tau_hat %>% qcut(n_calibration_quantiles))
  
  
  
  master <- master %>% left_join(calibration_df, by='sid') %>% mutate(blocknum = factor(blocknum))
  
  
  itt_vars_calibration <- c("d13andunder","d14","d15","d16","d17andover","dlearningdisabled","dfreelunch",
                            "dblack","dhispanic","dother","dgrade9","dgrade10","gpa_pre_zeros",
                            "numAs_pre","numBs_pre","numCs_pre","numDs_pre","numFs_pre","missing_gpa_pre",
                            "days_absent_pre_zeros","missing_attend_pre","mathxil_z_pre_np_zeros","readxil_z_pre_np_zeros",
                            "mathxil_z_pre_missing","readxil_z_pre_missing","oss_dis_pre_zeros","incidents_pre_zeros",
                            "any_arrests_pre","violent_pre","property_pre","drug_pre", "blocknum")
  
  # add 
  quantile_itts <- lm(paste0(outcome_of_interest, " ~ ",
                             "dmatch:calibration_quantiles + calibration_quantiles +",
                             paste(itt_vars_calibration, collapse='+')),
                      data=master) %>% broom::tidy() %>% filter(term %>% startsWith('dmatch:')) %>% 
    mutate( term = term %>% stringr::str_replace('dmatch:calibration_quantiles', "")) %>% 
    rename(calibration_quantiles = term)
  
  # get avg pte for each bucket
  avg_pte_df <- calibration_df %>% group_by(calibration_quantiles) %>% summarise(pte = mean(grf_tau_hat))
  
  
  
  
  # get the statistic from the linear test; this will be the plot subtitle
  calibration_test_statistic <- lm(paste0(outcome_of_interest, " ~ ",
                                          "dmatch + grf_tau_hat + grf_tau_hat*dmatch +",
                                          paste0(itt_vars_calibration, collapse="+")),
                                   data=master) %>% broom::tidy() %>% filter(term %in% "dmatch:grf_tau_hat")
  
  
  calibration_test_note <- paste0("PTE x Treatment interaction estimate is ",
                                  calibration_test_statistic$estimate %>% round(3), 
                                  " (", calibration_test_statistic$std.error %>% round(3), "), with a p-value of ",
                                  calibration_test_statistic$p.value %>% round(3))
  
  
  quantile_calibration_plot <- quantile_itts %>% left_join(avg_pte_df) %>%
    ggplot(aes(x=pte, y=estimate)) + 
    geom_point() +
    geom_smooth(method='lm', se=F, color='black', linetype = 'dashed') +
    geom_abline(intercept = 0, slope = 1) + 
    xlab("Average Quantile Individual PTE") + 
    ylab("ITT Estimate for Quantile") + 
    ggtitle(paste0("Naive Calibration Plot: ", outcome_of_interest, "(full sample)"),
            subtitle = calibration_test_note)
  
  
  
  
  
  return(quantile_calibration_plot)
  
}


make_isr_itt_quantile_tables <- function(tau_df,
                                         n_calibration_quantiles=2){
  qcut <- function(x, n) {
    cut(x, quantile(x, seq(0, 1, length = n + 1)), labels = seq_len(n),
        include.lowest = TRUE)
  }
  
  tau_df <- tau_df %>%
    mutate(calibration_quantiles = grf_tau_hat %>% qcut(n_calibration_quantiles))
  
  
  
  make_isr_quartile_table <- function(outcome_list, isr_datafile, raw=F){
    empty_tibble <- NULL
    
    # for a given outcome...
    for (outcome_label in names(outcome_list)){
      
      
      outcome_of_interest <- outcome_list[outcome_label]
      
      fitted_model_df <- lm(paste0(outcome_of_interest, " ~ ",
                                   "dmatch:calibration_quantiles + calibration_quantiles +",
                                   paste(itt_vars, collapse='+')),
                            data=isr_datafile,
                            weights = isr_weight) %>%
        lmtest::coeftest(vcov=sandwich::vcovHC(., type='HC1')) %>% 
        broom::tidy() %>% filter(term %>% startsWith('dmatch:')) %>% 
        mutate( term = term %>% stringr::str_replace('dmatch:calibration_quantiles', "")) %>% 
        rename(calibration_quantiles = term)
      
      
      
      if (raw==F){
        fitted_model_table_row <- fitted_model_df %>% select(-statistic) %>%
          mutate(estimate = round(estimate, 3),
                 std.error = round(std.error, 3),
                 stars = case_when(
                   p.value < 0.01 ~ '***',
                   p.value < 0.05 ~ "**",
                   p.value < 0.01 ~ "*",
                   TRUE ~ ""),
                 Estimate = paste0(estimate, stars, " (", std.error,")"),
                 calibration_quantiles = paste('Quantile', calibration_quantiles)) %>% 
          select(Estimate, calibration_quantiles) %>% 
          tidyr::pivot_wider(values_from = Estimate, names_from = calibration_quantiles) %>% 
          mutate(Question = outcome_label) %>% relocate(Question)
        
        
        
 
       linear_term_model <- lm(paste0(outcome_of_interest, " ~ ",
                                       "calibration_quantiles + dmatch:calibration_quantiles +",
                                       paste(itt_vars, collapse='+')),
                                data=isr_datafile %>% mutate(calibration_quantiles = as.numeric(calibration_quantiles)),
                                weights = isr_weight)  %>%
          lmtest::coeftest(vcov=sandwich::vcovHC(., type='HC1')) %>% 
          broom::tidy() %>% filter(term %>% startsWith('calibration_quantiles:')) %>% 
          select(-statistic) %>%
          mutate(estimate = round(estimate, 3),
                 std.error = round(std.error, 3),
                 stars = case_when(
                   p.value < 0.01 ~ '***',
                   p.value < 0.05 ~ "**",
                   p.value < 0.01 ~ "*",
                   TRUE ~ ""),
                 Estimate = paste0(estimate, stars, " (", std.error,")")) %>%
          select('Linear Quantile Tx Interaction Term' = Estimate)
        
        fitted_model_table_row <- fitted_model_table_row %>% bind_cols(linear_term_model)
        
        
      }
      
      if (raw==T){
        fitted_model_table_row <- fitted_model_df %>% select(-statistic) %>%
          mutate(calibration_quantiles = as.numeric(calibration_quantiles)) %>% 
          select(estimate, calibration_quantiles)
      }
      
      
      
      
      empty_tibble <- empty_tibble %>% bind_rows(fitted_model_table_row)}
    
    return(empty_tibble)
  }
  
  
  
  itt_vars <- c("d13andunder","d14","d15","d16","d17andover","dlearningdisabled","dfreelunch",
                "dblack","dhispanic","dother","dgrade9","dgrade10","gpa_pre_zeros",
                "numAs_pre","numBs_pre","numCs_pre","numDs_pre","numFs_pre","missing_gpa_pre",
                "days_absent_pre_zeros","missing_attend_pre","mathxil_z_pre_np_zeros","readxil_z_pre_np_zeros",
                "mathxil_z_pre_missing","readxil_z_pre_missing","oss_dis_pre_zeros","incidents_pre_zeros",
                "any_arrests_pre","violent_pre","property_pre","drug_pre", "blocknum")
  
  
  master <- load_master_dataset()
  
  master <- master %>% left_join(tau_df, by='sid') %>% mutate(blocknum = as_factor(blocknum))
  
  
  
  # adding isr
  library(readr)
  # get isr data
  isr_results_2014 <- read_csv("/export/projects/migrated-BAM/Match_AnalysisFiles_May2016/Results/isr_2014.csv")
  isr_results_2015 <- read_csv("/export/projects/migrated-BAM/Match_AnalysisFiles_May2016/Results/isr_2015.csv")
  source('isr_outcomes_and_labels_list.R')
  
  isr_outcome_list_2014 <- isr_outcomes_2014
  names(isr_outcome_list_2014) <- isr_outcome_labels_2014
  
  isr_outcome_list_2015 <- isr_outcomes_2015
  names(isr_outcome_list_2015) <- isr_outcome_labels_2015
  
  
  # weights:
  # -> weight post 1 for isr 2014
  # -> weight post 2 for isr 2015
  isr_2014_with_quantiles <- master %>% select(sid, calibration_quantiles, dmatch, all_of(itt_vars)) %>%
    left_join( isr_results_2014 %>% select(sid, isr_weight = weight_post1,
                                           all_of(isr_outcomes_2014)),
               by='sid') %>% 
    filter(!is.na(calibration_quantiles))
  
  
  isr_2015_with_quantiles <- master %>% select(sid, calibration_quantiles, dmatch, all_of(itt_vars)) %>%
    left_join( isr_results_2015 %>% select(sid, isr_weight = weight_post2,
                                           all_of(isr_outcomes_2015)),
               by='sid') %>% 
    filter(!is.na(calibration_quantiles))
  
  
  
  
  
  
  
  wave1_table <- make_isr_quartile_table(isr_outcome_list_2014, isr_2014_with_quantiles)
  wave2_table <- make_isr_quartile_table(isr_outcome_list_2015, isr_2015_with_quantiles)
  
  return(list(wave1_table, wave2_table))
}

make_tau_rank_order_plot <- function(tau_df, outcome_label){
  
  overall_n <- nrow(tau_df)
  
  tau_rank_plot_df <- tau_df %>% arrange(grf_tau_hat) %>% mutate(order=1:n())
  
  tau_rank_plot_summary_df <- tau_rank_plot_df %>% group_by(tau_quartile) %>%
    summarize(min_pred = min(grf_tau_hat),
              min_rank = min(order),
              avg_pred = mean(grf_tau_hat),
              #avg_pred = mean(avg_tau_hat),
              avg_rank = mean(order)) %>%
    mutate(q_name = paste0("Q", tau_quartile, " avg: ", round(avg_pred, 3)))
  
  
  tau_rank_plot <- tau_rank_plot_df %>% ggplot(aes(x=order, y=grf_tau_hat)) +
    geom_line() +
    #geom_line(aes(x=order, y=grf_tau_hat + 1.96*grf_sd), alpha = 0.2) +
    #geom_line(aes(x=order, y=grf_tau_hat - 1.96*grf_sd), alpha = 0.2) +
    geom_vline(data=tau_rank_plot_summary_df, #%>% filter(! tau_quartile %in% 1),
               aes(group=tau_quartile, xintercept = min_rank),
               linetype='dotted', color='gray20') + 
    geom_hline(yintercept = 0, linetype='dashed') +
    ylab(expression("Individual Predicted Treatment Effect"~hat(tau)^(-i)~(X[i]))) +
    xlab('Reverse-rank of sorted Predicted Treatment Effects') +
    theme_minimal() +
    ggtitle(paste0(outcome_label, ": Predicted Treatment Effect vs Rank, N=", overall_n)) +
    geom_text(data=tau_rank_plot_summary_df, aes(label=q_name, y=avg_pred+0.05, x=avg_rank))
  
  return(tau_rank_plot)
}

make_itt_quantile_tables <- function(tau_df,
                                     n_calibration_quantiles=2,
                                     outcome_list = append(outcome_and_label_list,
                                                           additional_table_outcome_label_list)){
  qcut <- function(x, n) {
    cut(x, quantile(x, seq(0, 1, length = n + 1)), labels = seq_len(n),
        include.lowest = TRUE)
  }
  
  tau_df <- tau_df %>%
    mutate(calibration_quantiles = grf_tau_hat %>% qcut(n_calibration_quantiles))
  
  
  
  make_quartile_itt_table <- function(outcome_list, datafile, raw=F){
    empty_tibble <- NULL
    
    # for a given outcome...
    for (outcome_label in names(outcome_list)){
      
      
      outcome_of_interest <- outcome_list[outcome_label]
      
      fitted_model_df <- lm(paste0(outcome_of_interest, " ~ ",
                                   "dmatch:calibration_quantiles + calibration_quantiles +",
                                   paste(itt_vars, collapse='+')),
                            data=datafile) %>%
        lmtest::coeftest(vcov=sandwich::vcovHC(., type='HC1')) %>% 
        broom::tidy() %>% filter(term %>% startsWith('dmatch:')) %>% 
        mutate( term = term %>% stringr::str_replace('dmatch:calibration_quantiles', "")) %>% 
        rename(calibration_quantiles = term)
      
      
      
      if (raw==F){
        fitted_model_table_row <- fitted_model_df %>% select(-statistic) %>%
          mutate(estimate = round(estimate, 3),
                 std.error = round(std.error, 3),
                 stars = case_when(
                   p.value < 0.01 ~ '***',
                   p.value < 0.05 ~ "**",
                   p.value < 0.01 ~ "*",
                   TRUE ~ ""),
                 Estimate = paste0(estimate, stars, " (", std.error,")"),
                 calibration_quantiles = paste('Quantile', calibration_quantiles)) %>% 
          select(Estimate, calibration_quantiles) %>% 
          tidyr::pivot_wider(values_from = Estimate, names_from = calibration_quantiles) %>% 
          mutate(Question = outcome_label) %>% relocate(Question)
        
        
        
        
        linear_term_model <- lm(paste0(outcome_of_interest, " ~ ",
                                       "calibration_quantiles + dmatch:calibration_quantiles +",
                                       paste(itt_vars, collapse='+')),
                                data=datafile %>% mutate(calibration_quantiles = as.numeric(calibration_quantiles)))  %>%
          lmtest::coeftest(vcov=sandwich::vcovHC(., type='HC1')) %>% 
          broom::tidy() %>% filter(term %>% startsWith('calibration_quantiles:')) %>% 
          select(-statistic) %>%
          mutate(estimate = round(estimate, 3),
                 std.error = round(std.error, 3),
                 stars = case_when(
                   p.value < 0.01 ~ '***',
                   p.value < 0.05 ~ "**",
                   p.value < 0.01 ~ "*",
                   TRUE ~ ""),
                 Estimate = paste0(estimate, stars, " (", std.error,")")) %>%
          select('Linear Quantile Tx Interaction Term' = Estimate)
        
        fitted_model_table_row <- fitted_model_table_row %>% bind_cols(linear_term_model)
        
        
      }
      
      if (raw==T){
        fitted_model_table_row <- fitted_model_df %>% select(-statistic) %>%
          mutate(calibration_quantiles = as.numeric(calibration_quantiles)) %>% 
          select(estimate, calibration_quantiles)
      }
      
      
      
      
      empty_tibble <- empty_tibble %>% bind_rows(fitted_model_table_row)}
    
    return(empty_tibble)
  }
  
  
  
  itt_vars <- c("d13andunder","d14","d15","d16","d17andover","dlearningdisabled","dfreelunch",
                "dblack","dhispanic","dother","dgrade9","dgrade10","gpa_pre_zeros",
                "numAs_pre","numBs_pre","numCs_pre","numDs_pre","numFs_pre","missing_gpa_pre",
                "days_absent_pre_zeros","missing_attend_pre","mathxil_z_pre_np_zeros","readxil_z_pre_np_zeros",
                "mathxil_z_pre_missing","readxil_z_pre_missing","oss_dis_pre_zeros","incidents_pre_zeros",
                "any_arrests_pre","violent_pre","property_pre","drug_pre", "blocknum")
  
  
  master <- load_master_dataset()
  
  master <- master %>% left_join(tau_df, by='sid') %>%
    mutate(blocknum = forcats::as_factor(blocknum)) %>% 
    filter(!is.na(calibration_quantiles))
  
  
  
  output_table <- make_quartile_itt_table(outcome_list, master)
  
  
  
  
  
  
  output_table
}

make_cross_pte_comparison_table <- function(one_iteration){
  cross_pte_comparisons <- NULL
  
  for (row in 1:nrow(bivariate_pte_comparisons)){
    outcome_1 <- bivariate_pte_comparisons[[row, 'outcome1']]
    outcome_2 <- bivariate_pte_comparisons[[row, 'outcome2']]
    
    
    outcome_1_label <- which(outcome_and_label_list == outcome_1) %>% names()
    outcome_2_label <- which(outcome_and_label_list == outcome_2) %>% names()
    
    
    # get pte's from both
    single_regression_row <- one_iteration[[outcome_1_label]]$tau_df_splitsample %>% select(sid, grf_tau_hat) %>% left_join(
      one_iteration[[outcome_2_label]]$tau_df_splitsample %>% select(sid, grf_tau_hat),
      by='sid') %>% lm(grf_tau_hat.y ~ grf_tau_hat.x, data=.) %>% 
      lmtest::coeftest(vcov=sandwich::vcovHC(., type="HC1")) %>% broom::tidy() %>% 
      filter(!term %in% '(Intercept)') %>% select(-statistic)
    
    single_regression_row[1,'term'] <- paste(outcome_1_label, "~", outcome_2_label)
    
    cross_pte_comparisons <- bind_rows(cross_pte_comparisons, single_regression_row)
  }
  cross_pte_comparisons %>% rename(Comparison = term)
}


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

