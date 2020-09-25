

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
  
  
  # -> add ventiles
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
  
  
  
  
  
  return(list(quantile_itts,
              calibration_test_statistic))
  
}

