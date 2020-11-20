# combining all predictions into one dataset
load_master_dataset <- function(x){0}
insertSource('/export/projects/migrated-BAM/Match_AnalysisFiles_May2016/Results/2020/grf/grf_globals.R',
             functions='load_master_dataset')
load_master_dataset <- load_master_dataset@.Data

master_pool <- load_master_dataset(addvars=F)

# get isr data
isr_results_2014 <- readr::read_csv("/export/projects/migrated-BAM/Match_AnalysisFiles_May2016/Results/isr_2014.csv")
isr_results_2015 <- readr::read_csv("/export/projects/migrated-BAM/Match_AnalysisFiles_May2016/Results/isr_2015.csv")
source('isr_outcomes_and_labels_list.R')


cuts <- 2


qcut <- function(x, n) {
  cut(x, quantile(x, seq(0, 1, length = n + 1)), labels = seq_len(n),
      include.lowest = TRUE)
}


# outcome_1 <- names(fullsample_models)[1]
# outcome_2 <- names(fullsample_models)[2]



make_qq_summary_table <- function(combined_df,
                                  new_avg_colname_1,
                                  new_avg_colname_2,
                                  outcome_1,
                                  outcome_2){
  
  
  
  
  
  
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
  
  all_table_vars <- c(outcomes_of_interest, table_baselines)#, all_isr_outcomes)
  all_table_var_labels <- c(outcomes_of_interest_labels, table_baselines_labels, all_isr_outcome_labels)
  
  # merging master dataset with ISR
  table_master_df <- combined_df %>%
    left_join(isr_results_2014 %>% select(sid, all_of(isr_outcomes_2014), weight_post1), by='sid') %>% 
    left_join(isr_results_2015 %>% select(sid, all_of(isr_outcomes_2015), weight_post2), by='sid') %>% 
    select(sid, 'QQ_group',
           !!sym(new_avg_colname_1),
           !!sym(new_avg_colname_2),
           all_of(all_table_vars),
           all_of(all_isr_outcomes), weight_post1, weight_post2)
  
  
  
  
  
  qq_summary_table <- list()
  
  for (qq_group in sort(unique(table_master_df$'QQ_group'))){
    
    # get N in qq_group and get avg tau hat 
    n_quant <- table_master_df %>% filter(QQ_group %in% qq_group) %>% nrow()
    avg_tau_hat_1 <- table_master_df %>% filter(QQ_group %in% qq_group) %>% pull(!!sym(new_avg_colname_1)) %>% mean(na.rm=T)
    avg_tau_hat_2 <- table_master_df %>% filter(QQ_group %in% qq_group) %>% pull(!!sym(new_avg_colname_2)) %>% mean(na.rm=T)
    
    # avg_tau_hat <- tau_avg_df %>% filter(QQ_group %in% qq_group) %>% pull(avg)
    
    # make column name
    col_name <- paste0("X/Y group ", qq_group)
    temp_quant_column <- c(avg_tau_hat_1, avg_tau_hat_2, n_quant)
    
    for (var in all_table_vars){
      mean_val <- table_master_df %>% filter(QQ_group %in% qq_group) %>% pull(var) %>% mean(na.rm=T)
      temp_quant_column <- c(temp_quant_column, mean_val)
    }
    
    for (var in isr_outcomes_2014){
      mean_val_values <- table_master_df %>% filter(QQ_group %in% qq_group) %>% pull(var) 
      mean_val_weights <- table_master_df %>% filter(QQ_group %in% qq_group) %>% pull(weight_post1)
      mean_val <- weighted.mean(mean_val_values, mean_val_weights, na.rm=T)
      temp_quant_column <- c(temp_quant_column, mean_val)
    }
    
    for (var in isr_outcomes_2015){
      mean_val_values <- table_master_df %>% filter(QQ_group %in% qq_group) %>% pull(var) 
      mean_val_weights <- table_master_df %>% filter(QQ_group %in% qq_group) %>% pull(weight_post2)
      mean_val <- weighted.mean(mean_val_values, mean_val_weights, na.rm=T)
      temp_quant_column <- c(temp_quant_column, mean_val)
    }
    
    qq_summary_table[[col_name]] <- temp_quant_column
  }
  
  
  qq_summary_table <- qq_summary_table %>% as_tibble()
  qq_summary_table <- bind_cols(
    tibble(Baseline = c(paste0("\\textit{Mean} $\\hat{\\tau}$ on outcome X: ", outcome_1),
                        paste0("\\textit{Mean} $\\hat{\\tau}$ on outcome Y:", outcome_2),
                        "\\textit{N}",
                        all_table_var_labels)),
    qq_summary_table)
  
  return(qq_summary_table)
}





make_cluster_scatterplot <- function(combined_df,
                                     new_avg_colname_1, new_avg_colname_2, outcome_1, outcome_2
                                     ){
  find_hull <- function(df, var1=new_avg_colname_1, var2=new_avg_colname_2){
    df[chull(df[,var1], df[,var2]), ]}
  
  # picking optimal number of clusters
  library(NbClust)
  best_clusters <- NbClust(data=combined_df[, c(new_avg_colname_1, new_avg_colname_2)],
                           min.nc = 2, max.nc = 5, method='kmeans', index='all')
  
  # 
  # kmeans_output <- combined_df[, c(new_avg_colname_1, new_avg_colname_2)] %>%
  #   kmeans(4, iter.max = 100)
  # cluster_df <- combined_df %>% bind_cols(clusters=as.factor(kmeans_output$cluster))
  cluster_df <- combined_df %>% bind_cols(clusters = as.factor(best_clusters$Best.partition))
  
  cluster_hulls <- plyr::ddply(cluster_df, 'clusters', find_hull)
  
  cluster_summary_df <- cluster_df %>% group_by(clusters) %>% 
    summarise(n=n(),
              avg_x = mean(!!sym(new_avg_colname_1)),
              avg_y = mean(!!sym(new_avg_colname_2)),
              x_y_tuple = paste0("(", round(avg_x,3), ", ", round(avg_y,3), "), N=", n))
  
  
  cluster_scatterplot <- cluster_df %>% ggplot(aes(x=!!sym(new_avg_colname_1),
                                                   y=!!sym(new_avg_colname_2),
                                                   color = clusters)) +
    geom_point(stroke=0, alpha=0.4) +
    geom_polygon(data=cluster_hulls, alpha=.1, linetype='dashed') +
    ggtitle(paste0("PTE Cluster Scatterplot: X=",
                   outcome_1, " vs Y=", outcome_2)) +
    geom_text(data=cluster_summary_df, aes(label=x_y_tuple, y=avg_y, x=avg_x, color=NULL))
  
  
  
  # make summary plot based on cluster group
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
  
  all_table_vars <- c(outcomes_of_interest, table_baselines)#, all_isr_outcomes)
  all_table_var_labels <- c(outcomes_of_interest_labels, table_baselines_labels, all_isr_outcome_labels)
  
  # merging master dataset with ISR
  table_master_df <- cluster_df %>%
    left_join(isr_results_2014 %>% select(sid, all_of(isr_outcomes_2014), weight_post1), by='sid') %>% 
    left_join(isr_results_2015 %>% select(sid, all_of(isr_outcomes_2015), weight_post2), by='sid') %>% 
    select(sid, clusters,
           !!sym(new_avg_colname_1),
           !!sym(new_avg_colname_2),
           all_of(all_table_vars),
           all_of(all_isr_outcomes), weight_post1, weight_post2)
  
  
  
  cluster_summary_table <- list()
  
  for (cluster_group in sort(unique(table_master_df$clusters))){
    
    # get N in cluster_group and get avg tau hat 
    n_quant <- table_master_df %>% filter(clusters %in% cluster_group) %>% nrow()
    avg_tau_hat_1 <- table_master_df %>% filter(clusters %in% cluster_group) %>% pull(!!sym(new_avg_colname_1)) %>% mean(na.rm=T)
    avg_tau_hat_2 <- table_master_df %>% filter(clusters %in% cluster_group) %>% pull(!!sym(new_avg_colname_2)) %>% mean(na.rm=T)
    
    # avg_tau_hat <- tau_avg_df %>% filter(clusters %in% cluster_group) %>% pull(avg)
    
    # make column name
    col_name <- paste0("cluster ", cluster_group)
    temp_quant_column <- c(avg_tau_hat_1, avg_tau_hat_2, n_quant)
    
    for (var in all_table_vars){
      mean_val <- table_master_df %>% filter(clusters %in% cluster_group) %>% pull(var) %>% mean(na.rm=T)
      temp_quant_column <- c(temp_quant_column, mean_val)
    }
    
    for (var in isr_outcomes_2014){
      mean_val_values <- table_master_df %>% filter(clusters %in% cluster_group) %>% pull(var) 
      mean_val_weights <- table_master_df %>% filter(clusters %in% cluster_group) %>% pull(weight_post1)
      mean_val <- weighted.mean(mean_val_values, mean_val_weights, na.rm=T)
      temp_quant_column <- c(temp_quant_column, mean_val)
    }
    
    for (var in isr_outcomes_2015){
      mean_val_values <- table_master_df %>% filter(clusters %in% cluster_group) %>% pull(var) 
      mean_val_weights <- table_master_df %>% filter(clusters %in% cluster_group) %>% pull(weight_post2)
      mean_val <- weighted.mean(mean_val_values, mean_val_weights, na.rm=T)
      temp_quant_column <- c(temp_quant_column, mean_val)
    }
    
    cluster_summary_table[[col_name]] <- temp_quant_column
  }
  
  
  cluster_summary_table <- cluster_summary_table %>% as_tibble()
  cluster_summary_table <- bind_cols(
    tibble(Baseline = c(paste0("\\textit{Mean} $\\hat{\\tau}$ on outcome X: ", outcome_1),
                        paste0("\\textit{Mean} $\\hat{\\tau}$ on outcome Y:", outcome_2),
                        "\\textit{N}",
                        all_table_var_labels)),
    cluster_summary_table)
  
  cluster_itt_table <- make_itt_group_tables(cluster_df, 'clusters')
  
  return(list(cluster_scatterplot, cluster_summary_table, cluster_itt_table))
}


make_itt_group_tables <- function(input_df,
                                  group_col,
                                  outcome_list = outcome_and_label_list){
  
  
  
  
  
  itt_vars <- c("d13andunder","d14","d15","d16","d17andover","dlearningdisabled","dfreelunch",
                "dblack","dhispanic","dother","dgrade9","dgrade10","gpa_pre_zeros",
                "numAs_pre","numBs_pre","numCs_pre","numDs_pre","numFs_pre","missing_gpa_pre",
                "days_absent_pre_zeros","missing_attend_pre","mathxil_z_pre_np_zeros","readxil_z_pre_np_zeros",
                "mathxil_z_pre_missing","readxil_z_pre_missing","oss_dis_pre_zeros","incidents_pre_zeros",
                "any_arrests_pre","violent_pre","property_pre","drug_pre", "blocknum")
  
  
  master <- load_master_dataset()
  
  master <- input_df %>% select(sid, group_col) %>% 
    left_join(master, by=c('sid')) %>%
    mutate(blocknum = forcats::as_factor(blocknum)) %>%
    filter(!is.na(!!sym(group_col))) %>% as_tibble()
  
  
  empty_tibble <- NULL
  
  # for a given outcome...
  for (outcome_label in names(outcome_list)){
    
    
    outcome_of_interest <- outcome_list[outcome_label]
    
    fitted_model_df <- lm(paste0(outcome_of_interest,
                                 " ~ dmatch:", group_col, " + ",
                                 group_col, " + ",
                                 paste(itt_vars, collapse='+')),
                          data=master) %>%
      lmtest::coeftest(vcov=sandwich::vcovHC(., type='HC1')) %>% 
      broom::tidy() %>% filter(term %>% startsWith('dmatch:')) %>% 
      mutate( term = term %>% stringr::str_replace(paste0('dmatch:', group_col), "")) %>% 
      rename(!!sym(group_col) := term)
    
    
    fitted_model_table_row <- fitted_model_df %>% select(-statistic) %>%
      mutate(estimate = round(estimate, 3),
             std.error = round(std.error, 3),
             stars = case_when(
               p.value < 0.01 ~ '***',
               p.value < 0.05 ~ "**",
               p.value < 0.01 ~ "*",
               TRUE ~ ""),
             Estimate = paste0(estimate, stars, " (", std.error,")"),
             !!sym(group_col) := paste('Group',  !!sym(group_col))) %>% 
      select(Estimate, !!sym(group_col)) %>% 
      tidyr::pivot_wider(values_from = Estimate, names_from = group_col) %>% 
      mutate(Outcome = outcome_label) %>% relocate(Outcome)
    
    
    
    empty_tibble <- empty_tibble %>% bind_rows(fitted_model_table_row)
  }
  
  return(empty_tibble) 
}


make_multi_pte_comparison_plots <- function(outcome_1,
                                         outcome_2){
  
  new_avg_colname_1 <- paste0(outcome_1, 
                              ", Estimate")
  
  new_se_colname_1 <- paste0(outcome_1, 
                             ", SE")
  
  new_quantile_colname_1 <- paste0(outcome_1, 
                                 ", PTE Quantile")
  
  new_avg_colname_2 <- paste0(outcome_2, 
                              ", Estimate")
  
  new_se_colname_2 <- paste0(outcome_2, 
                             ", SE")
  
  new_quantile_colname_2 <- paste0(outcome_2, 
                                 ", PTE Quantile")
  
  tau_df_1 <- fullsample_models[[outcome_1]]$tau_df %>%
    select(sid,
           !!new_avg_colname_1 := grf_tau_hat,
           !!new_se_colname_1 := grf_sd)
  
  tau_df_2 <- fullsample_models[[outcome_2]]$tau_df %>%
    select(sid,
           !!new_avg_colname_2 := grf_tau_hat,
           !!new_se_colname_2 := grf_sd)
  
  combined_df <- inner_join(tau_df_1, tau_df_2, by='sid') %>% 
    mutate(!!new_quantile_colname_1 := qcut(!!sym(new_avg_colname_1), cuts),
           !!new_quantile_colname_2 := qcut(!!sym(new_avg_colname_2), cuts),
           'QQ_group' = paste0(!!sym(new_quantile_colname_1), "x", !!sym(new_quantile_colname_2)))
    
  
  
  combined_df <- combined_df %>% left_join(master_pool, by='sid') %>% as_tibble()
  
  augmented_df <- combined_df %>%
    mutate(compliance_dummy = ifelse(dmatch==1, treat_post1, NA)) %>% 
    group_by(!!sym(new_quantile_colname_1),
             !!sym(new_quantile_colname_2)) %>% 
    summarise(Percent = n()/nrow(combined_df),
              Compliance = mean(compliance_dummy, na.rm=T),
              "% Black" = mean(dblack),
              "Average Age" = mean(age_pre),
              "Baseline Math Test Score" = mean(mathxil_z_pre_np, na.rm=T),
              "Missing baseline math test" = mean(mathxil_z_pre_missing),
              avg_x = mean(!!sym(new_avg_colname_1)),
              avg_y = mean(!!sym(new_avg_colname_2)),
              x_y_tuple = paste0("(", round(avg_x,3), ", ", round(avg_y,3), ")"))
  
  
  
  
  # scatterplot: find convex hull of each point set
  find_hull <- function(df, var1=new_avg_colname_1, var2=new_avg_colname_2){
    df[chull(df[,var1], df[,var2]), ]}
  
  hulls <- plyr::ddply(combined_df, "QQ_group", find_hull)
  
  pte_two_var_scatterplot <- combined_df %>% ggplot(aes(x=!!sym(new_avg_colname_1),
                                                        y=!!sym(new_avg_colname_2),
                                                        color = QQ_group)
  ) +
    geom_point(stroke=0, alpha=0.3) + geom_rug(alpha=.1, color='black') +
    geom_polygon(data=hulls, alpha=.1, linetype='dashed') +
    theme(legend.position = 'none') +
    ggtitle(paste0("PTE Scatterplot: X=",
                   outcome_1, " vs Y=", outcome_2)) +
    geom_text(data=augmented_df, aes(label=x_y_tuple, y=avg_y, x=avg_x, color=NULL))
  
  # adding in density plot version of the above
  # --> make a plot with the counts and percentages for each group
  qq_summary_df <- combined_df %>%
    group_by(QQ_group) %>%
    summarize(N = n(),
              !!new_avg_colname_1 := median(!!sym(new_avg_colname_1)),
              !!new_avg_colname_2 := median(!!sym(new_avg_colname_2))
    ) %>%
    ungroup() %>% 
    mutate(pct = N/sum(N),
           Percent = paste0(round(pct*100, 1), "%"),
           label = paste0(QQ_group, ": ", N, ", ", Percent))
  
  
  
  pte_two_var_densityplot <- combined_df %>% ggplot(aes(!!sym(new_avg_colname_2),
                                                        !!sym(new_avg_colname_1))) +
    #stat_density_2d(aes(fill = ..level..), geom = "polygon") +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      legend.position='none'
    ) +
    ggrepel::geom_text_repel(data=qq_summary_df, aes(label=label), color = 'white', size=3) + 
    geom_rug(alpha=.1, color='black') +
    geom_polygon(aes(color = QQ_group),
                 data=hulls, alpha=.1, linetype='dashed') + theme(legend.position = 'none') + ggtitle(paste0("PTE Densities: ",
                                                                                                             outcome_1, " vs ", outcome_2))
  
  
  qq_summary_table <- combined_df %>% make_qq_summary_table(new_avg_colname_1, new_avg_colname_2, outcome_1, outcome_2)
  
  # making ITT table on other outcomes
  combined_df$nugroupvar = combined_df$QQ_group
  qq_itt_table <- make_itt_group_tables(combined_df, 'nugroupvar')

  
  cluster_output <- combined_df %>% make_cluster_scatterplot(new_avg_colname_1, new_avg_colname_2, outcome_1, outcome_2)
  
  
  qq_table_name <- paste0("Outcome X: ", outcome_1,
                          ", Outcome Y: ", outcome_2,
                          " above/below median group summary table")
  
  qq_itt_table_name <- paste0("Outcome X: ", outcome_1,
                          ", Outcome Y: ", outcome_2,
                          " above/below median group ITT table")
  
  cluster_table_name <- paste0("Outcome X: ", outcome_1,
                          ", Outcome Y: ", outcome_2,
                          " cluster summary table")
  cluster_itt_table_name <- paste0("Outcome X: ", outcome_1,
                               ", Outcome Y: ", outcome_2,
                               " cluster ITT table")
  
  outlist <- list(pte_two_var_scatterplot, pte_two_var_densityplot, cluster_output[[1]])
  
  outlist[[qq_table_name]] <- qq_summary_table
  outlist[[qq_itt_table_name]] <- qq_itt_table
  outlist[[cluster_table_name]] <- cluster_output[[2]]
  outlist[[cluster_itt_table_name]] <- cluster_output[[3]]
  
  
  outlist
}



names(fullsample_models)

all_multi_pte_comparisons <- list(
  make_multi_pte_comparison_plots("Math test score (Z)", "Math GPA"),
  make_multi_pte_comparison_plots("Math test score (Z)", "Math Course Failures"),
  make_multi_pte_comparison_plots("Math test score (Z)", "Non-Math Core GPA (all non-math core courses)"),
  make_multi_pte_comparison_plots("Math GPA", "Non-Math Core GPA (all non-math core courses)"),
  make_multi_pte_comparison_plots("Math Course Failures", "Non-Math Core GPA (all non-math core courses)"))


