# combining all predictions into one dataset



cuts <- 2


qcut <- function(x, n) {
  cut(x, quantile(x, seq(0, 1, length = n + 1)), labels = seq_len(n),
      include.lowest = TRUE)
}


# outcome_1 <- names(all_fullsample_cf_models)[1]
# outcome_2 <- names(all_fullsample_cf_models)[2]



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
  
  tau_df_1 <- all_fullsample_cf_models[[outcome_1]]$tau_df %>%
    select(sid,
           !!new_avg_colname_1 := grf_tau_hat,
           !!new_se_colname_1 := grf_sd)
  
  tau_df_2 <- all_fullsample_cf_models[[outcome_2]]$tau_df %>%
    select(sid,
           !!new_avg_colname_2 := grf_tau_hat,
           !!new_se_colname_2 := grf_sd)
  
  combined_df <- inner_join(tau_df_1, tau_df_2, by='sid') %>% 
    mutate(!!new_quantile_colname_1 := qcut(!!sym(new_avg_colname_1), cuts),
           !!new_quantile_colname_2 := qcut(!!sym(new_avg_colname_2), cuts),
           'QQ Group' = paste0(!!sym(new_quantile_colname_1), "x", !!sym(new_quantile_colname_2)))
    
  
  
  combined_df <- combined_df %>% left_join(master_pool, by='sid')
  
  augmented_df <- combined_df %>%
    mutate(compliance_dummy = ifelse(dmatch==1, treat_post1, NA)) %>% 
    group_by(!!sym(new_quantile_colname_1),
             !!sym(new_quantile_colname_2)) %>% 
    summarise(Percent = n()/nrow(combined_df),
              Compliance = mean(compliance_dummy, na.rm=T),
              "% Black" = mean(dblack),
              "Average Age" = mean(age_pre),
              "Baseline Math Test Score" = mean(mathxil_z_pre_np, na.rm=T),
              "Missing baseline math test" = mean(mathxil_z_pre_missing))
  
  
  
  
  # scatterplot: find convex hull of each point set
  find_hull <- function(df, var1=new_avg_colname_1, var2=new_avg_colname_2){
    df[chull(df[,var1], df[,var2]), ]}
  
  hulls <- plyr::ddply(combined_df, "`QQ Group`", find_hull)
  
  pte_two_var_scatterplot <- combined_df %>% ggplot(aes(!!sym(new_avg_colname_2),
                                                        !!sym(new_avg_colname_1),
                                                        color = `QQ Group`)
                                                    ) +
    geom_point(stroke=0, alpha=0.3) + geom_rug(alpha=.1, color='black') +
    geom_polygon(data=hulls, alpha=.1, linetype='dashed') + theme(legend.position = 'none') + ggtitle(paste0("PTE Scatterplot: ",
                                                                                                             outcome_1, " vs ", outcome_2))
  
  # adding in density plot version of the above
  # --> make a plot with the counts and percentages for each group
  qq_summary_df <- combined_df %>%
    group_by(`QQ Group`) %>%
    summarize(N = n(),
              !!new_avg_colname_1 := median(!!sym(new_avg_colname_1)),
              !!new_avg_colname_2 := median(!!sym(new_avg_colname_2))
    ) %>%
    ungroup() %>% 
    mutate(pct = N/sum(N),
           Percent = paste0(round(pct*100, 1), "%"),
           label = paste0(`QQ Group`, ": ", N, ", ", Percent))
  
  
  
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
    geom_polygon(aes(color = `QQ Group`),
                 data=hulls, alpha=.1, linetype='dashed') + theme(legend.position = 'none') + ggtitle(paste0("PTE Densities: ",
                                                                                                             outcome_1, " vs ", outcome_2))
  
  
  
  outlist <- list(pte_two_var_scatterplot, pte_two_var_densityplot)
  
  outlist
}



names(all_fullsample_cf_models)

all_pte_quantile_plots <- list(
  make_multi_pte_comparison_plots("Math test score (Z)", "Math GPA"),
  make_multi_pte_comparison_plots("Math test score (Z)", "Math Course Failures"),
  make_multi_pte_comparison_plots("Math GPA", "Math Course Failures"),
  
  make_multi_pte_comparison_plots("Math test score (Z)", "Non-math GPA"),
  make_multi_pte_comparison_plots("Math GPA", "Non-math GPA"),
  make_multi_pte_comparison_plots("Math Course Failures", "Non-math GPA"),
  
  make_multi_pte_comparison_plots("Math test score (Z)", "Graduated on-time"),
  make_multi_pte_comparison_plots("Math GPA", "Graduated on-time"),
  make_multi_pte_comparison_plots("Math Course Failures", "Graduated on-time"),
  make_multi_pte_comparison_plots("Non-math GPA", "Graduated on-time"),
  
  make_multi_pte_comparison_plots("Math test score (Z)", "11th Grade Math GPA"),
  make_multi_pte_comparison_plots("Math GPA", "11th Grade Math GPA"),
  make_multi_pte_comparison_plots("Math Course Failures", "11th Grade Math GPA"),
  make_multi_pte_comparison_plots("Non-math GPA", "11th Grade Math GPA"),
  make_multi_pte_comparison_plots("11th Grade Math GPA", "Graduated on-time"),
  
  make_multi_pte_comparison_plots("11th Grade Math Test Score (Z)", "Math test score (Z)"),
  make_multi_pte_comparison_plots("11th Grade Math Test Score (Z)", "Graduated on-time"),
  make_multi_pte_comparison_plots("11th Grade Math Test Score (Z)", "11th Grade Math GPA")
  #
  
  )


