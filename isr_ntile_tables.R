
itt_vars <- c("d13andunder","d14","d15","d16","d17andover","dlearningdisabled","dfreelunch",
              "dblack","dhispanic","dother","dgrade9","dgrade10","gpa_pre_zeros",
              "numAs_pre","numBs_pre","numCs_pre","numDs_pre","numFs_pre","missing_gpa_pre",
              "days_absent_pre_zeros","missing_attend_pre","mathxil_z_pre_np_zeros","readxil_z_pre_np_zeros",
              "mathxil_z_pre_missing","readxil_z_pre_missing","oss_dis_pre_zeros","incidents_pre_zeros",
              "any_arrests_pre","violent_pre","property_pre","drug_pre", "blocknum")

n_calibration_quantiles = 4
forest <- final_forests_missingness[[2]]
qcut <- function(x, n) {
  cut(x, quantile(x, seq(0, 1, length = n + 1)), labels = seq_len(n),
      include.lowest = TRUE)
}

# for a given forest...
tau_df <- forest$tau_df %>% mutate(calibration_quantiles = grf_tau_hat %>% qcut(n_calibration_quantiles))

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
    }
    
    if (raw==T){
      fitted_model_table_row <- fitted_model_df %>% select(-statistic) %>%
        mutate(calibration_quantiles = as.numeric(calibration_quantiles)) %>% 
        select(estimate, calibration_quantiles)
    }
    
    
    
    
    empty_tibble <- empty_tibble %>% bind_rows(fitted_model_table_row)}
  
  return(empty_tibble)
}



make_isr_quartile_table(isr_outcome_list_2014, isr_2014_with_quantiles) %>% write_csv('mathgpa_wave1.csv')
make_isr_quartile_table(isr_outcome_list_2015, isr_2015_with_quantiles) %>% write_csv('mathgpa_wave2.csv')





raw_wave1_table <- make_isr_quartile_table(isr_outcome_list_2014, isr_2014_with_quantiles, raw=T)
raw_wave2_table <- make_isr_quartile_table(isr_outcome_list_2015, isr_2015_with_quantiles, raw=T)



raw_wave1_table %>%
  ggplot(aes(x=calibration_quantiles, y=estimate)) + geom_point(alpha=0.3) + geom_smooth(method='lm') + ggtitle('ISR Wave 1')

raw_wave2_table %>%
  ggplot(aes(x=calibration_quantiles, y=estimate)) + geom_point(alpha=0.3) + geom_smooth(method='lm') + ggtitle('ISR Wave 2')


# make labels and names a named list
# -> values are 


