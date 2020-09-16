

# we want each core to run ALL outcomes





for each iteration:
  randomly split sample (block/treatment stratfication)
  for each method:
    for each outcome:
      run method(outcome)  
  


seeds <- round(1000000*runif(n_splitsample_iterations))

# deduping and prepping data for split
input_dataset <- master_pool
input_controls <- controls_sans_missingness_dummies


duplicate_sids <- input_dataset %>% 
  group_by(sid) %>% filter(n()>1) %>% 
  distinct(sid) %>% pull(sid)

input_dataset <- input_dataset %>% filter(!sid %in% duplicate_sids)

input_dataset <- input_dataset %>%  drop_na(inv_prob_weight) %>%  mutate(row_id = row_number())





for (i in 1:n_iter){
  
  # set seed
  set.seed(i)
  
  estimation_df <- master_pool %>% sample_frac(splitsample_fraction)
  holdout_df <- dataset %>% anti_join(working_df, by='row_id') 
  
  # run cf for each outcome
  causal_forest <- make_causal_forest_for_each_outcome(input_dataset,
                                      input_controls,
                                      outcomes_of_interest,
                                      outcomes_of_interest_labels,
                                      flipped_outcomes)
  
  # run xrf for each outcome
  
  
}






# generate list of random seeds
