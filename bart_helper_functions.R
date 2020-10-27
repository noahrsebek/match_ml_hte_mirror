
make_average_effect_list <- function(tau.hat, Y.hat){
  ate.highest_25 <- avg_effect(Y = Y, Y.hat = Y.hat, W = W, W.hat = W.hat,
                               weights = sample_weights, predictions = tau.hat,
                               subset = tau_df$tau_quartile == 4)
  
  ate.bottom_25 <- avg_effect(Y = Y, Y.hat = Y.hat, W = W, W.hat = W.hat,
                              weights = sample_weights, predictions = tau.hat,
                              subset = tau_df$tau_quartile == 1)
  
  ate.quartile_2 <- avg_effect(Y = Y, Y.hat = Y.hat, W = W, W.hat = W.hat,
                               weights = sample_weights, predictions = tau.hat,
                               subset = tau_df$tau_quartile == 2)
  
  ate.quartile_3 <- avg_effect(Y = Y, Y.hat = Y.hat, W = W, W.hat = W.hat,
                               weights = sample_weights, predictions = tau.hat,
                               subset = tau_df$tau_quartile == 3)
  
  ate.bottom_75 <- avg_effect(Y = Y, Y.hat = Y.hat, W = W, W.hat = W.hat,
                              weights = sample_weights, predictions = tau.hat,
                              subset = tau_df$tau_quartile != 4)
  
  ate.high <- avg_effect(Y = Y, Y.hat = Y.hat, W = W, W.hat = W.hat,
                         weights = sample_weights, predictions = tau.hat,
                         subset = tau_df$tau_quartile %in% c(3,4))
  
  
  ate.low <- avg_effect(Y = Y, Y.hat = Y.hat, W = W, W.hat = W.hat,
                        weights = sample_weights, predictions = tau.hat,
                        subset = tau_df$tau_quartile %in% c(1,2))
  
  
  
  subsample_tau_avgs <- list('highest_quartile' = ate.highest_25,
                             'bottom_three_quartiles' = ate.bottom_75,
                             'above_median'=ate.high,
                             'below_median'=ate.low,
                             'quartile_1' = ate.bottom_25,
                             'quartile_2' = ate.quartile_2,
                             'quartile_3' = ate.quartile_3,
                             'quartile_4' = ate.highest_25)
  
  subsample_tau_avgs
}

make_calibration_table <- function(tau.hat, Y.hat){
  mean.pred <- mean(tau.hat)
  
  DF <- data.frame(target = unname(Y - Y.hat), 
                   mean.forest.prediction = unname(W - W.hat) * mean.pred,
                   differential.forest.prediction = unname(W - W.hat) * (tau.hat - mean.pred))
  
  # changed "observation.weight" to "sample_weights"
  best.linear.predictor <- lm(target ~ mean.forest.prediction + 
                                differential.forest.prediction + 0, weights = sample_weights, 
                              data = DF)
  
  
  blp.summary <- lmtest::coeftest(best.linear.predictor,
                                  vcov=sandwich::vcovHC(best.linear.predictor, type="HC1"))
  
  dimnames(blp.summary)[[2]][4] <- gsub("[|]", "", dimnames(blp.summary)[[2]][4])
  blp.summary[, 4] <- ifelse(blp.summary[, 3] < 0, 1 - blp.summary[, 
                                                                   4]/2, blp.summary[, 4]/2)
  forest_calibration <- blp.summary %>%
    broom::tidy() %>% 
    mutate(lower_CI = estimate - 1.96*std.error,
           upper_CI = estimate + 1.96*std.error) %>% 
    select(-statistic, std.error)
  
  forest_calibration
}


get_quantile_025 <- function(x){
  quantile(x, 0.025)
}
get_quantile_975 <- function(x){
  quantile(x, 0.975)
}




add_CIs_to_mashup_tables <- function(x){
  x %>% group_by_at(1) %>%
    summarize_all(list('SE' = sd,
                       'Low CI: 2.5%' = get_quantile_025,
                       'High CI: 97.5%' = get_quantile_975)) %>% ungroup
}  





