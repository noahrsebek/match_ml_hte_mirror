
# Match/Saga ML Treatment Heterogeneity Work


These scripts comprise the current state of the Match ML Tx Heterogeneity work. In brief, we explore 3 different ML methods for analyzing treatment heterogeneity: Bayesian Additive Regression Trees [`BART`], Generalized Random Forests [`GRF`], and X-learner (a meta-learner implemented with GRF) [`XRF`]

To explore heterogeneous treatment effects, we analyze these methods in the context of their full sample performance, as well as using the novel split-sample methods of [Chernozhukov et. al](https://arxiv.org/abs/1712.04802) to assure statistically-valid inference for each method. 



- `run_final_fullsample_models.R` runs the GRF models for each outcome using the full sample plugged into the GRF function (and stores statistics and tables and measures of interest in the `all_fullsample_cf_mdoels.Rds` file)

- `run_final_splitsample_models_nobart.R` runs the GRF models for the splitsample method (parrallelized), and saves each chunk of the runs in the `all_splitsample_models_final*.Rds` files, which are combined into `combined_final_grf_splitsample_models.Rds`

- `causal_forest_memo.Rmd` builds the (rather large) PDF file with all the tables and plots of interest




