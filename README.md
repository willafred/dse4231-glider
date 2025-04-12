# Comparative Study for the use of GLiDeR on Simulated Data (DSE4231)

The table below describes the folders and coding files

| Folder Name | Description|
|-------------|------------|
|old_plots/|Old plots|
|individual_model_results/|Results for individual models|
|plots/|Plots |
|processed_results/|Processed results|
|data_sim_fn.R|Generation of simulated data scenarios|
|glider_fn.R|GLiDeR functions|
|models_fn.R|All functions for models, and also includes the code to process the results to obtain MC ratios|
|visualisations.Rmd|Code to generate plots|

To Note:
1) Naming Convention of Files
overall_5var_scene10_50_simulation.rdata
- It is to note that Scenario 10 has a different number of predictors from the other scenarios
- scene10_50 indicates that the number of predictors in scenario 10 is 50
- In this case, 5var indicates that there are only 5 predictors for all scenarios (other than Scenario 10)

2) Method of saving files
- Given that each element has their own dataframe, we saved the dataframes as a list for easy access by the '$' command to call data name.
- For MC Ratios list, the object is called mc_ratio_var5_sim1000.rdata.
- There is also another object overall_var_5_simulations_results.rdata which refers to the RAW 1000 ate computations. For ease of access by '$' command, it is also saved as a list

3) MADR and backward selection has only 12 scenarios (missing scenario 10) instead of 13 given that Koch et al. (2018) did not do it as well
