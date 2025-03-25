# dse4231-glider

1) Naming convention: overall_5var_scene10_50_simulation.rdata
scene 10 has a diff number of predictors from the other scenarios
scene10 50 means it has 50 predictors

2) we saved the data as a list cos each element has their own dataframe, so you can access it by using $ to call data name. for mc ratios list (this image), the obj is called mc_ratio_var5_sim1000.rdata. there is also another obj overall_var_5_simulations_results.rdata which is the RAW 1000 ate computations so if yall need can use it also! this raw data is also saved as a list so u can access the model's results more easily

3) MADR and backward selection has only 12 scenarios (missing scenario 10) since paper didnt do also