# Code-for-River-silica-concentrations-rise-with-warming-temperatures-potential-geochemical-mechanisms
01_download_data.R is code to download silica data from the Gages dataset 
02_build_cq_master.R compiles silica and discharge data by site and common attributes from the gages dataset
03_smk_analysis.R is code to compute the seasonal mann kendall test over static and moving windows for all sites 
04_geology_assignment.R is code to extract geology by area for each watershed and sort into simplified groups
05_prism_tair_smk.R is code to extract air temperature timeseries by area for each watershed and compute trends
06_lag_analysis.R compares air temperature trends to silica trends
07_figure1.R code to plot figure 1
08_figure2.R code to plot figure 2
09_figure3.R code to plot figure 3

manuscript_models_post.zip is Crunchflow code and model outputs for all six reactive transport scenearios in a zipped folder. For ease, the six input files are also posted which correspond with Tabls S1-3: run1_T_cT.in, run2_T_dt_small.in, run3_T_dt_big.in, run4_K_cT.in, run5_K_dT_small.in, run6_K_dt.in 
