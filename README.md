README file for code used for "Random forest for survival data: which methods work best and under what conditions?" (2024) by Berkowitz et al

1. The first step is to run the code in the "Combined Data Simulation.R" file, which uses the .csv file "flc.csv". This file requires the 'survival' and  'flexsurv' packages.

2. Next, with all simulated data in your R environment, open "Full_Fact.R", which fits all the forest methods and computes the three error metrics. This file requires several packages listed at the top of the R document.

I ran the code on 1/3 of the FLCs at a time, separately for FLCs with: 1) No ExtraVars; 2) Noisy variables; and 3) Correlated variables.

See bottom of R file for how I assembled the results.

3. The "models_and_plots.R" reads in the simulation study results from the "Full_results_merged.csv" file. 
