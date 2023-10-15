## Getting results for AMR pathogen paper

The process includes several steps.
1. Run mcod_by_pathogen.py to create the mcod_by_pathogen.csv to be used in the gbd adjustments 
2. Getting adjustment values for several GBD causes using the gbd_fix_draws.py
3. Running save_combine_gbd_adjusted.py to combine all adjustments into one file
4. Getting H.pylori attributable stomach cancer values by running hpylori_stomach_cancer_draws.R and stomach_cancer_adj.py
5. Getting all the results with get_burden_estimates. 
6. Running checks for the results and outputting them 
