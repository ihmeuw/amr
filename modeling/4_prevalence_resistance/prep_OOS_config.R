rm(list = ls())
hdrive <- ADDRESS
setwd(hdrive)
library(tidyverse)

bugdrugs <- read.csv(FILEPATH)
bugdrugs <- bugdrugs[bugdrugs$NAME_exclude == 0 & bugdrugs$modeler == 'IHME',]
combos <- paste0(bugdrugs$pathogen, '-', bugdrugs$abx_class)

config_table <- c()
for(x in combos) {
  y <- substring(x,1,stringr::str_locate(x,"-")-1)[1]
  z <- substring(x,stringr::str_locate(x,"-")+1,)[1]
  description <- paste0(x)
  data <- read.csv(paste0(FILEPATH,y,'/', z,'/input_w_ho.csv'))
  for (holdout_num in 1:10) {
    ho_data <- data[data$master_fold_id != holdout_num,]
    ho_data$master_fold_id <- NULL
    write.csv(ho_data, paste0(FILEPATH,y,'/', z,'/holdout_',
                              holdout_num, '/ho_input.csv'), row.names = F)
    path_to_data <-  paste0(FILEPATH,y,'/', z,'/holdout_',
                            holdout_num, '/ho_input.csv')
    path_to_custom_stage_1 <-  paste0(FILEPATH,y,'/', z,'/holdout_',
                                      holdout_num, '/custom_stage1_df.csv') 
    
    location_set_id <- 35
    year_end <- 2018
    st_omega <- '1,1'
    st_zeta <- '0.075,0.001'
    st_lambda <- '0.4,0.4'
    prediction_age_group_ids <- 22
    prediction_sex_ids <- 3
    data_transform <- 'logit'
    decomp_step <- 'iterative'
    prediction_units <-'proportion resistance'
    modelable_entity_id <- 24920
    gpr_scale <- '20,20'
    gbd_round_id <- 6
    gpr_amp_method <- 'broken_stick'
    gpr_draws <- 1000
    holdout <- holdout_num
    
    year_start <- min(ho_data$year_id)
    by_loc_sum = group_by(ho_data, location_id) %>%
      summarize(
        num_dp = length(val)
      )
    density_cutoffs <- max(2, ceiling(quantile(by_loc_sum$num_dp, 0.25))[[1]])
    
    config <- c(modelable_entity_id,gbd_round_id,decomp_step,data_transform,prediction_units,st_lambda,st_omega,st_zeta,gpr_scale,path_to_data,path_to_custom_stage_1,location_set_id,year_start,year_end,prediction_age_group_ids,prediction_sex_ids,description,holdout,gpr_amp_method,gpr_draws,density_cutoffs)
    config_table <- rbind(config_table,config)
  }
}
  
  
  

config_table <- as.data.frame(config_table)
names(config_table) <- c('modelable_entity_id','gbd_round_id','decomp_step','data_transform','prediction_units','st_lambda','st_omega','st_zeta','gpr_scale','path_to_data','path_to_custom_stage_1','location_set_id','year_start','year_end','prediction_age_group_ids','prediction_sex_ids','description','holdout','gpr_amp_method','gpr_draws','density_cutoffs')
write.csv(config_table,FILEPATH,row.names = F)
