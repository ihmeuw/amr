rm(list = ls())
source(FILEPATH)
pacman::p_load(data.table, openxlsx, dplyr, ggplot2, gridExtra, argparse)
options(stringsAsFactors = FALSE)

location_md2 <- subset(get_location_metadata(location_set_id=35, gbd_round_id=6),select = c(location_id,level,super_region_id, region_id,ihme_loc_id,parent_id), droplevels = TRUE)

path_config <- FILEPATH
path_logs <- FILEPATH

config_file <- 'CONFIG_FILE_RUN'

config <- read.csv(paste0(path_config,'stgpr_config_preamp_', config_file, '.csv'))

stopifnot(unique(config$success) == 1)

for (i in 1:nrow(config)){
  amp <- model_load(config$run_id[i], 'amp_nsv')
  min_amp <- min(amp$st_amp)
  amp_factor <- (min_amp+0.05)/min_amp
  
  config$gpr_amp_factor[i] <- amp_factor 
}

finalconfig <-c()

central_root <- FILEPATH
setwd(central_root)
source(FILEPATH)
source(FILEPATH)

for (i in 1:nrow(config)){
  iter_temp <- config[i,]
  write.csv(iter_temp, paste0(path_config, 'temp_iterative_config.csv'),row.names=F) 
  run_id <- register_stgpr_model(path_to_config = paste0(path_config, 'temp_iterative_config.csv'))
  stgpr_sendoff(run_id, 
                'proj_amr', 
                log_path=path_logs)
  iter_temp$run_id <- run_id
  finalconfig<-rbind(finalconfig,iter_temp)
}

for(i in 1:nrow(finalconfig)) {
  finalconfig$success[i] <- check_run(finalconfig$run_id[i])
}

if (unique(finalconfig$success) == 1){
  write.csv(finalconfig, paste0(path_config,'stgpr_config_w_pt5_amp_', config_file, '.csv'), row.names = FALSE)
} else {
  reruns <- finalconfig[finalconfig$success != 1, !colnames(finalconfig) %in% "success"]
  
  rerunned <- c()
  for(i in 1:nrow(reruns)) {
    write.csv(temp_config,paste0(path_config,'temp_config.csv'),row.names=F) 
    run_id <- register_stgpr_model(path_to_config = paste0(path_config,'temp_config.csv'))
    stgpr_sendoff(run_id, 
                  'proj_amr', 
                  log_path=path_logs)
    temp_config$run_id <-run_id 
    rerunned<-rbind(rerunned,temp_config)
  }
  
  for(i in 1:nrow(rerunned)) {
    rerunned$success[i] <- check_run(rerunned$run_id[i])
  }
}

finalconfig <- finalconfig[finalconfig$success == 1,]
finalconfig <- rbind(finalconfig, rerunned)

write.csv(finalconfig, paste0(path_config,'stgpr_config_w_pt5_amp_', config_file, '.csv'), row.names = FALSE)
