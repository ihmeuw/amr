### Launching ST-GPR
rm(list = ls())
central_root <- ADDRESS
setwd(central_root)
source(FILEPATH)
source(FILEPATH)
source(FILEPATH)
source(FILEPATH)
source(FILEPATH)


path_config <- FILEPATH
path_logs <- FILEPATH

config_file <- 'CONFIG_FILE_RUN'

config  <- read.csv(paste0(path_config, 'stgpr_config_', config_file, '.csv'), stringsAsFactors = F)

finalconfig<-c()

for(i in 1:nrow(config)) {
  temp_config <- config[i,]
  write.csv(temp_config,paste0(path_config,'temp_config.csv'),row.names=F) 
  run_id <- register_stgpr_model(path_to_config = paste0(path_config,'temp_config.csv'))
  stgpr_sendoff(run_id, 
              'proj_amr', 
              log_path = path_logs)
  temp_config$run_id <-run_id 
  finalconfig<-rbind(finalconfig,temp_config)
}

finalconfig <- as.data.frame(finalconfig)

for(i in 1:nrow(finalconfig)) {
  finalconfig$success[i] <- check_run(finalconfig$run_id[i])
}

if (unique(finalconfig$success) == 1){
  write.csv(finalconfig, paste0(path_config,'stgpr_config_preamp_', config_file, '.csv'), row.names = FALSE)
} else {
  reruns <- finalconfig[finalconfig$success != 1, !colnames(finalconfig) %in% "success"]
  
  rerunned <- c()
  for(i in 1:nrow(reruns)) {
    temp_config <- reruns[i,]
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

write.csv(finalconfig, paste0(path_config,'stgpr_config_preamp_', config_file, '.csv'), row.names = FALSE)

