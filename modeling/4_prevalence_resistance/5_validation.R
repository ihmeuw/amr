rm(list = ls())
library(tidyverse)
library(readxl)
options(stringsAsFactors = FALSE)
repo <- FILEPATH
data_path <-FILEPATH

source(FILEPATH)
source(FILEPATH)

add_loc_info <- function(df){
  locs = get_location_metadata(location_set_id = 35, 
                               gbd_round_id = 6)[, c("location_id", "location_name", "region_name", "super_region_name")]
  df <- merge(df, locs, by = 'location_id', all.x = TRUE)
  return(df)
}


# generate confusion matrix to calculate accuracy
get_confusion_matrix <- function(dt) {
  dt[, pred_resistant := gpr_mean * sample_size]
  dt[, TP := pmin(resistant, pred_resistant)]
  dt[, TN := pmin(susceptible, sample_size - pred_resistant)]
  dt[, FP := pmax(0, pred_resistant - resistant)]
  dt[, FN := pmax(0, sample_size - pred_resistant - susceptible)]
  
  confusion_matrix <- dt[, lapply(.SD, sum), .SDcols = c("TP", "FP", "FN", "TN")]
  return(confusion_matrix)
}

get_bin_class_metrics <- function(cm) {
  cm <- copy(cm)
  cm[, sensitivity := TP / (TP + FN)] 
  cm[, specificity := TN / (TN + FP)]
  cm[, accuracy := (TP + TN) / (TP + TN + FP + FN)]
  return(cm)
}

pub.dir <- FILEPATH

if(!dir.exists(paste0(pub.dir, '/plots_insample_validity/'))){
  dir.create(paste0(pub.dir, '/plots_insample_validity/'), recursive = TRUE)
}

pathogen_meta <- read.csv(FILEPATH)
drug_meta <- read.csv(FILEPATH)
file_runs <- as.data.frame(read.csv(FILEPATH))
file_runs$combo <- paste0(file_runs$pathogen,"-",file_runs$abx_class)

location_md2 <- subset(get_location_metadata(location_set_id=35, gbd_round_id=6),select = c(location_id,level,super_region_id, region_id,ihme_loc_id,parent_id), droplevels = TRUE)

##### In-sample Validation Metrics #####
is_accuracy <- c()
is_summary <- c()

for (i in 1:nrow(file_runs)){
  gpr <- model_load(file_runs$resistance_run_id[i], 'gpr')
  moddata <- model_load(file_runs$resistance_run_id[i], 'data')
  
  if (file_runs$pathogen[i] == 'acinetobacter_baumanii' & file_runs$abx_class[i] == 'carbapenem'){
    moddata <- subset(moddata, sample_size < 100000)
  }
  
  preds <- merge(moddata, gpr, by = c('location_id', 'year_id', 'sex_id', 'age_group_id'))
  preds$pathogen <- file_runs$pathogen[i]
  preds$abx_class <- file_runs$abx_class[i]
  
  preds$gpr_var <- ((preds$gpr_mean-preds$gpr_upper)/1.96)^2
  preds$adj_var <- preds$gpr_var+preds$variance
  preds$gpr_lower_adj <- pmax(0, preds$gpr_mean - 1.96*sqrt(preds$adj_var))
  preds$gpr_upper_adj <- pmin(1, preds$gpr_mean + 1.96*sqrt(preds$adj_var))
  
  bugdrugsum <- preds %>% dplyr::group_by(pathogen, abx_class) %>%
    dplyr::summarize(
      bias = sum(gpr_mean-data)/length(data),
      rmse = sqrt(sum((gpr_mean-data)^2)/length(data)),
      coverage = length(data[data <= gpr_upper_adj & data >= gpr_lower_adj])/length(data),
      old_coverage = length(data[data <= gpr_upper & data >= gpr_lower])/length(data)
    )
  
  bugdrugsum$data_points <- nrow(moddata)
  bugdrugsum$n_locations <- length(unique(moddata$location_id))
  
  is_summary <- bind_rows(is_summary, bugdrugsum)
  
  preds <- as.data.table(preds)
  preds$resistant <- (preds$data*preds$sample_size)
  preds$susceptible <- preds$sample_size-preds$resistant
  cm <- get_confusion_matrix(preds)
  metrics <- get_bin_class_metrics(cm)
  metrics$pathogen <- file_runs$pathogen[i]
  metrics$abx_class <- file_runs$abx_class[i]
  
  is_accuracy <- bind_rows(is_accuracy, metrics)
}

is_summary <- merge(is_accuracy[, c('pathogen', 'abx_class', 'sensitivity', 'specificity', 'accuracy')],
                     is_summary, by = c('pathogen', 'abx_class'))

write.csv(is_summary, paste0(pub.dir, 'IS_validation_metrics'), row.names = F)

is_accuracy <- is_accuracy[is_accuracy$abx_class != 'prop_retreated',]
is_accuracy$sub_desig[grepl('retreated', is_accuracy$abx_class)] <- 'retreated'
is_accuracy$sub_desig[grepl('new', is_accuracy$abx_class)] <- 'new'
is_accuracy$abx_class <- gsub('_retreated', '', is_accuracy$abx_class)
is_accuracy$abx_class <- gsub('_new', '', is_accuracy$abx_class)
is_accuracy2 <- merge(is_accuracy, pathogen_meta[,c('pathogen', 'pathogen_name_long')], on = 'pathogen')
is_accuracy2 <- as.data.frame(is_accuracy2)
is_accuracy2 <- merge(is_accuracy2, drug_meta[,c('abx_class', 'abx_class_name_short')], on = 'abx_class')

is_accuracy2$pathogen <- as.character(is_accuracy2$pathogen_name_long)
is_accuracy2$pathogen[!is.na(is_accuracy2$sub_desig)] <- paste0(is_accuracy2$pathogen_name_long[!is.na(is_accuracy2$sub_desig)],
                                                                ' (', is_accuracy2$sub_desig[!is.na(is_accuracy2$sub_desig)], ')')
is_accuracy2$accuracy <- round(is_accuracy2$accuracy, 3)

##### Out-of-sample Validation Metrics #####
OOSconfig_base <- read.csv(FILEPATH)
OOSconfig_base <- OOSconfig_base[!OOSconfig_base$description %in% c(
  'escherichia_coli-sulfa', 'enterococcus_spp-fluoroquinolone'),]

OOSconfig_rerun <- read.csv(FILEPATH)
finalconfig <- bind_rows(OOSconfig_base, OOSconfig_rerun)
rm(OOSconfig_base, OOSconfig_rerun)

ooscombos <- unique(finalconfig$description)

oos_summary <- c()
oos_metrics <- c()
for (combo in ooscombos){
  bug <- str_split(combo, '-')[[1]][1]
  drug <- str_split(combo, '-')[[1]][2]
  input <- read.csv(paste0(FILEPATH, bug, '/', drug, '/input_w_ho.csv'))
  
  if (combo == 'acinetobacter_baumanii-carbapenem'){
    input <- input[input$sample_size < 100000, ]
  }
  
  allpreds <- c()
  
  for (holdout in unique(finalconfig$holdout)){
    holdouts <- input[input$master_fold_id == holdout,]
    gpr <- model_load(finalconfig$run_id[finalconfig$holdout == holdout & finalconfig$description == combo], 'gpr')
    pred_set <- merge(holdouts, gpr, by = c('location_id', 'year_id', 'sex_id'))
    stg1 <- model_load(finalconfig$run_id[finalconfig$holdout == holdout & finalconfig$description == combo], 'stage1')
    stg1$stage1 <- inv.logit(stg1$stage1)
    pred_set <- merge(pred_set, stg1[, c('location_id', 'year_id', 'sex_id', 'stage1')], by = c('location_id', 'year_id', 'sex_id'))
    pred_set$pathogen <- bug
    pred_set$abx_class <- drug
    allpreds <- bind_rows(allpreds, pred_set)
  }
  
  allpreds$gpr_var <- ((allpreds$gpr_mean-allpreds$gpr_upper)/1.96)^2
  allpreds$adj_var <- allpreds$gpr_var+allpreds$variance
  allpreds$gpr_lower_adj <- pmax(0, allpreds$gpr_mean - 1.96*sqrt(allpreds$adj_var))
  allpreds$gpr_upper_adj <- pmin(1, allpreds$gpr_mean + 1.96*sqrt(allpreds$adj_var))
  
  bugdrugsum <- allpreds %>% dplyr::group_by(pathogen, abx_class) %>%
    dplyr::summarize(
      bias = sum(gpr_mean-val)/length(val),
      bias_stg1 = sum(stage1-val)/length(val),
      rmse = sqrt(sum((gpr_mean-val)^2)/length(val)),
      rmse_stg1 = sqrt(sum((stage1-val)^2)/length(val)),
      coverage = length(val[val <= gpr_upper_adj & val >= gpr_lower_adj])/length(val),
      old_coverage = length(val[val <= gpr_upper & val >= gpr_lower])/length(val)
    )
  
  bugdrugsum$data_points <- nrow(input)
  bugdrugsum$n_locations <- length(unique(input$location_id))
  
  oos_summary <- bind_rows(oos_summary, bugdrugsum)
  
  allpreds <- as.data.table(allpreds)
  allpreds$resistant <- (allpreds$val*allpreds$sample_size)
  allpreds$susceptible <- allpreds$sample_size-allpreds$resistant
  cm <- get_confusion_matrix(allpreds)
  metrics <- get_bin_class_metrics(cm)
  metrics$pathogen <- bug
  metrics$abx_class <- drug
  oos_metrics <- bind_rows(oos_metrics, metrics)
}

oos_summary <- merge(oos_metrics[, c('pathogen', 'abx_class', 'sensitivity', 'specificity', 'accuracy')],
                     oos_summary, by = c('pathogen', 'abx_class'))

write.csv(oos_summary, paste0(pub.dir, 'OOS_validation_metrics.csv'), row.names = F)
