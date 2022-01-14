#####################################################
# In-sample validation for MRBRT relative risk models
#####################################################
library(dplyr)
library(readxl)
library(metafor)
library(ggplot2)
library(data.table)
library(mrbrt002, lib.loc = "FILENAME")
source("FILENAME")
source("FILENAME")
get_pval <- function(beta, beta_sd, one_sided = FALSE) {
  zscore <- abs(beta/beta_sd)
  if (one_sided) 1 - pnorm(zscore)
  if (!one_sided) (1 - pnorm(zscore))*2
}
# root mean square error
rmse <- function(m, o){
  sqrt(mean(abs(m - o)^2, na.rm = T))
}
# mean absolute error
mae <- function(m, o) {
  mean(abs(m - o), na.rm = T)
}

# Define paths and metadata
outpath <- 'FILEPATH'  
inpath <- 'FILEPATH'
pub.dir <- paste0(outpath, "FILEPATH")
username <- Sys.getenv("USER")
repo_path <- sprintf("FILEPATH", username)
abx_meta <- read.csv(paste0(repo_path,"FILEPATH"), stringsAsFactors = F) 
pathogen_meta <- read.csv(paste0(repo_path,"FILEPATH"), stringsAsFactors = F)
location_meta <- subset(get_location_metadata(location_set_id=35, gbd_round_id=6),select = c(location_id,level,super_region_name,super_region_id, region_id,ihme_loc_id,parent_id), droplevels = TRUE)
file_runs <- as.data.frame(read.csv(paste0(repo_path,"FILEPATH"), stringsAsFactors = F)) %>% filter(exclude==0)
file_runs$combo <- paste0(file_runs$pathogen,"-",file_runs$abx_class)

# define combinations to evaluate
file_runs<- file_runs[!grepl("new",file_runs$combo),]
file_runs<- file_runs[!grepl("ret",file_runs$combo),]
file_runs$abx_class[grepl("mdr",file_runs$abx_class)] <- 'sulfa'

#read in predicted data
summary_sterile <- read.csv(paste0(inpath,'FILENAME'), stringsAsFactors = FALSE)  %>% 
  mutate(sterile = 1, pred_var = ((mean_RR - exp(ub))/1.96)^2,
         pred_lb = exp(lb), pred_ub = exp(ub)) %>%
  select(sterile, pathogen, abx_class, mean_RR, pred_var, pred_ub, pred_lb) %>% 
  rename(predicted = mean_RR)

# Load observed data
  d2 <- read.csv(paste0(inpath,'FILENAME'), stringsAsFactors = FALSE) 

#Quantities to be used in the model
  d2$lnrr<-d2$yi
  d2$lnrr_se<-d2$sei
  d2$lnrr_var<-d2$vi
  d2$location_id[grep('REF1',d2$refvar)] <- 73

# Replace adjusted with crude estimate when the adjusted estimate is inplausible 
d2$yi[d2$sterile==0 & d2$pathogen == 'streptococcus_pneumoniae' & d2$refvar=='DATA' & d2$drug=='beta_lactamase_inhibitor'] <- d2$crude_lnrr[d2$sterile==0 & d2$pathogen == 'streptococcus_pneumoniae' & d2$refvar=='DATA' & d2$drug=='beta_lactamase_inhibitor']
d2$yi[d2$sterile==1 & d2$pathogen == 'streptococcus_pneumoniae' & d2$refvar=='DATA' & d2$drug=='methicillin'] <- d2$crude_lnrr[d2$sterile==1 & d2$pathogen == 'streptococcus_pneumoniae' & d2$refvar=='DATA' & d2$drug=='methicillin']
d2$yi[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='methicillin'] <- d2$crude_lnrr[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='methicillin']
d2$yi[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='sulfa'] <- d2$crude_lnrr[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='sulfa']
d2$yi[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='anti_pseudomonal_penicillin'] <- d2$crude_lnrr[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='anti_pseudomonal_penicillin']
d2$yi[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='fourth_gen_ceph'] <- d2$crude_lnrr[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='fourth_gen_ceph']
d2$yi[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='anti_pseudomonal_penicillin'] <- d2$crude_lnrr[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='anti_pseudomonal_penicillin']
d2$yi[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='carbapenem'] <- d2$crude_lnrr[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='carbapenem']
d2$yi[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='fluoroquinolone'] <- d2$crude_lnrr[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='fluoroquinolone']
d2$yi[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='fourth_gen_ceph'] <- d2$crude_lnrr[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='fourth_gen_ceph']
d2$yi[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='sulfa'] <- d2$crude_lnrr[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='sulfa']
d2$yi[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='third_gen_ceph'] <- d2$crude_lnrr[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='third_gen_ceph']
d2$yi[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='methicillin'] <- d2$crude_lnrr[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='methicillin']
d2$yi[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='fluoroquinolone'] <- d2$crude_lnrr[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='fluoroquinolone']
d2$yi[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='sulfa'] <- d2$crude_lnrr[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='sulfa']
d2$yi[d2$sterile==1 & d2$pathogen == 'coagulase_negative_staph' & d2$refvar=='DATA' & d2$drug=='sulfa'] <- d2$crude_lnrr[d2$sterile==1 & d2$pathogen == 'coagulase_negative_staph' & d2$refvar=='DATA' & d2$drug=='sulfa']
d2$yi[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='third_gen_ceph'] <- d2$crude_lnrr[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='third_gen_ceph']
d2$yi[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='aminopenicillin'] <- d2$crude_lnrr[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='aminopenicillin']
d2$yi[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='aminopenicillin'] <- d2$crude_lnrr[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='aminopenicillin']
d2$yi[d2$sterile==1 & d2$pathogen == 'staphylococcus_aureus' & d2$refvar=='DATA' & d2$drug=='aminoglycoside'] <- d2$crude_lnrr[d2$sterile==1 & d2$pathogen == 'staphylococcus_aureus' & d2$refvar=='DATA' & d2$drug=='aminoglycoside']
d2$sei[d2$sterile==0 & d2$pathogen == 'streptococcus_pneumoniae' & d2$refvar=='DATA' & d2$drug=='beta_lactamase_inhibitor'] <- d2$crude_lnrr_se[d2$sterile==0 & d2$pathogen == 'streptococcus_pneumoniae' & d2$refvar=='DATA' & d2$drug=='beta_lactamase_inhibitor']
d2$sei[d2$sterile==1 & d2$pathogen == 'streptococcus_pneumoniae' & d2$refvar=='DATA' & d2$drug=='methicillin'] <- d2$crude_lnrr_se[d2$sterile==1 & d2$pathogen == 'streptococcus_pneumoniae' & d2$refvar=='DATA' & d2$drug=='methicillin']
d2$sei[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='methicillin'] <- d2$crude_lnrr_se[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='methicillin']
d2$sei[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='sulfa'] <- d2$crude_lnrr_se[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='sulfa']
d2$sei[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='anti_pseudomonal_penicillin'] <- d2$crude_lnrr_se[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='anti_pseudomonal_penicillin']
d2$sei[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='fourth_gen_ceph'] <- d2$crude_lnrr_se[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='fourth_gen_ceph']
d2$sei[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='anti_pseudomonal_penicillin'] <- d2$crude_lnrr_se[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='anti_pseudomonal_penicillin']
d2$sei[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='carbapenem'] <- d2$crude_lnrr_se[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='carbapenem']
d2$sei[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='fluoroquinolone'] <- d2$crude_lnrr_se[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='fluoroquinolone']
d2$sei[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='fourth_gen_ceph'] <- d2$crude_lnrr_se[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='fourth_gen_ceph']
d2$sei[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='sulfa'] <- d2$crude_lnrr_se[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='sulfa']
d2$sei[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='third_gen_ceph'] <- d2$crude_lnrr_se[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='third_gen_ceph']
d2$sei[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='methicillin'] <- d2$crude_lnrr_se[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='methicillin']
d2$sei[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='fluoroquinolone'] <- d2$crude_lnrr_se[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='fluoroquinolone']
d2$sei[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='sulfa'] <- d2$crude_lnrr_se[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='sulfa']
d2$sei[d2$sterile==1 & d2$pathogen == 'coagulase_negative_staph' & d2$refvar=='DATA' & d2$drug=='sulfa'] <- d2$crude_lnrr_se[d2$sterile==1 & d2$pathogen == 'coagulase_negative_staph' & d2$refvar=='DATA' & d2$drug=='sulfa']
d2$sei[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='third_gen_ceph'] <- d2$crude_lnrr_se[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='third_gen_ceph']
d2$sei[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='aminopenicillin'] <- d2$crude_lnrr_se[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='aminopenicillin']
d2$sei[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='aminopenicillin'] <- d2$crude_lnrr_se[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='aminopenicillin']
d2$sei[d2$sterile==1 & d2$pathogen == 'staphylococcus_aureus' & d2$refvar=='DATA' & d2$drug=='aminoglycoside'] <- d2$crude_lnrr_se[d2$sterile==1 & d2$pathogen == 'staphylococcus_aureus' & d2$refvar=='DATA' & d2$drug=='aminoglycoside']
d2$vi[d2$sterile==0 & d2$pathogen == 'streptococcus_pneumoniae' & d2$refvar=='DATA' & d2$drug=='beta_lactamase_inhibitor'] <- d2$crude_lnrr_var[d2$sterile==0 & d2$pathogen == 'streptococcus_pneumoniae' & d2$refvar=='DATA' & d2$drug=='beta_lactamase_inhibitor']
d2$vi[d2$sterile==1 & d2$pathogen == 'streptococcus_pneumoniae' & d2$refvar=='DATA' & d2$drug=='methicillin'] <- d2$crude_lnrr_var[d2$sterile==1 & d2$pathogen == 'streptococcus_pneumoniae' & d2$refvar=='DATA' & d2$drug=='methicillin']
d2$vi[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='methicillin'] <- d2$crude_lnrr_var[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='methicillin']
d2$vi[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='sulfa'] <- d2$crude_lnrr_var[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='sulfa']
d2$vi[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='anti_pseudomonal_penicillin'] <- d2$crude_lnrr_var[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='anti_pseudomonal_penicillin']
d2$vi[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='fourth_gen_ceph'] <- d2$crude_lnrr_var[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='fourth_gen_ceph']
d2$vi[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='anti_pseudomonal_penicillin'] <- d2$crude_lnrr_var[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='anti_pseudomonal_penicillin']
d2$vi[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='carbapenem'] <- d2$crude_lnrr_var[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='carbapenem']
d2$vi[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='fluoroquinolone'] <- d2$crude_lnrr_var[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='fluoroquinolone']
d2$vi[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='fourth_gen_ceph'] <- d2$crude_lnrr_var[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='fourth_gen_ceph']
d2$vi[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='sulfa'] <- d2$crude_lnrr_var[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='sulfa']
d2$vi[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='third_gen_ceph'] <- d2$crude_lnrr_var[d2$sterile==1 & d2$pathogen == 'klebsiella_pneumoniae' & d2$refvar=='DATA' & d2$drug=='third_gen_ceph']
d2$vi[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='methicillin'] <- d2$crude_lnrr_var[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='methicillin']
d2$vi[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='fluoroquinolone'] <- d2$crude_lnrr_var[d2$sterile==0 & d2$pathogen == 'enterobacter_spp' & d2$refvar=='DATA' & d2$drug=='fluoroquinolone']
d2$vi[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='sulfa'] <- d2$crude_lnrr_var[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='sulfa']
d2$vi[d2$sterile==1 & d2$pathogen == 'coagulase_negative_staph' & d2$refvar=='DATA' & d2$drug=='sulfa'] <- d2$crude_lnrr_var[d2$sterile==1 & d2$pathogen == 'coagulase_negative_staph' & d2$refvar=='DATA' & d2$drug=='sulfa']
d2$vi[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='third_gen_ceph'] <- d2$crude_lnrr_var[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='third_gen_ceph']
d2$vi[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='aminopenicillin'] <- d2$crude_lnrr_var[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='aminopenicillin']
d2$vi[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='aminopenicillin'] <- d2$crude_lnrr_var[d2$sterile==1 & d2$pathogen == 'escherichia_coli' & d2$refvar=='DATA' & d2$drug=='aminopenicillin']
d2$vi[d2$sterile==1 & d2$pathogen == 'staphylococcus_aureus' & d2$refvar=='DATA' & d2$drug=='aminoglycoside'] <- d2$crude_lnrr_var[d2$sterile==1 & d2$pathogen == 'staphylococcus_aureus' & d2$refvar=='DATA' & d2$drug=='aminoglycoside']

# Merge observed to predicted and add metadata
  d2 <- d2 %>% left_join(summary_sterile, by = c('sterile','pathogen','drug' = 'abx_class')) %>% 
    left_join(location_meta,by='location_id') %>% left_join(file_runs, by = c('pathogen','drug' = 'abx_class')) %>%
    mutate(observed = exp(yi), obs_var = ((exp(yi) - exp(yi+(1.96*sei)))/1.96)^2, observed_lb = exp(yi - (1.96 * sei)), observed_ub = exp(yi + (1.96 * sei)))

  # Compute confidence intervals based on both observed and predicted variance
  d2$var_adj_is <- d2$obs_var + d2$pred_var
  d2$pred_lb_adj_is <- d2$predicted - 1.96*sqrt(d2$var_adj_is)
  d2$pred_ub_adj_is <- d2$predicted + 1.96*sqrt(d2$var_adj_is)
  
  # subset to combination and specimen of interest
  d2 <- d2[!is.na(d2$resistance_run_id),]
  d2 <- d2[d2$sterile==1,]
  totalref <- length(d2$refvar)
  
  # empty vector to append performance results 
  table_metrics<-c()

  # A loop to obtain evaluation metrics per antibiotic class 
  for(j in unique(d2$drug)) {
    data <- d2[d2$drug == paste0(j) & d2$sterile == 1,]
    data <- data[,c("location_id", "predicted", "pred_ub_adj_is", "pred_lb_adj_is","observed",'obs_var','super_region','drug','pathogen','combo','pred_ub','pred_lb')]
    data <- na.omit(data)

    # METRICS : RMSE
    rmse_is <- rmse(data$observed, !is.na(data$predicted))
    
    # MAE
    mae_is <- mae(data$observed, data$predicted)
    
    # coefficient of determination
    rss <- sum((data$observed - data$predicted) ^ 2)
    tss <- sum((data$observed - mean(data$observed)) ^ 2)
    rsq_is <- 1 - (rss/tss)
    
    # coverage
    data$covered <- ifelse((data$observed >= data$pred_lb_adj_is) & (data$observed <= data$pred_ub_adj_is),1,0)
    coverage_is <- nrow(data[data$covered == 1,])/length(data$covered)
        
    # all metrics
    is <- c(paste0(j),rsq_is,round(rmse_is,2),round(mae_is,2),round(coverage_is,3))
    table_metrics <- rbind(table_metrics,is)

    # a plot
    upper <- max(data$pred_ub)
    lower <- min(data$pred_lb)
    pal <- c("observed"="black","predicted"="#009999")
    isv_plot <- ggplot(data, aes(ymax = pred_ub, ymin = pred_lb)) + 
      scale_colour_manual(values = pal, limits = names(pal)) +
            geom_point(aes(y = observed, x = combo)) +geom_point(aes(y = predicted, x = combo,color="#009999")) +
      geom_errorbar(aes(x = combo, y = predicted), color = rep("#009999",length(data$pathogen)), width=0) +
      theme_bw() + labs(x = "", y = "Relative Risk") + coord_flip() + ylim(0,1.6) +
     labs(x = "Pathogen-drug combinations", y = "Relative risk and confidence interval", size = "Weight", title = (paste0("In-sample validation MRBRT - ",j)),
   subtitle = paste0("Coverage: ", round(coverage_is,2), ", RMSE: ", round(rmse_is,2), ", MAE: ", round(mae_is,2)))
    
    ggsave(paste0('FILENAME'), plot = isv_plot ,
           height = 18, width = 18, unit = 'cm')
    rm(isv_plot,data,is_rmse,is_r2,rsq,tss,rss,is,rsq_is,rmse_is,mae_is,coverage_is)
  }

table_metrics <- as.data.frame(table_metrics)
colnames(table_metrics) <- c('abx_class','R2_is','RMSE_is','MAE_is',"coverage_is")
table_metrics$sterile <- as.integer(1)
table_metrics <- left_join(table_metrics, abx_meta) 
write.csv(table_metrics[,c('abx_class_name_long','RMSE_is','MAE_is',"coverage_is")], paste0(pub.dir,"FILENAME"), row.names = F)



# Out of sample validation
# 4 random parts to hold out
d2$holdout <- '.'
set.seed(4321)
for(i in unique(d2$drug)) {
  references <- d2[d2$drug == paste0(i), c('pathogen','refvar','drug','location_id','total')]
  references$uid <- paste0(references$pathogen, references$refvar, references$location_id, references$total)
  abx_uid <- sample(references$uid,length(references$uid))
  fold_id <- cut(seq(1,length(abx_uid)),breaks=4,labels=FALSE)
  abx_uid <- as.data.table(cbind(abx_uid, fold_id))
  references <- references %>% left_join(abx_uid, by = c('uid' = 'abx_uid')) %>% dplyr::select(refvar,pathogen,drug,fold_id,location_id,total)
  d2 <- d2 %>% left_join(references, by = c('refvar','pathogen','drug','location_id','total'))
  d2$holdout <- ifelse(!is.na(d2$fold_id),d2$fold_id,d2$holdout)
  d2$fold_id <- NULL
  rm(references, abx_uid, fold_id)
  }

if (length(d2$refvar) == totalref) {print('references matched')} else {print('error in matching references')}

# initial vector of results
predicted_holdout<-c()
for(j in unique(d2$drug)) {
    for(k in 1:4) {
    print(paste0(j," k:",k))
    # initialise priors before holdouts
    overall_mean <- rma(data = d2[d2$drug == paste0(j) & d2$sterile == 1,], yi = yi , vi = vi)$beta
    print(paste0('prior: ',overall_mean))
    eff <- c('prior_beta_gaussian', overall_mean, 0.1)
    
    # subset data to exclude holdout and include sterile - abx_class  
    data <- d2[d2$drug == paste0(j) & d2$sterile == 1 & d2$holdout != paste0(k),]
    orgs<-paste(unique(data$pathogen),sep = ',')
    for (i in orgs) {
      data$new <- ifelse(data$pathogen == i,1,0)
      colnames(data)[colnames(data) == 'new'] <- paste0(i)
    }
  
    #this loop creates the priors matrix needed for the model
    priors <- c()
    for(i in orgs) {
      effs <- c(i,eff)
      priors <- rbind(priors, effs)
    }
    
    # Concatenate the definition of priors to enter the model
    definition_model <- paste0("LinearCovModel(\"",priors[,1], "\",",priors[,2]," = array(c(", as.double(priors[,3]), ",",as.double(priors[,4]),")))")
    model_definition <- paste(definition_model, collapse = ",\n")
    
    # MRBRT data
    dat1 <- MRData()
    dat1$load_df(
      data = data,  col_obs = "yi", col_obs_se = "sei",
      col_covs = as.list(orgs),
      col_study_id = "pathogen" )
    
    # Text to parse into the MRBRT model
    mrbrtschema <- paste0(
      "MRBRT(
  data = dat1,
  cov_models = list(LinearCovModel('intercept', use_re = TRUE, prior_beta_uniform = array(c(0,0))),",
      model_definition,
      "))"
    )
    
    # initialize MRBRT
    mod1 <- eval(parse(text = paste(mrbrtschema)))
    mod1$fit_model(inner_print_level = 5L, inner_max_iter = 1000L)
    samples <- mod1$sample_soln(sample_size = 100L)
    
    prediction <- c()
    for (x in 1:(length(orgs)+1)) {
      beta <- mod1$beta_soln[x] 
      beta_sd <- sd(samples[[1]][,x]) 
      vector <- c(mod1$cov_model_names[x],exp(beta),beta_sd,c(exp(beta-(1.96*beta_sd))),c(exp(beta+(1.96*beta_sd))))
      prediction <- rbind(prediction,vector)
      vector <-c()
    }
    prediction <- as.data.frame(prediction)
    names(prediction) <- c('pathogen','predicted_oos','pred_oos_se','pred_oos_lb','pred_oos_ub')
    prediction$drug <- paste0(j)
    prediction$sterile <- 1
    prediction$holdout <- paste0(k)
    predicted_holdout <-rbind(predicted_holdout,prediction)
  }
}  
predicted_holdout <- as.data.frame(predicted_holdout)
predicted_holdout$predicted_oos <- as.numeric(predicted_holdout$predicted_oos)    
predicted_holdout$pred_oos_lb <- as.numeric(predicted_holdout$pred_oos_lb)    
predicted_holdout$pred_oos_ub <- as.numeric(predicted_holdout$pred_oos_ub)    
predicted_holdout$holdout <- as.integer(predicted_holdout$holdout)
predicted_holdout <- predicted_holdout %>% mutate(pred_var_oos = ((predicted_oos - pred_oos_ub)/1.96)^2)

# merging predicted oos with data  
d2$holdout <- as.integer(d2$holdout)
d2 <- d2 %>% left_join (predicted_holdout, by = c('sterile','holdout', 'pathogen', 'drug')) 
if (length(d2$refvar) == totalref) {print('references matched')} else {print('error in matching references')}

# Compute confidence intervals based on both observed and predicted variance
d2$var_adj_oos <- d2$obs_var + d2$pred_var_oos
d2$pred_lb_adj_oos <- d2$predicted - 1.96*sqrt(d2$var_adj_oos)
d2$pred_ub_adj_oos <- d2$predicted + 1.96*sqrt(d2$var_adj_oos)

# metrics for out-of-sample prediction by antibiotic class
table_oos <- c()
for(j in unique(d2$drug)) {
  # subset data to exclude holdout and inluclude sterile - abx_class  
  data <- d2[d2$drug == paste0(j) & d2$sterile == 1,]
  data <- data[!is.na(data$predicted_oos),]
  
  #RMSE
  rmse_oos <- rmse(data$observed, data$predicted_oos)
  
  # coefficient of determination
  rss <- sum((data$observed - data$predicted_oos)^ 2)
  tss <- sum((data$observed - mean(data$observed) ^ 2))
  rsq <- 1 - (rss/tss)
  
  # MAE
  mae_oos <- mae(data$observed, data$predicted_oos)
  
  # coverage
  data$covered <- ifelse((data$observed >= data$pred_lb_adj_oos) & (data$observed <= data$pred_ub_adj_oos),1,0)
  coverage_oos <- nrow(data[data$covered == 1,])/length(data$covered)

  #table
  table <- c(paste0(j),round(rmse_oos,2),round(mae_oos,2),round(coverage_oos,3))
  table_oos <- rbind(table_oos, table)

  # plot
  upper <- max(data$pred_oos_ub)
  lower <- min(data$pred_oos_lb)
  pal <- c("observed"="black","predicted"="orange")

  oosv_plot <- ggplot(data, aes(ymax = pred_oos_ub, ymin = pred_oos_lb)) + 
    scale_colour_manual(values = pal, limits = names(pal)) +
    geom_point(aes(y = observed, x = combo)) +geom_point(aes(y = predicted_oos, x = combo,color="orange")) +
    geom_errorbar(aes(x = combo, y= predicted_oos), color = rep("orange",length(data$pathogen)), width=0) +
    theme_bw() + labs(x = "", y = "Relative Risk") + coord_flip() + ylim(lower,upper) +
    labs(x = "Pathogen-drug combinations", y = "Relative risk and confidence interval", size = "Weight",
           title = (paste0("Out-of-sample validation MRBRT - ",j)),
          subtitle = paste0("Coverage : ", round(coverage_oos,2), "RMSE : ", round(rmse_oos,2), ", MAE: ", round(mae_oos,2)))
  
  ggsave(paste0('FILENAME'), plot = oosv_plot ,
         height = 18, width = 18, unit = 'cm')
  rm(isv_plot,data,is_rmse,is_r2,rsq,tss,rss,is_rsq,rmse_is)
  rm(table,data,rmse_oos,rss,tss,rsq,oos_rsq,oosv_plot,coverage_oos,mae_oos,upper,lower)
}

table_oos <- as.data.frame(table_oos)
colnames(table_oos) <- c('abx_class','RMSE_oos','MAE_oos','coverage_oos')
table_oos$sterile <- as.integer(1)
validation <- table_metrics %>% left_join(table_oos, by = c('abx_class', 'sterile'))

# Add metadata and export summary table
validation[,c('sterile','abx_class','abx_class_name_short')] <- NULL
validation <- validation[order(validation$coverage_is,validation$abx_class_name_long),]
validation <- validation[order(validation$coverage_is,validation$abx_class_name_long),c('abx_class_name_long',"RMSE_is",'MAE_is','coverage_is',
              colnames(validation[grep("oos",colnames(validation))]))]
write.csv(validation,paste0('FILENAME'),row.names = F, na="")

