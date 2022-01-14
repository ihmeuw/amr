##################################################
## Project: Antimicrobial Resistance - CW Resistance Input Data
## Script purpose: Crosswalk tertiary resistance levels to
##                 non-tertiary prior to running Ensemble-STGPR
## Date: 8/5/2021
##################################################
##### R Initialization and Functions #####
hdrive <- ADDRESS
setwd(hdrive)
rm(list = ls())
library(dplyr)
library(ggplot2)
source(FILEPATH)
source(FILEPATH)
library(crosswalk, lib.loc = FILEPATH)
library(mrbrt001, lib.loc = FILEPATH)

add_loc_info <- function(df){
  locs = get_location_metadata(location_set_id = 35, 
                               gbd_round_id = 6)[, c("location_id", "location_name", "region_name", "super_region_name")]
  df <- merge(df, locs, by = 'location_id', all.x = TRUE)
  return(df)
}

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

set.seed(123)

modelpath <- FILEPATH

##### Read in Data #####
resdat <- fread(FILEPATH, stringsAsFactors = F)

combos <- read.csv(FILEPATH)
combos <- combos[combos$pathogen != 'neisseria_gonorrheae',]
resdat <- merge(resdat, combos[,c('pathogen', 'abx_class')], by = c('pathogen', 'abx_class'), all = FALSE)
SOURCE_1 <- resdat[resdat$source == 'SOURCE_1',]
resdat <- resdat[!resdat$source == 'SOURCE_1',]
combos <- combos[combos$modeler == 'MODELER',]
resdat <- merge(resdat, combos[,c('pathogen', 'abx_class')], by = c('pathogen', 'abx_class'), all = FALSE)

resdat2 <- read.csv(FILEPATH)
resdat2 <- resdat2[!is.na(resdat2$val),]
resdat <- bind_rows(resdat, resdat2)

tertdat <- resdat[resdat$source != 'ecdc_ears' | resdat$location_id != 95,]
tertdat <- bind_rows(tertdat, SOURCE_1)

tertdat <- add_loc_info(tertdat)
tertdat$tert_bin <- ifelse(tertdat$hospital_type %in% 'tertiary', 'tertiary', 'nontert-mixed')

htcounts <- tertdat[!is.na(tertdat$super_region_name),] %>% dplyr::group_by(super_region_name) %>%
  dplyr::summarize(
    non_tertiary = sum(cases[hospital_type == 'non-tertiary'])/sum(cases),
    tertiary = sum(cases[hospital_type == 'tertiary'])/sum(cases),
    mixed_unknown = sum(cases[hospital_type == 'm/u'])/sum(cases)
  )

tertdat <- tertdat %>% dplyr::group_by(pathogen, abx_class, source, tert_bin, super_region_name, year_id) %>%
  dplyr::summarize(
    resistant = sum(resistant),
    susceptible = sum(susceptible),
    cases = sum(cases)
  )

tertdat <- tertdat[tertdat$cases >= 5,]

tertdat$prev <- tertdat$resistant/tertdat$cases
tertdat$prev[tertdat$prev == 0] <- 0.001
tertdat$prev[tertdat$prev == 1] <- 0.999
tertdat$se <- sqrt((tertdat$prev*(1-tertdat$prev))/tertdat$cases)

tertdat$index <- seq.int(nrow(tertdat))

keepcols <- c('location_id', 'region_name', 'super_region_name', 'pathogen', 'abx_class',
              'year_id', 'source', 'hospital_type', 'prev', 'se', 'index', 'cases')

df_ref <- tertdat[tertdat$tert_bin == 'nontert-mixed', colnames(tertdat) %in% keepcols]
df_alt <- tertdat[tertdat$tert_bin == 'tertiary', colnames(tertdat) %in% keepcols]

df_matched <- merge(df_ref, df_alt, by = c('super_region_name', 'pathogen', 'abx_class'), suffixes = c('_ref', '_alt'), allow.cartesian=TRUE)

rm(df_ref, df_alt, SOURCE_1, resdat2, tertdat)

dat_diff <- as.data.frame(cbind(
  delta_transform(
    mean = df_matched$prev_alt, 
    sd = df_matched$se_alt,
    transformation = "linear_to_logit" ),
  delta_transform(
    mean = df_matched$prev_ref, 
    sd = df_matched$se_ref,
    transformation = "linear_to_logit")
))
names(dat_diff) <- c("mean_alt", "mean_se_alt", "mean_ref", "mean_se_ref")

df_matched[, c("logit_diff", "logit_diff_se")] <- calculate_diff(
  df = dat_diff, 
  alt_mean = "mean_alt", alt_sd = "mean_se_alt",
  ref_mean = "mean_ref", ref_sd = "mean_se_ref" )

df_matched$altvar <- 'tertiary'
df_matched$refvar <- 'nontert-mixed'

match5yr <- df_matched[abs(df_matched$year_id_alt - df_matched$year_id_ref) <= 5,]

match5yr$super_region_name <- case_when(match5yr$super_region_name == "Central Europe, Eastern Europe, and Central Asia" ~ "CEurCAsia",
                                        match5yr$super_region_name == "High-income" ~ "HI",
                                        match5yr$super_region_name == "Latin America and Caribbean" ~ "LatAm",
                                        match5yr$super_region_name == "North Africa and Middle East" ~ "NAfrME",
                                        match5yr$super_region_name == "South Asia" ~ "SAsia",
                                        match5yr$super_region_name == "Southeast Asia, East Asia, and Oceania" ~ "SEAsiaOce",
                                        match5yr$super_region_name == "Sub-Saharan Africa" ~ "SSAfr")


match5yr$bug_group <- case_when(match5yr$pathogen %in% c('enterococcus_faecalis', 'enterococcus_faecium', 
                                                         'enterococcus_spp', 'staphylococcus_aureus',
                                                         'group_b_strep', 'group_a_strep', 'streptococcus_pneumoniae') ~ 'strep_group',
                                match5yr$pathogen %in% c('citrobacter_spp', 'enterobacter_spp', 'escherichia_coli',
                                                         'haemophilus_influenzae', 'klebsiella_pneumoniae', 'proteus_spp', 'serratia_spp') ~ 'ecoli_group',
                                TRUE ~ 'acinetobacter_group')


match5yr$drug_group <- case_when(match5yr$abx_class %in% c('anti_pseudomonal_penicillin', 'beta_lactamase_inhibitor',
                                                           'carbapenem', 'fourth_gen_ceph', 'third_gen_ceph', 'aminopenicillin',
                                                           'methicillin', 'penicillin') ~ 'beta_lactam',
                                 TRUE ~ match5yr$abx_class)

match5yr$combo <- paste0(match5yr$pathogen, ' - ', match5yr$abx_class)

##############
# RUN MODELS
##############
glob_mod_combos <- match5yr %>% group_by(bug_group, drug_group) %>%
  summarize(
    n = length(logit_diff),
    nspr = length(unique(super_region_name))
  )

glob_mod_combos$mean_eff <- NA

for (i in 1:nrow(glob_mod_combos)){
  bug = glob_mod_combos$bug_group[i]
  drug = glob_mod_combos$drug_group[i]
  data <- match5yr[match5yr$bug_group == bug & match5yr$drug_group == drug,]
  
  mrdat <- MRData()
  mrdat$load_df(
    data = data,  col_obs = "logit_diff", col_obs_se = "logit_diff_se",
    col_study_id = "combo")
  
  combomod <- MRBRT(
    data = mrdat,
    cov_models = list(
      LinearCovModel("intercept", use_re = TRUE, prior_beta_uniform = array(c(0,Inf)))))
  
  combomod$fit_model(inner_print_level = 5L, inner_max_iter = 10000L)
  
  glob_mod_combos$mean_eff[i] <- combomod$summary()[[1]][[1]]
  
  py_save_object(object = combomod, filename = paste0(modelpath, "global/", bug, ":", drug, ".pkl"), pickle = "dill")
}

spr_mod_combos <- match5yr %>% group_by(bug_group, drug_group, super_region_name) %>%
  summarize(
    n = length(logit_diff)
  )

spr_mod_combos$mean_eff <- NA

for (i in 1:nrow(spr_mod_combos)){
  bug = spr_mod_combos$bug_group[i]
  drug = spr_mod_combos$drug_group[i]
  spr = spr_mod_combos$super_region_name[i]
  data <- match5yr[match5yr$bug_group == bug & match5yr$drug_group == drug & match5yr$super_region_name == spr,]
  
  mrdat <- MRData()
  mrdat$load_df(
    data = data,  col_obs = "logit_diff", col_obs_se = "logit_diff_se",
    col_study_id = "combo")
  
  combomod <- MRBRT(
    data = mrdat,
    cov_models = list(
      LinearCovModel("intercept", use_re = TRUE, prior_beta_uniform = array(c(0,Inf)))))
  
  combomod$fit_model(inner_print_level = 5L, inner_max_iter = 10000L)
  
  spr_mod_combos$mean_eff[i] <- combomod$summary()[[1]][[1]]
  
  if(!dir.exists(paste0(modelpath, spr, '/'))){
    dir.create(paste0(modelpath, spr, '/'), recursive = TRUE)
  }
  py_save_object(object = combomod, filename = paste0(modelpath, spr, "/", bug, ":", drug, ".pkl"), pickle = "dill")
}

##############
# GENERATE PREDICTIONS
##############
resdat <- add_loc_info(resdat)
resdat$super_region_name[resdat$country == 'Kosovo' & is.na(resdat$location_id)] <- 'Central Europe, Eastern Europe, and Central Asia'
resdat$super_region_name[resdat$country == 'Rwanda,Uganda' & is.na(resdat$location_id)] <- 'Sub-Saharan Africa'
resdat$super_region_name[resdat$country == 'Benin,Ghana' & is.na(resdat$location_id)] <- 'Sub-Saharan Africa'

tertprop <- resdat %>% group_by(super_region_name) %>%
  summarize(
    prop_n_tertiary = sum(cases[hospital_type == 'tertiary'])/sum(cases),
    prop_dp_tertiary = length(cases[hospital_type == 'tertiary' & cases >= 5])/length(cases[cases >= 5])
  )
tertprop <- tertprop[!is.na(tertprop$super_region_name),]
write.csv(tertprop, FILEPATH, row.names = F)

resdat$super_region_name <- case_when(resdat$super_region_name == "Central Europe, Eastern Europe, and Central Asia" ~ "CEurCAsia",
                                      resdat$super_region_name == "High-income" ~ "HI",
                                      resdat$super_region_name == "Latin America and Caribbean" ~ "LatAm",
                                      resdat$super_region_name == "North Africa and Middle East" ~ "NAfrME",
                                      resdat$super_region_name == "South Asia" ~ "SAsia",
                                      resdat$super_region_name == "Southeast Asia, East Asia, and Oceania" ~ "SEAsiaOce",
                                      resdat$super_region_name == "Sub-Saharan Africa" ~ "SSAfr")

resdat$prev <- resdat$resistant/resdat$cases
resdat$COprev <- resdat$prev
resdat$COprev[resdat$prev == 0] <- 0.001
resdat$COprev[resdat$prev == 1] <- 0.999
resdat$se <- sqrt((resdat$COprev*(1-resdat$COprev))/resdat$cases)

resdat <- resdat[!is.na(resdat$COprev),]

resdat <- as.data.table(resdat) %>% 
  cbind(delta_transform(resdat$COprev, resdat$se, "linear_to_logit"))

resdat$hospital_type <- ifelse(resdat$hospital_type == 'tertiary', 'tertiary', 'nontert-mixed')

resdat$bug_group <- case_when(resdat$pathogen %in% c('enterococcus_faecalis', 'enterococcus_faecium', 
                                                     'enterococcus_spp', 'staphylococcus_aureus',
                                                     'group_b_strep', 'group_a_strep', 'streptococcus_pneumoniae') ~ 'strep_group',
                              resdat$pathogen %in% c('citrobacter_spp', 'enterobacter_spp', 'escherichia_coli',
                                                     'haemophilus_influenzae', 'klebsiella_pneumoniae', 'proteus_spp', 'serratia_spp') ~ 'ecoli_group',
                              TRUE ~ 'acinetobacter_group')


resdat$drug_group <- case_when(resdat$abx_class %in% c('anti_pseudomonal_penicillin', 'beta_lactamase_inhibitor',
                                                       'carbapenem', 'fourth_gen_ceph', 'third_gen_ceph', 'aminopenicillin',
                                                       'methicillin', 'penicillin') ~ 'beta_lactam',
                               TRUE ~ resdat$abx_class)

resdat$index <- c(1:nrow(resdat))

#### Global Adjustment ####
globmods <- gsub('.pkl', '', list.files(paste0(modelpath, 'global/'), pattern = ':'))
tertadjglob <- as.data.frame(c())

for (mod in globmods){
  bug = strsplit(mod, ':')[[1]][1]
  drug = strsplit(mod, ':')[[1]][2]
  
  result = py_load_object(filename = paste0(modelpath, 'global/', mod, ".pkl"), pickle = "dill")
  
  tertdata = resdat[resdat$hospital_type == 'tertiary' & resdat$bug_group == bug & resdat$drug_group == drug,]
  tertpred <- MRData()
  
  tertpred$load_df(
    data = tertdata,
    col_covs=list('prev')
  )
  
  tertdata$adj <- result$predict(data = tertpred)
  tertdata$adj_logit <- tertdata$mean_logit - tertdata$adj
  tertdata$pred_prev <- logit2prob(tertdata$adj_logit)
  
  tertadjglob <- bind_rows(tertadjglob, tertdata)
}

adjdat <- bind_rows(resdat[resdat$hospital_type != 'tertiary',], tertadjglob)

adjdat$pred_prev[adjdat$hospital_type != 'tertiary'] <- adjdat$COprev[adjdat$hospital_type != 'tertiary']

#### Super Region Adjustment ####
spr_mod_combos_datarich <- spr_mod_combos[spr_mod_combos$n >= 250,]
tertadjspr <- as.data.frame(c())

for (i in 1:nrow(spr_mod_combos_datarich)){
  bug = spr_mod_combos_datarich$bug_group[i]
  drug = spr_mod_combos_datarich$drug_group[i]
  spr = spr_mod_combos_datarich$super_region_name[i]
  
  result = py_load_object(filename = paste0(modelpath, spr, '/', bug, ':', drug, ".pkl"), pickle = "dill")
  
  tertdata = resdat[resdat$hospital_type == 'tertiary' & resdat$bug_group == bug & resdat$drug_group == drug & resdat$super_region_name == spr,]
  
  tertpred <- MRData()
  tertpred$load_df(
    data = tertdata,
    col_covs=list('prev')
  )
  
  tertdata$adj <- result$predict(data = tertpred)
  tertdata$adj_logit <- tertdata$mean_logit - tertdata$adj
  tertdata$pred_prev <- logit2prob(tertdata$adj_logit)
  
  tertadjspr <- bind_rows(tertadjspr, tertdata)
}

adjdat <- merge(adjdat, tertadjspr[,c('pathogen', 'abx_class', 'year_id', 'source', 'hospital_type', 'super_region_name', 'cases', 'resistant', 'susceptible', 'index', 'pred_prev'),],
                  by = c('pathogen', 'abx_class', 'year_id', 'source', 'hospital_type', 'super_region_name', 'cases', 'resistant', 'susceptible', 'index'),
                  suffixes = c('_glob', '_spr'), all.x = TRUE)

adjdat$pred_prev_spr[is.na(adjdat$pred_prev_spr)] <- adjdat$pred_prev_glob[is.na(adjdat$pred_prev_spr)]


#### Adjusted Set 1 ####
adj1cols <- c(colnames(resdat2)[!colnames(resdat2) %in% 'n_resistant'], 'pred_prev_spr')
adjdat1 <- adjdat[!is.na(adjdat$country), ..adj1cols]

adjdat1$adjusted_resistant[adjdat1$resistant != 0 & adjdat1$hospital_type == 'tertiary'] <-
  adjdat1$pred_prev_spr[adjdat1$resistant != 0 & adjdat1$hospital_type == 'tertiary'] *
  adjdat1$cases[adjdat1$resistant != 0 & adjdat1$hospital_type == 'tertiary']

adjdat1$adjusted_resistant[adjdat1$resistant == 0 | adjdat1$hospital_type != 'tertiary'] <-
  adjdat1$resistant[adjdat1$resistant == 0 | adjdat1$hospital_type != 'tertiary']

adjdat1 <- adjdat1[, pred_prev_spr:= NULL]

write.csv(adjdat1, FILEPATH, row.names = F)

table(adjdat1$resistant != adjdat1$adjusted_resistant)

#### Adjusted Set 2 ####
adjdat2 <- adjdat[is.na(adjdat$country),]

adjdat2$resistant[adjdat2$resistant != 0 & adjdat2$hospital_type == 'tertiary'] <-
  adjdat2$pred_prev_spr[adjdat2$resistant != 0 & adjdat2$hospital_type == 'tertiary'] *
  adjdat2$cases[adjdat2$resistant != 0 & adjdat2$hospital_type == 'tertiary']

adjdat2$susceptible <- adjdat2$cases - adjdat2$resistant

keepcols2 <- c('source', 'location_id', 'year_id', 'pathogen', 'abx_class', 'hospital_type', 'cases',
                  'resistant', 'susceptible')

adjdat2 <- adjdat2[, ..keepcols2]

aggadj2 <- adjdat2 %>% dplyr::group_by(source, location_id, year_id, pathogen, abx_class) %>%
  dplyr::summarize(
    cases = sum(cases),
    resistant = sum(resistant),
    susceptible = sum(susceptible),
    nsub_dp = length(hospital_type)
  )

aggadj2 <- aggadj2[, !colnames(aggadj2) %in% 'nsub_dp']

write.csv(aggadj2, FILEPATH, row.names = F)
