rm(list = ls())
hdrive <- ADDRESS
setwd(hdrive)
library(ggplot2)
library(data.table)
library(gbm)
library(xgboost)
library(mgcv)
library(tidyverse)
library(glmnet)
library(matrixStats)
library(quadprog)
library(gtools)
library(nnet)
library(randomForest)
library(stats)
library(lme4)
library(Cubist, lib.loc = "R_packages")
library(caret, lib.loc = "R_packages")
source(FILEPATH)
source(FILEPATH)

#~~~~~~~~~~~~~#
# i. Setup ####
#~~~~~~~~~~~~~#

oos <- F
n_holdouts <- 10
oos_method <- 'country'

data_path <- FILEPATH
abx <- read.csv(FILEPATH) %>% filter(year>1999) %>%
  dplyr::mutate(as.vector(scale(ddd_per_1000, center = TRUE, scale = TRUE))) %>%
  dplyr::mutate(as.vector(scale(J01A, center = TRUE, scale = TRUE))) %>%
  dplyr::mutate(as.vector(scale(J01B, center = TRUE, scale = TRUE))) %>%
  dplyr::mutate(as.vector(scale(J01C, center = TRUE, scale = TRUE))) %>%
  dplyr::mutate(as.vector(scale(J01D, center = TRUE, scale = TRUE))) %>%
  dplyr::mutate(as.vector(scale(J01E, center = TRUE, scale = TRUE))) %>%
  dplyr::mutate(as.vector(scale(J01F, center = TRUE, scale = TRUE))) %>%
  dplyr::mutate(as.vector(scale(J01G, center = TRUE, scale = TRUE))) %>%
  dplyr::mutate(as.vector(scale(J01M, center = TRUE, scale = TRUE)))
abx <- cbind(abx$loc_id,abx$year,abx[,12:20])
colnames(abx) <-c('location_id','year_id','ddd_per_1000','j01a','j01b','j01c','j01d','j01e','j01f','j01g','j01m')

temperature <- get_covariate_estimates(covariate_id = 71, gbd_round_id = 6 , year_id = 2000:2018, decomp_step="iterative") %>%
  dplyr::select(location_id, year_id, mean_value) %>%
  dplyr::mutate(as.vector(scale(mean_value, center = TRUE, scale = TRUE))) 
colnames(temperature) <-c('location_id','year_id','raw_mean_temperature','mean_temperature')
temperature <- temperature %>% dplyr::select(location_id,year_id,mean_temperature)

pigs <- get_covariate_estimates(covariate_id = 100, gbd_round_id = 6 , year_id = 2000:2018, decomp_step="iterative") %>%
  dplyr::select(location_id, year_id, mean_value) %>%
  dplyr::mutate(as.vector(scale(mean_value, center = TRUE, scale = TRUE))) 
colnames(pigs) <-c('location_id','year_id','raw_pigs_pc','pigs_pc')
pigs <- pigs %>% dplyr::select(location_id,year_id,pigs_pc)

haqi <- get_covariate_estimates(covariate_id = 1099, gbd_round_id = 6 , year_id = 2000:2018, decomp_step="iterative") %>%
  dplyr::select(location_id, year_id, mean_value)%>%
  dplyr::mutate(as.vector(scale(mean_value, center = TRUE, scale = TRUE))) 
colnames(haqi) <-c('location_id','year_id','haqi_u','haqi')
haqi <- haqi %>% dplyr::select(location_id,year_id,haqi)

location_md <- subset(get_location_metadata(location_set_id=35, gbd_round_id=6),select = c(location_id, location_name,super_region_name,super_region_id,region_id,region_name,level,ihme_loc_id,parent_id), droplevels = TRUE)

initial_data <- fread(FILEPATH, stringsAsFactors = F)
initial_data <- initial_data[!initial_data$source %in% c('SOURCE_1', 'SOURCE_2')]
initial_data <- initial_data[!(initial_data$source == 'SOURCE_3' & initial_data$year_id >=2016)]

fqgono <- read.csv(FILEPATH)
initial_data <- bind_rows(initial_data, fqgono)
rm(fqgono)

SOURCE_4 <- read.csv(FILEPATH)
initial_data <- bind_rows(initial_data, SOURCE_4)
rm(SOURCE_4)

ir <- read.csv(FILEPATH, stringsAsFactors = F) %>% filter(include == 1) %>% dplyr::select(abx_class, pathogen, atc)
ir$atc[ir$abx_class == 'vancomycin'] <- 'j01d'
ir$combo <- paste0(ir$pathogen,"-",ir$abx_class)
ir <- ir %>% add_row(abx_class = "fluoroquinolone", pathogen = "neisseria_gonorrheae",
                     atc = 'J01M', combo = 'neisseria_gonorrheae-fluoroquinolone')

bugdrugs <- read.csv(FILEPATH)

ir <- merge(ir, bugdrugs[bugdrugs$modeler == 'MODELER', c('pathogen', 'abx_class')])

combos <- unique(ir$combo)
primary_combos <- c('escherichia_coli-fluoroquinolone','escherichia_coli-third_gen_ceph','klebsiella_pneumoniae-carbapenem','klebsiella_pneumoniae-third_gen_ceph','staphylococcus_aureus-methicillin','streptococcus_pneumoniae-penicillin','non_typhoidal_salmonellae-fluoroquinolone','shigella_spp-fluoroquinolone')

small_combos <- c('acinetobacter_baumanii-methicillin','acinetobacter_baumanii-macrolide','acinetobacter_baumanii-penicillin','acinetobacter_baumanii-vancomycin','klebsiella_pneumoniae-vancomycin','streptococcus_pneumoniae-vancomycin','enterobacter_spp-vancomycin','haemophilus_influenzae-beta_lactamase_inhibitor','enterococcus_faecalis-fourth_gen_ceph','group_a_strep-sulfa','group_b_strep-sulfa','group_a_strep-fluoroquinolone','haemophilus_influenzae-sulfa','morganella_spp-anti_pseudomonal_penicillin','morganella_spp-sulfa','morganella_spp-aminoglycoside','morganella_spp-macrolide','enterococcus_faecium-sulfa','providencia_spp-sulfa')
glm_stacker <- c() 

combos <- combos[!combos %in% c(primary_combos)]
combos <- combos[!combos %in% c(small_combos)]

organism <- c('acinetobacter_baumanii','pseudomonas_aeruginosa','enterobacter_spp','enterococcus_faecalis','enterococcus_faecium','enterococcus_spp','group_a_strep','group_b_strep','haemophilus_influenzae','escherichia_coli','klebsiella_pneumoniae','staphylococcus_aureus','streptococcus_pneumoniae','serratia_spp','neisseria_meningitidis','citrobacter_spp','proteus_spp','providencia_spp','listeria','morganella_spp','legionella_spp') #'moraxella_spp','non_typhoidal_salmonellae','salmonella_paratyphi','salmonella_typhi','salmonella_typhi_paratyphi') #'shigella_spp','mycoplasma')
antibiotic_group<-c('third_gen_ceph','carbapenem','fluoroquinolone','penicillin','aminopenicillin','beta_lactamase_inhibitor','anti_pseudomonal_penicillin','methicillin','vancomycin','fourth_gen_ceph','sulfa','aminoglycoside','macrolide')

list <- c()
for(i in organism) {
  for(j in antibiotic_group) {
    potential <- paste0(i,"-",j)
    list <- c(list,potential)
    } 
  }
list <-list[list %in% c(combos)]

gonoFQoneoff <- FALSE
if (gonoFQoneoff){
  initial_data <- read.csv(FILEPATH)
  list <- "neisseria_gonorrheae-fluoroquinolone"
  ir <- ir %>% add_row(abx_class = "fluoroquinolone", pathogen = "neisseria_gonorrheae",
                 atc = 'J01M', combo = 'neisseria_gonorrheae-fluoroquinolone')
}

vancorerun <- FALSE
if (vancorerun){
  list <- "staphylococcus_aureus-vancomycin"
}

strep3GC <- FALSE
if (strep3GC){
  list <- "streptococcus_pneumoniae-third_gen_ceph"
}

acinetobacter_fix <- FALSE
if (acinetobacter_fix){
  list <- c("acinetobacter_baumanii-aminoglycoside", "acinetobacter_baumanii-fluoroquinolone")
}

if (oos){
  args <- commandArgs(trailingOnly = T)
  list <- as.character(args[1])
}

for(x in list) {
  y <- substring(x,1,stringr::str_locate(x,"-")-1)[1]
  z <- substring(x,stringr::str_locate(x,"-")+1,)[1]
  
outputdir <-  paste0(FILEPATH,y,'/', z,'/') 
dir.create(outputdir, showWarnings = F, recursive = T)
setwd(outputdir)
mydata <- initial_data

mydata <- dplyr::rename(mydata, sample_size = cases)
mydata <- mydata[mydata$sample_size>=5 & mydata$pathogen == paste0(y) & mydata$abx_class == paste0(z),] 

mydata <- mydata %>% dplyr::rename(nid = source) %>% mutate(val = resistant/sample_size)

mydata <- mydata %>% left_join(location_md, by = 'location_id')
mydata$location_id[mydata$level == 4] <- mydata$parent_id[mydata$level == 4]
mydata$location_id[mydata$location_id == 44767] <- 95
mydata[,c('parent_id','level','super_region_name','region_name','location_name','super_region_id','region_id','ihme_loc_id')] <- NULL  

mydata <- group_by(mydata, location_id, year_id, nid, pathogen, abx_class) %>%
  summarize(
    resistant = sum(resistant),
    susceptible = sum(susceptible),
    sample_size = sum(sample_size)
  )
mydata$val <- mydata$resistant/mydata$sample_size

covs <- left_join(haqi,temperature, by = c('location_id','year_id')) %>%
  left_join(pigs, by = c('location_id','year_id')) %>%
  left_join(abx, by=c('location_id','year_id')) %>% left_join(location_md, by = 'location_id')

abx_for_combo <- tolower(unique(ir$atc[ir$abx_class==paste0(z)]))
covs_to_include <-  c("haqi",
                      "mean_temperature",
                      "pigs_pc","ddd_per_1000", abx_for_combo) 

covs <- covs[colnames(covs) %in% covs_to_include | colnames(covs) %in% c('location_id','year_id','super_region_id','region_id','ihme_loc_id')]
covs <- data.table(covs)
covs <- na.omit(covs)
mydata <- merge(mydata, covs, by = c('location_id', 'year_id'))
mydata <- data.table(mydata)
response = cbind(successes = mydata$resistant,
                 failures = mydata$susceptible)

text <- paste0("haqi + mean_temperature + pigs_pc + ",abx_for_combo)
textmodel <-paste0("glmer(response ~ 1 +",text,"+ (1 |super_region_id / region_id / ihme_loc_id ) , data = mydata, family = 'binomial')")
model1 <- eval(parse(text = (textmodel)))

covs$ihme_loc_id <- as.character(covs$ihme_loc_id) 
covs$pred <- predict(model1, newdata = covs, type = 'response', allow.new.levels = TRUE)
summary(covs$pred)
other<- mydata[,c('ihme_loc_id', 'year_id', 'val', 'sample_size')]
covs <- merge(covs, other, by = c('ihme_loc_id', 'year_id'), all.x = T, all.y = T)
covs <- data.table(covs)
covs$upper_bound <-  NULL
covs$lower_bound <-  NULL
mydata$upper_bound <-  NULL
mydata$lower_bound <-  NULL
mydata$is_outlier <- 0

value_deviation <-2
MADs <-  covs[,.(upper_bound = pred + value_deviation*mad(pred[!is.na(val)], val[!is.na(val)]),
                 lower_bound = pred - value_deviation*mad(pred[!is.na(val)],val[!is.na(val)])),
              by = c('ihme_loc_id','location_id')]

MADs <- MADs[,c('lower_bound','upper_bound')]
covs <- cbind(covs, MADs)
covs$upper_bound[covs$upper_bound>1] <- 1
covs$lower_bound[covs$lower_bound<0] <- 0
covs <- covs[!is.na(covs$super_region_id),]

MADs <- covs[,.(ihme_loc_id, year_id, lower_bound, upper_bound)]
MADs <-  unique(MADs)
mydata <- merge(mydata, MADs, by = c('ihme_loc_id', 'year_id'))
mydata$is_outlier[mydata$val<mydata$lower_bound |mydata$val>mydata$upper_bound] <- 1
mydata$is_outlier[mydata$val>mydata$lower_bound & mydata$val<mydata$upper_bound] <- 0
mydata$age_group_id <- 22
mydata$sex_id<- 3
mydata$measure_id<-18

if (x == 'streptococcus_pneumoniae-third_gen_ceph'){
  mydata$is_outlier[mydata$ihme_loc_id == 'NGA' & mydata$nid == 'SOURCE_5' & mydata$year_id == 2017] <- 1
}

if (x == 'enterococcus_spp-vancomycin'){
  mydata$is_outlier[mydata$ihme_loc_id == 'IRN' & mydata$nid == 'SOURCE_6' & mydata$val > 0.9] <- 1
  mydata$is_outlier[mydata$ihme_loc_id == 'IRN' & mydata$nid == 'SOURCE_6' & mydata$val <= 0.9] <- 0
}
if (x == 'group_b_strep-fluoroquinolone'){
  mydata$is_outlier[mydata$ihme_loc_id == 'IND' & mydata$nid == 'SOURCE_6' & mydata$val > 0.9] <- 1
}
if (x == 'group_b_strep-macrolide'){
  mydata$is_outlier[mydata$ihme_loc_id == 'SAU' & mydata$nid == 'SOURCE_6' & mydata$val > 0.9] <- 1
}

if(x %in% c('staphylococcus_aureus-vancomycin', 'group_b_strep-penicillin')){
  mydata$val[mydata$val == 0] <- 0.005
} else {
  mydata$val[mydata$val==0]<-0.02
}
mydata$val[mydata$val==1]<-0.98

mydata$val[mydata$val > 0.98 & mydata$resistant%%1!=0] <- 0.98

mydata$variance <- ((mydata$val)*(1-mydata$val))/mydata$sample_size
mydata$variance[mydata$variance == 0] <- 0.0001
mydata <- mydata[,c("location_id","year_id","nid","super_region_id","region_id","ihme_loc_id","age_group_id","sex_id","measure_id","is_outlier",'resistant',"sample_size","val","variance")]
write.csv(mydata,paste0(outputdir,'input.csv'), row.names = F)
mydata <- mydata[mydata$is_outlier == 0,]

if (oos) {
  if(oos_method == 'random'){
    mydata <- mydata[sample(nrow(mydata)),]
    mydata[,master_fold_id := cut(seq(1,nrow(mydata)),breaks=n_holdouts,labels=FALSE)]
  }
  
  if(oos_method == 'country'){
    country <- unique(mydata[, location_id])
    country <- country[sample(length(country))]
    master_fold_id <- cut(seq(1,length(country)),breaks=n_holdouts,labels=FALSE)
    folds <- as.data.table(cbind(country, master_fold_id))
    folds <- rename(folds, location_id = country)
    
    mydata <- merge(mydata, folds, by = c('location_id'))
    mydata$master_fold_id <- as.numeric(mydata$master_fold_id)
    rm(country, master_fold_id)
  }
}

write.csv(mydata, paste0(outputdir,'input_w_ho.csv'), row.names = F)

if (oos){
  mydata <- mydata[,c("location_id","year_id","nid",'resistant',"sample_size","val","master_fold_id")]
} else {
  mydata <- mydata[,c("location_id","year_id","nid",'resistant',"sample_size","val")]
}

child_models <- c('gam', 'ridge', 'lasso', 'enet', 'nnet', 'rf', 'cubist')  

if(x %in% c(glm_stacker) ){
  stacker <- 'GLM'
} else {
  stacker <- 'RWM'
}

family <- 'binomial'

transformation <- NULL
centre_scale <- FALSE
include_year <-  TRUE

covs <- left_join(haqi,temperature, by = c('location_id','year_id')) %>%
  left_join(pigs, by = c('location_id','year_id')) %>%
  left_join(abx, by=c('location_id','year_id'))

abx_for_combo <- tolower(unique(ir$atc[ir$abx_class==paste0(z)]))
covs_to_include <-  c("haqi",
                      "mean_temperature",
                      "pigs_pc",abx_for_combo) 

p <- 'val'       
n <- 'resistant'         
d <- 'sample_size'  
w <- NULL            

min_year <- min(mydata$year_id)
max_year <- 2018

colnames(mydata)[colnames(mydata)==d] <- 'd' 
if(!is.null(p)) {colnames(mydata)[colnames(mydata)==p] <- 'p'} 
if(!is.null(n)) {colnames(mydata)[colnames(mydata)==n] <- 'n'} 

if(is.null(n) &!is.null(p)&!is.null(d)){mydata$n <- mydata$p*mydata$d}

if(is.null(p) &!is.null(n)&!is.null(d)){mydata$p <- mydata$n/mydata$d}
mydata$p[mydata$p==0]<-0.02
mydata$p[mydata$p==1]<-0.98

if(is.null(w)){
  mydata$w <- 1
} else {
  colnames(mydata)[colnames(mydata)==w] <- 'w' 
}

if(is.null(transformation)){
} else if(transformation == 'log'){
  if(family == 'binomial'){
    mydata$n <- log(mydata$n)
    mydata$d <- log(mydata$d)
    mydata$p <- log(mydata$p)
  } else if(family == 'gaussian')
    mydata$n <- log(mydata$n)
} else if(transformation == 'logit'){  #ln(p/1-p)
  if(family == 'binomial'){
    mydata$p <- log(mydata$p/(1-mydata$p))
  } else if(family == 'gaussian'){
    message('should not be using logit transformation with gaussian data')
  }
}

covs <- covs[colnames(covs) %in% covs_to_include | colnames(covs)=='location_id' | colnames(covs) =='year_id']
covs <- data.table(covs)
covs <- na.omit(covs)

if(centre_scale == TRUE){
  covs <- data.frame(covs)
  covs[colnames(covs) %in% covs_to_include] <- data.frame(scale(covs[colnames(covs) %in% covs_to_include]))
  covs$year <- scale(covs$year_id)
  covs <-  data.table(covs)
}

mydata <- merge(mydata, covs, by = c('location_id', 'year_id'))
mydata <- data.table(mydata)

if(family == 'binomial'){
  mydata    <- na.omit(mydata, c('n', 'd', 'p', names(covs)))
}

if(family == 'gaussian'){
  mydata    <- na.omit(mydata, c('n', names(covs)))
}

if (oos){
  loopval <- n_holdouts
  master_data <- mydata
  staticoutput <- paste0(FILEPATH,y,'/', z,'/')
  staticcovs <- covs
  staticchildren <- child_models
} else {
  loopval <- 1
}

holdout_method <- 'random'

for(h in 1:loopval){
  if (oos){
    message(paste0('Running holdout ', h))
    mydata <- master_data[master_data$master_fold_id !=h,]
    mydata$master_fold_id <- NULL
    covs <- staticcovs
    child_models <- staticchildren
  }
 
  if(holdout_method == 'random'){
    mydata <- mydata[sample(nrow(mydata)),]
    mydata[,fold_id := cut(seq(1,nrow(mydata)),breaks=5,labels=FALSE)]
  }
  
  if(holdout_method == 'country'){
    country <- unique(mydata[, country])
    country <- country[sample(length(country))]
    fold_id <- cut(seq(1,length(country)),breaks=5,labels=FALSE)
    folds <- as.data.table(cbind(country, fold_id))
    
    mydata <- merge(mydata, folds, by = c('country'))
    mydata$fold_id <- as.numeric(mydata$fold_id)
    rm(country, fold_id)
  }

  mydata <- mydata[mydata$year_id >= min_year & mydata$year_id <= max_year,]
  covs <- covs[covs$year_id >= min_year & covs$year_id <= max_year,]
  
  mydata[, a_rowid := seq(1:nrow(mydata))]
  
  if(include_year == TRUE){
    covs_to_include <- c('year_id', covs_to_include)
  }
  
  if (oos){
    outputdir <- paste0(staticoutput, '/holdout_', h, '/')
    dir.create(outputdir, showWarnings = F)
  }
  
#~~~~~~~~~~~~~~~~~~~~~#
# Fit child models ####
#~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~#
# 1. XGBoost ####
#~~~~~~~~~~~~~~~#

if('xgboost' %in% child_models){
  dir.create(paste0(outputdir, '/xgboost'), showWarnings = F)
  
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  xg_grid <- expand.grid(nrounds = c(50, 100, 200),
                         max_depth = c(4, 6, 8),
                         eta = (3:6) / 100,
                         colsample_bytree = .5,
                         min_child_weight = 1,
                         subsample = 1,
                         gamma = 0)
  
  train_control <- trainControl(selectionFunction = "oneSE",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]))
  
  xg_fit <- caret::train(form,
                  data = mydata,
                  trControl = train_control,
                  verbose = T,
                  tuneGrid = xg_grid,
                  metric = "RMSE",
                  method = "xgbTree",
                  objective = if(family == 'binomial'){"reg:logistic"}else if(family == 'gaussian'){"reg:linear"}else{message('Family of model not compatiable')},
                  weights = mydata$w)
  
  saveRDS(xg_fit, paste0(outputdir, "/xgboost/xg_fit.RDS"))
  
  write.csv(xg_fit$bestTune, paste0(outputdir, 'xgboost/xgboost_best_tune_.csv'))
  xg_best_tune <- xg_fit$bestTune
  
  xg_grid_final <- expand.grid(nrounds = xg_best_tune$nrounds,
                               max_depth = xg_best_tune$max_depth,
                               eta = xg_best_tune$eta,
                               colsample_bytree = .5,
                               min_child_weight = 1,
                               subsample = 1,
                               gamma = 0)
  
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  message("Fitting xgboost on final tuned hyperparameters")
  xg_fit_final <- train(form,
                        data = mydata,
                        trControl = train_control_final,
                        verbose = F,
                        tuneGrid = xg_grid_final,
                        metric = "RMSE",
                        method = "xgbTree",
                        objective = if(family == 'binomial'){"reg:logistic"}else if(family == 'gaussian'){"reg:linear"}else{message('Family of model not compatiable')},
                        weights = mydata$w)
  
  cov_plot <-
    ggplot(varImp(xg_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/xgboost/_covariate_importance.png'),
         plot = cov_plot)
  
  mydata[, 'xgboost_cv_pred'   := arrange(xg_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'xgboost_full_pred' := predict(xg_fit_final, mydata)]
  
  xg_fit_final$model_name <- "xgboost"
  saveRDS(xg_fit_final, paste0(outputdir, '/xgboost/full_xgboost.RDS'))
  
  covs[, 'xgboost' := predict(xg_fit_final, covs)]
  
  rm(form, cov_plot, train_control, train_control_final, xg_best_tune, xg_fit, xg_fit_final, xg_grid, xg_grid_final)
}

#~~~~~~~~~~~#
# 2. GAM ####
#~~~~~~~~~~~#
if('gam' %in% child_models){
  dir.create(paste0(outputdir, '/gam/'), showWarnings = F)
  
  if(family == 'binomial'){
    response <- cbind(sucesses = mydata$n, 
                      failures = mydata$d - mydata$n)
  } else if (family == 'gaussian'){
    response <- mydata$n
  }
  
  gam_formula <- paste0('response ~ 1+ s(', paste(covs_to_include, collapse = ", bs = 'ts', k = 3) + s("), ", bs = 'ts', k = 3)")
  gam_formula <- as.formula(gam_formula)
  
  full_gam = mgcv::gam(gam_formula, 
                       data = mydata, 
                       family = if(family =='binomial'){'quasibinomial'}else if(family == 'gaussian'){'gaussian'}, 
                       weights = mydata$w, 
                       control = list(nthreads = 2))
  full_gam$model_name = 'GAM'
  
  mydata[,'gam_full_pred' := predict(full_gam, mydata, type = 'response')]
  
  for(i in 1:5){
    if(family == 'binomial'){
      response <- cbind(successes = mydata$n[mydata$fold_id!=i], 
                        failures = mydata$d[mydata$fold_id!=i] - mydata$n[mydata$fold_id!=i])
    } else if (family == 'gaussian'){
      response <- mydata$n[mydata$fold_id!=i] 
    }      
    
    baby_gam = mgcv::gam(gam_formula, 
                         data = mydata[mydata$fold_id!=i], 
                         family = if(family =='binomial'){'quasibinomial'}else if(family == 'gaussian'){'gaussian'}, 
                         weights = mydata$weight[mydata$fold_id!=i], 
                         control = list(nthreads = 2))
    
    mydata[fold_id==i, 'gam_cv_pred' := predict(baby_gam, mydata[fold_id==i,],type = 'response')] 
    
  }
  
  saveRDS(full_gam, paste0(outputdir, '/gam/full_gam.RDS'))
  
  covs[,'gam' := predict(full_gam, covs, type = 'response')]
  
  pdf(paste0(outputdir, '/gam/plots2.pdf'))
  gam.check(full_gam)
  dev.off()
  
  rm(baby_gam, full_gam, gam_formula, response)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Penalised regression (E-net/Ridge/Lasso) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
if('enet' %in% child_models | 'ridge' %in% child_models | 'lasso' %in% child_models){
  
  dir.create(paste0(outputdir, '/glmnet'),showWarnings = F)
  
  if(family == 'binomial'){
    response <- cbind(failures   = mydata$d - mydata$n, 
                      successes = mydata$n)
  }else if(family == 'gaussian'){
    response <- mydata$n
  }
  
  vars <- as.matrix(mydata[, covs_to_include, with = F])
  colnames(vars) <- covs_to_include
  
  cv_lambda0 = cv.glmnet(x = vars , y= response, family = family, alpha = 0, weights = mydata$w, nfolds = 5, foldid = mydata$fold_id)
  cv_lambda0.25 = cv.glmnet(x = vars , y= response, family = family, alpha = 0.25, weights = mydata$w, nfolds = 5, foldid = mydata$fold_id)
  cv_lambda0.5 = cv.glmnet(x = vars , y= response, family = family, alpha = 0.5, weights = mydata$w, nfolds = 5, foldid = mydata$fold_id)
  cv_lambda0.75 = cv.glmnet(x = vars , y= response, family = family, alpha = 0.75, weights = mydata$w, nfolds = 5, foldid = mydata$fold_id)
  cv_lambda1 = cv.glmnet(x = vars , y= response, family = family, alpha = 1, weights = mydata$w, nfolds = 5, foldid = mydata$fold_id)
  
  pdf(paste0(outputdir, '/glmnet/parameter_selection.pdf'))
  par(mfrow=c(3,2))
  plot(cv_lambda0)
  plot(cv_lambda0.25)
  plot(cv_lambda0.5)
  plot(cv_lambda0.75)
  plot(cv_lambda1)
  plot(log(cv_lambda0$lambda),cv_lambda0$cvm,pch=19,col="red",xlab="log(Lambda)",ylab=cv_lambda0$name)
  points(log(cv_lambda0.25$lambda),cv_lambda0.25$cvm,pch=19,col="pink")
  points(log(cv_lambda0.5$lambda),cv_lambda0.5$cvm,pch=19,col="blue")
  points(log(cv_lambda0.75$lambda),cv_lambda0.75$cvm,pch=19,col="yellow")
  points(log(cv_lambda1$lambda),cv_lambda1$cvm,pch=19,col="green")
  legend("bottomright",legend=c("alpha= 1","alpha= .75", "alpha= .5", "alpha= .25","alpha 0"),pch=19,col=c("green","yellow","blue","pink","red"))
  dev.off()
  
  if('ridge' %in% child_models){full_ridge = glmnet(x = vars , y= response, family = family, alpha = 0, weights = mydata$w)}
  if('enet' %in% child_models){full_enet = glmnet(x = vars , y= response, family = family, alpha = 0.5, weights = mydata$w)}
  if('lasso' %in% child_models){full_lasso = glmnet(x = vars , y= response, family = family, alpha = 1, weights = mydata$w)}
  
  if('ridge' %in% child_models){mydata[,'ridge_full_pred' := predict(full_ridge,newx = vars, s = cv_lambda0$lambda.1se, type = 'response')]}
  if('lasso' %in% child_models){mydata[,'lasso_full_pred' := predict(full_lasso,newx = vars, s = cv_lambda1$lambda.1se, type = 'response')]}
  if('enet' %in% child_models){mydata[,'enet_full_pred' := predict(full_enet,newx = vars, s = cv_lambda0.5$lambda.1se, type = 'response')]}
  
  for(i in 1:5){
    if(family == 'binomial'){
      response <- cbind(failures = mydata$d[mydata$fold_id!=i] - mydata$n[mydata$fold_id!=i], 
                        successes = mydata$n[mydata$fold_id!=i])
    } else if(family == 'gaussian'){
      response <- mydata$n[mydata$fold_id!=i] 
    }
    
    vars <- as.matrix(mydata[fold_id != i, covs_to_include, with = F])
    colnames(vars) <- covs_to_include
    
    if('ridge' %in% child_models){baby_ridge = glmnet(x = vars , y= response, family = family, lambda = cv_lambda0$lambda.1se, alpha = 0, weights = mydata$w[mydata$fold_id!=i])}
    if('lasso' %in% child_models){baby_lasso = glmnet(x = vars , y= response, family = family, lambda = cv_lambda1$lambda.1se, alpha = 1, weights = mydata$w[mydata$fold_id!=i])}
    if('enet' %in% child_models){baby_enet = glmnet(x = vars , y= response, family = family, lambda = cv_lambda0.5$lambda.1se, alpha = 0.5, weights = mydata$w[mydata$fold_id!=i])}
    
    new_vars <- as.matrix(mydata[fold_id == i, covs_to_include, with = F])
    
    if('ridge' %in% child_models){mydata[fold_id==i,'ridge_cv_pred' := predict(baby_ridge,newx = new_vars, s = cv_lambda0$lambda.1se, type = 'response')]}
    if('lasso' %in% child_models){mydata[fold_id==i,'lasso_cv_pred' := predict(baby_lasso,newx = new_vars, s = cv_lambda1$lambda.1se, type = 'response')]}
    if('enet' %in% child_models){mydata[fold_id==i,'enet_cv_pred' := predict(baby_enet,newx = new_vars, s = cv_lambda0.5$lambda.1se, type = 'response')]}
  }
  
  if('ridge' %in% child_models){saveRDS(cv_lambda0, paste0(outputdir, '/glmnet/full_ridge.rds'))}
  if('enet' %in% child_models){saveRDS(cv_lambda0.5, paste0(outputdir, '/glmnet/full_enet.rds'))}
  if('lasso' %in% child_models){saveRDS(cv_lambda1, paste0(outputdir, '/glmnet/full_lasso.rds'))}
  
  all_names <- names(covs) 
  new_covs <- as.matrix(covs)
  names(new_covs) <- all_names
  if('ridge' %in% child_models){covs[,'ridge' := predict(full_ridge,newx = new_covs[,rownames(full_ridge$beta)], s = cv_lambda0$lambda.1se, type = 'response')]}
  if('enet' %in% child_models){covs[,'enet' := predict(full_enet,newx = new_covs[,rownames(full_enet$beta)], s = cv_lambda0.5$lambda.1se, type = 'response')]}
  if('lasso' %in% child_models){covs[,'lasso' := predict(full_lasso,newx = new_covs[,rownames(full_lasso$beta)], s = cv_lambda1$lambda.1se, type = 'response')]}
  
  rm(cv_lambda1, cv_lambda0.5, cv_lambda0, cv_lambda0.25, cv_lambda0.75, full_lasso, full_enet, full_ridge, baby_lasso, baby_ridge, baby_enet, new_vars, response, vars, i, new_covs, all_names)    
}

#~~~~~~~~~~~~~~~~~~~~~#
# 4. Random forest ####
#~~~~~~~~~~~~~~~~~~~~~#
if('rf' %in% child_models){
  dir.create(paste0(outputdir, '/rf'), showWarnings = F)
  
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <- expand.grid(.mtry=c(1:length(covs_to_include)))
  
  rf_fit <- train(form,
                  data = mydata,
                  trControl = train_control,
                  verbose = T,
                  tuneGrid = tunegrid,
                  metric = "RMSE",
                  method = "rf",
                  weights = mydata$w)
  
  saveRDS(rf_fit, paste0(outputdir, "/rf/rf_fit.RDS"))
  png(paste0(outputdir, '/rf/rf_fit.png'))
  plot(rf_fit)  
  dev.off()
  
  mtry_tune <- rf_fit$bestTune$mtry
  write.csv(mtry_tune, paste0(outputdir, '/rf/rf_params.csv'), row.names = F)
  tunegrid_final <- expand.grid(.mtry=mtry_tune)
  
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  rf_fit_final <- train(form,
                        data = mydata,
                        trControl = train_control_final,
                        verbose = T,
                        tuneGrid = tunegrid_final,
                        metric = "RMSE",
                        method = "rf",
                        importance=T,
                        weights = mydata$w)
  
  mydata[, 'rf_cv_pred'   := arrange(rf_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'rf_full_pred' := predict(rf_fit_final, mydata)]
  
  rf_fit_final$model_name <- "rf"
  saveRDS(rf_fit_final, paste0(outputdir, '/rf/full_rf.RDS'))
  
  cov_plot <-
    ggplot(varImp(rf_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/rf/rf_covariate_importance.png'),
         plot = cov_plot)
  
  covs[, 'rf' := predict(rf_fit_final, covs)]
  covs$rf <- ifelse(covs$rf < 0, 0, covs$rf)
  
  rm(form, cov_plot, train_control, train_control_final, rf_fit, rf_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Neural networks ####
#~~~~~~~~~~~~~~~~~~~~~~~#
if('nnet' %in% child_models){
  dir.create(paste0(outputdir, '/nnet'), showWarnings = F)
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
    linoutval = FALSE
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
    linoutval = TRUE
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <- expand.grid(.decay = c(1, 0.5, 0.1, 0.01, 0.001, 0.0001, 0.00001), .size = c(4, 5, 6, 7, 8, 9))
  
  nn_fit <- train(form,
                  data = mydata,
                  trControl = train_control,
                  verbose = T,
                  tuneGrid = tunegrid,
                  metric = "RMSE",
                  method = "nnet",
                  linout = linoutval,
                  maxit = 1000,
                  weights = mydata$w)
  
  saveRDS(nn_fit, paste0(outputdir, "/nnet/nn_fit.RDS"))
  png(paste0(outputdir, '/nnet/nn_fit.png'))
  plot(nn_fit)  
  dev.off()
  
  write.csv(nn_fit$bestTune, paste0(outputdir, '/nnet/nnet_best_tune.csv'))
  
  tunegrid_final <- expand.grid(.decay=nn_fit$bestTune$decay, .size=nn_fit$bestTune$size)
  
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
    nn_fit_final <- train(form,
                        data = mydata,
                        trControl = train_control_final,
                        verbose = T,
                        tuneGrid = tunegrid_final,
                        metric = "RMSE",
                        method = "nnet",
                        linout = linoutval,
                        maxit = 1000,
                        weights = mydata$w)
  
  mydata[, 'nnet_cv_pred'   := arrange(nn_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'nnet_full_pred' := predict(nn_fit_final, mydata)]
  
  nn_fit_final$model_name <- "nn"
  saveRDS(nn_fit_final, paste0(outputdir, '/nnet/full_nn.RDS'))
  
  cov_plot <-
    ggplot(varImp(nn_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/nnet/nn_covariate_importance.png'),
         plot = cov_plot)
  
  covs[, 'nnet' := predict(nn_fit_final, covs)]
  
  rm(nn_fit, nn_fit_final, train_control, train_control_final, tunegrid, tunegrid_final, cov_plot)
}

#~~~~~~~~~~~~~~~~~~~~#
# 6. Cubist model ####
#~~~~~~~~~~~~~~~~~~~~#
if('cubist' %in% child_models){
  dir.create(paste0(outputdir, '/cubist'), showWarnings = F)
  
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(covs_to_include, collapse = " + ")))
  } else if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(covs_to_include, collapse = " + ")))
  }
  
  train_control <- trainControl(selectionFunction = "best",
                                method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                index = list(mydata$a_rowid[mydata$fold_id!=1],
                                             mydata$a_rowid[mydata$fold_id!=2],
                                             mydata$a_rowid[mydata$fold_id!=3],
                                             mydata$a_rowid[mydata$fold_id!=4],
                                             mydata$a_rowid[mydata$fold_id!=5]),
                                indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                               mydata$a_rowid[mydata$fold_id==2],
                                               mydata$a_rowid[mydata$fold_id==3],
                                               mydata$a_rowid[mydata$fold_id==4],
                                               mydata$a_rowid[mydata$fold_id==5]),
                                search = 'grid')
  
  tunegrid <- expand.grid(.committees = c(seq(1, 40, 5)), 
                          .neighbors  = c(0, 3, 6, 9))
  
  cubist_fit <- train(form,
                      data = mydata,
                      trControl = train_control,
                      verbose = T,
                      tuneGrid = tunegrid,
                      metric = "RMSE",
                      method = "cubist",
                      control = Cubist::cubistControl(),
                      weights = mydata$w)
  
  saveRDS(cubist_fit, paste0(outputdir, "/cubist/cubist_fit.RDS"))
  png(paste0(outputdir, '/cubist/cubist_fit.png'))
  plot(cubist_fit)  
  dev.off()
  
  write.csv(cubist_fit$bestTune, paste0(outputdir, '/cubist/cubist_best_tune.csv'))
  
  tunegrid_final <- expand.grid(.committees=cubist_fit$bestTune$committees, .neighbors=cubist_fit$bestTune$neighbors)
  
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  
  cubist_fit_final <- train(form,
                            data = mydata,
                            trControl = train_control_final,
                            verbose = T,
                            tuneGrid = tunegrid_final,
                            metric = "RMSE",
                            method = "cubist",
                            weights = mydata$w)
  
  mydata[, 'cubist_cv_pred'   := arrange(cubist_fit_final$pred, rowIndex)[,"pred"]]
  mydata[, 'cubist_full_pred' := predict(cubist_fit_final, mydata)]
  
  mydata$cubist_cv_pred <- ifelse(mydata$cubist_cv_pred > 1, 1, mydata$cubist_cv_pred)
  mydata$cubist_full_pred <- ifelse(mydata$cubist_full_pred > 1, 1, mydata$cubist_full_pred)
  
  cubist_fit_final$model_name <- "cubist"
  saveRDS(cubist_fit_final, paste0(outputdir, '/cubist/full_cubist.RDS'))
  
  cov_plot <-
    ggplot(varImp(cubist_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(outputdir, '/cubist/cubist_covariate_importance.png'),
         plot = cov_plot)
  
  covs[, 'cubist' := predict(cubist_fit_final, covs)]
  covs$cubist <- ifelse(covs$cubist > 1, 1, covs$cubist)
  
  rm(form, cov_plot, train_control, train_control_final, cubist_fit, cubist_fit_final, tunegrid, tunegrid_final)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Print out correlations between data and predictions ####
# this is to aid selection of child models               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

table_r2 <-c()
for(i in 1:length(child_models)){
  if(family == 'binomial'){
    r2 <- c(child_models[i], formatC(cor(mydata$p, mydata[,get(paste0(child_models[i], '_cv_pred'))])^2), format = 'f', digits = 6)
    table_r2 <- rbind(table_r2,r2)
    message(paste0(child_models[i], ' correlation: ', round(cor(mydata$p, mydata[,get(paste0(child_models[i], '_cv_pred'))])^2,2)))
  }
  if(family == 'gaussian'){
    r2 <- c(child_models[i], formatC(cor(mydata$n, mydata[,get(paste0(child_models[i], '_cv_pred'))])^2), format = 'f', digits = 6)
    table_r2 <- rbind(table_r2,r2)
    message(paste0(child_models[i], ' correlation: ', round(cor(mydata$n, mydata[,get(paste0(child_models[i], '_cv_pred'))])^2,2)))
  }
}
colnames(table_r2) <- c('model','r2','format','digits')
table_r2 <- table_r2[,c('model','r2')]
table_r2 <- data.table(table_r2)
table_r2$r2 <- as.numeric(table_r2$r2)
table_r2 <- table_r2[order(-table_r2$r2),]
write.csv(table_r2, paste0(outputdir, '/table_r2.csv'), row.names = F)
table_r2 <- na.omit(table_r2)
table_r2 <- table_r2[order(-table_r2$r2),]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Check the correlation of the stackers ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

child_models <- unique(table_r2$model)
final_child_models <- unique(table_r2$model)
table_cor <- c()
for(i in 1:length(child_models)){
  for(j in 1:length(child_models)){
    if(i==j){
    } else{
      test1 <- sd(mydata[,get(paste0(child_models[i], '_cv_pred'))])
      test2 <- sd(mydata[,get(paste0(child_models[j], '_cv_pred'))])
      if(test1 > 0 & test2 > 0) {
        cor1 <- c(child_models[i],child_models[j],formatC(cor(mydata[,get(paste0(child_models[i], '_cv_pred'))],mydata[,get(paste0(child_models[j], '_cv_pred'))])^2, format = 'f', digits = 6))
        delete <- ifelse(cor1[[3]] < 0.8,"",ifelse(table_r2$r2[table_r2$model == child_models[i]] < table_r2$r2[table_r2$model == child_models[j]],child_models[i],child_models[j]))
        final_child_models <- final_child_models[!final_child_models == delete]
        row1 <- c(cor1,delete)
        table_cor <- rbind(table_cor,row1)
        rm(cor1, delete, row1)
        if(cor(mydata[,get(paste0(child_models[i], '_cv_pred'))],mydata[,get(paste0(child_models[j], '_cv_pred'))])^2>0.80){message(paste0(child_models[i],  ' and ', child_models[j], ' correlated, remove one'))
        }
      }
    }
  }
}

write.csv(table_cor, paste0(outputdir, '/table_correlation.csv'), row.names = F)
write.csv(final_child_models, paste0(outputdir, '/table_child_models.csv'), row.names = F)

final_child_columns <- c(paste(paste0(final_child_models,'_cv_pred'),sep = ","), paste(paste0(final_child_models,'_full_pred'),sep=","))
child_models <- final_child_models

if(centre_scale == TRUE){
  mydata$year <-  NULL
  covs$year <- NULL
}
write.csv(mydata, paste0(outputdir, '/fitted_child_models.csv'), row.names = F)
write.csv(covs, paste0(outputdir, '/child_model_preds.csv'), row.names = F)

mydata <- read.csv(paste0(outputdir, '/fitted_child_models.csv'), stringsAsFactors =  F)
covs <- read.csv(paste0(outputdir, '/child_model_preds.csv'), stringsAsFactors = F)
mydata <- data.table(mydata)
covs <- data.table(covs)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Combined the child model estimates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
headcols <- colnames(mydata)[!colnames(mydata) %in% colnames(mydata)[grep('pred',colnames(mydata))]]
headcols <- c(headcols,final_child_columns)
stackers <- subset(mydata, select = headcols, drop.levels=TRUE)
stackers <- data.frame(stackers)
X <- as.matrix(stackers[colnames(stackers)[(grep('cv_pred', colnames(stackers)))]])
Y = if(family == 'binomial'){stackers$p}else if(family == 'gaussian'){stackers$n}

C <- data.frame(covs)
C <- as.matrix(C[c(child_models)])


if(stacker == 'CWM'){
  s <- solve.QP( 
    t(X) %*% X, t(Y) %*% X, 
    cbind(  
      matrix(1, nr=length(child_models), nc=1),
      diag(length(child_models)),
      -diag(length(child_models))
    ),
    c(1, 
      rep(0.000001, length(child_models)),
      rep(-1, length(child_models))), 
    meq = 1
  )
  
  mydata$stacked_preds <- rowWeightedMeans(X, w = s$solution) 
  
  covs$cv_custom_stage_1 <- rowWeightedMeans(C, w = s$solution)
}

if(stacker == 'RWM'){
  r2 <- rep(NA, length(child_models))
  for(i in 1:length(child_models)){
    if(family == 'binomial'){
      r2[i] <- round(cor(mydata$p, mydata[,get(paste0(child_models[i], '_cv_pred'))])^2,2)
    }
    if(family == 'gaussian'){
      r2[i] <- round(cor(mydata$n, mydata[,get(paste0(child_models[i], '_cv_pred'))])^2,2)
    }
  }
  
  total <-  sum(r2)
  weights <- r2/total
  mydata$stacked_preds <- rowWeightedMeans(X, w = weights)   
  covs$cv_custom_stage_1 <- rowWeightedMeans(C, w = weights)
}

if(stacker == 'GBM'){
  if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(colnames(stackers)[(grep('cv_pred', colnames(stackers)))], collapse = " + ")))
  }
  if(family == 'binomial'){
    form <- as.formula(paste0('logit(p) ~ ', paste(colnames(stackers)[(grep('cv_pred', colnames(stackers)))], collapse = " + ")))
  }
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  model_gbm<- 
    train(form, data = mydata, method='gbm', trControl=train_control_final, tuneLength=3)
  
  mydata[, 'stacked_preds'   := inv.logit(predict(model_gbm, mydata))]
  covs <- data.frame(covs)
  colnames(covs)[colnames(covs) %in% child_models] <- colnames(stackers)[(grep('cv_pred', colnames(stackers)))]
  covs$cv_custom_stage_1 <- inv.logit(predict(model_gbm, covs[colnames(covs)[(grep('cv_pred', colnames(covs)))]])) 
  covs <- data.table(covs)
  
  jpeg(paste0(outputdir, 'stacker_cov_imporatance.jpeg'),
       height = 10, width = 10, units = 'cm', res = 150)
  ggplot(varImp(model_gbm, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  
  dev.off()
}
if(stacker == 'GLM'){
  if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(colnames(stackers)[(grep('cv_pred', colnames(stackers)))], collapse = " + ")))
  }
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(colnames(stackers)[(grep('cv_pred', colnames(stackers)))], collapse = " + ")))
  }
  
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  model_glm<- 
    train(form, data = mydata, method='glm', trControl=train_control_final, tuneLength=3)
  
  mydata[, 'stacked_preds'   := predict(model_glm, mydata)]
  
  covs <- data.frame(covs)
  colnames(covs)[colnames(covs) %in% child_models] <- colnames(stackers)[(grep('cv_pred', colnames(stackers)))]
  covs$cv_custom_stage_1 <- predict(model_glm, covs[colnames(covs)[(grep('cv_pred', colnames(covs)))]]) 
  covs <- data.table(covs)
}

if(stacker == 'nnet'){
  if(family == 'gaussian'){
    form <- as.formula(paste0('n ~ ', paste(colnames(stackers)[(grep('cv_pred', colnames(stackers)))], collapse = " + ")))
  }
  if(family == 'binomial'){
    form <- as.formula(paste0('p ~ ', paste(colnames(stackers)[(grep('cv_pred', colnames(stackers)))], collapse = " + ")))
  }  
  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final",
                                      index = list(mydata$a_rowid[mydata$fold_id!=1],
                                                   mydata$a_rowid[mydata$fold_id!=2],
                                                   mydata$a_rowid[mydata$fold_id!=3],
                                                   mydata$a_rowid[mydata$fold_id!=4],
                                                   mydata$a_rowid[mydata$fold_id!=5]),
                                      indexOut =list(mydata$a_rowid[mydata$fold_id==1],
                                                     mydata$a_rowid[mydata$fold_id==2],
                                                     mydata$a_rowid[mydata$fold_id==3],
                                                     mydata$a_rowid[mydata$fold_id==4],
                                                     mydata$a_rowid[mydata$fold_id==5]))
  
  model_nnet<- 
    train(form, data = mydata, method='nnet', trControl=train_control_final, tuneLength=3)
  
  mydata[, 'stacked_preds'   := predict(model_nnet, mydata)]
  
  covs <- data.frame(covs)
  colnames(covs)[colnames(covs) %in% child_models] <- colnames(stackers)[(grep('cv_pred', colnames(stackers)))]
  covs$cv_custom_stage_1 <- predict(model_gbm, covs[colnames(covs)[(grep('cv_pred', colnames(covs)))]]) 
  covs <- data.table(covs)
}



stg1 <- covs[, .(location_id, year_id,
                 age_group_id = rep(22, length(covs$location_id)),
                 sex_id = rep(3, length(covs$location_id)),
                 cv_custom_stage_1)]

mydata[, colnames(mydata)[grep('^cv_', colnames(mydata))] := NULL]

max(stg1$cv_custom_stage_1, na.rm = T)

stg1$cv_custom_stage_1[stg1$cv_custom_stage_1 == 0] <- 0.0001
stg1$cv_custom_stage_1[stg1$cv_custom_stage_1 == 1] <- 0.9999

write.csv(mydata, paste0(outputdir, '/fitted_stackers.csv'), row.names = F)
write.csv(stg1, paste0(outputdir, '/custom_stage1_df.csv'), row.names = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Calculate metrics for all models ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
child_model_metrics <- data.frame(child_models)
child_model_metrics$r_sq <- NA

mydata <- data.frame(mydata)
for(i in 1:length(child_models)){
  cm <- paste0(child_models, '_cv_pred')[i]
  pred <- mydata[c(cm)]
  pred <-  unlist(pred)
  child_model_metrics$rmse[child_model_metrics$child_models == child_models[i]] <- round(RMSE(pred, if(family == 'binomial'){mydata$p}else if(family == 'gaussian'){mydata$n}),4)
  child_model_metrics$r_sq[child_model_metrics$child_models == child_models[i]] <- round(cor(pred, if(family == 'binomial'){mydata$p}else if(family == 'gaussian'){mydata$n})^2,2)
}

child_model_metrics$child_models <- as.character(child_model_metrics$child_models)
child_model_metrics <- rbind(child_model_metrics, c('Stackers', NA, NA) )
child_model_metrics$rmse[child_model_metrics$child_models=='Stackers'] <- round(RMSE(mydata$stacked_preds, if(family == 'binomial'){mydata$p}else if(family == 'gaussian'){mydata$n}),4)
child_model_metrics$r_sq[child_model_metrics$child_models=='Stackers'] <- round(cor(mydata$stacked_preds, if(family == 'binomial'){mydata$p}else if(family == 'gaussian'){mydata$n})^2,2)

write.csv(child_model_metrics, paste0(outputdir, '/national_stacker_metrics.csv'), row.names = F)

############################
# Config template for stgpr
###########################
if (!oos){
  config_table <- c()
  for(x in list) {
    y <- substring(x,1,stringr::str_locate(x,"-")-1)[1]
    z <- substring(x,stringr::str_locate(x,"-")+1,)[1]
    description <- paste0(x)
    path_to_data <-  paste0(FILEPATH,y,'/', z,'/input.csv') 
    path_to_custom_stage_1 <-  paste0(FILEPATH,y,'/', z,'/custom_stage1_df.csv') 
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
    
    dat <- read.csv(path_to_data, stringsAsFactors = FALSE)
    year_start <- min(dat$year_id[dat$is_outlier==0])
    by_loc_sum = group_by(dat, location_id) %>%
      summarize(
        num_dp = length(val[is_outlier==0])
      )
    density_cutoffs <- max(2, ceiling(quantile(by_loc_sum$num_dp, 0.25))[[1]])
    
    config <- c(modelable_entity_id,gbd_round_id,decomp_step,data_transform,prediction_units,st_lambda,st_omega,st_zeta,gpr_scale,path_to_data,path_to_custom_stage_1,location_set_id,year_start,year_end,prediction_age_group_ids,prediction_sex_ids,description,gpr_amp_method,gpr_draws,density_cutoffs)
    config_table <- rbind(config_table,config)
  }
  
  config_table <- as.data.frame(config_table)
  names(config_table) <- c('modelable_entity_id','gbd_round_id','decomp_step','data_transform','prediction_units','st_lambda','st_omega','st_zeta','gpr_scale','path_to_data','path_to_custom_stage_1','location_set_id','year_start','year_end','prediction_age_group_ids','prediction_sex_ids','description','gpr_amp_method','gpr_draws','density_cutoffs')
  write.csv(config_table,FILEPATH,row.names = F)
}
