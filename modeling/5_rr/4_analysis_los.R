#######
# This script 1) appends the information on relative length of stay for resistant infections
# from literature and primary data,
# 2) sets up model function for MRBRT and
# 3) outputs estimates and draws 
#########
# 0) Define antibiotic classes and pathogens

organism <- c('acinetobacter_baumanii','citrobacter_spp','enterobacter_spp','enterococcus_faecalis','enterococcus_faecium','enterococcus_spp','escherichia_coli','group_a_strep','group_b_strep','haemophilus_influenzae','klebsiella_pneumoniae','moraxella_spp','morganella_spp','neisseria_meningitidis','non_typhoidal_salmonellae','proteus_spp','providencia_spp','pseudomonas_aeruginosa','pseudomonas_spp','neisseria_gonorrheae','shigella_spp','mycobacterium_tuberculosis','salmonella_paratyphi','salmonella_typhi','salmonella_typhi_paratyphi','serratia_spp','staphylococcus_aureus','streptococcus_pneumoniae','mycoplasma','listeria','legionella_spp')
antibiotic_group<-c('third_gen_ceph','carbapenem','fluoroquinolone','penicillin','aminopenicillin','beta_lactamase_inhibitor','anti_pseudomonal_penicillin','methicillin','vancomycin','fourth_gen_ceph','sulfa','aminoglycoside','macrolide')

# 1) Define paths, load packages, functions, data
library(dplyr)
library(data.table)
library(mrbrt002, lib.loc = "FILEPATH")
library(readxl)
library(metafor)
library(msm)
source("FILENAME")
source("FILENAME")

get_pval <- function(beta, beta_sd, one_sided = FALSE) {
  zscore <- abs(beta/beta_sd)
  if (one_sided) 1 - pnorm(zscore)
  if (!one_sided) (1 - pnorm(zscore))*2
}

#location metadata
location_md <- subset(get_location_metadata(location_set_id=35, gbd_round_id=6),select = c(location_id, location_name), droplevels = TRUE)

# define paths
data_path <- 'FILEPATH'
model_results <- 'FILEPATH'
output <- 'FILEPATH'

#load data
DATA_A <-as.data.frame(readxl::read_xlsx("FILENAME"))
DATA_B <- read.csv("FILENAME", stringsAsFactors = F) %>% rename(author=refvar,yi = coeff_adjusted_days, sei = se_adjusted_days, cases_r=admit_r, cases_s = admit_s) 
DATA_C <- read.csv("FILENAME", stringsAsFactors = F) %>% rename(author=refvar,yi = coeff_adjusted_days, sei = se_adjusted_days, cases_r=admit_r, cases_s = admit_s) 
DATA_D <- read.csv("FILENAME", stringsAsFactors = F) %>% rename(author=refvar,yi = coeff_adjusted_days, sei = se_adjusted_days, cases_r=admit_r, cases_s = admit_s) 

# For some studies, median and interquartile range inform mean and standard deviation 
DATA_A$sd_r[(is.na(DATA_A$sd_r) & !is.na(DATA_A$half_IQR_r))] <- DATA_A$half_IQR_r[(is.na(DATA_A$sd_r) & !is.na(DATA_A$half_IQR_r))]
DATA_A$sd_s[(is.na(DATA_A$sd_s) & !is.na(DATA_A$half_IQR_s))] <- DATA_A$half_IQR_s[(is.na(DATA_A$sd_s) & !is.na(DATA_A$half_IQR_s))]
DATA_A$mean_LOS_r[(is.na(DATA_A$mean_LOS_r) & !is.na(DATA_A$median_LOS_r))] <- DATA_A$median_LOS_r[(is.na(DATA_A$mean_LOS_r) & !is.na(DATA_A$median_LOS_r))]
DATA_A$mean_LOS_s[(is.na(DATA_A$mean_LOS_s) & !is.na(DATA_A$median_LOS_s))] <- DATA_A$median_LOS_s[(is.na(DATA_A$mean_LOS_s) & !is.na(DATA_A$median_LOS_s))]
# Omit missing data
DATA_A <- DATA_A[!is.na(DATA_A$sd_r),]
DATA_A <- DATA_A[,c('pathogen','abx_class','author','location','year_start','year_end','sample_size','cases_s','cases_r','mean_LOS_r','sd_r','mean_LOS_s','sd_s')]

# No correlation between resistant vs susceptible observations is assumed
DATA_A$ri<-0
# Use metafor to calculate the raw mean difference and uncertainty 
DATA_A <- escalc(data=DATA_A,measure="ROMC",m1i=mean_LOS_r,sd1i=sd_r,m2i=mean_LOS_s,sd2i = sd_s,ni=sample_size,ri= ri, append = TRUE)
DATA_A$sei <- sqrt(DATA_A$vi)
DATA_A <- DATA_A[!is.na(DATA_A$sei),]
DATA_A <- DATA_A %>% left_join(location_md, by = c('location_name'))
DATA_A$measure2 <- 'unadjusted'
DATA_A <- DATA_A[,c('pathogen','abx_class','author','location_id','sterile','cases_s','cases_r','mean_LOS_r','sd_r','mean_LOS_s','sd_s','measure2','yi','sei','vi')]

# append together all sources of data
data_or <- rbind(DATA_B, DATA_C, DATA_D)
data_or <- data_or[,c('pathogen','abx_class','author','location_id','sterile','cases_s','cases_r','mean_LOS_r','sd_r','mean_LOS_s','sd_s','measure2','yi','sei')]
data_or$vi <- data_or$sei^2
data_or <- rbind(data_or,DATA_A)
#Identify adjusted vs unadjusted
data_or$prior_days_unadjusted <- ifelse(data_or$measure2 == 'unadjusted',1,0)

# create sample size from resistant and susceptible cases
data_or$sample_size <- data_or$cases_r + data_or$cases_s

# correct naming of pathogens
data_or$pathogen[data_or$pathogen == 'acinetobacter_baumannii'] <- 'acinetobacter_baumanii'

#2) Sets up model function 
# outputs estimates and draws 
# A loop that goes through the antibiotic classes 
antibiotic_group <- c('third_gen_ceph','carbapenem','penicillin','methicillin','vancomycin','fluoroquinolone','aminopenicillin','beta_lactamase_inhibitor','anti_pseudomonal_penicillin','sulfa','aminoglycoside','macrolide','fourth_gen_ceph','isoniazid','rifampicin')

prior_draws <- c()
table_priors <- c()
for (j in antibiotic_group) { 
    data <- data_or[((data_or$abx_class == paste0(j)) & (data_or$sterile == 1)),]
    orgs<-paste(unique(data$pathogen),sep = ',')
    for (i in orgs) {
      data$new <- ifelse(data$pathogen == i,1,0)
      colnames(data)[colnames(data) == 'new'] <- paste0(i)
    }

    # initialize MRBRT for prior
    dat_prior <- MRData()
    dat_prior$load_df(
      data = data,  col_obs = "yi", col_obs_se = "sei",
      col_study_id = "author" )
    
    mrbrtprior <- MRBRT(
      data = dat_prior,
      cov_models = list(LinearCovModel('intercept', use_re = TRUE)))
    
    mrbrtprior$fit_model(inner_print_level = 5L, inner_max_iter = 1000L)
    
    # Create samples for prior draws
    n_samples <- 1000L
    prior_samples <- mrbrtprior$sample_soln(sample_size = n_samples)
    prior_samples <- as.data.frame(exp(prior_samples[[1]]))
    prior_samples <- t(prior_samples)
    prior_samples <- as.data.frame(prior_samples)
    colnames(prior_samples) <- paste(lapply(0:999, function(x) paste0("draw_", x)))
    prior_samples$abx_class <- paste0(j)
    prior_draws <- rbind(prior_draws,prior_samples)
    overall_mean <- mrbrtprior$beta_soln[1] 
    
    eff <- c('prior_beta_gaussian', overall_mean, 0.1)
    report <- c(1,j,overall_mean,0.1) 
    table_priors <- rbind(table_priors,report) 
    rm(report)
    
    # priors' matrix 
    priors <- c()
    for(i in orgs) {
      effs <- c(i,eff)
      priors <- rbind(priors, effs)
    }
    
    # concatenate priors 
    definition_model <- paste0("LinearCovModel(\"",priors[,1], "\",",priors[,2]," = array(c(", as.double(priors[,3]), ",",as.double(priors[,4]),")))")
    model_definition <- paste(definition_model, collapse = ",\n")
    

    # 2) Covariate selection using CovFinder 
    # adjusted estimates for 5 antibiotic classes 
    if (paste0(j) %in% c('third_gen_ceph','carbapenem','penicillin','methicillin','vancomycin','macrolide','beta_lactamase_inhibitor','anti_pseudomonal_penicillin')) {
        candidate_covs <- c('prior_days_unadjusted') 
        preselected_covs <- orgs 
        cov_names <- c(candidate_covs,preselected_covs)
        
        # load data
        mrdata <- MRData()
        mrdata$load_df(
          data = data, 
          col_obs = 'yi', col_obs_se = 'sei', 
          #the randome effect will be for each bug drug-combination
          col_study_id = 'pathogen', col_covs = as.list(cov_names))
        # run CovFinder
        covfinder <- CovFinder(
          data = mrdata, 
          covs = as.list(candidate_covs),
          pre_selected_covs = as.list(preselected_covs), 
          normalized_covs = TRUE,  num_samples = 6000L,   power_range = list(-4, 4),   power_step_size = 1,   laplace_threshold = 1e-5)
        covfinder$select_covs(verbose = TRUE)
        
        selected_covs <- covfinder$selected_covs
        if (length(orgs) != length(selected_covs)) {
        intercept_and_bias <- paste0("LinearCovModel('intercept', use_re = TRUE, prior_beta_uniform = array(c(0,0))), LinearCovModel(\'", 
                                     selected_covs[!selected_covs %in% preselected_covs], "'),")
        } else {
          intercept_and_bias <- "LinearCovModel('intercept', use_re = TRUE, prior_beta_uniform = array(c(0,0))),"
        }
    } else {
    
    # For unadjusted combinations we only include intercept and priors into the model
      intercept_and_bias <- "LinearCovModel('intercept', use_re = TRUE, prior_beta_uniform = array(c(0,0))),"
      selected_covs <- orgs
      cov_names <- orgs
      }
    
    # We load the data onto MRBRT
    dat1 <- MRData()
    dat1$load_df(
      data = data,  col_obs = "yi", col_obs_se = "sei",
      col_covs = as.list(cov_names),
      col_study_id = "pathogen" )
    
    # Text to parse into the MRBRT model
    mrbrtschema <- paste0(
      "MRBRT(
        data = dat1,
        cov_models = list(", intercept_and_bias,
        model_definition,
        "))"
    )
    
    # initialize MRBRT
    mod1 <- eval(parse(text = paste(mrbrtschema)))
    mod1$fit_model(inner_print_level = 5L, inner_max_iter = 1000L)
    
    # Create samples and draws to export
    n_samples <- 1000L
    samples <- mod1$sample_soln(sample_size = n_samples)
    draws <- as.data.frame(exp(samples[[1]]))
    names(draws) <- mod1$cov_model_names
    draws <- t(draws)
    draws <- as.data.frame(draws)
    colnames(draws) <- paste(lapply(0:999, function(x) paste0("draw_", x)))
    write.csv(draws,paste0(output,'FILEPATH','_',j,'FILENAME'), row.names = TRUE)
    rm(draws)
    
    # summary table to export
    table <- c()
    vector <-c()
    for (x in 1:(length(selected_covs)+1)) {
      beta <- mod1$beta_soln[x] 
      beta_sd <- sd(samples[[1]][,x]) 
      pval <- get_pval(beta,beta_sd)
      vector <- c(mod1$cov_model_names[x],exp(beta),ifelse(pval<0.05,"*",""),beta,beta_sd,c(beta-(1.96*beta_sd)),c(beta+(1.96*beta_sd)),pval)
      table <- rbind(table,vector)
      vector <-c()
    }
    table <- as.data.frame(table)
    names(table) <- c('pathogen','mean_RR','pval05','mean_ln(RR)','SE','lb','ub','pval')
    summary <- data %>% group_by(pathogen) %>% summarize (data_points = n(), sample_size = sum(sample_size)) 
    table <- right_join(summary,table,by = 'pathogen')
    table$abx_class<-paste0(j)
    
    #I export summary table into the summary folder
    write.csv(table,paste0(output,'FILEPATH/',j,'FILENAME'), row.names = FALSE)
    rm(table,samples,mod1)
}


# 3) Write and compile summary and draws estimates

write.csv(table_priors,paste0(output,'FILENAME'), row.names = FALSE)
write.csv(prior_draws,paste0(output,'FILENAME'), row.names = TRUE)

#appending draws together for loading into AMRCompiler
draws_sterile <- c() 
for (j in antibiotic_group) {
  test<-read.csv(paste0(output,'FILEPATH',j,'FILENAME'), stringsAsFactors = FALSE)
  test$abx_class <- paste0(j)
  draws_sterile <- rbind(draws_sterile,test)
}

names(draws_sterile)[names(draws_sterile) == 'X'] <- 'pathogen'
draws_sterile <- draws_sterile[!draws_sterile$pathogen %in% c('intercept','prior_days_unadjusted'),]

priordraws <- as.data.frame(read.csv(paste0(output,'FILENAME'), stringsAsFactors = FALSE ))

# append antibiotic class draws for combinations with no estimate
faecalis <- priordraws[priordraws$abx_class == 'aminoglycoside',]
faecalis$pathogen = 'enterococcus_faecalis'
faecalis$X <- NULL
groupb <- priordraws[priordraws$abx_class == 'fluoroquinolone',]
groupb$pathogen <- 'group_b_strep'
groupb$X <- NULL

draws_sterile<-rbind(draws_sterile,faecalis,groupb)
write.csv(draws_sterile,paste0(output,'FILENAME'), row.names = FALSE)
