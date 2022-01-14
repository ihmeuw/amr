# MRBRT relative risk model 
#0) Define antibiotic classes and pathogens
#1) Define paths, packages, functions, data
#2) Sets up model function and outputs 
##################################################
rm(list = ls())

# 0) Defines antibiotic class and pathogens

organism <- c('acinetobacter_baumanii','citrobacter_spp','enterobacter_spp','enterococcus_faecalis','enterococcus_faecium','enterococcus_spp','escherichia_coli','group_a_strep','group_b_strep','haemophilus_influenzae','klebsiella_pneumoniae','moraxella_spp','morganella_spp','neisseria_meningitidis','non_typhoidal_salmonellae','proteus_spp','providencia_spp','pseudomonas_aeruginosa','pseudomonas_spp','neisseria_gonorrheae','shigella_spp','mycobacterium_tuberculosis','salmonella_paratyphi','salmonella_typhi','salmonella_typhi_paratyphi','serratia_spp','staphylococcus_aureus','streptococcus_pneumoniae','mycoplasma','listeria','legionella_spp')
antibiotic_group<-c('third_gen_ceph','carbapenem','fluoroquinolone','penicillin','aminopenicillin','beta_lactamase_inhibitor','anti_pseudomonal_penicillin','methicillin','vancomycin','fourth_gen_ceph','sulfa','aminoglycoside','macrolide')

# define whether the model will use priors per antibiotic class
with_priors <- TRUE

# 1) Define paths, load packages, functions, data
require(dplyr,readxl,metafor,data.table)
library(dplyr)
library(readxl)
library(metafor)
library(data.table)
library(mrbrt001, lib.loc = "FILEPATH")
library(msm, lib.loc = "FILEPATH")
source("FILENAME")
source("FILENAME")

get_pval <- function(beta, beta_sd, one_sided = FALSE) {
  zscore <- abs(beta/beta_sd)
  if (one_sided) 1 - pnorm(zscore)
  if (!one_sided) (1 - pnorm(zscore))*2
}
ir <- read.csv('FILENAME', stringsAsFactors = FALSE) %>% rename(drug=abx_class) %>% select(pathogen,drug,include)
outpath <- 'FILEPATH'
inpath <- 'FILEPATH'

# Load data
d2 <- read.csv('FILENAME', stringsAsFactors = FALSE) 
d2<-d2%>%left_join(ir, by.x = c('pathogen', 'drug'))
d2$include[is.na(d2$include)] <- 1
d2 <- d2[!d2$include==0,]

#Quantities to be used in the model
d2$lnrr<-d2$yi
d2$lnrr_se<-d2$sei
d2$lnrr_var<-d2$vi

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

#2) A loop sets up model function and saves estimates and draws per antibiotic class
table_priors <- c()
for (s in 0:1) {
  for (j in antibiotic_group) {
    data <- d2[d2$drug == paste0(j) & d2$sterile == paste0(s),]
    orgs<-paste(unique(data$pathogen),sep = ',')
    for (i in orgs) {
      data$new <- ifelse(data$pathogen == i,1,0)
      colnames(data)[colnames(data) == 'new'] <- paste0(i)
    }
    
    overall_mean <- rma(yi,vi,data=data)$beta
    eff <- c('prior_beta_gaussian', overall_mean, 0.1)
    report <- c(s,j,overall_mean,0.1) 
    table_priors <- rbind(table_priors,report) 
    rm(report)
    
    priors <- c()
    for(i in orgs) {
      effs <- c(i,eff)
      priors <- rbind(priors, effs)
    }

    # Concatenate the definition of priors to be entered into the model
    definition_model <- paste0("LinearCovModel(\"",priors[,1], "\",",priors[,2]," = array(c(", as.double(priors[,3]), ",",as.double(priors[,4]),")))")
    model_definition <- paste(definition_model, collapse = ",\n")
    
  # Initiate MRBRT data
    dat1 <- MRData()
    dat1$load_df(
      data = data,  col_obs = "yi", col_obs_se = "sei",
      col_covs = as.list(orgs),
      col_study_id = "pathogen" )
    
  # Text to parse into model
  mrbrtschema <- paste0(
      "MRBRT(
  data = dat1,
  cov_models = list(LinearCovModel('intercept', use_re = TRUE, prior_beta_uniform = array(c(0,0))),",
      model_definition,
      "))"
    )
    
  # initiate MRBRT model
    mod1 <- eval(parse(text = paste(mrbrtschema)))
    mod1$fit_model(inner_print_level = 5L, inner_max_iter = 1000L)
    
    n_samples <- 1000L
    samples <- mod1$sample_soln(sample_size = n_samples)
    
  # Get draws from samples and save estimates
    draws <- as.data.frame(exp(samples[[1]]))
    names(draws) <- mod1$cov_model_names
    draws <- t(draws)
    write.csv(draws,'FILENAME', row.names = TRUE)
    rm(draws)
    
  # Summary table to export
    table <- c()
    vector <-c()
    for (x in 1:(length(unique(data$pathogen))+1)) {
      beta <- mod1$beta_soln[x] 
      beta_sd <- sd(samples[[1]][,x]) 
      pval <- get_pval(beta,beta_sd)
      vector <- c(mod1$cov_model_names[x],exp(beta),ifelse(pval<0.05,"*",""),beta,beta_sd,c(beta-(1.96*beta_sd)),c(beta+(1.96*beta_sd)),pval)
      table <- rbind(table,vector)
      vector <-c()
    }
    table <- as.data.frame(table)
    names(table) <- c('pathogen','mean_RR','pval05','mean_ln(RR)','SE','lb','ub','pval')
    summary <- data %>% group_by(pathogen) %>% summarize (data_points = n(), studies = n_distinct(unique(refvar)), sample_size = sum(total)) 
    table <- left_join(summary,table,by = 'pathogen')
    table$abx_class<-paste0(j)
    
    write.csv(table,'FILENAME', row.names = FALSE)
    rm(table,samples,mod1)
  }
}


# appending draws and summary table
draws_sterile <- c() 
for (j in antibiotic_group) {
  test<-read.csv('FILENAME', stringsAsFactors = FALSE)
  test$abx_class <- paste0(j)
  draws_sterile <- rbind(draws_sterile,test)
}
names(draws_sterile)[names(draws_sterile) == 'X'] <- 'pathogen'
draws_sterile <- draws_sterile[!draws_sterile$pathogen == 'intercept',]
write.csv(draws_sterile,'FILENAME')

summary_sterile <- c()
test<-c()
for (j in antibiotic_group) {
  test<-read.csv('FILENAME', stringsAsFactors = FALSE)
  test$abx_class <- paste0(j)
  summary_sterile <- rbind(summary_sterile,test)
}

write.csv(summary_sterile,'FILENAME')