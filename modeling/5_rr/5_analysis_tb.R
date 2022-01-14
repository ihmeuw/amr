#####
# Mycobacterium tuberculosis has a unique set of antibiotics 
# a) XDR , b) mono-INH and c) mono-RIF
# and analysis for these mentioned is carried out in this script
######

## LOAD FUNCTIONS
library(ggplot2)
library(msm)
library(readxl)
library(dplyr)
library(data.table)
library(mrbrt001, lib.loc = "FILEPATH")
library(msm, lib.loc = "FILEPATH")
source("FILENAME")
source("FILENAME")
source("FILENAME"))
get_pval <- function(beta, beta_sd, one_sided = FALSE) {
  zscore <- abs(beta/beta_sd)
  if (one_sided) 1 - pnorm(zscore)
  if (!one_sided) (1 - pnorm(zscore))*2
}

#A) XDR compared to non-MDR and susceptible MTB

#1)read in data and compute mean RR and SE
xdr2<-read.csv("FILENAME", stringsAsFactors = F, na = "") %>% filter(drug=="XDR")
xdr2$refvar2<-paste0(xdr2$refvar,xdr2$location_id)
xdr2 <- xdr2[xdr2$died_r >0 & xdr2$died_s >0,]
xdr2$rr <- (xdr2$died_r/xdr2$admit_r)/(xdr2$died_s/xdr2$admit_s)
xdr2$lnrr <- log((xdr2$died_r/xdr2$admit_r)/(xdr2$died_s/xdr2$admit_s))
xdr2$lnrr_var <- (1/xdr2$died_r) - (1/xdr2$admit_r) + (1/xdr2$died_s) - (1/xdr2$admit_s)
xdr2$lnrr_se <- sqrt(xdr2$lnrr_var)
xdr2$lnrr_upper <- (xdr2$lnrr+1.96*xdr2$lnrr_se)
xdr2$lnrr_lower <- (xdr2$lnrr-1.96*xdr2$lnrr_se)
xdr2$rr_upper <- exp(xdr2$lnrr_upper)
xdr2$rr_lower <- exp(xdr2$lnrr_lower)

# GRAPHICS
fit1 <- run_mr_brt(output_dir  = 'FILEPATH',model_label = "xdr",data        = xdr2,mean_var    = "lnrr",se_var      = "lnrr_se",
  method      = "trim_maxL",study_id    = "refvar2",trim_pct    = 0.15,overwrite_previous = TRUE)
plot_mr_brt(fit1, continuous_vars = "intercept", dose_vars = "intercept")

ggplot(xdr2, aes(ymax = rr_upper, ymin = rr_lower)) + 
  geom_point(aes(y = rr, x = refvar)) + geom_errorbar(aes(x = refvar), width=0) +
  theme_bw() + labs(x = "", y = "Relative Risk") + coord_flip() + ylim(-1,5) #+
forest(rma(data=xdr2,yi=lnrr, sei=lnrr_se,slab=paste0(xdr2$refvar,"_",xdr2$location_name)),cex=.75, atransf = exp)

#2) Sets up MRBRT function and outputs estimates and draws

antibiotic_group<-c('xdr')
data <- xdr2[xdr2$drug == 'XDR',]
  #Definition 
          definition_model <- paste0("LinearCovModel('",as.list('intercept'), "')")
         model_definition <- paste(definition_model, collapse = ",\n")

  # MRBRT Data
  dat1 <- MRData()
  dat1$load_df(
    data = data,  col_obs = "lnrr", col_obs_se = "lnrr_se",
    col_study_id = "refvar2" )
  
  # Text to evaluate MRBRT model
  mrbrtschema <- paste0(
    "MRBRT(data = dat1,cov_models = list(LinearCovModel('intercept', use_re = TRUE)))"
  )
  
  # initialize MRBRT
  mod1 <- eval(parse(text = paste(mrbrtschema)))
  mod1$fit_model(inner_print_level = 5L, inner_max_iter = 1000L)
  
  # Create samples and export draws
  n_samples <- 1000L
  samples <- mod1$sample_soln(sample_size = n_samples)
  draws <- as.data.table(exp(samples[[1]]))
  draws <- t(draws)  
  draws <- as.data.frame(draws)
  colnames(draws) <- paste(lapply(0:999, function(x) paste0("draw_", x)))
  draws$pathogen <- 'mycobacterium_tuberculosis'
  draws$abx_class <- 'xdr'
  write.csv(draws,'FILENAME', row.names = TRUE)
  rm(draws)
  
    beta <- mod1$beta_soln[1] 
    beta_sd <- sd(samples[[1]][,1]) 
    pval <- get_pval(beta,beta_sd)
    vector <- c(mod1$cov_model_names[x],exp(beta),ifelse(pval<0.05,"*",""),beta,beta_sd,exp(d+(1.96*c)),c(beta+(1.96*beta_sd)),pval)
  table <- as.data.frame(vector)
  names(table) <- c('pathogen','mean_RR','pval05','mean_ln(RR)','SE','lb','ub','pval')
  summary <- data %>% group_by(pathogen) %>% summarize (data_points = n(), sample_size = sum(sample_size)) 
  table <- left_join(summary,table,by = 'pathogen')
  table$abx_class<-paste0('xdr')
  write.csv(table,"FILENAME", row.names = F)

  
#B) mono-INH 

#1) read in data and compute mean RR and SE
library(tidyverse)
inh<-read.csv("FILENAME", stringsAsFactors = F, na = "") %>% filter(grepl("isoniazid",drug))
inh <- inh[grepl("pan",inh$drug)|(inh$refvar=="REF1")|(inh$refvar=="REF2")|(inh$refvar=="REF3"),]
inh <- inh[inh$died_r >0 & inh$died_s >0,]
inh$rr <- (inh$died_r/inh$admit_r)/(inh$died_s/inh$admit_s)
inh$lnrr <- log((inh$died_r/inh$admit_r)/(inh$died_s/inh$admit_s))
inh$lnrr_var <- (1/inh$died_r) - (1/inh$admit_r) + (1/inh$died_s) - (1/inh$admit_s)
inh$lnrr_se <- sqrt(inh$lnrr_var)

# 2) Plots
ggplot(inh, aes(ymax = rr_upper, ymin = rr_lower)) + geom_point(aes(y = rr, x = refvar2)) + geom_errorbar(aes(x = refvar2), width=0) + theme_bw() + labs(x = "", y = "Relative Risk") + coord_flip() + ylim(-2,6) 
forest(rma(data=inh,yi=lnrr, sei=lnrr_se,slab=paste0(inh$refvar2,(inh$died_r+inh$died_s))),cex=.75, atransf = exp)
fit1 <- run_mr_brt(
  output_dir  = 'FILEPATH',model_label = "inh_pansusceptible", data  = inh,mean_var = "lnrr",  se_var = "lnrr_se", study_id = "refvar",
  method      = "trim_maxL",  trim_pct    = 0.15,overwrite_previous = TRUE)
plot_mr_brt(fit1, continuous_vars = "intercept", dose_vars = "intercept")


#3) Sets up MRBRT function and outputs
#Definition 
definition_model <- paste0("LinearCovModel('",as.list('intercept'), "')")
model_definition <- paste(definition_model, collapse = ",\n")

# We load the data onto MRBRT
dat1 <- MRData()
dat1$load_df(
  data = inh,  col_obs = "lnrr", col_obs_se = "lnrr_se",
  col_study_id = "refvar" )

# Text to parse into the MRBRT model
mrbrtschema <- paste0(
  "MRBRT(
data = dat1,
cov_models = list(LinearCovModel('intercept', use_re = TRUE)))"
)

# initialize MRBRT
mod1 <- eval(parse(text = paste(mrbrtschema)))
mod1$fit_model(inner_print_level = 5L, inner_max_iter = 1000L)

#Create samples and export draws
n_samples <- 1000L
samples <- mod1$sample_soln(sample_size = n_samples)
draws <- as.data.table(exp(samples[[1]]))
draws <- t(draws)
draws <- as.data.frame(draws)
colnames(draws) <- paste(lapply(0:999, function(x) paste0("draw_", x)))
draws$pathogen <- 'mycobacterium_tuberculosis'
draws$abx_class <- 'isoniazid'
write.csv(draws,'FILENAME', row.names = TRUE)
rm(draws)

beta <- mod1$beta_soln[1] 
beta_sd <- sd(samples[[1]][,1]) 
pval <- get_pval(beta,beta_sd)
vector <- c('mycobacterium_tuberculosis','inh',sum(inh$sample_size),exp(beta),ifelse(pval<0.05,"*",""),beta,beta_sd,exp(beta+(1.96*beta_sd)),c(beta+(1.96*beta_sd)),pval)
write.csv(vector,"FILENAME", row.names = F)
print(vector)


#C) mono-RIF

rif<-read.csv("FILENAME", stringsAsFactors = F, na = "") %>% filter(grepl("rifampicin",drug))
table(rif$refvar[grepl("pan",rif$drug)],rif$drug[grepl("pan",rif$drug)])
rif <- rif[grepl("pan",rif$drug)|(rif$refvar=="REF1")|(rif$refvar=="REF2"),]
rif <- rif[rif$died_r >0 & rif$died_s >0,]
rif$rr <- (rif$died_r/rif$admit_r)/(rif$died_s/rif$admit_s)
rif$lnrr <- log((rif$died_r/rif$admit_r)/(rif$died_s/rif$admit_s))
rif$lnrr_var <- (1/rif$died_r) - (1/rif$admit_r) + (1/rif$died_s) - (1/rif$admit_s)
rif$lnrr_se <- sqrt(rif$lnrr_var)

rif <-rif[!is.na(rif$lnrr_se),]
rif$lnrr <- log(rif$rr)
rif$refvar2<-paste0(rif$refvar,rif$location)

#2 plots
ggplot(rif, aes(ymax = rr_upper, ymin = rr_lower)) + geom_point(aes(y = rr, x = refvar2)) + geom_errorbar(aes(x = refvar2), width=0) + theme_bw() + labs(x = "", y = "Relative Risk") + coord_flip() + ylim(-2,6) 
forest(rma(data=rif,yi=lnrr, sei=lnrr_se,slab=paste0(rif$refvar,"_",rif$location,"_",(rif$died_r+rif$died_s))),cex=.75, atransf = exp)
fit1 <- run_mr_brt(  output_dir  = 'FILEPATH', model_label = "rif", data        = rif, mean_var    = "lnrr",  se_var      = "lnrr_se",  study_id    = "refvar",
  method      = "trim_maxL",  trim_pct    = 0.15, overwrite_previous = TRUE)
plot_mr_brt(fit1, continuous_vars = "intercept", dose_vars = "intercept")


#3) Sets up MRBRT model and outputs
definition_model <- paste0("LinearCovModel('",as.list('intercept'), "')")
model_definition <- paste(definition_model, collapse = ",\n")

# We load the data onto MRBRT
dat1 <- MRData()
dat1$load_df(
  data = rif,  col_obs = "lnrr", col_obs_se = "lnrr_se",
  col_study_id = "refvar" )

# Text to parse into the MRBRT model
mrbrtschema <- paste0(
  "MRBRT(
data = dat1,
cov_models = list(LinearCovModel('intercept', use_re = TRUE)))"
)

# initialize MRBRT
mod1 <- eval(parse(text = paste(mrbrtschema)))
mod1$fit_model(inner_print_level = 5L, inner_max_iter = 1000L)

#Create samples and export draws
n_samples <- 1000L
samples <- mod1$sample_soln(sample_size = n_samples)
draws <- as.data.table(exp(samples[[1]]))
draws <- t(draws)
draws <- as.data.frame(draws)
colnames(draws) <- paste(lapply(0:999, function(x) paste0("draw_", x)))
draws$pathogen <- 'mycobacterium_tuberculosis'
draws$abx_class <- 'rifampicin'
write.csv(draws,'FILENAME', row.names = TRUE)

beta <- mod1$beta_soln[1] 
beta_sd <- sd(samples[[1]][,1]) 
pval <- get_pval(beta,beta_sd)
vector <- c('mycobacterium_tuberculosis','rif',sum(rif$sample_size),exp(beta),ifelse(pval<0.05,"*",""),beta,beta_sd,exp(beta+(1.96*beta_sd)),c(beta+(1.96*beta_sd)),pval)
write.csv(vector,"FILENAME", row.names = F)
print(vector)
