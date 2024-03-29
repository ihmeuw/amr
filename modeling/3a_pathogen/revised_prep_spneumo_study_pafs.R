## Clean data, calculate study-level Strep pneumo PAF in the absence of vaccine ##
########################################################################################################
## Prep ##
########################################################################################################
pacman::p_load(data.table, openxlsx, reshape2, boot, ggplot2, plyr, matrixStats, msm, metafor)
invisible(sapply(list.files("FILEPATH", full.names = T), source))
library(crosswalk, lib.loc = "FILEPATH")
source("FILEPATH")
source("FILEPATH")
source("FILEPATH")
out_dir <- "FILEPATH"
model_name <- "GBD19_stepwise_revised"
out_dir <- paste0(out_dir, model_name, "/")

vtype <- "PCV7"

locs <- read.csv("FILEPATH")

age_info <- read.csv("FILEPATH")

pcv_cov_pull <- data.frame(get_covariate_estimates(covariate_id=210, location_id="all", year_id=1990:2019, decomp_step="step4", gbd_round_id = 6))
setnames(pcv_cov_pull, c("mean_value","lower_value","upper_value"), c("pcv_cov","pcv_cov_lower","pcv_cov_upper"))
pcv_cov <- merge(pcv_cov_pull, locs[,c("location_id","super_region_name","region_name")], by="location_id")

sero_cov <- read.csv("FILEPATH")

pcv_cov$region <- with(pcv_cov, ifelse(region_name=="Oceania","Oceania",
                           ifelse(region_name=="High-income North America","North America",
                           ifelse(region_name=="Western Europe","Europe",
                           ifelse(region_name=="Central Europe", "Europe",
                           ifelse(region_name=="Eastern Europe","Europe",
                           ifelse(region_name=="Australasia","North America",
                           ifelse(region_name=="High-income Asia Pacific","Asia",
                           ifelse(super_region_name=="Sub-Saharan Africa","Africa",
                           ifelse(super_region_name=="North Africa and Middle East","Africa",
                           ifelse(super_region_name=="Southeast Asia, East Asia, and Oceania","Asia",
                           ifelse(super_region_name=="South Asia","Asia",
                           ifelse(super_region_name=="Central Europe, Eastern Europe, and Central Asia","Asia","LAC"))
                           )))))))))))

pcv_cov <- join(pcv_cov, sero_cov, by=c("region"))

pcv_cov$pcv_vt_cov <- pcv_cov$pcv_cov * pcv_cov$covmean

pcv_cov$pcv_cov_std <- (pcv_cov$pcv_cov_upper - pcv_cov$pcv_cov_lower) / 2 / qnorm(0.975)
pcv_cov$vtcov_std <- (pcv_cov$covupper - pcv_cov$covlower) / 2 / qnorm(0.975)
pcv_cov$pcv_vt_cov_std <- with(pcv_cov, sqrt(vtcov_std^2 * pcv_cov^2 + pcv_cov_std^2 * covmean^2 + vtcov_std^2 * pcv_cov_std^2))

pcv_data <- read.xlsx("FILEPATH")

w <- grep("ve_", colnames(pcv_data))
pcv_data[w] <- lapply( pcv_data[w], function(x) as.numeric(as.character(x)) )
pcv_data[w] <- lapply( pcv_data[w], function(x) x/100 )

pcv_data$age_cat <- ifelse(pcv_data$age_start <5, "Under-5", "Over-5")

pcv_data$ve_pp_vt_std <- (pcv_data$ve_pp_vt_upper / 100 - pcv_data$ve_pp_vt_lower / 100) / qnorm(0.975) / 2
pcv_data$ve_all_pneumo_std <- (pcv_data$ve_all_pneumo_upper / 100 - pcv_data$ve_all_pneumo_lower / 100) / qnorm(0.975) / 2
pcv_data <- as.data.table(pcv_data)

rma_ve_pp_vt <- rma(yi = ve_pp_vt, sei=ve_pp_vt_std, data=pcv_data[!is.na(ve_pp_vt) & age_start < 5 & Study.type == "RCT",], 
                    slab=pcv_data[!is.na(ve_pp_vt) & age_start < 5 & Study.type == "RCT",]$Author)
sink("FILEPATH")
print(summary(rma_ve_pp_vt))
sink()
pdf("FILEPATH")
p <- forest(rma_ve_pp_vt, cex = 1.5)
print(p)
dev.off()
ve_pp_vt_pool <- as.numeric(rma_ve_pp_vt$b)
ve_pp_vt_pool_se <- as.numeric(rma_ve_pp_vt$se)
ve_pp_vt_pool_lower <- as.numeric(rma_ve_pp_vt$ci.lb)
ve_pp_vt_pool_upper <- as.numeric(rma_ve_pp_vt$ci.ub)
ve_pp_vt_df <- data.frame(mean = ve_pp_vt_pool, se = ve_pp_vt_pool_se)
write.csv(ve_pp_vt_df, "FILEPATH", row.names = F)

dt_tmp <- copy(pcv_data[!is.na(ve_pp_vt) & age_start > 5 & Study.type == "RCT",])
if (nrow(dt_tmp)>1) stop("Meta-analysis required for adults as well")
ve_pp_vt_adult <- as.numeric(dt_tmp$ve_pp_vt)
ve_pp_vt_adult_se <- as.numeric(dt_tmp$ve_pp_vt_std)
ve_pp_vt_adult_lower <- as.numeric(dt_tmp$ve_pp_vt_lower)
ve_pp_vt_adult_upper <- as.numeric(dt_tmp$ve_pp_vt_upper)

pcv_data <- subset(pcv_data, !is.na(ve_all_pneumo))

setnames(pcv_data, "vaccine_type", "vtype")
pcv_data$vtype <- ifelse(pcv_data$vtype=="PCV23","PCV13", pcv_data$vtype)
pcv_data$vtype <- ifelse(pcv_data$vtype=="PCV7/PCV13","PCV13", pcv_data$vtype)

pcv_data$year_id <- ifelse(is.na(pcv_data$year_cov), pcv_data$year_end, pcv_data$year_cov)
loc_meta <- get_location_metadata(location_set_id = 35, gbd_round_id = 6, decomp_step = "step4")
pcv_data <- merge(pcv_data, loc_meta[,.(location_id, ihme_loc_id)], by="ihme_loc_id")
pcv_data <- join(pcv_data, pcv_cov[,c("location_id","year_id","region","pcv_cov","pcv_cov_lower","pcv_cov_upper","pcv_vt_cov","pcv_vt_cov_std","vtype", "covmean", "covlower", "covupper")], by=c("location_id","year_id","vtype"))

pcv_data$raw_coverage_lower <- ifelse(is.na(pcv_data$reported_cov_lower), pcv_data$pcv_cov_lower, pcv_data$reported_cov_lower)
pcv_data$raw_coverage_upper <- ifelse(is.na(pcv_data$reported_cov_lower), pcv_data$pcv_cov_upper, pcv_data$reported_cov_upper)
pcv_data$pcv_vt_cov <- pcv_data$raw_coverage * pcv_data$covmean

pcv_data$raw_coverage_std <- (pcv_data$raw_coverage_upper - pcv_data$raw_coverage_lower) / 2 / qnorm(0.975)
pcv_data$vtcov_std <- (pcv_data$covupper - pcv_data$covlower) / 2 / qnorm(0.975)
pcv_data$pcv_vt_cov_std <- with(pcv_data, sqrt(vtcov_std^2 * raw_coverage^2 + raw_coverage_std^2 * covmean^2 + vtcov_std^2 * raw_coverage_std^2))

pcv_data$covstd <- (pcv_data$covupper / 100 - pcv_data$covlower / 100) / qnorm(0.975) / 2
pcv_data$study_coverage <- ifelse(pcv_data$Study.type == "Before-after", pcv_data$pcv_vt_cov, pcv_data$covmean)
pcv_data$study_coverage_std <- ifelse(pcv_data$Study.type == "Before-after", pcv_data$pcv_vt_cov_std, pcv_data$covstd)

pcv_data$ve_pp_vt_pooled <- ifelse(pcv_data$age_start < 5, ve_pp_vt_pool, ve_pp_vt_adult)
pcv_data$ve_pp_vt_pooled_se <- ifelse(pcv_data$age_start < 5, ve_pp_vt_pool_se, ve_pp_vt_adult_se)
pcv_data$PAF <- pcv_data$ve_all_pneumo/(pcv_data$ve_pp_vt_pooled*pcv_data$study_coverage)
pcv_paf <- as.data.table(pcv_data)[,.(Author, nid, Study.type, ihme_loc_id, year_start, year_end, age_start, age_end, endpoint_group, ve_all_pneumo, ve_pp_vt_pooled, study_coverage, PAF)]
pcv_vax_coverage <- as.data.table(pcv_data)[,.(Author, nid, Study.type, ihme_loc_id, year_start, year_end, age_start, age_end, endpoint_group, pcv_cov, reported_cov, raw_coverage, covmean, study_coverage)]
write.xlsx(pcv_paf, "FILEPATH")
write.xlsx(pcv_vax_coverage, "FILEPATH")

pcv_data <- pcv_data[!is.na(pcv_data$ve_all_pneumo_std)]

logit_cov_df <- as.data.frame(delta_transform(pcv_data$study_coverage, pcv_data$study_coverage_std, transformation = "linear_to_logit"))
pcv_data$logit_cov <- logit_cov_df$mean_logit
pcv_data$logit_cov_se <- logit_cov_df$sd_logit

log_vevt_df <- as.data.frame(delta_transform(pcv_data$ve_pp_vt_pooled, pcv_data$ve_pp_vt_pooled_se, transformation = "linear_to_log"))
pcv_data$ln_ve_vt <- log_vevt_df$mean_log
pcv_data$ln_ve_vt_se <- log_vevt_df$sd_log

log_veall_df <- as.data.frame(delta_transform(pcv_data$ve_all_pneumo, pcv_data$ve_all_pneumo_std, transformation = "linear_to_log"))
pcv_data$ln_ve_all <- log_veall_df$mean_log
pcv_data$ln_ve_all_se <- log_veall_df$sd_log

write.csv(pcv_data, "FILEPATH", row.names=F)
pcv_data <- data.frame(pcv_data)


pcv_df <- data.frame(nid=pcv_data$nid)
pcv_df_log <- data.frame(nid=pcv_data$nid)
for(i in 1:1000){
  ve_all <- exp(rnorm(n=length(pcv_data$nid), mean=pcv_data$ln_ve_all, sd=pcv_data$ln_ve_all_se))
  ve_vt <- exp(rnorm(n=length(pcv_data$nid), mean=pcv_data$ln_ve_vt, sd=pcv_data$ln_ve_vt_se))
  cov <- inv.logit(rnorm(n=length(pcv_data$nid), mean=pcv_data$logit_cov, sd=pcv_data$logit_cov_se))
  
  paf <- ve_all/(ve_vt*cov)
  
  pcv_df[,paste0("draw_",i)] <- paf
  pcv_df_log[,paste0("draw_",i)] <- log(paf)
}

pcv_data$mean_paf <- rowMeans(pcv_df[,2:1001])
draws_dt <- as.data.table(pcv_df[,2:1001])
names(draws_dt) <- paste0("paf_draw_", 0:999)

pcv_data$std_paf <- apply(pcv_df[,2:1001], 1, sd)
pcv_data$paf_lower <- apply(pcv_df[,2:1001], 1, function(x) quantile(x, 0.025, na.rm=T))
pcv_data$paf_upper <- apply(pcv_df[,2:1001], 1, function(x) quantile(x, 0.975, na.rm=T))

pcv_data$mean_log_paf <- rowMeans(pcv_df_log[,2:1001])
pcv_data$std_log_paf <- apply(pcv_df_log[,2:1001], 1, sd)
rma_paf <- rma(mean_paf, sei=std_paf, data=pcv_data, method="REML", slab=Author)
forest(rma_paf)

pcv_paf2 <- as.data.table(pcv_data)[,.(Author, nid, Study.type, ihme_loc_id, year_start, year_end, age_start, age_end, endpoint_group, ve_all_pneumo, ve_pp_vt_pooled, study_coverage, PAF, mean_paf, paf_lower, paf_upper, std_paf)]
pcv_paf2 <- cbind(pcv_paf2, draws_dt)
write.xlsx(pcv_paf2, "FILEPATH")

############################################################################################
## Prep a file that can be run in MR-BRT ##
############################################################################################
pcv_data <- data.frame(pcv_data)
pcv_data$age <- floor((pcv_data$age_end + pcv_data$age_start) / 2)
pcv_data$before_after <-ifelse(pcv_data$Study.type=="RCT",0,1)
setnames(pcv_data, "Study.type", "study_type")
mrbrt <- pcv_data[,c("nid","age","meas_value","meas_stdev","mean_paf","std_paf","mean_log_paf","std_log_paf","before_after","age_start","age_end","age_cat","ihme_loc_id","study_type","vtype",
                 "ve_all_pneumo","ve_pp_vt","year_start","pcv_cov","pcv_vt_cov", "Author")]

##########################################################
## Run the age curve in MR-BRT ##

mrbrt$age <- (mrbrt$age_end + mrbrt$age_start) / 2
mrbrt <- as.data.table(mrbrt)
mrbrt <- mrbrt[mrbrt$age_cat=="Under-5" | mrbrt$Author == "Bonten" | mrbrt$Author == "Ochoa-Gondar"]

mrbrt_df <- as.data.frame(delta_transform(mrbrt$mean_paf, mrbrt$std_paf, transformation = "linear_to_log"))
mrbrt$log_value <- mrbrt_df$mean_log
mrbrt$log_meas_se <- mrbrt_df$sd_log
mrbrt$age_scaled <- mrbrt$age / 100
mrbrt <- mrbrt[mean_paf < 1]
mrbrt$weight <- 1/mrbrt$log_meas_se^2
write.csv(mrbrt, "FILEPATH", row.names=F)
library(mrbrt001, lib.loc = "FILEPATH")
dat3 <- MRData()
dat3$load_df(
  data = mrbrt, col_obs = "log_value", col_obs_se = "log_meas_se",
  col_covs = list("age_scaled", "before_after"),
  col_study_id = "nid")

mod5 <- MRBRT(
  data = dat3,
  cov_models = list(
    LinearCovModel("intercept", use_re = T), 
    LinearCovModel("before_after", use_re = F),
    LinearCovModel(
      alt_cov = "age_scaled",
      use_spline = TRUE,
      spline_knots = array(c(min(mrbrt$age_scaled), max(mrbrt$age_scaled))),
      spline_degree = 2L,
      spline_knots_type = 'frequency',
      spline_r_linear = TRUE,
      spline_l_linear = FALSE,
      prior_spline_monotonicity = 'decreasing'
    )),
  inlier_pct = 0.9
)



py_save_object(object = mod5, filename = "FILEPATH", pickle = "dill")

mod5$fit_model(inner_print_level = 5L, inner_max_iter = 500L, outer_max_iter = 100L)

df_pred <- data.frame(expand.grid(intercept=1, age_scaled = seq(0.01,1,.01), before_after = c(0,1)))

dat_pred1 <- MRData()

dat_pred1$load_df(
  data = df_pred, 
  col_covs = list("age_scaled", "before_after")
)

n_samples <- 1000L
samples2 <- mod5$sample_soln(sample_size = n_samples)

draws3 <- mod5$create_draws(
  data = dat_pred1,
  beta_samples = samples2[[1]],
  gamma_samples = samples2[[2]],
  random_study = FALSE )

df_pred$pred3 <- mod5$predict(data = dat_pred1, sort_by_data_id = TRUE)
df_pred$pred3_lo <- apply(draws3, 1, function(x) quantile(x, 0.025))
df_pred$pred3_hi <- apply(draws3, 1, function(x) quantile(x, 0.975))
df_pred$mean_lin <- exp(df_pred$pred3)
df_pred$mean_lo_lin <- exp(df_pred$pred3_lo)
df_pred$mean_hi_lin <- exp(df_pred$pred3_hi)
draws3 <- as.data.frame(draws3)
names(draws3) <- paste("draw", 0:(n_samples-1), sep = "_")
df_pred <- cbind(df_pred, draws3)

pdf("FILEPATH")
p <- ggplot(df_pred) + geom_line(aes(x=100*age_scaled, y=mean_lin, col=factor(before_after))) +
geom_ribbon(aes(x=100*age_scaled, ymin=mean_lo_lin, ymax= mean_hi_lin, fill=factor(before_after)), alpha=0.25) + ggtitle("Strep pneumo attributable fraction") + theme_bw() +
xlab("Age mid") + ylab("Mean attributable fraction") + geom_hline(yintercept=0) + guides(fill=F) +
geom_point(data=mrbrt, aes(x=age, y=mean_paf, col=factor(before_after),  size = weight), alpha=0.75) +
scale_color_manual("Study type", values = c("#30A9DE","#E53A40"), label=c("RCT","Before-after")) + scale_fill_manual("Study type", values = c("#30A9DE","#E53A40"), label=c("RCT","Before-after")) + 
scale_y_continuous(limits = c(-0.5,1.5))
print(p)
dev.off()

write.csv(df_pred, "FILEPATH", row.names=F)