## Pull in the results from MR-BRT (1000 draws) of the age-specific estimate of S. pneumoniae PAF
## in the absence of the vaccine. Take those results and calculate the age-year-location draws
## based on vaccine coverage for the final S pneumo PAF draws.
########################################################################################################
## Prep ##
########################################################################################################
library(metafor)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(boot)
library(openxlsx)
library(gridExtra)
library(crosswalk, lib.loc = "FILEPATH")
invisible(sapply(list.files("FILEPATH", full.names = T), source))

# Directories
in_dir <- "FILEPATH"
out_dir <- "FILEPATH"

# Age metadata
age_meta <- read.csv("FILEPATH")
age_meta <- subset(age_meta, age_pull==1)
age_meta$age_mid <- (age_meta$age_start + age_meta$age_end)/2
loc_meta <- get_location_metadata(location_set_id = 35, gbd_round_id = 7, decomp_step = "iterative")
locs <- fread("FILEPATH")
loc_loop <- locs$location_id

years <- c(1990, 1995, 2000, 2005, 2010, 2015, 2017, 2019)

pcv_data <- read.csv("FILEPATH")
pcv_data <- as.data.table(pcv_data)
pcv_data$age_mid <- pcv_data$age_scaled*100
pcv_data <- pcv_data[before_after == 0]
pcv_data[, paste0("draw_",0:999) := exp(.SD), .SDcols = paste0("draw_",0:999)] 

DT <- data.table(pcv_data, key = c("age_mid"))
tm <- data.table(age_meta, key = c("age_mid"))

paf_age <- DT[tm, roll='nearest']
setnames(paf_age, paste0("draw_",0:999), paste0("age_",1:1000))

pcv_cov_pull <- data.frame(get_covariate_estimates(covariate_id=210, location_id=loc_loop, year_id=1990:2019, decomp_step="step4", gbd_round_id = 6))
setnames(pcv_cov_pull, c("mean_value","lower_value","upper_value"), c("pcv_cov","pcv_cov_lower","pcv_cov_upper"))
pcv_cov <- merge(pcv_cov_pull, loc_meta[,c("location_id","super_region_name","region_name")], by="location_id")

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

pcv_cov <- subset(pcv_cov, vtype=="PCV13")
pcv_cov <- subset(pcv_cov, !is.na(covmean))

pcv_cov <- subset(pcv_cov, year_id %in% years)

pcv_cov_zeros <- subset(pcv_cov, pcv_vt_cov == 0)
pcv_cov <- subset(pcv_cov, pcv_vt_cov != 0)

logit_cov_df <- as.data.frame(delta_transform(pcv_cov$pcv_vt_cov, pcv_cov$pcv_vt_cov_std, transformation = "linear_to_logit"))
pcv_cov$logit_mean <- logit_cov_df$mean_logit
pcv_cov$logit_std <- logit_cov_df$sd_logit
pcv_cov <- rbind.fill(pcv_cov, pcv_cov_zeros)

for(i in 1:1000){
  pcv_cov[,paste0("vtcov_",i)] <- inv.logit(rnorm(length(pcv_cov$covariate_id), pcv_cov$logit_mean, pcv_cov$logit_std))
  pcv_cov[pcv_cov$pcv_vt_cov==0,paste0("vtcov_",i)] <- 0
}

ve_optimal <- read.csv("FILEPATH")
for(i in 1:1000){
  ve_optimal[,paste0("ve_optimal_",i)] <- rnorm(1, ve_optimal$mean, ve_optimal$se)
}

n <- 1
for(l in loc_loop){
  print(paste0("On location ", n, " of ",length(loc_loop)))
  df <- expand.grid(age_group_id = age_meta$age_group_id, year_id=years, location_id=l, sex_id=c(1,2))
  pcv_df <- expand.grid(age_group_id = age_meta$age_group_id, year_id=years, location_id=l, sex_id=c(1,2))
  pcv_df <- join(pcv_df, pcv_cov, by=c("location_id","year_id"))
  pcv_df <- join(pcv_df, paf_age, by=c("age_group_id"))
  # adjust down the VE_optimal based on the relative VE in adults as compared to children in child vaccine studies
  # This could be done using the Grijalva study where the VE for children 2-4 was ~74% and for adults 65+ VE was ~20%
  # so multiply by correction VE adult/VE child Where VE_adult/VE_child = .2/.74.
  pcv_df <- as.data.table(pcv_df)
  pcv_df[age_mid>5, paste0("vtcov_",1:1000) := .SD * (.2/.74), .SDcols = paste0("vtcov_",1:1000)] 
  pcv_df <- as.data.frame(pcv_df)
  
  for(i in 1:1000){
    base <- pcv_df[,paste0("age_",i)]
    vt <- pcv_df[,paste0("vtcov_",i)]
    v_effectiveness <- ve_optimal[,paste0("ve_optimal_",i)]
    paf <- (base * (1-vt*v_effectiveness)) / (1-base*vt*v_effectiveness)
    df[,paste0("draw_",i)] <- paf
  
    df[df$age_group_id < 4, paste0("draw_",i)] <- 0

  }
  
  df$cause_id <- 322
  df$rei_id <- 188
  df$modelable_entity_id <- 1263
  setnames(df, paste0("draw_",1:1000), paste0("draw_",0:999))
  
  fwrite(df, "FILEPATH")
  fwrite(df, "FILEPATH")
  
  n <- n + 1
}

for(l in loc_loop){
  for (measure in c("yll", "yld")){
    df <- fread("FILEPATH")
    draws <- df[,grep("draw", names(df)), with = F]
    if (any(draws > 1) + any(draws < -1) > 0) {
      print(paste("invalid values for", l, measure)) 
    }
  }
}

gbd19 <- get_outputs(topic = "rei",
                     rei_id = 188,
                     cause_id = 322,
                     measure_id = c(3,4),
                     metric_id = 2,
                     gbd_round_id = 6,
                     decomp_step = "step5",
                     sex_id = c(1,2),
                     age_group_id = "all",
                     year_id = 2019,
                     location_id = loc_loop)
locs <- read.csv("FILEPATH")
locs <- subset(locs, is_estimate==1 & most_detailed==1 & level < 4)
loc_loop <- locs$location_id
for (measure in c("yll", "yld")){
  pdf("FILEPATH", width = 8, height = 5)
  for(l in loc_loop){
    print(l)
    my_loc_name <- subset(locs, location_id == l)$location_name
    df <- fread("FILEPATH")
    if (measure == "yll") m_id <- 4 else m_id <- 3
    df <- df[year_id == 2019]
    df$mean_paf <- rowMeans(df[,paste0("draw_",0:999)])
    df_19 <- gbd19[location_id == l & measure_id == m_id]
    setnames(df_19, "val", "mean_paf")
    df_19$analysis <- "GBD 2019 Final"
    df$analysis <- "AMR Custom"
    merged <- rbind(df[,.(age_group_id, location_id, sex_id, year_id, mean_paf, analysis)], 
                    df_19[,.(age_group_id, location_id, sex_id, year_id, mean_paf, analysis)])
    merged <- merge(merged, data.table(age_group_id = age_meta$age_group_id, age_mid = age_meta$age_mid), by = "age_group_id")
    g <- ggplot(merged, aes(x = age_mid, y = mean_paf, color = analysis, shape = as.factor(sex_id)), alpha = 0.5) + geom_point() + xlab("Age") +
      ylab("PAF") + ggtitle(paste(my_loc_name, measure)) 
    print(g)
  }
  dev.off()
}

locs <- loc_meta[location_id %in% loc_loop]
measure <- "yll"
df_list <- list()
df_list <- lapply(loc_loop, function(l){
  df <- fread("FILEPATH")    
  df$mean_paf <- rowMeans(df[,paste0("draw_",0:999), with = F])
  return(df)
})
df_all <- rbindlist(df_list)
df_all <- merge(df_all, select(locs, c("location_id", "super_region_name")), by = "location_id")
df_all <- data.table(df_all)
df_all[, mean_paf_year:=mean(mean_paf), by = c("year_id", "super_region_name")]
df_all[, mean_paf_age:=mean(mean_paf), by = c("year_id", "super_region_name", "age_group_id")]
df_all <- merge(df_all, age_meta, by="age_group_id")
gbd_all <- get_outputs(topic = "rei",
                       rei_id = 188,
                       cause_id = 322,
                       measure_id = c(3,4),
                       metric_id = 2,
                       gbd_round_id = 6,
                       decomp_step = "step5",
                       sex_id = 3,
                       age_group_id = c(22,unique(age_meta$age_group_id)),
                       year_id = "all",
                       location_id = loc_meta[level == 1]$location_id)


pdf("FILEPATH", width = 12, height = 8)
g1 <- ggplot() + geom_point(data = df_all, aes(x = year_id, y = mean_paf_year, color = super_region_name)) +
  ggtitle("AMR Custom PAF, Averaged All Ages") + theme_bw() + theme(legend.position= "bottom") 
g1a <- ggplot() + geom_point(data = gbd_all[age_group_id == 22], aes(x = year_id, y = val, color = location_name)) +
  ggtitle("GBD 2019 PAF, All Ages") + theme_bw() + theme(legend.position= "bottom")
g2 <- ggplot() + geom_point(data = df_all, aes(x = year_id, y = mean_paf_age, color = super_region_name)) +
  ggtitle("AMR Custom PAF, By Age") + facet_wrap(~age_group_name)+ theme_bw() + theme(legend.position= "bottom")
g2a <- ggplot() + geom_point(data = gbd_all, aes(x = year_id, y = val, color = location_name)) +
  ggtitle("GBD 2019 PAF, By Age") + facet_wrap(~age_group_name)+ theme_bw() + theme(legend.position= "bottom")
grid.arrange(g1, g1a, nrow = 1); print(g2); print(g2a); dev.off()

gbd_all <- merge(gbd_all, age_meta, by = "age_group_id")
pdf(paste0("FILEPATH"), width = 12, height = 8)
g1 <- ggplot() + geom_point(data = df_all[year_id == 2019], aes(x = age_mid, y = mean_paf_age, color = super_region_name)) +
  ggtitle("AMR Custom PAF, Super Region Average") + theme_bw() + theme(legend.position= "bottom") + scale_y_continuous(limits=(c(0,0.8)))
g1a <- ggplot() + geom_point(data = gbd_all[year_id == 2019], aes(x = age_mid, y = val, color = location_name)) +
  ggtitle("GBD 2019 PAF, Super Region Average") + theme_bw() + theme(legend.position= "bottom") + scale_y_continuous(limits=(c(0,0.8)))
print(g1); print(g1a); dev.off()

  pdf("FILEPATH", width = 8, height = 5)
  for(l in loc_loop){
    print(l)
    my_loc_name <- subset(locs, location_id == l)$location_name
    df <- fread("FILEPATH")
    if (measure == "yll") m_id <- 4 else m_id <- 3
    df <- df[year_id == 2019]
    df$mean_paf <- rowMeans(df[,paste0("draw_",0:999)])
    df_19 <- gbd19[location_id == l & measure_id == m_id]
    setnames(df_19, "val", "mean_paf")
    df_19$analysis <- "GBD 2019 Final"
    df$analysis <- "AMR Custom"
    merged <- rbind(df[,.(age_group_id, location_id, sex_id, year_id, mean_paf, analysis)],
                    df_19[,.(age_group_id, location_id, sex_id, year_id, mean_paf, analysis)])
    merged <- merge(merged, data.table(age_group_id = age_meta$age_group_id, age_mid = age_meta$age_mid), by = "age_group_id")
    g <- ggplot(merged, aes(x = age_mid, y = mean_paf, color = analysis, shape = as.factor(sex_id)), alpha = 0.5) + geom_point() + xlab("Age") +
      ylab("PAF") + ggtitle(paste(my_loc_name, measure))
    print(g)
  }
  dev.off()