##############
# Aim: To obtain days elapsed between specimen date and hospital discharge
# Adjusted for days between hospital admission and specimen 
# It first loads and standardizes primary data to produce 
# adjusted length of infection estimates by pathogen and abx_class
##########################
library(lmtest)
library(tidyverse)
library(sandwich)
library(lubridate)
library(dplyr)
library(stringr)
library(data.table)
library(sqldf)
source("FILENAME")
source("FILENAME")

# defining username, paths and maps 
username <- Sys.getenv("USER")
maps_path <- sprintf("FILEPATH", username)
output <- 'FILEPATH'
antibiotic_map <- read.csv(paste0(maps_path,'FILENAME'),stringsAsFactors = FALSE)
pathogen_map <- read.csv(paste0(maps_path,'FILENAME'),stringsAsFactors = FALSE)
specimen_map <- read.csv(paste0(maps_path,'FILENAME'),stringsAsFactors = FALSE)
maps_path <- 'FILEPATH'
level_amr_map <- read.csv(paste0(maps_path,'icd_map_l2.csv'),stringsAsFactors = FALSE)
level_amr <- level_amr_map %>% distinct(acause,level_amr)
level_amr_map <- level_amr_map %>% distinct(cause_code,level_amr)
cause_id <- get_ids("cause") %>% left_join(level_amr) %>% distinct()
ucause_map <- read.csv(paste0(maps_path,'ucause_and_infectious_syndrome_map_non_icd.csv'),stringsAsFactors = FALSE)
ucause_map <- ucause_map %>% left_join(cause_id) %>% distinct(cause, .keep_all = TRUE)
ucause_map$level_amr[is.na(ucause_map$level_amr)] <- ucause_map$acause[is.na(ucause_map$level_amr)]

locs <- subset(get_location_metadata(location_set_id=35, gbd_round_id=6),select = c(location_id, location_name,super_region_name,super_region_id,region_id,region_name), droplevels = TRUE)
loc<-get_ids("location", return_all_columns=TRUE)

organism <- c('acinetobacter_baumanii','citrobacter_spp','enterobacter_spp','enterococcus_faecalis','enterococcus_faecium','enterococcus_spp','escherichia_coli','group_b_strep','group_a_strep','haemophilus_influenzae','klebsiella_pneumoniae','moraxella_spp','morganella_spp','neisseria_meningitidis','non_typhoidal_salmonellae','proteus_spp','providencia_spp','pseudomonas_aeruginosa','pseudomonas_spp','neisseria_gonorrheae','shigella_spp','mycobacterium_tuberculosis','salmonella_paratyphi','salmonella_typhi','salmonella_typhi_paratyphi','serratia_spp','staphylococcus_aureus','streptococcus_pneumoniae','mycoplasma','listeria','legionella_spp')
antibiotic_group<-c('third_gen_ceph','carbapenem','fluoroquinolone','penicillin','aminopenicillin','beta_lactamase_inhibitor','anti_pseudomonal_penicillin','methicillin','vancomycin','fourth_gen_ceph','sulfa','aminoglycoside','macrolide','nalidixic_acid')

# DATA_a
#####
dates_map <- fread('FILENAME')
DATA_A <- fread("FILENAME", stringsAsFactors=F) %>% 
            dplyr::select(patient_id,specimen_id,collected_date_id,discharge_date_id) 
DATA_A <- unique(DATA_A)

DATA_A <- DATA_A %>% left_join(dates_map , by = c("collected_date_id" = "date_id"))
DATA_A$specdt <- as.Date(DATA_A$date, "%m/%d/%Y")
DATA_A$date <- NULL
DATA_A <- DATA_A %>% left_join(dates_map , by = c("discharge_date_id" = "date_id"))
DATA_A$disdt <- as.Date(DATA_A$date, "%m/%d/%Y")
DATA_A$date <- NULL

DATA_A <- DATA_A %>% left_join(dates_map , by = c("admit_date_id" = "date_id"))
DATA_A$admdt <- as.Date(DATA_A$date, "%m/%d/%Y")
DATA_A$date <- NULL

DATA_A$los <- as.numeric(difftime(as.Date(DATA_A$disdt),as.Date(DATA_A$specdt), units = "days"))
DATA_A$prior_days <- as.numeric(difftime(as.Date(DATA_A$specdt),as.Date(DATA_A$admdt), units = "days"))

DATA_A$sp <- paste0('0',DATA_A$specimen_id)
DATA_A$sample_id <- paste0(DATA_A$patient_id, ".", DATA_A$sp)
DATA_A <- DATA_A[!is.na(DATA_A$patient_id),c('sample_id','los','prior_days','patient_id')]

#Loads sterile specimen data
DATA_A <- fread("FILENAME") 
DATA_A <- DATA_A %>% dplyr::select(year_id,sample_id,deaths, cases,raw_specimen,raw_pathogen,raw_antibiotic, resistance, hosp, cause,age_group_id) 


DATA_A <- DATA_A %>% left_join(DATA_A)
rm(dates_map)
year_id <- round(mean(unique(DATA_A$year_id)))
DATA_A$year_id <- year_id
DATA_A <- DATA_A %>% left_join(antibiotic_map)
DATA_A <- DATA_A %>% left_join(pathogen_map) 
DATA_A <- DATA_A %>% left_join(specimen_map)
DATA_A <- DATA_A %>% left_join(level_amr_map, by =  c('cause' = 'cause_code'))

DATA_A$resistant[DATA_A$resistance == 'resistant'] <- 1
DATA_A$resistant[DATA_A$resistance == 'sensitive'] <- 0
DATA_A <- DATA_A %>% filter(!is.na(los) & !is.na(abx_class) & !is.na(pathogen) & !is.na(resistant) & sterile == 1) 

# Groups and pivots data
DATA_A <- DATA_A %>%  group_by(sample_id,year_id,cause,sterile,pathogen,abx_class,age_group_id,level_amr,hosp) %>% 
  summarize(resistance = max(resistant), days = max(los), days_prior = max(prior_days))

# Obtain adjusted estimates in a loop for each combination
DATA_A$combo <- paste0(DATA_A$pathogen,"-",DATA_A$abx_class)

adjusted_coeff <- c()
all_coeff <- c()
for (s in c(1:1)) {
  for (i in organism) {
    for (j in antibiotic_group) {
      text_for_data <- paste0("DATA_A[DATA_A$pathogen == '",i,"' & DATA_A$abx_class == '",j,"' & DATA_A$sterile == ",s," & !is.na(DATA_A$resistance) & !is.na(DATA_A$days),]")
      data <- eval(parse(text = paste(text_for_data)))
      data <- data[!is.na(data$abx_class),]
      total <- data %>% dplyr::group_by(abx_class, resistance) %>% summarize(admit = n(), los = mean(days,na.rm = TRUE), sd = sd(days,na.rm = TRUE))
      if(!is.na(total[1,3]) & !is.na(total[2,3]) & !is.na(total[2,5]) & !is.na(total[1,5])) {
        model1 <- glm(data = data, formula = days ~ resistance, family = 'poisson'(link='log'))
        model2 <- glm(data = data, formula = days ~ resistance + days_prior , family = 'poisson'(link='log')) 
        if (data$combo %in% c('legionella-nalidixic_acid','streptococcus_pneumoniae-aminoglycoside','staphylococcus_aureus-fourth_gen_ceph','salmonella_typhi-nalidixic_acid','salmonella_typhi-aminoglycoside','salmonella_typhi-methicillin','salmonella_typhi-anti_pseudomonal_penicillin','salmonella_typhi-beta_lactamase_inhibitor','salmonella_paratyphi-aminopenicillin','salmonella_paratyphi-fluoroquinolone','shigella_spp-aminoglycoside','shigella_spp-beta_lactamase_inhibitor','proteus_spp-macrolide','klebsiella_pneumoniae-nalidixic_acid','klebsiella_pneumoniae-macrolide','klebsiella_pneumoniae-vancomycin','klebsiella_pneumoniae-penicillin','escherichia_coli-vancomycin','escherichia_coli-macrolide','escherichia_coli-penicillin','enterococcus_faecium-third_gen_ceph','enterococcus_faecium-carbapenem','enterococcus_faecium-anti_pseudomonal_penicillin','enterococcus_faecium-methicillin','enterococcus_faecium-fourth_gen_ceph')) {
          model3 <- model2
        }
        else {
          model3 <- glm(data = data, formula = days ~ resistance + days_prior + as.factor(age_group_id) + as.factor(level_amr), family = 'poisson'(link='log'))
        }
        r1 <- c('DATA_A', i, j, s, LOCATION_ID,'LOCATION_NAME','model1',summary.glm(model1)$coefficients[2,1:2],'model2',summary.glm(model2)$coefficients[2,1:2],'model3',summary.glm(model3)$coefficients[2,1:2],total[1,3],total[1,4],total[1,5],total[2,3],total[2,4],total[2,5])
        adjusted_coeff<- rbind(adjusted_coeff,r1)
        r11 <- c('DATA_A', i, j, s, LOCATION_ID,'LOCATION_NAME','model2',summary.glm(model2)$coefficients[2:3,1:2],'model3',summary.glm(model3)$coefficients[2:3,1:2],total[1,3],total[1,4],total[1,5],total[2,3],total[2,4],total[2,5])
        all_coeff <- rbind(all_coeff,r11)
      }
    }
  }
}
print("writing estimates")
colnames(adjusted_coeff) <- c('refvar','pathogen','abx_class','sterile','location_id','location','measure1','coeff_unadjusted','se_unadjusted','measure2','coeff_adjusted_days','se_adjusted_days','measure3','coeff_adjusted_cause','se_adjusted_cause','admit_s','mean_LOS_s','sd_s','admit_r','mean_LOS_r','sd_r')
adjusted_coeff$combo <- paste0(adjusted_coeff$pathogen,"-",adjusted_coeff$abx_class)
adjusted_coeff$measure3[adjusted_coeff$combo %in% c('legionella-nalidixic_acid','streptococcus_pneumoniae-aminoglycoside','staphylococcus_aureus-fourth_gen_ceph','salmonella_typhi-nalidixic_acid','salmonella_typhi-aminoglycoside','salmonella_typhi-methicillin','salmonella_typhi-anti_pseudomonal_penicillin','salmonella_typhi-beta_lactamase_inhibitor','salmonella_paratyphi-aminopenicillin','salmonella_paratyphi-fluoroquinolone','shigella_spp-aminoglycoside','shigella_spp-beta_lactamase_inhibitor','proteus_spp-macrolide','klebsiella_pneumoniae-nalidixic_acid','klebsiella_pneumoniae-macrolide','klebsiella_pneumoniae-vancomycin','klebsiella_pneumoniae-penicillin','escherichia_coli-vancomycin','escherichia_coli-macrolide','escherichia_coli-penicillin','enterococcus_faecium-third_gen_ceph','enterococcus_faecium-carbapenem','enterococcus_faecium-anti_pseudomonal_penicillin','enterococcus_faecium-methicillin','enterococcus_faecium-fourth_gen_ceph')] <- 'unadjusted'
write.csv(adjusted_coeff, paste0(output,'FILENAME'), row.names = F)


# DATA_B
# Reads laboratory data 
DATA_B <- read.csv('FILENAME', stringsAsFactors = FALSE) %>% 
  mutate(date = as.Date(Collection_date)) %>% select(ID1,ID2,date)
DATA_B<-as.data.table(unique(DATA_B))

# Reads hospital data
data_b <-read.csv('FILENAME', stringsAsFactors = FALSE) %>% 
  filter(ID3 == 1 & ID2 == 1 & !is.na(DIAGNOSIS)) %>% 
  select(DIAGNOSIS,admdt,disddt,ID) 
data_b$admdt <- as.Date(data_b$disdate)
data_b$disdt <- as.Date(data_b$disdate)
data_b$date<-data_b$admdt
data_b<-as.data.table(data_b)

#link dates to obtain difference
Data_B <- sqldf("select t.*, q.admdt, q.disdt, min(abs(t.date - q.admdt)), min(q.disdt - t.date) from DATA_B t left join data_b q on t.ID = q.ID and (t.date - q.admdt) > -2 and (q.disdt - t.date) > -1 group by t.rowid")
Data_B <- Data_B[!is.na(Data_B$admdt),]
Data_B$specdt <- Data_B$date
Data_B$los <- interval(Data_B$specdt,Data_B$disdt)/duration(num = 1, units = 'days') 
Data_B$prior_days <- as.numeric(difftime(as.Date(Data_B$specdt),as.Date(Data_B$admdt), units = "days"))
Data_B$ID <- paste(Data_B$ID1,Data_B$ID2,sep = "-")

#merge with standardized data 
DataB <- read.csv('FILENAME', stringsAsFactors = FALSE) 
DataB <- DataB %>% left_join (Data_B, by = "ID")

DataB$resistant[DataB$resistance == 'resistant'] <- 1
DataB$resistant[DataB$resistance == 'sensitive'|DataB$resistance == 'susceptible'] <- 0

DataB <- DataB %>% left_join(level_amr_map[,c('cause_code','level_amr')], by=c("cause" = "cause_code"))
DataB$cause[is.na(DataB$level_amr)] <- substring(DataB$cause[is.na(DataB$level_amr)],1,3)
DataB <- DataB %>% left_join(level_amr_map[,c('cause_code','level_amr')], by=c("cause" = "cause_code")) 
DataB$level_amr.x <- NULL
colnames(DataB)[colnames(DataB) == 'level_amr.y'] <- 'level_amr' 
DataB <- DataB %>% left_join(pathogen_map) %>% left_join(antibiotic_map) %>% left_join(specimen_map)

DataB$year_id <- round(mean(unique(DataB$year_id)))
DataB <- DataB %>% filter(!is.na(los) & !is.na(abx_class) & !is.na(pathogen) & !is.na(resistant)) 

# Group data for analysis
DataB <- DataB %>%  dplyr::group_by(sample_id,year_id,location_id,sterile,pathogen,abx_class,age_group_id,level_amr,hosp) %>% 
  summarize(resistance = max(resistant),days = max(los), days_prior=max(prior_days))

# define combinations present in data
antibiotic_group_DataB <-antibiotic_group[antibiotic_group %in% unique(DataB$abx_class)]
organism_DataB <- organism[organism %in% unique(DataB$pathogen)]
DataB$combo <- paste0(DataB$pathogen,"-",DataB$abx_class)

# A Loop obtains adjusted coefficients
adjusted_coeff <- c()
all_coeff <- c()
for (s in c(0:1)) {
  for (i in organism_DataB) {
    for (j in antibiotic_group_DataB) {
      text_for_data <- paste0("DataB[DataB$pathogen == '",i,"' & DataB$abx_class == '",j,"' & DataB$sterile == ",s," & !is.na(DataB$resistance) & !is.na(DataB$days),]")
      data <- eval(parse(text = paste(text_for_data)))
      data <- data[!is.na(data$abx_class),]
      text_for_table<- paste0("data %>% dplyr::group_by(abx_class, resistance) %>% summarize(admit = n(), los = mean(days,na.rm = TRUE), sd = sd(days,na.rm = TRUE))")
      total <- eval(parse(text = paste(text_for_table)))
      if(!is.na(total[1,3]) & !is.na(total[2,3]) & !is.na(total[2,5]) & !is.na(total[1,5])) {
        model1 <- glm(data = data, formula = days ~ resistance, family = 'poisson'(link='log'))
        model2 <- glm(data = data, formula = days ~ resistance + days_prior , family = 'poisson'(link='log')) 
        model2 <- coeftest(model2, vcov = vcovHC(model2, type = "HC1"))
        model3 <- glm(data = data, formula = days ~ resistance + days_prior + as.factor(level_amr) + as.factor(age_group_id), family = 'poisson'(link='log'))
        model3 <- coeftest(model3, vcov = vcovHC(model3, type = "HC1"))
      r2 <- c('DataB', i, j, s, LOCATION_ID,'LOCATION_NAME','model1',summary.glm(model1)$coefficients[2,1:2],'model2',model2[2,1:2],'model3',model3[2,1:2],total[1,3],total[1,4],total[1,5],total[2,3],total[2,4],total[2,5])
      adjusted_coeff<- rbind(adjusted_coeff,r2)
      r22 <- c('DataB', i, j, s, LOCATION_ID,'LOCATION_NAME','model2',model2[2:3,1:2],'model3',model3[2:4,1:2],total[1,3],total[1,4],total[1,5],total[2,3],total[2,4],total[2,5])
      all_coeff <- rbind(all_coeff,r22)
      }
    }
  }
}
colnames(adjusted_coeff) <- c('refvar','pathogen','abx_class','sterile','location_id','location','measure1','coeff_unadjusted','se_unadjusted','measure2','coeff_adjusted_days','se_adjusted_days','measure3','coeff_adjusted_cause','se_adjusted_cause','admit_s','mean_LOS_s','sd_s','admit_r','mean_LOS_r','sd_r')
adjusted_coeff <- as.data.frame(adjusted_coeff)
write.csv(adjusted_coeff, paste0(output,'FILENAME'), row.names = F)

# DATA_C
####
# load hospital data
Data_c1 <-read.csv('FILENAME', stringsAsFactors = FALSE) 
Data_c2 <-read.csv('FILENAME', stringsAsFactors = FALSE)  
Data_c3 <-read.csv('FILENAME', stringsAsFactors = FALSE)  

Data_c_dates <- rbind(Data_c1, Data_c2, Data_c3)
rm(Data_c1,Data_c2,Data_c3)
Data_c_dates <- Data_c_dates %>% select(ID, TYPE, specdt,admdt,disdt) 
Data_c_dates$admdt <- as.Date(Data_c_dates$admdt, "%d/%m/%Y")
Data_c_dates$disdt <- as.Date(Data_c_dates$disdt, "%d/%m/%Y")
Data_c_dates$specdt <- as.Date(Data_c_dates$specdt, "%d/%m/%Y")
Data_c_dates$los <- interval(Data_c_dates$specdt,Data_c_dates$disdt)/duration(num = 1, units = 'days') 
Data_c_dates$prior_days <- as.numeric(difftime(as.Date(Data_c_dates$specdt),as.Date(Data_c_dates$admdt), units = "days"))
Data_c_dates$sample_id <- paste0('Tha-', Data_c_dates$ID, "-", str_to_lower(Data_c_dates$TYPE), "-", Data_c_dates$specdt)
Data_c_dates <- Data_c_dates %>% distinct()

# merge with standardized data to map for sterile specimen
DATA_c <- read.csv('FILENAME', stringsAsFactors = FALSE) 
DATA_c$resistant[DATA_c$resistance == 'resistant'] <- 1
DATA_c$resistant[DATA_c$resistance == 'susceptible'] <- 0

# merge to maps
DATA_c <- DATA_c %>% left_join(pathogen_map) %>% left_join(antibiotic_map) %>% left_join(specimen_map)
DATA_c$year_id <- round(mean(unique(DATA_c$year_id)))
DATA_c <- DATA_c %>% left_join(level_amr_map, by = c('cause' = 'cause_code'))
DATA_c$cause[is.na(DATA_c$level_amr)] <- substring(DATA_c$cause[is.na(DATA_c$level_amr)],1,3)
DATA_c <- DATA_c %>% left_join(level_amr_map[,c('cause_code','level_amr')], by=c("cause" = "cause_code")) 
DATA_c$level_amr.x <- NULL
colnames(DATA_c)[colnames(DATA_c) == 'level_amr.y'] <- 'level_amr' 

# merge to dates
DATA_c <- DATA_c %>% left_join (Data_c_dates, by = "sample_id")

# discard nulls
DATA_c <- DATA_c %>% filter(!is.na(los) & !is.na(abx_class) & !is.na(pathogen) & !is.na(resistant)) 

# reshape for analysis
DATA_c <- DATA_c %>%  dplyr::group_by(sample_id,year_id,location_id,sterile,pathogen,abx_class,age_group_id,level_amr) %>% 
  summarize(resistance = max(resistant),days = max(los), days_prior=max(prior_days))

# Define combinations present in data
antibiotic_group_DATA_c <-antibiotic_group[antibiotic_group %in% unique(DATA_c$abx_class)]
organism_DATA_c <- organism[organism %in% unique(DATA_c$pathogen)]
organism_DATA_c <- unique(DATA_c$pathogen)[!unique(DATA_c$pathogen) %in% c("gram_negative_other","gram_negative_rod", "stenotrophomonas_maltophilia", "streptococcus_unsp", "contaminant", "burkholderia_mallei", "burkholderia_pseudomallei")]
DATA_c$combo <- paste0(DATA_c$pathogen,"-",DATA_c$abx_class)


#A Loop obtains adjusted coefficients
adjusted_coeff_DATA_c <- c()
all_coeff_DATA_c <- c()
for (s in c(1:1)) {
  for (i in organism_DATA_c) {
    for (j in antibiotic_group_DATA_c) {
      text_for_data <- paste0("DATA_c[DATA_c$pathogen == '",i,"' & DATA_c$abx_class == '",j,"' & DATA_c$sterile == ",s," & !is.na(DATA_c$resistance) & !is.na(DATA_c$days),]")
      data <- eval(parse(text = paste(text_for_data)))
      data <- data[!is.na(data$abx_class),]
      total <- data %>% dplyr::group_by(abx_class, resistance) %>% summarize(admit = n(), los = mean(days,na.rm = TRUE), sd = sd(days,na.rm = TRUE))
      if(!is.na(total[1,3]) & !is.na(total[2,3]) & !is.na(total[2,5]) & !is.na(total[1,5])) {
        model1 <- glm(data = data, formula = days ~ resistance, family = 'poisson'(link='log'))
        model1 <- coeftest(model1, vcov = vcovHC(model1, type = "HC1"))
        if (data$combo %in% c('acinetobacter_baumanii-beta_lactamase_inhibitor','haemophilus_influenzae-macrolide','aeromonas_spp-anti_pseudomonal_penicillin')) {
          model2 <- rbind(model1,model1)
          model3 <- model2
          } else {
          model2 <- glm(data = data, formula = days ~ resistance + days_prior , family = 'poisson'(link='log')) 
          model2 <- coeftest(model2, vcov = vcovHC(model2, type = "HC1"))
            if (data$combo %in% c('group_b_strep-sulfa','group_b_strep-macrolide')) {
              model3 <- rbind(model2,model2)
            } else {
            model3 <- glm(data = data, formula = days ~ resistance + days_prior + as.factor(level_amr) + as.factor(age_group_id), family = 'poisson'(link='log'))
            model3 <- coeftest(model3, vcov = vcovHC(model3, type = "HC1"))
            }
        }
        r3 <- c('DATA_c', i, j, s, LOCATION_ID,'LOCATION_NAME','unadjusted',model1[2,1:2],'adjusted',model2[2,1:2],'model3',model3[2,1:2],total[1,3],total[1,4],total[1,5],total[2,3],total[2,4],total[2,5])
        adjusted_coeff_DATA_c <- rbind(adjusted_coeff_DATA_c,r3)
        r33 <- c('DATA_c', i, j, s, LOCATION_ID,'LOCATION_NAME','model2',model2[2:3,1:2],'model3',model3[2:4,1:2],total[1,3],total[1,4],total[1,5],total[2,3],total[2,4],total[2,5])
        all_coeff_DATA_c <- rbind(all_coeff_DATA_c,r33)
      }
    }
  }
}
adjusted_coeff_DATA_c <- as.data.frame(adjusted_coeff_DATA_c)
colnames(adjusted_coeff_DATA_c) <- c('refvar','pathogen','abx_class','sterile','location_id','location','measure1','coeff_unadjusted','se_unadjusted','measure2','coeff_adjusted_days','se_adjusted_days','measure3','coeff_adjusted_cause','se_adjusted_cause','admit_s','mean_LOS_s','sd_s','admit_r','mean_LOS_r','sd_r')
adjusted_coeff_DATA_c$combo <- paste0(adjusted_coeff_DATA_c$pathogen,"-",adjusted_coeff_DATA_c$abx_class)
adjusted_coeff_DATA_c$measure2[adjusted_coeff_DATA_c$combo %in% c('acinetobacter_baumanii-beta_lactamase_inhibitor','haemophilus_influenzae-macrolide','aeromonas_spp-anti_pseudomonal_penicillin')] <- 'unadjusted'
adjusted_coeff_DATA_c$measure3[adjusted_coeff_DATA_c$combo %in% c('group_b_strep-sulfa','group_b_strep-macrolide','acinetobacter_baumanii-beta_lactamase_inhibitor','haemophilus_influenzae-macrolide','aeromonas_spp-anti_pseudomonal_penicillin')] <- 'unadjusted'
adjusted_coeff_DATA_c$combo <- NULL
write.csv(adjusted_coeff_DATA_c, paste0(output,'FILENAME'), row.names = F)
