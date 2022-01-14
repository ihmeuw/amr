## Script to estimate relative risk of dying from resistant drug bug combinations ##
# from primary data (hospital or laboratory with discharge information)
#0) Defines antibiotic class and pathogens
#1) Defines paths, packages, functions, data
#2) Outputs relative risk adjusted for confounders (diagnosis, age, hospital onset) when possible  
##################################################
rm(list = ls())

###0) Define antibiotic class and pathogens

library(tidyverse)
library(lmtest)
library(sandwich)

organism <- c('acinetobacter_baumanii','citrobacter_spp','enterobacter_spp','enterococcus_faecalis','enterococcus_faecium','enterococcus_spp','escherichia_coli','group_a_strep','group_b_strep','haemophilus_influenzae','klebsiella_pneumoniae','moraxella_spp','morganella_spp','neisseria_meningitidis','non_typhoidal_salmonellae','proteus_spp','providencia_spp','pseudomonas_aeruginosa','pseudomonas_spp','neisseria_gonorrheae','shigella_spp','mycobacterium_tuberculosis','salmonella_paratyphi','salmonella_typhi','salmonella_typhi_paratyphi','serratia_spp','staphylococcus_aureus','streptococcus_pneumoniae','mycoplasma','listeria','legionella_spp')
antibiotic_group<-c('third_gen_ceph','carbapenem','fluoroquinolone','penicillin','aminopenicillin','beta_lactamase_inhibitor','anti_pseudomonal_penicillin','methicillin','vancomycin','fourth_gen_ceph','sulfa','aminoglycoside','macrolide')

source("FILENAME")
source("FILENAME")
maps_path <- 'FILEPATH'
antibiotic_map <- read.csv(paste0(maps_path,'FILENAME'),stringsAsFactors = FALSE)
pathogen_map <- read.csv(paste0(maps_path,'FILENAME'),stringsAsFactors = FALSE)
specimen_map <- read.csv(paste0(maps_path,'FILENAME'),stringsAsFactors = FALSE)
level_amr_map <- read.csv(paste0(maps_path,'FILENAME'),stringsAsFactors = FALSE)
level_amr <- level_amr_map %>% distinct(acause,level_amr)
cause_id <- get_ids("cause") %>% left_join(level_amr) %>% distinct()
ucause_map <- read.csv(paste0(maps_path,'FILENAME'),stringsAsFactors = FALSE)
ucause_map <- ucause_map %>% left_join(cause_id) %>% distinct(cause, .keep_all = TRUE)
ucause_map$level_amr[is.na(ucause_map$level_amr)] <- ucause_map$acause[is.na(ucause_map$level_amr)]

locs <- subset(get_location_metadata(location_set_id=1, gbd_round_id=6),select = c(location_id, location_name,super_region_name,super_region_id,region_id,region_name), droplevels = TRUE)
loc<-get_ids("location", return_all_columns=TRUE)

# Data A
DATA_A <- read.csv('FILENAME', stringsAsFactors = FALSE)
DATA_A <- DATA_A %>% 
  rename(raw_pathogen = pathogen, third_gen_ceph = X_3rd_gen_ceph_r, fourth_gen_ceph = X_4th_gen_ceph_r)
names(DATA_A) <- gsub("_r","",names(DATA_A))
names(DATA_A)[names(DATA_A) == 'beta_lactamase_inhib'] <- 'beta_lactamase_inhibitor'
DATA_A$raw_pathogen <- str_to_lower(DATA_A$raw_pathogen)
DATA_A$pathogen <-  gsub(" ","_",DATA_A$raw_pathogen)
DATA_A$pathogen[grep('various',DATA_A$pathogen)] <-  'non_typhoidal_salmonellae'
DATA_A$pathogen[grep('enterobacter',DATA_A$pathogen)] <-  'enterobacter_spp'
DATA_A$pathogen[grep('morganella',DATA_A$pathogen)] <-  'morganella_spp'
DATA_A$pathogen[grep('proteus',DATA_A$pathogen)] <-  'proteus_spp'
DATA_A$pathogen[grep('shigella',DATA_A$pathogen)] <-  'shigella_spp'
DATA_A$pathogen[grep('providencia',DATA_A$pathogen)] <-  'providencia_spp'
DATA_A$pathogen[grep('serratia',DATA_A$pathogen)] <-  'serratia_spp'
DATA_A$pathogen[grep('citrobacter',DATA_A$pathogen)] <-  'citrobacter_spp'
DATA_A$sterile[DATA_A$specimen %in% c('CSV: Blood','Respiratory','intraabdominal','GU: Urine')] <- 1
DATA_A$sterile[is.na(DATA_A$sterile)] <- 0
DATA_A <- DATA_A[DATA_A$pathogen %in% organism,]

#Loop to obtain coefficients of resistance adjusted for underlying cause, age, and hospital acquired infection, and sample size
adjusted_coeff <- c()
for (s in c(0:1)) {
for (i in organism) {
  for (j in antibiotic_group) {
    text_for_data <- paste0("DATA_A[DATA_A$pathogen == '",i,"' & !is.na(DATA_A$",j,") & DATA_A$sterile == ",s,",]")
    data <- eval(parse(text = paste(text_for_data)))
    text_for_table<- paste0("data %>% group_by(", j, ") %>% summarize(admit = sum(!is.na(Died)), died = sum(Died,na.rm = TRUE))")
    total <- eval(parse(text = paste(text_for_table)))
    if(!is.na(total[1,3] > 0 & total[2,3] > 0) & total[1,3] > 0 & total[2,3] > 0) {
        m <- paste0("glm(data = data, formula = Died ~", j , " + level_amr + age_spec + hai, family = 'poisson'(link='log'))")
        model <- eval(parse(text = paste(m)))
        model <- coeftest(model, vcov = vcovHC(model, type = "HC1"))
        r1 <- c('DATA_A', i, j, s, LOCATIONID,'LOCATION_NAME',model[2,1:2],total[1,1],total[1,2],total[1,3],total[2,1],total[2,2],total[2,3]) 
        adjusted_coeff <- rbind(adjusted_coeff,r1)
      }
    }
  }
}
rm(DATA_A)

# Data B
DATA_B <- read.csv('FILENAME', stringsAsFactors = FALSE) %>%
rename(raw_pathogen = pathogen, third_gen_ceph = X_3rd_gen_ceph_r, fourth_gen_ceph = X_4th_gen_ceph_r)
names(DATA_B) <- gsub("_r","",names(DATA_B))
names(DATA_B)[names(DATA_B) == 'beta_lactamase_inhib'] <- 'beta_lactamase_inhibitor'
DATA_B$raw_pathogen <- str_to_lower(DATA_B$raw_pathogen)
DATA_B$pathogen <-  gsub(" ","_",DATA_B$raw_pathogen)
DATA_B$sterile <- 1

s <- 1
for (i in organism) {
  for (j in antibiotic_group) {
    text_for_data <- paste0("DATA_B[ox$pathogen == '",i,"' & !is.na(DATA_B$",j,") & DATA_B$sterile == ",s,",]")
      data <- eval(parse(text = paste(text_for_data)))
      text_for_table<- paste0("data %>% group_by(", j, ") %>% summarize(admit = sum(!is.na(Died)), died = sum(Died,na.rm = TRUE))")
      total <- eval(parse(text = paste(text_for_table)))
      if(!is.na(total[1,3] > 0 & total[2,3] > 0) & total[1,3] > 0 & total[2,3] > 0) {
      m <- paste0("glm(data = data, formula = Died ~", j , "+ level_amr1 + age_spec + hai, family = 'poisson'(link='log'))")
      model <- eval(parse(text = paste(m)))
      model <- coeftest(model, vcov = vcovHC(model, type = "HC1"))
      r2 <- c('DATA_B', i, j, s, LOCATION_ID,'LOCATION_NAME',model[2,1:2],total[1,1],total[1,2],total[1,3],total[2,1],total[2,2],total[2,3]) 
      adjusted_coeff<- rbind(adjusted_coeff,r2)
    }
  }
}
rm(DATA_B)

# Data C
DATA_C <- read.csv('FILENAME', stringsAsFactors = FALSE) %>%
  rename(raw_pathogen = isolate1, third_gen_ceph = X3rd.gen.ceph_r, fourth_gen_ceph = X4th.gen.ceph_r)
names(DATA_C) <- gsub("_r","",names(DATA_C))
names(DATA_C)[names(DATA_C) == 'beta.latcam.beta.lactamase.inhibitor'] <- 'beta_lactamase_inhibitor'
names(DATA_C)[names(DATA_C) == 'anti.pseudomonal.penicillin'] <- 'anti_pseudomonal_penicillin'
DATA_C$raw_pathogen <- str_to_lower(DATA_C$raw_pathogen)
DATA_C <- left_join(DATA_C,pathogen_map)
DATA_C$sterile <- 1

# Small sample size does not allow confounder-adjustments for 'unadjusted' organisms
org_adjusted <- c('escherichia_coli','klebsiella_pneumoniae','non_typhoidal_salmonellae','coagulase_negative_staph','staphylococcus_aureus','streptococcus_pneumoniae')
org_unadjusted <- organism[!organism %in% org_adjusted]

s <- 1
for (i in org_adjusted) {
  for (j in antibiotic_group) {
    text_for_data <- paste0("DATA_C[DATA_C$pathogen == '",i,"' & !is.na(DATA_C$",j,") & DATA_C$sterile == ",s,",]")
    data <- eval(parse(text = paste(text_for_data)))
    text_for_table<- paste0("data %>% group_by(", j, ") %>% summarize(admit = sum(!is.na(Died)), died = sum(Died,na.rm = TRUE))")
    total <- eval(parse(text = paste(text_for_table)))
    if(!is.na(total[1,3] > 0 & total[2,3] > 0) & total[1,3] > 0 & total[2,3] > 0) {
      m <- paste0("glm(data = data, formula = Died ~", j , "+ acause , family = 'poisson'(link='log'))")
      model <- eval(parse(text = paste(m)))
      model <- coeftest(model, vcov = vcovHC(model, type = "HC1"))
      r3 <- c('DATA_C', i, j, s, LOCATION_ID,'LOCATION_NAME',model[2,1:2],total[1,1],total[1,2],total[1,3],total[2,1],total[2,2],total[2,3]) 
      adjusted_coeff<- rbind(adjusted_coeff,r3)
    }
  }
}

s <- 1
for (i in org_unadjusted) {
  for (j in antibiotic_group) {
    text_for_data <- paste0("DATA_C[DATA_C$pathogen == '",i,"' & !is.na(DATA_C$",j,") & DATA_C$sterile == ",s,",]")
    data <- eval(parse(text = paste(text_for_data)))
    text_for_table<- paste0("data %>% group_by(", j, ") %>% summarize(admit = sum(!is.na(Died)), died = sum(Died,na.rm = TRUE))")
    total <- eval(parse(text = paste(text_for_table)))
    if(!is.na(total[1,3] > 0 & total[2,3] > 0) & total[1,3] > 0 & total[2,3] > 0) {
      m <- paste0("glm(data = data, formula = Died ~", j , ", family = 'poisson'(link='log'))")
      model <- eval(parse(text = paste(m)))
      model <- coeftest(model, vcov = vcovHC(model, type = "HC1"))
      r3 <- c('DATA_C', i, j, s, LOCATION_ID,'LOCATION_NAME',model[2,1:2],total[1,1],total[1,2],total[1,3],total[2,1],total[2,2],total[2,3]) 
      adjusted_coeff<- rbind(adjusted_coeff,r3)
    }
  }
}
rm(DATA_C)

# Data D
DATA_D <- read.csv('FILENAME', stringsAsFactors = FALSE)
DATA_D$resistant[DATA_D$resistance == 'resistant'] <- 1
DATA_D$resistant[DATA_D$resistance == 'susceptible'] <- 0
DATA_D <- DATA_D %>% left_join(pathogen_map) %>% left_join(antibiotic_map) %>% left_join(ucause_map) %>% left_join(specimen_map)
antibiotic_group_D <- paste(unique(DATA_D$abx_class[DATA_D$abx_class %in% antibiotic_group],sep =','))
DATA_D <- DATA_D %>% filter(!is.na(deaths) & !is.na(abx_class) & !is.na(pathogen) & !is.na(resistant)) 
DATA_D$aux <- substring(DATA_D$sample_id,6)
aux2 <- str_locate(DATA_D$aux,"-")[,1]
DATA_D$admid <- substring(DATA_D$aux,1,aux2-1) 
rm(aux2)

# Reshape data to wide form
DATA_D <- DATA_D %>%  group_by(admid,location_id,sterile,pathogen,level_amr,abx_class) %>% 
  summarize(resistance = max(resistant),Died = max(deaths)) 
DATA_D <- DATA_D %>% pivot_wider(id_cols = c(admid,sterile,location_id,pathogen,level_amr,Died), names_from = abx_class, values_from = resistance)

# Small sample size does not allow confounder-adjustments for 'unadjusted' organisms
org_adjusted <- c('escherichia_coli','coagulase_negative_staph','staphylococcus_aureus','salmonella_typhi')
org_unadjusted <- organism[!organism %in% org_adjusted ]

s <- 1
for (i in org_adjusted) {
  for (j in antibiotic_group_D) {
    text_for_data <- paste0("DATA_D[DATA_D$pathogen == '",i,"' & !is.na(DATA_D$",j,") & DATA_D$sterile == ",s,",]")
    data <- eval(parse(text = paste(text_for_data)))
    text_for_table<- paste0("data %>% group_by(", j, ") %>% summarize(admit = sum(!is.na(Died)), died = sum(Died,na.rm = TRUE))")
    total <- eval(parse(text = paste(text_for_table)))
    if(!is.na(total[1,3] > 0 & total[2,3] > 0) & total[1,3] > 0 & total[2,3] > 0) {
      m <- paste0("glm(data = data, formula = Died ~", j , "+ level_amr , family = 'poisson'(link='log'))")
      model <- eval(parse(text = paste(m)))
      model <- coeftest(model, vcov = vcovHC(model, type = "HC1"))
      r4 <- c('DATA_D', i, j, s, LOCATION_ID,'LOCATION_NAME',model[2,1:2],total[1,1],total[1,2],total[1,3],total[2,1],total[2,2],total[2,3]) 
      adjusted_coeff<- rbind(adjusted_coeff,r4)
    }
  }
}

s <- 1
for (i in org_unadjusted) {
  for (j in antibiotic_group_D) {
    text_for_data <- paste0("DATA_D[DATA_D$pathogen == '",i,"' & !is.na(DATA_D$",j,") & DATA_D$sterile == ",s,",]")
    data <- eval(parse(text = paste(text_for_data)))
    text_for_table<- paste0("data %>% group_by(", j, ") %>% summarize(admit = sum(!is.na(Died)), died = sum(Died,na.rm = TRUE))")
    total <- eval(parse(text = paste(text_for_table)))
    if(!is.na(total[1,3] > 0 & total[2,3] > 0) & total[1,3] > 0 & total[2,3] > 0) {
      m <- paste0("glm(data = data, formula = Died ~", j , ", family = 'poisson'(link='log'))")
      model <- eval(parse(text = paste(m)))
      model <- coeftest(model, vcov = vcovHC(model, type = "HC1"))
      r4 <- c('DATA_D', i, j, s, LOCATION_ID,'LOCATION_NAME',model[2,1:2],total[1,1],total[1,2],total[1,3],total[2,1],total[2,2],total[2,3]) 
      adjusted_coeff<- rbind(adjusted_coeff,r4)
    }
  }
}


# Data E
DATA_E <- read.csv('FILENAME', stringsAsFactors = FALSE)
yr_E <- round(mean(unique(DATA_E$year_id)))
DATA_E$resistant[DATA_E$resistance == 'resistant'] <- 1
DATA_E$resistant[DATA_E$resistance == 'susceptible'] <- 0
DATA_E$admid <- substring(DATA_E$sample_id,1,17)
DATA_E$raw_pathogen[DATA_E$raw_pathogen == 'p. aeruginosa'] <- 'p. aeruginosa '
DATA_E <- left_join(DATA_E, level_amr_map, by = c('cause' = 'icd_code_with_decimal'))
DATA_E <- DATA_E %>% left_join(antibiotic_map) %>% left_join(pathogen_map) 
DATA_E <- DATA_E %>% filter(!is.na(deaths) & !is.na(abx_class) & !is.na(pathogen) & !is.na(resistant)) 

# Reshape data to wide form
DATA_E <- DATA_E %>%  group_by(admid,location_id,pathogen,level_amr,abx_class) %>% 
  summarize(resistance = max(resistant),Died = max(deaths)) 

DATA_E <- DATA_E %>% pivot_wider(id_cols = c(admid,location_id,pathogen,level_amr,Died), names_from = abx_class, values_from = resistance)

# Small sample size does not allow confounder-adjustments for 'unadjusted' organisms
antibiotic_group_E <- c(paste(names(DATA_E)[names(DATA_E) %in% antibiotic_group], sep = ","))
DATA_E <- as.data.frame(DATA_E)
for (i in organism) {
  for (j in antibiotic_group_E) {
    text_for_data <- paste0("DATA_E[DATA_E$pathogen == '",i,"' & !is.na(DATA_E$",j,"),]")
    data <- eval(parse(text = paste(text_for_data)))
    text_for_table<- paste0("data %>% group_by(", j, ") %>% summarize(admit = sum(!is.na(Died)), died = sum(Died,na.rm = TRUE))")
    total <- eval(parse(text = paste(text_for_table)))
    if(!is.na(total[1,3] > 0 & total[2,3] > 0) & total[1,3] > 0 & total[2,3] > 0) {
      m <- paste0("glm(data = data, formula = Died ~", j , "+ level_amr , family = 'poisson'(link='log'))")
      model <- eval(parse(text = paste(m)))
      model <- coeftest(model, vcov = vcovHC(model, type = "HC1"))
      r6 <- c('DATA_E', i, j,1,LOCATION_ID,'LOCATION_NAME',model[2,1:2],total[1,1],total[1,2],total[1,3],total[2,1],total[2,2],total[2,3]) 
      adjusted_coeff<- rbind(adjusted_coeff,r6)
    }
  }
}


# Data F
DATA_F <- read.csv('FILENAME'), stringsAsFactors = FALSE)
yr_F <- round(mean(unique(DATA_F$year_id)))
DATA_F$resistant[DATA_F$resistance == 'resistant'] <- 1
DATA_F$resistant[DATA_F$resistance == 'susceptible'] <- 0
DATA_F <- DATA_F %>% left_join(pathogen_map) %>% left_join(antibiotic_map) %>% left_join(ucause_map) %>% left_join(specimen_map)
DATA_F$admid <- substring(ctmrf$sample_id,1,17)
DATA_F <- DATA_F %>% filter(!is.na(deaths) & !is.na(abx_class) & !is.na(pathogen) & !is.na(resistant)) 

# Reshape data to wide form
DATA_F <- DATA_F %>%  group_by(admid,location_id,level_amr,sterile,pathogen,abx_class) %>% 
  summarize(resistance = max(resistant),Died = max(deaths)) 

DATA_F <- DATA_F %>% pivot_wider(id_cols = c(admid,location_id,sterile,pathogen,level_amr,Died), names_from = abx_class, values_from = resistance)

# Small sample size does not allow confounder-adjustments for 'unadjusted' organisms
antibiotic_group_F <- c(paste(names(DATA_F)[names(DATA_F) %in% antibiotic_group], sep = ","))
DATA_F <- as.data.frame(DATA_F)

for (s in c(0:1)) {
for (i in organism) {
  for (j in antibiotic_group_F) {
    text_for_data <- paste0("DATA_F[DATA_F$pathogen == '",i,"' & !is.na(DATA_F$",j,") & DATA_F$sterile == ",s,",]")
    data <- eval(parse(text = paste(text_for_data)))
    text_for_table<- paste0("data %>% group_by(", j, ") %>% summarize(admit = sum(!is.na(Died)), died = sum(Died,na.rm = TRUE))")
    total <- eval(parse(text = paste(text_for_table)))
    if(!is.na(total[1,3] > 0 & total[2,3] > 0) & total[1,3] > 0 & total[2,3] > 0) {
      m <- paste0("glm(data = data, formula = Died ~", j , "+ level_amr , family = 'poisson'(link='log'))")
      model <- eval(parse(text = paste(m)))
      model <- coeftest(model, vcov = vcovHC(model, type = "HC1"))
      r7 <- c('DATA_F', i, j, s,LOCATION_ID,'LOCATION_NAME',model[2,1:2],total[1,1],total[1,2],total[1,3],total[2,1],total[2,2],total[2,3]) 
      adjusted_coeff<- rbind(adjusted_coeff,r7)
      }
    }
  }
}

# Data G
DATA_G <- read.csv("FILENAME"), stringsAsFactors = FALSE)
yr_G <- round(mean(unique(DATA_G$year_id)))
DATA_G$aux <- substring(DATA_G$sample_id,10,)
aux2 <- str_locate(DATA_G$aux,"-")[,1]
DATA_G$Anon_BabyID <- substring(DATA_G$aux,1,aux2-1) 
rm(aux2)

DATA_G <- DATA_G %>% left_join(antibiotic_map) %>% left_join(pathogen_map) 
DATA_G$resistant[DATA_G$resistance == 'resistant'] <- 1
DATA_G$resistant[DATA_G$resistance == 'sensitive'] <- 0
DATA_G <- DATA_G %>% filter(!is.na(deaths) & !is.na(abx_class) & !is.na(pathogen) & !is.na(resistant)) 

DATA_G <- DATA_G %>%  group_by(Anon_BabyID,location_id,pathogen,abx_class) %>% 
  summarize(resistance = max(resistant),Died = max(deaths)) 

DATA_G <- DATA_G %>% pivot_wider(id_cols = c(Anon_BabyID,location_id,pathogen,Died), names_from = abx_class, values_from = resistance)

antibiotic_group_G <- c(paste(names(DATA_G)[names(DATA_G) %in% antibiotic_group], sep = ","))
DATA_G <- as.data.frame(DATA_G)

for (i in organism) {
  for (j in antibiotic_group_G) {
    text_for_data <- paste0("DATA_G[DATA_G$pathogen == '",i,"' & !is.na(DATA_G$",j,"),]")
    data <- eval(parse(text = paste(text_for_data)))
    text_for_table<- paste0("data %>% group_by(", j, ") %>% summarize(admit = sum(!is.na(Died)), died = sum(Died,na.rm = TRUE))")
    total <- eval(parse(text = paste(text_for_table)))
    if(!is.na(total[1,3] > 0 & total[2,3] > 0) & total[1,3] > 0 & total[2,3] > 0) {
      m <- paste0("glm(data = data, formula = Died ~", j , ", family = 'poisson'(link='log'))")
      model <- eval(parse(text = paste(m)))
      model <- coeftest(model, vcov = vcovHC(model, type = "HC1"))
      r8 <- c('DATA_G', i, j, 1, LOCATION_ID,'LOCATION_NAME',model[2,1:2],total[1,1],total[1,2],total[1,3],total[2,1],total[2,2],total[2,3]) 
      adjusted_coeff<- rbind(adjusted_coeff,r8)
    }
  }
}


# Data H
DATA_H <- read.csv("FILENAME", stringsAsFactors = FALSE) %>% 
      filter(!is.na(deaths))
yr_H <- round(mean(unique(DATA_H$year_id)))
DATA_H$aux <- str_locate(DATA_H$sample_id,"-")[,1]
DATA_H$admid <- substring(DATA_H$sample_id,DATA_H$aux+1,)
DATA_H$aux <- NULL

DATA_H <- DATA_H %>% left_join(antibiotic_map) %>% left_join(pathogen_map) %>% left_join(specimen_map)
DATA_H$sterile[DATA_H$raw_specimen == 'synovial fluid'] <- 1
DATA_H$resistant[DATA_H$resistance == 'resistant'] <- 1
DATA_H$resistant[DATA_H$resistance == 'sensitive'] <- 0
DATA_H <- DATA_H %>% filter(!is.na(deaths) & !is.na(abx_class) & !is.na(pathogen) & !is.na(resistant)) 

DATA_H <- DATA_H %>%  group_by(admid,location_id,sterile,pathogen,abx_class) %>% 
  summarize(resistance = max(resistant),Died = max(deaths)) 

DATA_H <- DATA_H %>% pivot_wider(id_cols = c(admid,location_id,sterile,pathogen,Died), names_from = abx_class, values_from = resistance)

antibiotic_group_H <- c(paste(names(DATA_H)[names(DATA_H) %in% antibiotic_group], sep = ","))
DATA_H <- as.data.frame(DATA_H)

for (s in 0:1) {
  for (i in organism) {
    for (j in antibiotic_group_H) {
    text_for_data <- paste0("DATA_H[DATA_H$pathogen == '",i,"' & !is.na(DATA_H$",j,") & DATA_H$sterile == ",s,",]")
    data <- eval(parse(text = paste(text_for_data)))
    text_for_table<- paste0("data %>% group_by(", j, ") %>% summarize(admit = sum(!is.na(Died)), died = sum(Died,na.rm = TRUE))")
    total <- eval(parse(text = paste(text_for_table)))
    if(!is.na(total[1,3] > 0 & total[2,3] > 0) & total[1,3] > 0 & total[2,3] > 0) {
      m <- paste0("glm(data = data, formula = Died ~", j , ", family = 'poisson'(link='log'))")
      model <- eval(parse(text = paste(m)))
      model <- coeftest(model, vcov = vcovHC(model, type = "HC1"))
      r9 <- c('DATA_H', i, j, s, LOCATION_ID,'LOCATION_NAME',model[2,1:2],total[1,1],total[1,2],total[1,3],total[2,1],total[2,2],total[2,3]) 
      adjusted_coeff<- rbind(adjusted_coeff,r9)
      }
    }
  }
}


# Data I
DATA_I <- read.csv('FILENAME', stringsAsFactors = FALSE)
org_unadjusted_I <- organism[organism %in% c('shigella_spp','salmonella_typhi','non_typhoidal_salmonellae','streptocuccus_pneumoniae','enterobacter_spp','enterococcus_faecalis','enterococcus_faecium','citrobacter_spp','acinetobacter_spp','proteus_spp','moraxella_spp','listeria','pseudomonas_aeruginosa','group_b_strep','group_a_strep','klebsiella_spp','haemophilus_influenzae','listeria','neisseria_gonorrheae')]
org_adjusted_I <- organism[!organism %in% org_unadjusted_I]
DATA_I <- DATA_I %>% filter(!is.na(Died) & !is.na(pathogen)) 
DATA_I <- as.data.frame(DATA_I)
antibiotic_group_I <- c(paste(names(DATA_I)[names(DATA_I) %in% antibiotic_group], sep = ","))

for (i in org_adjusted_I) {
  for (j in antibiotic_group_I) {
    text_for_data <- paste0("DATA_I[DATA_I$pathogen == '",i,"' & !is.na(DATA_I$",j,"),]")
    data <- eval(parse(text = paste(text_for_data)))
    text_for_table<- paste0("data %>% group_by(", j, ") %>% summarize(admit = sum(!is.na(Died)), died = sum(Died,na.rm = TRUE))")
    total <- eval(parse(text = paste(text_for_table)))
    if(!is.na(total[1,3] > 0 & total[2,3] > 0) & total[1,3] > 0 & total[2,3] > 0) {
      m <- paste0("glm(data = data, formula = Died ~", j , "+ level_amr + age_spec , family = 'poisson'(link='log'))")
      model <- eval(parse(text = paste(m)))
      model <- coeftest(model, vcov = vcovHC(model, type = "HC1"))
      r10 <- c('DATA_I', i, j, 1, LOCATION_ID,'LOCATION_NAME',model[2,1:2],total[1,1],total[1,2],total[1,3],total[2,1],total[2,2],total[2,3]) 
      adjusted_coeff<- rbind(adjusted_coeff,r10)
    }
  }
}

for (i in org_unadjusted_I) {
  for (j in antibiotic_group_I) {
    text_for_data <- paste0("DATA_I[DATA_I$pathogen == '",i,"' & !is.na(DATA_IS$",j,"),]")
    data <- eval(parse(text = paste(text_for_data)))
    text_for_table<- paste0("data %>% group_by(", j, ") %>% summarize(admit = sum(!is.na(Died)), died = sum(Died,na.rm = TRUE))")
    total <- eval(parse(text = paste(text_for_table)))
    if(!is.na(total[1,3] > 0 & total[2,3] > 0) & total[1,3] > 0 & total[2,3] > 0) {
      m <- paste0("glm(data = data, formula = Died ~", j , " , family = 'poisson'(link='log'))")
      model <- eval(parse(text = paste(m)))
      model <- coeftest(model, vcov = vcovHC(model, type = "HC1"))
      r10 <- c('DATA_I', i, j, 1, LOCATION_ID,'LOCATION_NAME',model[2,1:2],total[1,1],total[1,2],total[1,3],total[2,1],total[2,2],total[2,3]) 
      adjusted_coeff<- rbind(adjusted_coeff,r10)
    }
  }
}


# Data J
DATA_J <- read.csv("FILENAME", stringsAsFactors = FALSE)
DATA_J <-DATA_J %>% left_join(antibiotic_map) %>% left_join(pathogen_map) %>% left_join(specimen_map) %>% left_join(ucause_map)
yr_J <- round(mean(unique(DATA_J$year_id)))
DATA_J$admid <- substring(DATA_J$sample_id,1,7)
DATA_J$resistant[DATA_J$resistance == 'resistant'] <- 1
DATA_J$resistant[DATA_J$resistance == 'sensitive'] <- 0
DATA_J <- DATA_J %>% filter(!is.na(deaths) & !is.na(abx_class) & !is.na(pathogen) & !is.na(resistant)) 

DATA_J <- DATA_J %>%  group_by(admid,location_id,level_amr,sterile,pathogen,abx_class) %>% 
  summarize(resistance = max(resistant),Died = max(deaths)) 

DATA_J <- DATA_J %>% pivot_wider(id_cols = c(admid,location_id,level_amr,sterile,pathogen,Died), names_from = abx_class, values_from = resistance)

antibiotic_group_J <- c(paste(names(DATA_J)[names(DATA_J) %in% antibiotic_group], sep = ","))
DATA_J <- as.data.frame(DATA_J)

for (s in 0:1) {
  for (i in organism) {
    for (j in antibiotic_group_J) {
      text_for_data <- paste0("DATA_J[DATA_J$pathogen == '",i,"' & !is.na(DATA_J$",j,") & DATA_J$sterile == ",s,",]")
      data <- eval(parse(text = paste(text_for_data)))
      text_for_table<- paste0("data %>% group_by(", j, ") %>% summarize(admit = sum(!is.na(Died)), died = sum(Died,na.rm = TRUE))")
      total <- eval(parse(text = paste(text_for_table)))
      if(!is.na(total[1,3] > 0 & total[2,3] > 0) & total[1,3] > 0 & total[2,3] > 0) {
        m <- paste0("glm(data = data, formula = Died ~", j , "+ level_amr, family = 'poisson'(link='log'))")
        model <- eval(parse(text = paste(m)))
        model <- coeftest(model, vcov = vcovHC(model, type = "HC1"))
        r11 <- c('DATA_J', i, j, s, LOCATION_ID,'LOCATION_NAME',model[2,1:2],total[1,1],total[1,2],total[1,3],total[2,1],total[2,2],total[2,3]) 
        adjusted_coeff<- rbind(adjusted_coeff,r11)
      }
    }
  }
}

# Data K
DATA_K <- read.csv('FILENAME', stringsAsFactors = FALSE)
DATA_K <- DATA_K %>% left_join(antibiotic_map) %>% left_join(pathogen_map) %>% left_join(ucause_map) %>% left_join(locs)
yr_K <- round(mean(unique(DATA_K$year_id)))
DATA_K$admid <- substring(DATA_K$sample_id,6,)

DATA_K$resistant[DATA_K$resistance == 'resistant'] <- 1
DATA_K$resistant[DATA_K$resistance == 'sensitive'] <- 0
DATA_K <- DATA_K %>% filter(!is.na(deaths) & !is.na(abx_class) & !is.na(pathogen) & !is.na(resistant)) 

DATA_K <- DATA_K %>%  group_by(admid,super_region_id,super_region_name,level_amr,pathogen,abx_class) %>% 
  summarize(resistance = max(resistant),Died = max(deaths)) 

DATA_K <- DATA_K %>% pivot_wider(id_cols = c(admid,super_region_id,super_region_name,level_amr,pathogen,Died), names_from = abx_class, values_from = resistance)

antibiotic_group_K <- c(paste(names(DATA_K)[names(DATA_K) %in% antibiotic_group], sep = ","))
DATA_K <- as.data.frame(DATA_K)
adjusted_coeff<-c()
for (r in unique(DATA_K$super_region_id)) {
  for (i in organism) {
    for (j in antibiotic_group_K) {
      text_for_data <- paste0("DATA_K[DATA_K$pathogen == '",i,"' & !is.na(DATA_K$",j,") & DATA_K$super_region_id == ", r,",]")
      data <- eval(parse(text = paste(text_for_data)))
      text_for_table<- paste0("data %>% group_by(", j, ") %>% summarize(admit = sum(!is.na(Died)), died = sum(Died,na.rm = TRUE))")
      total <- eval(parse(text = paste(text_for_table)))
      if(!is.na(total[1,3] > 0 & total[2,3] > 0) & total[1,3] > 0 & total[2,3] > 0) {
        m <- paste0("glm(data = data, formula = Died ~", j , " + level_amr, family = 'poisson'(link='log'))")
        model <- eval(parse(text = paste(m)))
        model <- coeftest(model, vcov = vcovHC(model, type = "HC1"))
        r12 <- c('DATA_K', i, j, 1, r,'',model[2,1:2],total[1,1],total[1,2],total[1,3],total[2,1],total[2,2],total[2,3]) 
        adjusted_coeff <- rbind(adjusted_coeff,r12)
      }
    }
  }
}

# Data L
DATA_L_deaaths <- read.csv("FILENAME"), stringsAsFactors = FALSE)
DATA_L_deaths <- DATA_L_deaths %>% select('ID',  'Outcome') %>% filter(Outcome < 2)
DATA_L <- read.csv("FILENAME"), stringsAsFactors = FALSE)
DATA_L$raw_antibiotic <- gsub('resistance to ','',DATA_L$raw_antibiotic)
yr_L <- round(mean(unique(DATA_L$year_id)))
DATA_L$aux <- substring(DATA_L$sample_id,7,)
aux2 <- str_locate(DATA_L$aux,"_")[,1]
DATA_L$ID <- substring(DATA_L$aux,1,aux2-1) 
DATA_L$ID <- as.numeric(DATA_L$ID) 
DATA_L$aux <- NULL
rm(aux2)
DATA_L <- left_join(DATA_L,DATA_L_deaths, by = "ID")
DATA_L$Died <- DATA_L$IOutcome

DATA_L <- DATA_L%>%left_join(locs,by='location_id')
DATA_L$super_region_id[DATA_L$location_id %in% c(496,513,514,515,53614)] <- 4
DATA_L$super_region_id[DATA_L$location_id %in% c(53660,53662,53663,53667,53674,44955)] <- 31
DATA_L$super_region_id[DATA_L$location_id == 35497] <- 64
DATA_L$super_region_id[DATA_L$location_id == 44870] <- 137
unique(DATA_L$location_id[is.na(DATA_L$super_region_id)])
DATA_L$super_region_id[is.na(DATA_L$super_region_id)] <- 158

DATA_L <- DATA_L %>% left_join(antibiotic_map) %>% left_join(pathogen_map) %>% left_join(specimen_map) %>% left_join(ucause_map)
DATA_L$resistant[DATA_L$resistance == 'resistant'] <- 1
DATA_L$resistant[DATA_L$resistance == 'sensitive'] <- 0
DATA_L <- DATA_L %>% filter(!is.na(Died) & !is.na(abx_class) & !is.na(pathogen) & !is.na(resistant)) 

DATA_L <- DATA_L %>%  group_by(sample_id,ID,super_region_id,super_region_name,level_amr,sterile,pathogen,abx_class) %>% 
  summarize(resistance = max(resistant),Died = max(Died)) 

antibiotic_group_L <-paste(names(DATA_L)[names(DATA_L) %in% antibiotic_group], sep = ',')
org_unadjusted_L <- organism[organism %in% c('shigella_spp','salmonella_paratyphi','pseudomonas_spp','salmonella_typhi','non_typhoidal_salmonellae','neisseria_meningitidis','moraxella_spp','listeria','legionella_spp','group_b_strep','group_a_strep','klebsiella_spp','haemophilus_influenzae','neisseria_gonorrheae')]
org_adjusted_L <- organism[!organism %in% org_unadjusted_L]

for (s in 0:1) {
  for (r in unique(DATA_L$super_region_id)) {
    for (i in org_adjusted_L) {
      for (j in antibiotic_group_L) {
        text_for_data <- paste0("DATA_L[DATA_L$pathogen == '",i,"' & !is.na(DATA_L$", j,") & DATA_L$sterile == ", s,") & DATA_L$super_region_id == ", r,",]")
            data <- eval(parse(text = paste(text_for_data)))
            text_for_table<- paste0("data %>% group_by(", j, ") %>% summarize(admit = sum(!is.na(Died)), died = sum(Died,na.rm = TRUE))")
            total <- eval(parse(text = paste(text_for_table)))
            if(!is.na(total[1,3] > 0 & total[2,3] > 0) & total[1,3] > 0 & total[2,3] > 0) {
              m <- paste0("glm(data = data, formula = Died ~", j , " + level_amr, family = 'poisson'(link='log'))")
              model <- eval(parse(text = paste(m)))
              model <- coeftest(model, vcov = vcovHC(model, type = "HC1"))
              r13 <- c('DATA_L', i, j, s, r,'',model[2,1:2],total[1,1],total[1,2],total[1,3],total[2,1],total[2,2],total[2,3]) 
              adjusted_coeff <- rbind(adjusted_coeff,r13)
            }
          }
        }
      }
}

for (s in 0:1) {
  for (i in org_unadjusted_L) {
    for (j in antibiotic_group_L) {
        text_for_data <- paste0("DATA_L[DATA_L$pathogen == '",i,"' & !is.na(DATA_L$", j,") & DATA_L$sterile == ", s,",]")
        data <- eval(parse(text = paste(text_for_data)))
        text_for_table<- paste0("data %>% group_by(", j, ") %>% summarize(admit = sum(!is.na(Died)), died = sum(Died,na.rm = TRUE))")
        total <- eval(parse(text = paste(text_for_table)))
        if(!is.na(total[1,3] > 0 & total[2,3] > 0) & total[1,3] > 0 & total[2,3] > 0) {
          m <- paste0("glm(data = data, formula = Died ~", j , ", family = 'poisson'(link='log'))")
          model <- eval(parse(text = paste(m)))
          model <- coeftest(model, vcov = vcovHC(model, type = "HC1"))
          r13 <- c('DATA_L', i, j, s, LOCATION_ID,'LOCATION_NAME',model[2,1:2],total[1,1],total[1,2],total[1,3],total[2,1],total[2,2],total[2,3]) 
          adjusted_coeff <- rbind(adjusted_coeff,r13)
        }
      }
    }
  }


# Appending adjusted estimates 
adjusted_coeff <- as.data.frame(adjusted_coeff, stringsAsFactors = FALSE)
colnames(adjusted_coeff) <- c('refvar', 'pathogen','abx_class','sterile','location_id','location_name','lnrr','lnrr_se','susc','admit_s','died_s','res','admit_r','died_r')
adjusted_coeff$year_id <- 2016
adjusted_coeff$year_id[adjusted_coeff$location_name %in% c('LOCATION_NAME')] <- 2018
adjusted_coeff$year_id[adjusted_coeff$location_name %in% c('LOCATION_NAME')] <- 2006
adjusted_coeff$year_id[adjusted_coeff$location_name %in% c('LOCATION_NAME')] <- 2010
adjusted_coeff$year_id[adjusted_coeff$location_name %in% c('LOCATION_NAME')] <- 2009
adjusted_coeff$year_id[adjusted_coeff$refvar=='DATA_A'] = yr_A
adjusted_coeff$year_id[adjusted_coeff$refvar=='DATA_B'] = yr_B
adjusted_coeff$year_id[adjusted_coeff$refvar=='DATA_C'] = yr_C
adjusted_coeff$year_id[adjusted_coeff$refvar=='DATA_D'] = yr_D
adjusted_coeff$year_id[adjusted_coeff$refvar=='DATA_E'] = yr_E
#adjusted_coeff$location_id[adjusted_coeff$location_id == ''] <- LOCATION_ID
adjusted_coeff$location_name <- NULL
adjusted_coeff$susc <- NULL
adjusted_coeff$res <- NULL
adjusted_coeff$lnrr <- as.double(adjusted_coeff$lnrr)
adjusted_coeff$admit_r <- as.double(adjusted_coeff$admit_r)
adjusted_coeff$died_r <- as.double(adjusted_coeff$died_r)
adjusted_coeff$admit_s <- as.double(adjusted_coeff$admit_s)
adjusted_coeff$died_s <- as.double(adjusted_coeff$died_s)
adjusted_coeff$rr <- exp(adjusted_coeff$lnrr)
adjusted_coeff$lnrr_var <- as.double(adjusted_coeff$lnrr_se)^2
adjusted_coeff$crude_lnrr <- log((adjusted_coeff$died_r/adjusted_coeff$admit_r)/(adjusted_coeff$died_s/adjusted_coeff$admit_s))
adjusted_coeff$crude_lnrr_var <- (1/adjusted_coeff$died_r) - (1/adjusted_coeff$admit_r) + (1/adjusted_coeff$died_s) - (1/adjusted_coeff$admit_s)
adjusted_coeff$crude_lnrr_se <- sqrt(adjusted_coeff$crude_lnrr_var)

# Data used for unadjusted estimates 
# Data M
DATA_M <- read.csv('FILENAME', stringsAsFactors = FALSE)
DATA_M$sterile <- 1
DATA_M <- DATA_M %>% left_join(antibiotic_map) %>% left_join(pathogen_map) %>% left_join(locs,by='location_id')
DATA_M$abx_class[DATA_M$raw_antibiotic == 'carbapenems (imipenem, meropenem, doripenem)'] <- 'carbapenem'

organism_M_small <-c('moraxella_spp','enterococcus_faecium','legionella_spp','neisseria_meningitidis','non_typhoidal_salmonellae','group_a_strep','salmonella_typhi_paratyphi','haemophilus_influenzae')
organism_M <-paste(unique(DATA_M$pathogen[DATA_M$pathogen %in% organism]), sep = ',')
organism_M <-organism_M[!organism_M %in% organism_M_small]

crude_coeff <- c()
for (i in organism_M) {
    data <- DATA_M[DATA_M$pathogen %in% paste0(i),]
      data <- data %>% filter(location_id != "" & deaths != "" & abx_class != "") %>%
        select(deaths,location_id,sterile,abx_class,resistance,year_id,pathogen) %>% 
        group_by(location_id,sterile,abx_class,resistance,pathogen) %>% 
        summarize(admit=n(),died = sum(deaths), year_id=round(mean(year_id))) 
      data <- as.data.frame(data)

      data$resistance[data$resistance == 'resistant'] <- 'r'
      data$resistance[data$resistance == 'susceptible'] <- 's'
      data <- reshape(data,direction = "wide",v.names = c('admit','died'), timevar = 'resistance', idvar = c('location_id','sterile','pathogen','abx_class','year_id'),sep='_')
      data <- na.omit(data)
      data$refvar <- 'DATA_M'
      crude_coeff<- rbind(crude_coeff,data)
      data <- c()
    }

# Small sample size only allows for a super region estimate
organism_M_small <-c('haemophilus_influenzae')
for (i in organism_M_small) {
  data <- DATA_M[DATA_M$pathogen %in% paste0(i),]
  data <- data %>% filter(super_region_id != "" & deaths != "" & abx_class != "") %>%
    select(deaths,super_region_id,sterile,abx_class,resistance,year_id,pathogen) %>% 
    group_by(super_region_id,sterile,abx_class,resistance,pathogen) %>% 
    summarize(admit=n(),died = sum(deaths), year_id=round(mean(year_id))) 
  data <- as.data.frame(data)
  
  data$resistance[data$resistance == 'resistant'] <- 'r'
  data$resistance[data$resistance == 'susceptible'] <- 's'
  data <- reshape(data,direction = "wide",v.names = c('admit','died'), timevar = 'resistance', idvar = c('super_region_id','sterile','pathogen','abx_class','year_id'),sep='_')
  names(data)[names(data) == 'super_region_id'] <- 'location_id'
  data <- na.omit(data)
  data$refvar <- 'DATA_M'
  crude_coeff<- rbind(crude_coeff,data)
  data <- c()
}


# Data N
DATA_N <- read.csv('FILENAME', stringsAsFactors = FALSE)
DATA_N$sterile <- 1
DATA_N <- DATA_N %>% left_join(antibiotic_map) %>% left_join(ucause_map) 

DATA_N$pathogen <- 'PATHOGEN'

for (i in c('streptococcus_pneumoniae')) {
  data <- DATA_N[DATA_N$pathogen %in% paste0(i),]
  data <- data %>% filter(location_id != "" & deaths != "" & abx_class != "" & resistance != "unknown") %>%
    select(deaths,location_id,sterile,abx_class,resistance,year_id,pathogen) %>% 
    group_by(location_id,sterile,abx_class,resistance,pathogen) %>% 
    summarize(admit=n(),died = sum(deaths), year_id=round(mean(year_id))) 
  data <- as.data.frame(data)
  
  data$resistance[data$resistance == 'resistant'] <- 'r'
  data$resistance[data$resistance == 'susceptible'] <- 's'
  data <- reshape(data,direction = "wide",v.names = c('admit','died'), timevar = 'resistance', idvar = c('location_id','sterile','pathogen','abx_class','year_id'),sep='_')
  data <- na.omit(data) 
  data$refvar <- 'DATA_N'
  crude_coeff<- rbind(crude_coeff,data)
  data <- c()
}


# Data O
DATA_O <- read.csv('FILENAME', stringsAsFactors = FALSE)
DATA_O$sterile <- 1
DATA_O <- DATA_O %>% left_join(antibiotic_map) %>% left_join(pathogen_map)
organism_O <-paste(unique(DATA_O$pathogen[DATA_O$pathogen %in% organism]), sep = ',')
org_O_no_res <- c('PATHOGEN')
organism_O <-organism_O[!organism_O %in% org_O_no_res]

for (i in organism_O) {
  data <- DATA_O[DATA_O$pathogen %in% paste0(i),]
  data <- data %>% filter(location_id != "" & deaths != "" & abx_class != "") %>%
    select(deaths,location_id,sterile,abx_class,resistance,year_id,pathogen) %>% 
    group_by(location_id,sterile,abx_class,resistance,pathogen) %>% 
    summarize(admit=n(),died = sum(deaths), year_id=round(mean(year_id))) 
  data <- as.data.frame(data)
  
  data$resistance[data$resistance == 'resistant'] <- 'r'
  data$resistance[data$resistance == 'susceptible'] <- 's'
  data <- reshape(data,direction = "wide",v.names = c('admit','died'), timevar = 'resistance', idvar = c('location_id','sterile','pathogen','abx_class','year_id'),sep='_')
  data$refvar <- 'DATA_O'
  data <- na.omit(data) 
  crude_coeff<- rbind(crude_coeff,data)
  data <- c()
}

# Data P
DATA_P <- read.csv("FILENAME", stringsAsFactors = FALSE)
DATA_P <- DATA_P %>%  left_join(specimen_map) %>% left_join(pathogen_map) %>% left_join(antibiotic_map)
organism_P <-paste(unique(DATA_P$pathogen[DATA_P$pathogen %in% organism]), sep = ',')

for (i in organism_P) {
  data <- DATA_P[DATA_P$pathogen %in% paste0(i),]
  data <- data %>% filter(location_id != "" & deaths != "" & abx_class != "") %>%
    select(deaths,cases,location_id,sterile,abx_class,resistance,pathogen) %>% 
    group_by(location_id,sterile,abx_class,resistance,pathogen) %>% 
    summarize(admit=sum(cases),died = sum(deaths), year_id=round(mean(year_id))) 
  data <- as.data.frame(data)
  
  data$resistance[data$resistance == 'resistant'] <- 'r'
  data$resistance[data$resistance == 'susceptible'] <- 's'
  data <- reshape(data,direction = "wide",v.names = c('admit','died'), timevar = 'resistance', idvar = c('location_id','sterile','pathogen','abx_class','year_id'),sep='_')
  data$refvar <- 'DATA_P'
  data <- na.omit(data) 
  crude_coeff<- rbind(crude_coeff,data)
  data <- c()
}


crude_coeff$location_id[is.na(crude_coeff$location_id)] <- crude_coeff$super_region_id
crude_coeff$super_region_id<-NULL
crude_coeff <- crude_coeff[crude_coeff$died_r >0 & crude_coeff$died_s >0,]
crude_coeff$rr <- (crude_coeff$died_r/crude_coeff$admit_r)/(crude_coeff$died_s/crude_coeff$admit_s)
crude_coeff$lnrr <- log((crude_coeff$died_r/crude_coeff$admit_r)/(crude_coeff$died_s/crude_coeff$admit_s))
crude_coeff$lnrr_var <- (1/crude_coeff$died_r) - (1/crude_coeff$admit_r) + (1/crude_coeff$died_s) - (1/crude_coeff$admit_s)
crude_coeff$lnrr_se <- sqrt(crude_coeff$lnrr_var)
crude_coeff$crude_lnrr <- log((crude_coeff$died_r/crude_coeff$admit_r)/(crude_coeff$died_s/crude_coeff$admit_s))
crude_coeff$crude_lnrr_var <- (1/crude_coeff$died_r) - (1/crude_coeff$admit_r) + (1/crude_coeff$died_s) - (1/crude_coeff$admit_s)
crude_coeff$crude_lnrr_se <- sqrt(crude_coeff$lnrr_var)
results <- rbind(adjusted_coeff,crude_coeff)
results$sample_size <- results$admit_s + results$admit_r
repo <- 'FILEPATH'
data_path <-'FILEPATH'

library(data.table)
fwrite(results, file = paste0(data_path,'FILENAME'),row.names = FALSE)
results <- read.csv(paste0(data_path,'FILENAME'), stringsAsFactors = FALSE)

##LOADING LITERATURE REVIEW
source("FILENAME")
locs <- subset(get_location_metadata(location_set_id=1, gbd_round_id=5),select = c(location_id, location_name), droplevels = TRUE)
first_litrev <- read.csv(paste0(data_path,'FILENAME'), stringsAsFactors = FALSE)[1:10]
names(first_litrev) <- c('sterile','pathogen','refvar','year_id','abx_class','location_name','sample_size','lnrr','lnrr_var','lnrr_se')
first_litrev$rr <- exp(first_litrev$lnrr)
for (i in c('admit_r','died_r','admit_s','died_s','specimen')) {
  first_litrev$new <- 'NA' 
  names(first_litrev)[names(first_litrev)=='new'] <- paste0(i) 
}
first_litrev$crude_lnrr <- first_litrev$lnrr
first_litrev$crude_lnrr_var <- first_litrev$lnrr_var
first_litrev$crude_lnrr_se <- first_litrev$lnrr_se

second_litrev <- read.csv(paste0(data_path,'FILENAME'), stringsAsFactors = FALSE)  
second_litrev$sterile <- ifelse(second_litrev$specimen == 'BSI',1,0)
second_litrev <- second_litrev[second_litrev$died_r >0 & second_litrev$died_s >0,]
second_litrev$rr <- (second_litrev$died_r/second_litrev$admit_r)/(second_litrev$died_s/second_litrev$admit_s)
second_litrev$lnrr <- log((second_litrev$died_r/second_litrev$admit_r)/(second_litrev$died_s/second_litrev$admit_s))
second_litrev$lnrr_var <- (1/second_litrev$died_r) - (1/second_litrev$admit_r) + (1/second_litrev$died_s) - (1/second_litrev$admit_s)
second_litrev$lnrr_se <- sqrt(second_litrev$lnrr_var)
second_litrev$crude_lnrr <- log((second_litrev$died_r/second_litrev$admit_r)/(second_litrev$died_s/second_litrev$admit_s))
second_litrev$crude_lnrr_var <- (1/second_litrev$died_r) - (1/second_litrev$admit_r) + (1/second_litrev$died_s) - (1/second_litrev$admit_s)
second_litrev$crude_lnrr_se <- sqrt(second_litrev$lnrr_var)

litrev <- rbind(first_litrev,second_litrev) %>% left_join(locs)
litrev$location_name <- NULL
litrev$specimen <- NULL

all_data <- rbind(litrev,results)
names(all_data)[names(all_data)=='abx_class']='drug'
names(all_data)[names(all_data)=='sample_size']='total'

d2<-all_data[!is.na(all_data$refvar),]
d2<-d2[!is.na(d2$refvar),]
d2<-d2[!is.na(d2$vi),]
d2<-d2[!d2$refvar=='NA',]
unique(d2$refvar[d2$vi == 0])
d2<-d2[!d2$sei == 0,]
d2<-d2[!(d2$pathogen == 'group_b_strep' & d2$refvar == 'DATA_D' & d2$drug == 'macrolide'),]
d2<-d2[!(d2$pathogen == 'group_b_strep' & d2$refvar == 'DATA_D' & d2$drug == 'penicillin'),]
d2$pathogen[grep('baum',d2$pathogen)] <- 'acinetobacter_baumanii'

#Save data for modelling
write.csv(d2, paste0(data_path,'FILENAME'), row.names = FALSE)
