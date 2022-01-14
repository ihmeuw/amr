rm(list = ls())

library(tidyverse)

source(FILEPATH)
maps_path <- FILEPATH
antibiotic_map <- read.csv(paste0(maps_path,'antibiotic_map.csv'),stringsAsFactors = FALSE)
pathogen_map <- read.csv(paste0(maps_path,'pathogen_map.csv'),stringsAsFactors = FALSE)
specimen_map <- read.csv(paste0(FILEPATH),stringsAsFactors = FALSE)

location_md <- subset(get_location_metadata(location_set_id=35, gbd_round_id=6),select = c(location_id, location_name,super_region_name,super_region_id,region_id,region_name), droplevels = TRUE)

SOURCE_1 <- read.csv(FILEPATH, stringsAsFactors = FALSE) 
yr_SOURCE_1 <- round(mean(unique(SOURCE_1$year_id)))

SOURCE_1 <- SOURCE_1 %>% left_join(antibiotic_map) %>% left_join(pathogen_map) %>% left_join(specimen_map)
SOURCE_1$pathogen[grep('k. pneumoniae',SOURCE_1$raw_pathogen)] <- 'klebsiella_pneumoniae'
SOURCE_1$pathogen[grep('salmonella grupo',SOURCE_1$raw_pathogen)] <- 'non_typhoidal_salmonellae'
SOURCE_1$pathogen[grep('salmonella choleraesuis',SOURCE_1$raw_pathogen)] <- 'non_typhoidal_salmonellae'
SOURCE_1$pathogen[grep('salmonella enteritidis',SOURCE_1$raw_pathogen)] <- 'non_typhoidal_salmonellae'
SOURCE_1$pathogen[grep('salmonella_sp',SOURCE_1$raw_pathogen)] <- 'non_typhoidal_salmonellae'
SOURCE_1$pathogen[grep('k. oxytoca',SOURCE_1$raw_pathogen)] <-  'klebsiella_spp'
SOURCE_1$pathogen[grep('klebsiella sp.',SOURCE_1$raw_pathogen)] <-  'klebsiella_spp'
SOURCE_1$pathogen[grep('iwoffii',SOURCE_1$raw_pathogen)] <-  'acinetobacter_spp'
SOURCE_1$pathogen[grep('aumannii',SOURCE_1$raw_pathogen)] <-  'acinetobacter_baumanii'
SOURCE_1$pathogen[grep('accineto',SOURCE_1$raw_pathogen)] <-  'acinetobacter_spp'
SOURCE_1$pathogen[grep('acinetobacter cloacae',SOURCE_1$raw_pathogen)] <-  'acinetobacter_spp'
SOURCE_1$pathogen[grep('pseudomonas vulgaris',SOURCE_1$raw_pathogen)] <-  'pseudomonas_spp'
SOURCE_1$pathogen[grep('pseudomonas fluorescens',SOURCE_1$raw_pathogen)] <-  'pseudomonas_spp'
SOURCE_1$pathogen[grep('stutzeri',SOURCE_1$raw_pathogen)] <-  'pseudomonas_spp'

SOURCE_1 <- SOURCE_1 %>% filter(!is.na(abx_class) & !is.na(pathogen) & !is.na(resistance)) 
SOURCE_1 <- SOURCE_1[SOURCE_1$sterile==1,]
SOURCE_1$res <- ifelse(SOURCE_1$resistance == 'resistant',1,0)
SOURCE_1$sus <- ifelse(SOURCE_1$resistance == 'susceptible',1,0)
SOURCE_1$cas <- 1
SOURCE_1 <- subset(SOURCE_1,select = c('location_id','year_id','abx_class','pathogen','res','sus','cas', 'hospital_type'))
SOURCE_1 <- SOURCE_1 %>% dplyr::group_by(location_id,year_id,abx_class,pathogen,hospital_type) %>% 
  dplyr::summarize(cases = sum(cas), resistant = sum(res), susceptible = sum(sus)) 
SOURCE_1$source <- 'Peru' 
SOURCE_1 <- as.data.frame(SOURCE_1)

SOURCE_2a <- read.csv(FILEPATH, stringsAsFactors = FALSE) 
SOURCE_2a$pathogen[SOURCE_2a$pathogen == 'acinetobacter_spp'] <- 'acinetobacter_baumanii'
SOURCE_2a$hospital_type <- 'm/u'
SOURCE_2b <- read.csv(FILEPATH, stringsAsFactors = FALSE) 
SOURCE_2b$pathogen[SOURCE_2b$pathogen == 'acinetobacter_spp'] <- 'acinetobacter_baumanii'
SOURCE_2b$hospital_type <- 'm/u'

SOURCE_3 <- read.csv(FILEPATH, stringsAsFactors = FALSE) 
SOURCE_3$hospital_type <- 'm/u'

SOURCE_4a <- read.csv(FILEPATH, stringsAsFactors = FALSE) 
SOURCE_4a$hospital_type <- 'm/u'
SOURCE_4b <- read.csv(FILEPATH, stringsAsFactors = FALSE) 
SOURCE_4b$hospital_type <- 'm/u'

SOURCE_5 <- read.csv(FILEPATH, stringsAsFactors = FALSE) 
SOURCE_5$location_name[SOURCE_5$location_name=='Bolivia'] <- "Bolivia (Plurinational State of)"
SOURCE_5$location_name[SOURCE_5$location_name=='Venezuela'] <- "Venezuela (Bolivarian Republic of)"
SOURCE_5 <- SOURCE_5 %>% left_join(location_md)
summary(is.na(SOURCE_5))
SOURCE_5[,c('location_name','super_region_name','super_region_id','region_id','region_name')] <-NULL
SOURCE_5$hospital_type <- 'm/u'

SOURCE_6 <- read.csv(FILEPATH, stringsAsFactors = FALSE)
SOURCE_6$susceptible <- SOURCE_6$cases - SOURCE_6$resistant
SOURCE_6$location_name <-NULL
SOURCE_6$hospital_type <- 'm/u'

SOURCE_7 <- read.csv(FILEPATH, stringsAsFactors = FALSE)
SOURCE_7$susceptible <- SOURCE_7$cases - SOURCE_7$resistant
SOURCE_7$location_name <-NULL
SOURCE_7$hospital_type <- 'm/u'

SOURCE_8 <- read.csv(FILEPATH, stringsAsFactors = FALSE) %>%  dplyr::rename(location_name = country)
SOURCE_8$location_name[SOURCE_8$location_name=='Swaziland'] <- "Eswatini"
SOURCE_8 <- left_join(SOURCE_8,location_md)
SOURCE_8$source <- 'SOURCE_8'
SOURCE_8$susceptible <- SOURCE_8$cases - SOURCE_8$resistant
SOURCE_8[,c('location_name','super_region_name','super_region_id','region_id','region_name')] <-NULL
SOURCE_8$hospital_type <- 'm/u'

SOURCE_9 <- read.csv(FILEPATH, stringsAsFactors = FALSE) %>% dplyr::rename(year_id = year)
SOURCE_9$susceptible <- SOURCE_9$cases - SOURCE_9$resistant
SOURCE_9$location_id <- 67
SOURCE_9$v<-NULL
SOURCE_9 <- na.omit(SOURCE_9)
SOURCE_9$hospital_type <- 'm/u'

SOURCE_10 <- read.csv(FILEPATH, stringsAsFactors = FALSE) 
SOURCE_10$susceptible <- SOURCE_10$cases - SOURCE_10$resistant
SOURCE_10$location_name<-NULL
SOURCE_10$hospital_type <- 'm/u'

SOURCE_11 <- read.csv(FILEPATH, stringsAsFactors = FALSE) %>% dplyr::filter(year_id %in% c(2012,2013,2014,2015))
SOURCE_11$susceptible <- SOURCE_11$cases - SOURCE_11$resistant
SOURCE_11$location_name <- NULL
SOURCE_11$source <- "SOURCE_11"
SOURCE_11$hospital_type <- 'm/u'

SOURCE_12 <- read.csv(FILEPATH, stringsAsFactors = FALSE) 
SOURCE_12$susceptible <- SOURCE_12$cases - SOURCE_12$resistant
SOURCE_12$location_name <-NULL
SOURCE_12$age_group_id <- NULL
SOURCE_12$sex_id <- NULL
SOURCE_12$specimen <- NULL
SOURCE_12$hospital_type <- 'm/u'

SOURCE_13 <- read.csv(FILEPATH, stringsAsFactors = FALSE) %>% filter(pathogen == 'acinetobacter_spp')
SOURCE_13$pathogen <- 'acinetobacter_baumanii'
SOURCE_13$res <- ifelse(SOURCE_13$resistance == 'resistant',SOURCE_13$cases,0)
SOURCE_13$sus <- ifelse(SOURCE_13$resistance == 'susceptible',SOURCE_13$cases,0)
SOURCE_13 <- subset(SOURCE_13,select = c('location_id','year_id','abx_class','pathogen','res','sus','cases'))
SOURCE_13 <- SOURCE_13 %>% dplyr::group_by(location_id,year_id,abx_class,pathogen) %>% 
  dplyr::summarize(resistant = sum(res), susceptible = sum(sus), case = sum(cases)) 
SOURCE_13$cases <- SOURCE_13$resistant + SOURCE_13$susceptible 
SOURCE_13$case <- NULL
SOURCE_13$source <- 'SOURCE' 
SOURCE_13 <- as.data.frame(SOURCE_13)
SOURCE_13$hospital_type <- 'm/u'

extra_data <- rbind(SOURCE_1,SOURCE_2a,SOURCE_2b,SOURCE_3,SOURCE_4a,SOURCE_4b,SOURCE_5,SOURCE_6,SOURCE_7,SOURCE_8,SOURCE_9,SOURCE_10,SOURCE_11,SOURCE_12,SOURCE_13)
initial <- read.csv(FILEPATH, stringsAsFactors = F)

SOURCES <- c('SOURCE_1','SOURCE_2a','SOURCE_2b','SOURCE_3','SOURCE_4a','SOURCE_4b','SOURCE_5','SOURCE_6','SOURCE_7','SOURCE_8','SOURCE_9','SOURCE_10','SOURCE_11','SOURCE_12','SOURCE_13')
initial <- initial[!initial$source %in% c(SOURCES), !colnames(initial) %in% 'index']

initial_data <- rbind(initial,extra_data)
write.csv(initial_data,FILEPATH, row.names = F)