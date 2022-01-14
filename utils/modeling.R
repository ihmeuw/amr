get_crosswalk_filters <- function(crosswalk_number){
  
  df <- fread(paste0(H_ROOT, "FILEPATH"))
  filters <- df[xwalk_number==crosswalk_number]
  if (nrow(filters)!=1){
    stop("you need more specific crosswalk filters")
    quit(save="yes")
  }
  
  reference <- filters[, reference]
  alternates <- filters[, alternates]
  pathogen <- filters[, pathogen]
  
  alternates <- trimws(strsplit(alternates, ";")[[1]])
  return(list(alternates=alternates, reference=reference, pathogen=pathogen))
  
  
}

get_case_definitions <- function(pathogen_int, syndrome, timestamp, model_label){
  
  df <- fread(paste0(H_ROOT, "FILEPATH"))
  filters <- df[pathogen_integer==pathogen_int & infectious_syndrome==syndrome & date==timestamp & model==model_label]
  if (nrow(filters)!=1){
    stop("you need more specific crosswalk filters")
    quit(save="yes")
  }
  reference <- filters[, reference]
  alternates <- filters[, alternates]
  print(alternates)
  
  alternates <- trimws(strsplit(alternates, ";")[[1]])
  return(list(alternates=alternates, reference=reference))
}

deltamethod_logit <- function(dt, mean_col, se_col, transformed_col, to_normal = FALSE){
  dt <- copy(dt)
  if (!to_normal) {
    dt[, eval(transformed_col) := sqrt((1/(get(mean_col) - get(mean_col)^2))^2 * get(se_col)^2)]
  } else {
    dt[, eval(transformed_col) := sqrt((exp(get(mean_col))/ (1 + exp(get(mean_col)))^2)^2 * get(se_col)^2)]
  }
  
  return(dt[])
}

print_log_message <- function(message){
  print(paste(message, ":", Sys.time()))
}