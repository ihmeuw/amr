Sys.setenv(MKL_VERBOSE=0)

.libPaths(c(.libPaths(), "FILEPATH"))
library(data.table)
library(boot)
library(tidyr)
library(arm)
library(merTools)
library(utils)
library(MASS)
library(dplyr)

source(paste0(Sys.getenv("HOME"), "FILEPATH"))

args <- commandArgs(trailingOnly = T)
description <- as.character(args[1])
int_cause <- as.character(args[2])
dir <- paste0(as.character(args[3]), "/")

task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))
if (is.na(task_id)) {
    year_id <- as.integer(args[4])
    cause_id <- as.integer(args[5])
} else {
    array_df <- fread(paste0(dir, "predictions_array_job_template.csv"))
    year_id <- as.integer(array_df[task_id, year_id])
    cause_id <- as.integer(array_df[task_id, cause_id])
    print(paste0("task id ", task_id))
}

print(paste0("year id ", year_id))
print(paste0("cause id ", cause_id))
print(paste0("model description ", description))

get_prediction_interval <- function(model, template, by_age=FALSE) {
    if ((nrow(template) == 0) & (by_age)) {
        return(data.frame())
    } else {
        template <- convert_factor_variables(template)
        print(paste0(Sys.time(), " Getting predictions for 1000 simluations"))
        df_pred <- predictInterval(merMod=model, newdata=template, level=0.95, n.sims=1000, stat="mean",
                                   type="linear.prediction", returnSims=TRUE, include.resid.var=FALSE,
                                   seed=69, which="fixed")
        draws <- attr(df_pred, "sim.results")

        lvl1_cause <- as.character(unique(template$level_1))
        lvl2_cause <- as.character(unique(template$level_2))
        stopifnot(length(lvl1_cause) == 1)
        stopifnot(length(lvl2_cause) == 1)
        cause_effects <- ranef(model)
        lvl1_effect <- cause_effects$level_1[lvl1_cause, ]
        lvl2_effect <- cause_effects[['level_2:level_1']][paste(lvl2_cause, lvl1_cause, sep = ":"), ]
        if (!any(is.na(c(lvl1_effect, lvl2_effect)))) {
            draws <- draws + lvl1_effect + lvl2_effect
        } else if (!is.na(lvl1_effect)) {
            print("No conditional level 2 effect available, using only level 1")
            draws <- draws + lvl1_effect
        } else if (!is.na(lvl2_effect)) {
            stop("You have a null level 1 effect but a not null level 2, how is this possible?")
        } else {
            print("No conditional random effects available, using only fixed effects")
        }

        draws <- inv.logit(draws)
        template <- cbind(template, draws %>% as.data.table %>% setnames(paste0("draw_", 0:999)))
        return(template)
    }
}

get_prediction_ci <- function(model) {
    set.seed(69)
    draws <- MASS::mvrnorm(n=1000, mu=coef(model), Sigma=vcov(model))
    targets <- fread(paste0(dir, "cause_list.csv"))
    template <- data.frame()
    for (cause_id in targets$keep_causes) {
        df <- fread(paste0(dir, year_id, "/", cause_id, "_template.csv"))
        template <- rbind(template, df)
    }
    template <- convert_factor_variables(template)
    df_pred <- model.matrix(formula(model)[-2], data = template)
    y_hat <- (df_pred %>% as.matrix) %*% (draws %>% as.matrix %>% t)
    y_hat <- inv.logit(y_hat)
    cbind(template, y_hat %>% as.data.table %>% setnames(paste0("draw_", 0:999)))
}

set_sepsis_fraction_to_1 <- function(predictions, dir, year_id, cause_id) {
    if (is.element(cause_id, c(368, 383))) {
        write.csv(predictions, paste0(dir, year_id, "/", cause_id, "_diagnostic.csv"), row.names=F)
        predictions[, paste0("draw_", 0:999)] = 1
    }
    return(predictions)
}

if (grepl("simple_glm", description)) {
    model <- readRDS(paste0(dir, int_cause, "_model.rds"))
    predictions <- get_prediction_ci(model)
    print(paste0(Sys.time(), " Saving output"))
    for (cause in unique(predictions$cause_id)) {
        print(paste0(Sys.time(), " Saving output for cause_id ", cause))
        write.csv(predictions[predictions$cause_id == cause], paste0(dir, year_id, "/", cause, ".csv"), row.names=F)
    }
} else {
    if (grepl("by_age", description)) {
        sub_dirs <- list.files(dir)
        age_group_ids <- sub_dirs[nchar(sub_dirs) < 4]
        predictions <- data.frame()
        for (age in age_group_ids) {
            print(paste0("Working on age_group_id: ", age))
            age_dir <- paste0(dir, age, '/')
            model <- readRDS(paste0(age_dir, int_cause, "_model.rds"))
            template <- fread(paste0(age_dir, year_id, '/', cause_id, '_template.csv'))
            if (dim(template)[1] != 0) {
                age_predictions <- get_prediction_interval(model, template, by_age=TRUE)
                predictions <- rbind(predictions, age_predictions, fill=TRUE)
            }
        }
    } else {
        model <- readRDS(paste0(dir, int_cause, "_model.rds"))
        template <- fread(paste0(dir, year_id, '/', cause_id, '_template.csv'))
        predictions <- get_prediction_interval(model, template)
    }

    if ("detailed_age_group_id" %in% colnames(predictions)) { 
        predictions$detailed_age_group_id <- as.factor(predictions$detailed_age_group_id)
        predictions[is.na(detailed_age_group_id), detailed_age_group_id := age_group_id]
        
        predictions <- predictions %>%
            dplyr::rename(
                agg_age_group_id = age_group_id,
                age_group_id = detailed_age_group_id
                )
        predictions <- subset(predictions, select = -c(agg_age_group_id) )
    }

    print(paste0(Sys.time(), " Saving output"))
    if ((int_cause == "sepsis") & grepl("mortality", dir)) {predictions <- set_sepsis_fraction_to_1(predictions, dir, year_id, cause_id)}
    write.csv(predictions, paste0(dir, year_id, "/", cause_id, ".csv"), row.names=F)
}