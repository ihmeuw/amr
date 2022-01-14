rm(list=ls())
library(ggplot2)
library(RColorBrewer)
library(plotly)
library(data.table)
library(argparse)

convert_to_factor <- function(df, factor_cols) {
  # Convert relevant columns to factors
  factor_cols <- factor_cols[factor_cols %in% names(df)]
  if (length(factor_cols) > 0) {
    df[, (factor_cols) := lapply(.SD, factor), .SDcols = factor_cols]
  }
  return(df)
}

compute_rmse <- function(obs, pred) {
  return(sqrt(mean((obs - pred)^2)))
}

compute_coef_deter <- function(obs, pred) {
  # Compute coefficient of determination (R^2)
  rss <- sum((pred - obs) ^ 2)
  tss <- sum((obs - mean(obs)) ^ 2)
  return(1 - rss / tss)
}

compute_bias <- function(obs, pred) {
  return(mean(pred - obs))
}

make_plots <- function(data_orig, resids, prop_resids, predictions,
                       continuous_covs) {
  plots <- list()
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  set.seed(42)

  # Input data
  for (cov in continuous_covs) {
    plots[[cov]] <- ggplot(data_orig, aes_string(cov, fill = "source")) +
      geom_histogram() +
      scale_fill_manual(values = sample(getPalette(length(unique(data_orig$source)))))
  }

  # Model fit - separate for out of sample and in-sample
  resids[, fit_type := ifelse(holdout == 'no_holdout', 'in_sample', 'out_of_sample')]
  resids[, log_ratio_pred := log_ratio - resid]
  rmse_table <- data.table()
  for (ft in c("in_sample", "out_of_sample")) {
    plot_resids <- resids[fit_type == ft]
    plots[[paste0('resid_vs_pred_', ft)]] <- ggplot(
      plot_resids, aes(x = log_ratio_pred, y = resid)) +
      geom_hex() +
      geom_abline(slope = 0, intercept = 0)

    plots[[paste0('resid_vs_pred_by_pathogen_', ft)]] <- ggplot(
      plot_resids, aes(x = log_ratio_pred, y = resid)) +
      geom_hex() +
      geom_abline(slope = 0, intercept = 0) +
      facet_grid(rows = vars(pathogen_x), cols = vars(pathogen_y))

    plots[[paste0('pred_vs_obs_', ft)]] <- ggplot(
      plot_resids, aes(x = log_ratio, y = log_ratio_pred)) +
      geom_hex() +
      geom_abline(slope = 1, intercept = 0)

    rmse <- compute_rmse(exp(plot_resids$log_ratio), exp(plot_resids$log_ratio_pred))
    log_rmse <- compute_rmse(plot_resids$log_ratio, plot_resids$log_ratio_pred)
    rsq <- compute_coef_deter(plot_resids$log_ratio, plot_resids$log_ratio_pred)
    bias <- compute_bias(plot_resids$log_ratio, plot_resids$log_ratio_pred)
    rmse_table <- rbind(rmse_table, data.table(
      "fit_type" = ft, "unit" = "ratio",
      "rmse" = rmse, "log_rmse" = log_rmse,
      "rsq" = rsq, "bias" = bias))
  }
  # Model fit in proportion space
  prop_resids[, fit_type := ifelse(holdout == 'no_holdout', 'in_sample', 'out_of_sample')]
  prop_resids[, resid := prop - prop_pred]
  for (ft in c("in_sample", "out_of_sample")) {
    plot_resids <- prop_resids[fit_type == ft]

    plots[[sprintf("resid_vs_pred_%s_%s", ft, "prop")]] <- ggplot(
      plot_resids, aes(x = prop_pred, y = resid)) +
      geom_hex() +
      geom_abline(slope = 0, intercept = 0)

    plots[[sprintf("pred_vs_obs_%s_%s", ft, "prop")]] <- ggplot(
      plot_resids, aes(x = prop, y = prop_pred)) +
      geom_hex() +
      geom_abline(slope = 1, intercept = 0)

    rmse <- compute_rmse(plot_resids$prop, plot_resids$prop_pred)
    rsq <- compute_coef_deter(plot_resids$prop, plot_resids$prop_pred)
    bias <- compute_bias(plot_resids$prop, plot_resids$prop_pred)
    rmse_table <- rbind(rmse_table, data.table(
      "fit_type" = ft, "unit" = "proportion",
      "rmse" = rmse, "rsq" = rsq, "bias" = bias), fill = TRUE)
  }
  plots[['rmse']] <- rmse_table

  # Distributions
  ages <- fread("FILEPATH")
  ages[, age_group_id := factor(age_group_id)]
  predictions[, agg_age_group_id := factor(agg_age_group_id)]
  predictions <- merge(
    predictions, ages[, c("age_group_id", "age_group_name")],
    by.x = c("agg_age_group_id"), by.y = c("age_group_id"),
    all.x = TRUE, all.y = FALSE
  )
  stopifnot(!any(is.na(predictions$age_group_name)))
  colors <- sample(getPalette(length(unique(predictions$pathogen))))
  for (hosp_val in unique(predictions$hosp)) {
    for (prop in c("prop_deaths", "prop_cases")) {
      if (hosp_val != 'all') {
        plot_title <- paste0(prop, hosp_val)
      } else {
        plot_title <- prop
      }
      plots[[plot_title]] <- ggplot(
          predictions[hosp == hosp_val,],
          aes_string(fill = "pathogen", y = prop, x = "super_region_name")) + 
        geom_bar(position = "stack", stat = "summary", fun = mean) +
        scale_fill_manual(values = colors) +
        facet_wrap(~ age_group_name) +
        theme(axis.text.x = element_text(angle = 90))
    }
  }

  for (age_group in unique(predictions$age_group_name)) {
    plot_df <- predictions[age_group_name == age_group & hosp != 'unknown', ]
    plot_df <- rbindlist(
      list(
        plot_df[, .(pathogen, super_region_name, hosp, prop_deaths)],
        plot_df[, .(pathogen, super_region_name, hosp, prop_cases)]
      ), fill = TRUE
    )
    plot_df[, metric := "cases"][is.na(prop_cases), metric := 'deaths']
    plot_df[, value := prop_cases][is.na(prop_cases), value := prop_deaths]
    plots[[paste0("prop_", age_group)]] <- ggplot(
      plot_df, aes_string(fill = "pathogen", y = "value", x = "super_region_name")) +
      geom_bar(position = "stack", stat = "summary", fun = mean) +
      scale_fill_manual(values = colors) +
      facet_grid(rows = vars(hosp), cols = vars(metric)) +
      theme(axis.text.x = element_text(angle = 90))
  }
  return(plots)
}

save_plots <- function(model_dir, plots) {
  for (plot_name in names(plots)) {
    print(sprintf("Working on plot %s", plot_name))
    if (grepl("resid_vs_pred|pred_vs_obs", plot_name)) {
      if (grepl("resid_vs_pred_by_pathogen", plot_name)) {
        ggsave(
          paste0(plot_name, ".pdf"),
          plot = plots[[plot_name]],
          device = "pdf",
          path = model_dir,
          height = 20,
          width = 20,
          units = "in"
        )
      } else {
        ggsave(
          paste0(plot_name, ".pdf"),
          plot = plots[[plot_name]],
          device = "pdf",
          path = model_dir
        )
      }
    } else if (plot_name == 'rmse') {
      write.csv(
        plots[[plot_name]], "FILEPATH",
        row.names = FALSE
      )
    } else {
      fig <- ggplotly(plots[[plot_name]])
      htmlwidgets::saveWidget(fig, "FILEPATH")
    }
  }
}

main <- function(model_dir, infectious_syndrome, encode_cols,
                 cov_cols) {
  continuous_covs <- cov_cols[!(cov_cols %in% encode_cols)]

  data_orig <- fread("FILEPATH")
  resids <- fread("FILEPATH")
  prop_resids <- fread("FILEPATH")
  preds <- fread("FILEPATH")

  data_orig <- convert_to_factor(data_orig, encode_cols)
  resids <- convert_to_factor(resids, encode_cols)
  prop_resids <- convert_to_factor(prop_resids, encode_cols)
  preds <- convert_to_factor(preds, encode_cols)

  plots <- make_plots(
    data_orig, resids, prop_resids, preds, continuous_covs)
  save_plots(model_dir, plots)
}

parser <- ArgumentParser(description = "Arguments for plotting")
parser$add_argument('model_dir', type = "character")
parser$add_argument('infectious_syndrome', type = "character")
parser$add_argument('--encode_cols', type = "character", nargs="*", default=c())
parser$add_argument('--cov_cols', type = "character", nargs="*", default=c())
args <- parser$parse_args()
if (!grepl("/$", args$model_dir)) {
  args$model_dir <- paste0(args$model_dir, "/")
}
main(
  args$model_dir,
  args$infectious_syndrome,
  args$encode_cols,
  args$cov_cols
)
