convert_factor_variables <- function(df) {
    df[, c("age_group_id", "sex_id") := lapply(.SD, as.factor), .SDcols = c("age_group_id", "sex_id")]
    if ("cause_id" %in% colnames(df)) {
        df$cause_id <- factor(df$cause_id)
    } else if ("level_2" %in% colnames(df)) {
        df[, c("level_2", "level_1") := lapply(.SD, as.factor), .SDcols = c("level_2", "level_1")]
    }
    return(df)
}

get_covariates <- function(int_cause) {
  covariates_df <- fread(paste0("FILEPATH1", "FILEPATH2"))
  int_cause_str <- int_cause
  covariates <- covariates_df[covariates_df$int_cause == int_cause_str, covariates]
  covariates <- unlist(strsplit(covariates, ", "))
  return(covariates)
}

run_model <- function(data, formula, int_cause, data_dir, fit_all = FALSE, save=FALSE, fold="None") {
  print(formula)

  if (any(grepl("\\(1 \\| ", format(formula)))) {
    print(paste0(Sys.time(), " Running mixed effect binomial for ", int_cause))
    model <- glmer(
    formula = formula, data = data, family = binomial(link = "logit"),
    verbose = 1, control = glmerControl(optimizer = c("bobyqa", "Nelder_Mead"),
      optCtrl = list(maxfun = 2e5, verbose = 1))
    )
    converged <- check_model_converged(model)

  } else {
    print(paste0(Sys.time(), " Running simple binomial for ", int_cause))
    model <- glm(formula = formula, data = data, family = binomial(link = "logit"))
  }
if (save == TRUE){
  print(paste0(Sys.time(), " Saving model output"))
  print(summary(model))
  if (fold != "None"){
    saveRDS(model, paste0(data_dir, "/", int_cause, "_model_fold_", fold, ".rds"))
  }else{
    saveRDS(model, paste0(data_dir, "/", int_cause, "_model.rds"))
  }
}
return(model)
}

check_model_converged <- function(model) {
  if (isSingular(model, tol=1e-05)) {check1 <- FALSE} else {check1 <- TRUE}
  devfun <- update(model, devFunOnly = TRUE)
  if (isLMM(model)) {
    pars <- getME(model, "theta")
  } else {
    pars <- getME(model, c("theta", "fixef"))
  }
  if (require("numDeriv")) {
    cat("hess:\n")
    print(hess <- hessian(devfun, unlist(pars)))
    cat("grad:\n")
    print(grad <- grad(devfun, unlist(pars)))
    cat("scaled gradient:\n")
  }
  dd <- model@optinfo$derivs
  if (with(dd, max(abs(solve(Hessian, gradient))) < 2e-3)) {
    check2 <- TRUE
  } else {
    check2 <- FALSE
  }
  if (check2 & check1) {converged <- TRUE} else {converged <- FALSE}
  print(paste0("Model passed convergence checks? ", converged))
  return(converged)
}

get_rmse <- function(df) {
  df[, residual_squared := (pred_frac - obs_fraction)^2]
  rmse <- sqrt(mean(df$residual_squared))
  return(rmse)
}

save_pre_model_diagnostics <- function(df, dir, covariates) {
  if ("level_1" %in% colnames(df)) {
    for (level in unique(df$level_1)) {
        level_data <- subset(df, level_1 == level)
        level_data[, logit_fraction := logit(obs_fraction)]
        ggplot(level_data, aes(x=logit_fraction, color=sex_id)) + geom_histogram()
        ggsave(paste0(dir, "/hist_", level, ".pdf"), width=11, height=8.5)
    }
  }
  for (covariate in covariates) {
    ggplot(data, aes_string(x=covariate, fill="iso3", color="iso3")) +
      geom_histogram(position="identity", alpha=0.5) + theme_classic(base_size = 15) +
      labs(x = covariate, y = "Count") +
      scale_y_continuous(labels = scales::comma)
    ggsave(paste0(dir, "/hist_", covariate, ".pdf"), width = 11, height = 8.5)
  }

  ggplot(df, aes(x = iso3, y = obs_fraction, fill = sex_id)) + geom_boxplot()
  ggsave(paste0(dir, "/boxplot_iso3.pdf"), width = 11, height = 8.5)
  ggplot(df, aes(x = level_1, y = obs_fraction, fill = sex_id)) + geom_boxplot()
  ggsave(paste0(dir, "/boxplot_level1.pdf"), width = 11, height = 8.5)
}

make_predictions <- function(dt, mod) {
  dt$age_group_id <- as.factor(dt$age_group_id)
  dt$sex_id <- as.factor(dt$sex_id)
  dt$pred_frac <- predict(mod, newdata = dt, type = "response", allow.new.levels=TRUE)
  return(dt)
}

get_roc_auc <- function(dt) {
  thresholds <- 0:50 / 50
  results <- as.list(rep(NA, length(thresholds)))
  counter <- 1
  for (i in thresholds) {
    dt[, (c("TP", "TN", "FP", "FN")) := 0]
    dt[pred_frac >= i, TP := successes]
    dt[pred_frac >= i, FP := failures]
    dt[pred_frac < i, FN := successes]
    dt[pred_frac < i, TN := failures]
    result <- dt[, lapply(.SD, sum), .SDcols = c("TP", "FP", "FN", "TN")]
    result[, threshold := i]
    results[[counter]] <- result
    counter <- counter + 1
  }
  results <- rbindlist(results)
  results[, tpr := TP / (TP + FN)]
  results[, fpr := FP / (FP + TN)]
  roc_curve <- ggplot(results, aes(fpr, tpr)) + geom_line() + geom_point() +
    ggtitle("Receiver operating characteristic (ROC) curve")
  auc_score <- AUC(results$fpr, results$tpr, method = "trapezoid")
  return(list(roc_curve = roc_curve, auc_score = auc_score))
}

get_confusion_matrix <- function(dt) {
  dt[, pred_successes := pred_frac * deaths]
  dt[, TP := pmin(successes, pred_successes)]
  dt[, TN := pmin(failures, deaths - pred_successes)]
  dt[, FP := pmax(0, pred_successes - successes)]
  dt[, FN := pmax(0, deaths - pred_successes - failures)]
  
  confusion_matrix <- dt[, lapply(.SD, sum), .SDcols = c("TP", "FP", "FN", "TN")]
  return(confusion_matrix)
}

get_bin_class_metrics <- function(cm) {
  cm <- copy(cm)
  cm[, sensitivity := TP / (TP + FN)]
  cm[, specificity := TN / (TN + FP)]
  cm[, accuracy := (TP + TN) / (TP + TN + FP + FN)]
  return(cm)
}

get_regularized_residuals <- function(dt) {
  dt[, reg_res := (obs_fraction - pred_frac) * sqrt(deaths) / sqrt(pred_frac * (1 - pred_frac))]
  return(dt)
}

save_betas <- function(model, df, dir, int_cause, fold){
  df <- make_predictions(df, model)

  se <- sqrt(diag(vcov(model)))
  tab <- data.table(cbind(beta = fixef(model), lower = fixef(model) - 1.96 * se, upper = fixef(model) + 1.96 * se))
  exp_tab <- setnames(exp(tab), c("exp_beta", "exp_lower", "exp_upper"))
  covs <- cbind(covariates = names(fixef(model)))
  write.csv(cbind(covs, tab, exp_tab), paste0("FILEPATH", "_model_coefficients.csv"), row.names = FALSE)
}

save_diagnostics <- function(dt, dir, int_cause) {
  print(paste0(Sys.time(), " Saving diagnostics"))

  
  rmse <- get_rmse(dt)
  res <- get_roc_auc(dt)
  cm <- get_confusion_matrix(dt)
  cm <- get_bin_class_metrics(cm)
  cm[, auc_score := res$auc_score]
  dt <- get_regularized_residuals(dt)
  reg_res_plot <- ggplot(dt, aes(pred_frac, reg_res)) +
    geom_hex() +
    ggtitle("Regularized residuals versus predicted fraction")
  
  out_dir <- dir
  ggsave(
    sprintf("FILENAME", int_cause),
    plot = res$roc_curve,
    device = "pdf",
    path = out_dir
  )
  file_name <- sprintf(
    "FILENAME", int_cause)
  write.csv(cm, file = paste0(out_dir, "/", file_name), row.names = FALSE)
  ggsave(
    sprintf("FILENAME", int_cause),
    plot = reg_res_plot,
    device = "pdf",
    path = out_dir
  )

  file_name_rmse <- sprintf(
    "IS_%s_RMSE.csv", int_cause)
  write.csv(rmse, file = paste0(out_dir, "/", file_name_rmse), row.names = FALSE)
}