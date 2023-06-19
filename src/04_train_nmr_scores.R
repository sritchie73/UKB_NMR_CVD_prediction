library(data.table)
library(foreach)
library(survival)
library(glmnet)
registerDoMC(10) # recommend requesting more cores than this in sbatch for memory reasons

# Setup array task information
tasklist <- expand.grid(prediction_cv_testfold = 0:4, endpoint = c("CHD", "Stroke"), sex=c("Male", "Female"))
setDT(tasklist)
this_task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
this_task <- tasklist[this_task]

this_test_fold <- this_task$prediction_cv_testfold
this_endpoint <- this_task$endpoint
this_sex <- this_task$sex
message(sprintf("Training NMR score for %s in %ss in prediction test fold %s\n", this_endpoint, this_sex, this_test_fold))

# Load imputed data
dat <- fread("data/imputed/analysis_cohort.txt")

# set up factor column for smoking
dat[, smoking := factor(smoking, levels=c(FALSE, TRUE))]

# Set up model matrix formula
mf <- as.formula(sprintf("~ %s + %s", 
  "age*smoking + age*sbp + age*tchol + age*hdl", # SCORE2 formula
  paste(sprintf("age*%s", ukbnmr::nmr_info[Type == "Non-derived", Biomarker]), collapse=" + ")
))

# Set up SCORE2 terms we want to not apply the lasso penalty to
score2_terms <- c("age", "smokingTRUE", "sbp", "tchol", "hdl", "age:smokingTRUE", "age:sbp", "age:tchol", "age:hdl")

# Make output directory
out_dir <- sprintf("analyses/nmr_score_training/test_fold_%s/%s/%s", this_test_fold, this_endpoint, this_sex)
system(sprintf("mkdir -p %s", out_dir))

# Extract subset of training data
this_dat <- dat[sex == this_sex & prediction_cv_foldid != this_test_fold]

# Extract endpoint
if (this_endpoint == "CHD") {
  survdat <- this_dat[,.(followup=incident_followup, event=cvd_primarily_chd)]
} else if (this_endpoint == "Stroke") {
  survdat <- this_dat[,.(followup=incident_followup, event=cvd_primarily_stroke)]
}

# Create model matrix of predictor terms
inmat <- model.matrix(mf, this_dat)

# Create vector of penalties: 1 = apply elasticnet, 0 = always include
penalties <- rep(1, ncol(inmat))
names(penalties) <- colnames(inmat)
penalties[match(score2_terms, colnames(inmat))] <- 0

# Setup list of alpha mixing parameters to search across (controls balance of ridge vs. lasso)
alphas <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1) # 0 = ridge, 1 = lasso

# Run elastinet with the given alpha
cv.coxnet.list <- foreach(this_alpha = alphas, .inorder=TRUE) %do% {
  # Fit elastinet in 10-fold cross-validation with the given alpha penalty
  cv.coxnet <- cv.glmnet(inmat, Surv(survdat$followup, survdat$event), family="cox",
                         alpha = this_alpha, penalty.factor = penalties, 
                         foldid = this_dat$elasticnet_cv_foldid + 1, # Errors where foldid = 0
                         trace.it=1, parallel=TRUE)

  # Plot
  pdf(width=7.2, height=5, file=sprintf("%s/glmnet_alpha_%s.pdf", out_dir, this_alpha))
  plot(cv.coxnet)
  dev.off()

  # Return
  return(cv.coxnet)
} 
names(cv.coxnet.list) <- alphas

# Save
saveRDS(cv.coxnet.list, file=sprintf("%s/cv_coxnet_list.rds", out_dir))

# build table of measure of best fits (partial likelihood deviance) for each lambda and alpha
fitdt <- foreach(idx = seq_along(alphas), .combine=rbind) %do% {
data.table(
  alpha=alphas[idx], 
  lambda=cv.coxnet.list[[idx]][["lambda"]],
  fit_metric=cv.coxnet.list[[idx]][["name"]],
  fit_mean=cv.coxnet.list[[idx]][["cvm"]], 
  fit_sd=cv.coxnet.list[[idx]][["cvsd"]],
  fit_mean_minus_sd=cv.coxnet.list[[idx]][["cvlo"]], 
  fit_mean_plus_sd=cv.coxnet.list[[idx]][["cvup"]],
  nonzero=cv.coxnet.list[[idx]][["nzero"]],
  lambda.min=cv.coxnet.list[[idx]][["lambda.min"]],
  lambda.1se=cv.coxnet.list[[idx]][["lambda.1se"]])
}
fwrite(fitdt, sep="\t", quote=FALSE, file=sprintf("%s/cv_coxnet_list_all_fits.txt", out_dir))

# Build table containing information about the best fits in each case:
bestfit <- fitdt[lambda == lambda.min | lambda == lambda.1se]
bestfit[lambda == lambda.min, model := "lambda.min"]
bestfit[lambda == lambda.1se, model := "lambda.1se"]
bestfit <- bestfit[, .(alpha, lambda, model, nonzero, fit_metric, fit_mean, fit_sd, fit_mean_minus_sd, fit_mean_plus_sd)]
fwrite(bestfit, sep="\t", quote=FALSE, file=sprintf("%s/cv_coxnet_best_fits.txt", out_dir))

# Choose best lambda.min and lambda.1se across models
bestbestfit <- bestfit[,.SD[which.min(fit_mean)],by=model]
fwrite(bestbestfit,  sep="\t", quote=FALSE, file=sprintf("%s/cv_coxnet_best_best_fits.txt", out_dir))

# Extract non-zero coefficients for best best fits
model_info <- data.table(prediction_cv_testfold=this_test_fold, endpoint=this_endpoint, sex=this_sex, samples=survdat[,.N], cases=survdat[,sum(event)])
active_coef <- foreach(mIdx = bestbestfit[,.I], .combine=rbind) %do% {
  this_alpha <- bestbestfit[mIdx, alpha]
  this_lambda <- bestbestfit[mIdx, lambda]
  this_lambda_type <- bestbestfit[mIdx, model]
  this_cv_coxnet <- cv.coxnet.list[[as.character(this_alpha)]]
  active <- coef(this_cv_coxnet, s=this_cv_coxnet[[this_lambda_type]])
  active <- as.matrix(active)
  active <- as.data.table(active, keep.rownames=TRUE)
  setnames(active, c("coef", "beta"))
  active <- active[beta != 0]
  active[, coef_type := ifelse(coef %in% score2_terms, "SCORE2", "NMR")]
  active <- cbind(model_info, alpha=this_alpha, lambda=this_lambda, lambda.fit=this_lambda_type, active)
  return(active)
}
fwrite(active_coef, sep="\t", quote=FALSE, file=sprintf("%s/best_best_fits_coefficients.txt", out_dir))
