library(data.table)
library(foreach)
library(caret)
library(survival)
library(glmnet)
registerDoMC(10) # recommend requesting more cores than this in sbatch for memory reasons
source("src/utils/SCORE2.R")

# Setup array task information
tasklist <- expand.grid(prediction_cv_testfold = 1:5, endpoint = c("CAD", "Stroke"), sex=c("Male", "Female")) # 20 tasks
setDT(tasklist)
this_task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
this_task <- tasklist[this_task]

this_test_fold <- this_task$prediction_cv_testfold
this_endpoint <- this_task$endpoint
this_sex <- this_task$sex
message(sprintf("Training NMR score for %s in %ss in prediction test fold %s\n", this_endpoint, this_sex, this_test_fold))

# Make output directory
out_dir <- sprintf("analyses/nmr_score_training/test_fold_%s/%s/%s", this_test_fold, this_endpoint, this_sex)
system(sprintf("mkdir -p %s", out_dir))

# Load imputed data
dat <- fread("data/imputed/analysis_cohort.txt")

# Extract subset of training data
dat <- dat[sex == this_sex & prediction_cv_foldid != this_test_fold]

# Compute SCORE2 linear predictor to be held constant in model training
dat[, SCORE2_LP := score2(sex, age, smoking, sbp, tchol, hdl, type="linear predictor")]

# Extract columns required for model training
if (this_endpoint == "CAD") {
  dat <- dat[!(prevalent_cad) | is.na(prevalent_cad), .(eid, age, SCORE2_LP, event=incident_cad, followup=incident_cad_followup)]
} else if (this_endpoint == "Stroke") {
  dat <- dat[!(prevalent_stroke) | is.na(prevalent_stroke), .(eid, age, SCORE2_LP, event=incident_stroke, followup=incident_stroke_followup)]
}

# Allocate samples to 10-folds for cross-validation, balancing folds by cases status
dat[, elasticnet_cv_foldid := createFolds(event, k=10, list=FALSE)]

# Write out foldid allocation
fwrite(dat[,.(eid, elasticnet_cv_foldid)], sep="\t", quote=FALSE, file=sprintf("%s/cross_validation_fold_allocation.txt", out_dir))

# Add in standardized NMR biomarker terms
nmr_scaled <- fread("data/standardised/non_derived_nmr.txt")
dat <- dat[nmr_scaled, on = .(eid), nomatch=0]

# Convert age into 5-year age group centered at 60 years to match SCORE2 interaction terms
dat[, age := (age - 60)/5]

# Set up model matrix formula
mf <- as.formula(sprintf("~0 + %s", 
  paste(sprintf("age*%s", ukbnmr::nmr_info[Type == "Non-derived", Biomarker]), collapse=" + ")
))

# Create model matrix of predictor terms
inmat <- model.matrix(mf, dat)

# Drop the age column created by the model matrix (due to the age interactions) -
# age is captured in the SCORE2 linear predictor offset
inmat <- inmat[, -which(colnames(inmat) == "age")]

# Setup list of alpha mixing parameters to search across (controls balance of ridge vs. lasso)
alphas <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1) # 0 = ridge, 1 = lasso

# Run elastinet with the given alpha
cv.coxnet.list <- foreach(this_alpha = alphas, .inorder=TRUE) %do% {
  # Fit elastinet in 10-fold cross-validation with the given alpha penalty
  cv.coxnet <- cv.glmnet(inmat, Surv(dat$followup, dat$event), family="cox",
                         alpha = this_alpha, offset=dat$SCORE2_LP,
                         foldid = dat$elasticnet_cv_foldid,
                         trace.it=1, parallel=TRUE, standardize=FALSE)

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

# Curate model information
model_info <- data.table(prediction_cv_testfold=this_test_fold, endpoint=this_endpoint, sex=this_sex, samples=dat[,.N], cases=dat[,sum(event)])

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
fitdt <- cbind(model_info, fitdt)
fwrite(fitdt, sep="\t", quote=FALSE, file=sprintf("%s/cv_coxnet_list_all_fits.txt", out_dir))

# Build table containing information about the best fits in each case:
bestfit <- fitdt[lambda == lambda.min | lambda == lambda.1se]
bestfit[lambda == lambda.min, model := "lambda.min"]
bestfit[lambda == lambda.1se, model := "lambda.1se"]
bestfit <- bestfit[, .(alpha, lambda, model, nonzero, fit_metric, fit_mean, fit_sd, fit_mean_minus_sd, fit_mean_plus_sd)]
bestfit <- cbind(model_info, bestfit)
fwrite(bestfit, sep="\t", quote=FALSE, file=sprintf("%s/cv_coxnet_best_fits.txt", out_dir))

# Choose best lambda.min and lambda.1se across models
bestbestfit <- bestfit[,.SD[which.min(fit_mean)],by=model]
bestbestfit <- cbind(model_info, bestbestfit)
fwrite(bestbestfit,  sep="\t", quote=FALSE, file=sprintf("%s/cv_coxnet_best_best_fits.txt", out_dir))

# Extract non-zero coefficients for best best fits
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
  active <- cbind(model_info, alpha=this_alpha, lambda=this_lambda, lambda.fit=this_lambda_type, active)
  return(active)
}
fwrite(active_coef, sep="\t", quote=FALSE, file=sprintf("%s/best_best_fits_coefficients.txt", out_dir))
