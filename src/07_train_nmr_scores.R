library(data.table)
library(foreach)
library(caret)
library(survival)
library(glmnet)
registerDoMC(10) # recommend requesting more cores than this in sbatch for memory reasons

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

# Load NMR info
nmr_info <- fread("data/ukb/NMR_metabolomics/biomarker_information.txt")

# Extract subset of training data
dat <- dat[sex == this_sex & cvd_prediction_foldid != this_test_fold]

# Extract columns required for model training
if (this_endpoint == "CAD") {
  dat <- dat[!(prevalent_cad) | is.na(prevalent_cad), .(eid, age, SCORE2_excl_UKB, event=incident_cad, followup=incident_cad_followup)]
} else if (this_endpoint == "Stroke") {
  dat <- dat[!(prevalent_stroke) | is.na(prevalent_stroke), .(eid, age, SCORE2_excl_UKB, event=incident_stroke, followup=incident_stroke_followup)]
}

# Allocate samples to 10-folds for cross-validation, balancing folds by cases status
dat[, elasticnet_cv_foldid := createFolds(event, k=10, list=FALSE)]

# Write out foldid allocation
fwrite(dat[,.(eid, elasticnet_cv_foldid)], sep="\t", quote=FALSE, file=sprintf("%s/cross_validation_fold_allocation.txt", out_dir))

# Add in standardized NMR biomarker terms
nmr_scaled <- fread("data/standardised/nmr_concentrations.txt")
dat <- dat[nmr_scaled, on = .(eid), nomatch=0]

# Convert age into 5-year age group centered at 60 years to match SCORE2 interaction terms
dat[, age := (age - 60)/5]

# Setup list of alpha mixing parameters to search across (controls balance of ridge vs. lasso)
alphas <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1) # 0 = ridge, 1 = lasso

# We're going to train two different scores: 
# (1) Using the subset of 21 clinically accredited biomarkers (excluding 3 composite and 9 ratios, and HDL, LDL, and total cholesterol)
# (2) Using the subset of 106 non-derived biomarkers

training_runs <- foreach(this_type = c("non-derived", "clinical")) %do% {
  # Get list of candidate biomarkers
  if (this_type == "non-derived") {
    biomarkers <- nmr_info[Type == "Non-derived" & Biomarker != "Clinical_LDL_C", Biomarker]
  } else if (this_type == "clinical") {
    biomarkers <- c("VLDL_C", "Total_TG", "ApoB", "ApoA1", 
      "Omega_3", "Omega_6", "MUFA", "SFA", "DHA", "Ala", "Gly", "His", "Ile", "Leu", "Val", "Phe", "Tyr",
      "Glucose", "Creatinine", "Albumin", "GlycA")
  }

  # Set up model matrix formula
  mf <- as.formula(sprintf("~0 + %s", 
    paste(sprintf("age*%s", biomarkers), collapse=" + ")
  ))

  # Create model matrix of predictor terms
  inmat <- model.matrix(mf, dat)

  # Drop the age column created by the model matrix (due to the age interactions) -
  # age is captured in the SCORE2 linear predictor offset
  inmat <- inmat[, -which(colnames(inmat) == "age")]

  # Run elastinet with the given alpha
  cv.coxnet.list <- foreach(this_alpha = alphas, .inorder=TRUE) %do% {
    # Fit elastinet in 10-fold cross-validation with the given alpha penalty
    cv.coxnet <- cv.glmnet(inmat, Surv(dat$followup, dat$event), family="cox",
                           alpha = this_alpha, offset=dat$SCORE2_excl_UKB,
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
  return(cv.coxnet.list)
}
names(training_runs) <- c("non-derived", "clinical")

# Save
saveRDS(training_runs, file=sprintf("%s/training_runs.rds", out_dir))

# Curate model information
model_info <- data.table(prediction_cv_testfold=this_test_fold, endpoint=this_endpoint, sex=this_sex, samples=dat[,.N], cases=dat[,sum(event)])

# build table of measure of best fits (partial likelihood deviance) for each lambda and alpha
fitdt <- foreach(this_type = c("non-derived", "clinical"), .combine=rbind) %:% 
  foreach(idx = seq_along(alphas), .combine=rbind) %do% {
    data.table(
      type=this_type, alpha=alphas[idx], 
      lambda=training_runs[[this_type]][[idx]][["lambda"]],
      fit_metric=training_runs[[this_type]][[idx]][["name"]],
      fit_mean=training_runs[[this_type]][[idx]][["cvm"]], 
      fit_sd=training_runs[[this_type]][[idx]][["cvsd"]],
      fit_mean_minus_sd=training_runs[[this_type]][[idx]][["cvlo"]], 
      fit_mean_plus_sd=training_runs[[this_type]][[idx]][["cvup"]],
      nonzero=training_runs[[this_type]][[idx]][["nzero"]],
      lambda.min=training_runs[[this_type]][[idx]][["lambda.min"]],
      lambda.1se=training_runs[[this_type]][[idx]][["lambda.1se"]])
}
fitdt <- cbind(model_info, fitdt)
fwrite(fitdt, sep="\t", quote=FALSE, file=sprintf("%s/cv_coxnet_all_fits.txt", out_dir))

# Build table containing information about the best fits in each case:
bestfit <- rbind(idcol="lambda.fit",
  "lambda.min"=fitdt[lambda == lambda.min],
  "lambda.1se"=fitdt[lambda == lambda.1se]
) 
bestfit[, c("lambda.min", "lambda.1se") := NULL]
fwrite(bestfit, sep="\t", quote=FALSE, file=sprintf("%s/cv_coxnet_best_fits.txt", out_dir))

# Choose best lambda.min and lambda.1se across models
bestbestfit <- bestfit[,.SD[which.min(fit_mean)],by=.(type, lambda.fit)]
fwrite(bestbestfit,  sep="\t", quote=FALSE, file=sprintf("%s/cv_coxnet_best_best_fits.txt", out_dir))

# Extract non-zero coefficients for best best fits
active_coef <- foreach(mIdx = bestbestfit[,.I], .combine=rbind) %do% {
  this_score_type <- bestbestfit[mIdx, type]
  this_alpha <- bestbestfit[mIdx, alpha]
  this_lambda <- bestbestfit[mIdx, lambda]
  this_lambda_type <- bestbestfit[mIdx, lambda.fit]
  this_cv_coxnet <- training_runs[[this_score_type]][[as.character(this_alpha)]]
  active <- coef(this_cv_coxnet, s=this_cv_coxnet[[this_lambda_type]])
  active <- as.matrix(active)
  active <- as.data.table(active, keep.rownames=TRUE)
  setnames(active, c("coef", "beta"))
  active <- active[beta != 0]
  active <- cbind(model_info, type=this_score_type, alpha=this_alpha, lambda=this_lambda, lambda.fit=this_lambda_type, active)
  return(active)
}
fwrite(active_coef, sep="\t", quote=FALSE, file=sprintf("%s/best_best_fits_coefficients.txt", out_dir))
