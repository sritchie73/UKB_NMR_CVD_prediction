library(data.table)
library(ukbnmr)
library(foreach)
library(survival)
library(caret)
library(glmnet)
source("src/utils/factor_by_size.R")

# Make output directory
system("mkdir -p analyses/train/")

# Load preprocessed training dataset
train <- fread("data/processed/training/processed_training_data.txt")

# Code factors, using non-risk/lower-risk group as reference
train[, sex := factor(sex, levels=c("Female", "Male"))]
train[, diabetes := factor(diabetes, levels=c("FALSE", "TRUE"))]
train[, smoking := factor(smoking, levels=c("FALSE", "TRUE"))]
train[, family_history_cvd := factor(family_history_cvd, levels=c("FALSE", "TRUE"))]
train[, assessment_centre := factor_by_size(assessment_centre)]
train[, earliest_hospital_nation := factor_by_size(earliest_hospital_nation)]
train[, latest_hospital_nation := factor_by_size(latest_hospital_nation)]

# Balance split by case/control status, type, and sex. Any more factors
# and the groups get too small to be meaningfully useful
train[, foldgrp := paste(incident_cvd, cvd_primarily_stroke, cvd_is_fatal, cvd_is_primary_cause, sex)]
train[, foldid := createFolds(foldgrp, k=10, list=FALSE)]
fwrite(train[,.(eid, foldid)], sep="\t", file="analyses/train/training_folds.txt")

# Load biomarker information sheets
bio_info <- fread("data/ukb/biomarkers/output/biomarker_info.txt")
nmr_info <- ukbnmr::nmr_info

# Curate list of variables that define the different models
conv_rf <- c("age", "sex", "sbp", "diabetes", "smoking", "family_history_cvd")
conv_rf_lipids <- c(conv_rf, c("tchol", "hdl"))
nmr <- nmr_info$Biomarker
covariates <- c("assessment_centre", "earliest_hospital_nation", "latest_hospital_nation")
clin_add <- setdiff(intersect(names(train), bio_info$var), c("tchol", "hdl", "ldl"))

# Train with lasso regression a model allowing NMR biomarkers to replace Total and HDL cholesterol
# alongside conventional risk factors (which remain fixed in the model, i.e. always selected as 
# important features)
active <- foreach(model = "nmr", .combine=rbind) %do% {
    # Extract event
    survdat <- train[, .(incident_followup, incident_cvd)]
    setnames(survdat, c("followup", "event"))

    # select predictor columns
    features <- c(covariates, conv_rf, nmr)
    inmat <- model.matrix(~ 0 + ., train[, .SD, .SDcols=features])

    # Create vector of penalties: 1 = apply elasticnet, 0 = always include
    penalties <- rep(1, ncol(inmat))
    names(penalties) <- colnames(inmat)
    penalties[apply(sapply(conv_rf, function(x) { names(penalties) %like% x }), 1, any)] <- 0 

    # Fit cox regression lasso for feature selection
    cv.coxnet <- cv.glmnet(inmat, Surv(survdat$followup, survdat$event), family="cox",
                           alpha = 1, penalty.factor = penalties, foldid = train$foldid,
                           trace.it=1, parallel=TRUE, standardize=FALSE)
  
    # Save
    saveRDS(cv.coxnet, file=sprintf("analyses/train/cox_lasso_feature_selection_%s.rds", model))

    # Plot
    pdf(width=12, height=5, file=sprintf("analyses/train/cox_lasso_feature_selection_%s.pdf", model))
    plot(cv.coxnet)
    dev.off()

    # Extract non-zero coefficients
    foreach(ss = c("lambda.min", "lambda.1se"), .combine=rbind) %do% {
      active <- coef(cv.coxnet, s=cv.coxnet[[ss]])
      active <- as.matrix(active)
      active <- as.data.table(active, keep.rownames=TRUE)
      setnames(active, c("coef", "beta"))
      active <- active[beta != 0]
      active[, coef_type := fcase(
        coef %in% names(penalties[penalties == 0]), "conventional",
        coef %in% nmr, "NMR",
        default = "Dataset-specific covariate")]
      active <- cbind(data.table(endpoint = "CVD", samples = survdat[,.N], cases = survdat[, sum(event)], model = model, lambda.fit = ss), active)
      return(active)
    }
}
fwrite(active, sep="\t", quote=FALSE, file="analyses/train/cox_lasso_features.txt")

