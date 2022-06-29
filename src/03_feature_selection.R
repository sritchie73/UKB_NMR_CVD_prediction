library(data.table)
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

# Split training data into 10-folds for elasticnet cross-validation
# Balance split by case/control status, sex, prevalent diabetes, smoking status, recruitment centre, and nation
##  train[, foldgrp := paste(incident_cvd, cvd_is_primary_cause, cvd_is_fatal, cvd_primarily_chd, cvd_primarily_stroke,
##                         earliest_hospital_nation, latest_hospital_nation, assessment_centre, sex, diabetes, smoking, family_history_cvd,
##                         is.na(CAD_metaGRS), genetic_white_british, nmr_excess_miss, biochemistry_excess_miss)]

# Balance split by case/control status, type, and sex. Any more factors (e.g. like above)
# and the groups get too small to be meaningfully useful
train[, foldgrp := paste(incident_cvd, cvd_primarily_stroke, cvd_is_fatal, cvd_is_primary_cause, sex)]
train[, foldid := createFolds(foldgrp, k=10, list=FALSE)]
fwrite(train[,.(eid, foldid)], sep="\t", file="analyses/train/training_folds.txt")

# Load biomarker information sheets
bio_info <- fread("data/ukb/biomarkers/output/biomarker_info.txt")
nmr_info <- fread("data/ukb/NMR_metabolomics/output/biomarker_information.txt")

# Curate list of variables that define the different models
conv_rf <- c("age", "sex", "sbp", "diabetes", "smoking", "family_history_cvd")
conv_rf_lipids <- c(conv_rf, c("tchol", "hdl"))
nmr <- nmr_info$Biomarker
covariates <- c("assessment_centre", "earliest_hospital_nation", "latest_hospital_nation")
clin_add <- setdiff(intersect(names(train), bio_info$var), c("tchol", "hdl", "ldl"))

# We want to train two models with lasso regression:
#
# (1) Conventional risk factors + NMR biomarkers, allowing NMR to replace Total and HDL cholesterol
# (2) Conventional risk factors + clinical chemistry biomarkers not measured on the NMR platform, so
#     we can assess potential for improvement e.g. with mass-spec or proteomics technologies as they
#     come down in cost.
#
active <- foreach(model = c("nmr", "assays"), .combine=rbind) %do% {
    # Extract the samples who can be used for glmnet
    if (model == "nmr") {
      this_dat <- train[!(nmr_excess_miss)]
    } else if (model == "assays") {
      this_dat <- train[!(biochemistry_excess_miss) & !is.na(hdl) & !is.na(tchol)]
    }

    # Extract event
    survdat <- this_dat[, .(incident_followup, incident_cvd)]
    setnames(survdat, c("followup", "event"))

    # select predictor columns
    if (model == "nmr") {
      features <- c(covariates, conv_rf, nmr)
    } else if (model == "assays") {
      features <- c(covariates, conv_rf_lipids, clin_add)
    }
    inmat <- model.matrix(~ 0 + ., this_dat[, .SD, .SDcols=features])

    # Create vector of penalties: 1 = apply elasticnet, 0 = always include
    penalties <- rep(1, ncol(inmat))
    names(penalties) <- colnames(inmat)
    if (model == "nmr") {
      penalties[apply(sapply(conv_rf, function(x) { names(penalties) %like% x }), 1, any)] <- 0 
    } else if (model == "assays") {
      penalties[apply(sapply(conv_rf_lipids, function(x) { names(penalties) %like% x }), 1, any)] <- 0
    }

    # Fit cox regression lasso for feature selection
    cv.coxnet <- cv.glmnet(inmat, Surv(survdat$followup, survdat$event), family="cox",
                           alpha = 1, penalty.factor = penalties, foldid = this_dat$foldid,
                           trace.it=1, parallel=TRUE, standardize=FALSE)

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
        coef %in% clin_add, "Assays",
        coef %in% nmr, "NMR",
        default = "Dataset-specific covariate")]
      active <- cbind(data.table(endpoint = "CVD", samples = survdat[,.N], cases = survdat[, sum(event)], model = model, lambda.fit = ss), active)
      return(active)
    }
}
fwrite(active, sep="\t", quote=FALSE, file="analyses/train/cox_lasso_features.txt")

