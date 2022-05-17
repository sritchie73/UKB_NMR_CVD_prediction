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
train[, foldgrp := paste(incident_cvd, cvd_is_primary_cause, cvd_is_fatal, cvd_primarily_chd, cvd_primarily_stroke,
                         earliest_hospital_nation, latest_hospital_nation, assessment_centre, sex, diabetes, smoking, family_history_cvd)]
train[, foldid := createFolds(foldgrp, k=10, list=FALSE)]
fwrite(train[,.(eid, foldid)], sep="\t", file="analyses/train/training_folds.txt")

# Load biomarker information sheets
bio_info <- fread("data/ukb/biomarkers/output/biomarker_info.txt")
nmr_info <- fread("data/ukb/NMR_metabolomics/output/biomarker_information.txt")

# Curate list of variables that define the different models
standard_rf <- c("age", "sex", "tchol", "hdl", "sbp", "diabetes", "smoking", "family_history_cvd")
covariates <- c("assessment_centre", "earliest_hospital_nation", "latest_hospital_nation")
night <- nmr_info$Biomarker
blood_clin <- setdiff(intersect(bio_info$var, names(train)), standard_rf)
pgs <- c("CAD_metaGRS", "Stroke_metaGRS")

# Run lasso regression to identify (1) clinical chemistry biomarkers that significantly add
# to CVD prediction above and beyond standard clinical risk factors, (2) nightingale 
# biomarkers that add to CVD prediction above and beyond clinical risk factors, and 
# (3) nightingale and clinical chemistry biomarkers in conjunction, and (4-6) the above also 
# with CAD and Stroke PGS added.
active <- foreach(with_pgs = c(TRUE, FALSE), .combine=rbind) %:%
  foreach(model = c("blood", "nightingale", "blood+nightingale"), .combine=rbind) %do% {
    # Extract event
    survdat <- train[, .(incident_followup, incident_cvd)]
    setnames(survdat, c("followup", "event"))

    # select columns
    features <- c(covariates, standard_rf)
    if (with_pgs) features <- c(features, pgs)
    if (grepl("blood", model)) features <- c(features, blood_clin)
    if (grepl("nightingale", model)) features <- c(features, night)
    inmat <- model.matrix(~ 0 + ., train[, .SD, .SDcols=features])

    # Create vector of penalties: 1 = apply elasticnet, 0 = always include
    penalties <- rep(0, ncol(inmat))
    names(penalties) <- colnames(inmat)
    penalties[intersect(c(blood_clin, night), names(penalties))] <- 1 # always apply elasticnet to non-conventional risk factor biomarkers
    for (cov in covariates) {
      penalties[names(penalties) %like% cov] <- 1 # and too follow-up time related covaraties
    }

    # Fit cox regression lasso for feature selection
    cv.coxnet <- cv.glmnet(inmat, Surv(survdat$followup, survdat$event), family="cox",
                           foldid = train$foldid, alpha = 1, penalty.factor = penalties,
                           trace.it=1, parallel=TRUE, standardize=FALSE)

    # Plot
    pdf(width=12, height=5, file=sprintf("analyses/train/cox_lasso_feature_selection_%s%s.pdf",
        model, ifelse(with_pgs, paste0("_", "PGS"), "")))
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
        coef %in% pgs, "pgs",
        coef %in% names(penalties[penalties == 0]), "conventional",
        coef %in% blood_clin, "blood",
        coef %in% nmr_info$Biomarker, "nightingale",
        default = "Dataset-specific covariate")]
      active <- cbind(data.table(endpoint = "CVD", model = model, PGS = with_pgs, lambda.fit = ss), active)
      return(active)
    }
}
fwrite(active, sep="\t", quote=FALSE, file="analyses/train/cox_lasso_features.txt")

