library(data.table)
library(foreach)
library(survival)

# Make output directory
system("mkdir -p analyses/CVD_weight_training", wait=TRUE)

# Load discovery data
dat <- fread("data/cleaned/analysis_cohort.txt")

# Load and add predicted NMR scores in the discovery data
nmr_scores <- fread("analyses/nmr_score_training/aggregate_test_non_derived_NMR_scores.txt")
dat <- dat[nmr_scores, on = .(eid)]

# Load FDR-significant clinical chemistry biomarkers
bio_res <- fread("analyses/univariate/cindices_sensitivity_analysis.txt")
bio_res <- bio_res[cohort == "pooled" & sex == "Sex-stratified" & endpoint == "ASCVD" & score == "SCORE2" & model_type == "Clinical biochemistry assay" & deltaC.fdr < 0.05]

# Fit multivariable models
models <- c("NMR scores", "Biochemistry", "PRS", "NMR scores + PRS", "Biochemistry + PRS")
cvd_weights <- foreach(this_sex = c("Male", "Female"), .combine=rbind) %:%
  foreach(this_model = models, .combine=rbind) %do% {
    # Build model formula
    mf <- "Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2_excl_UKB)"
    if (this_model %in% c("NMR scores", "NMR scores + PRS"))
      mf <- paste(mf, "+ scale(CAD_NMR_score) + scale(Stroke_NMR_score)")
    if (this_model %in% c("Biochemistry", "Biochemistry + PRS")) 
      mf <- paste(mf, "+", bio_res[, paste(sprintf("scale(%s)", biomarker), collapse=" + ")])
    if (this_model %in% c("PRS", "NMR scores + PRS", "Biochemistry + PRS"))
      mf <- paste(mf, "+ scale(CAD_metaGRS) + scale(Stroke_metaGRS)")

    # Fit survival model
    cx <- coxph(as.formula(mf), data=dat[sex == this_sex])

    # Extra hazard ratios
    cf <- coef(summary(cx))
    ci <- confint(cx)

    HR <- cf[,2]
    HR.SE <- cf[,3]
    HR.L95 <- exp(ci[,1])
    HR.U95 <- exp(ci[,2])
    HR.pval <- cf[,5]

    # Get mean and sd information on variables so we can apply HRs in replication cohort
    vars <- gsub("scale\\(", "", gsub("\\)", "", names(HR)))
    var_dat <- dat[sex == this_sex][as.integer(rownames(cx$y)), .SD, .SDcols=vars]
    var_means <- apply(var_dat, 2, mean)
    var_sds <- apply(var_dat, 2, sd)

    # Extract relevant info and return
    data.table(
      model = paste("SCORE2 +", this_model), sex = this_sex, samples=cx$n, events=cx$nevent,
      variable_name = NA_real_, variable_col = vars, mean = var_means, sd = var_sds, 
      HR, HR.SE, HR.L95, HR.U95, HR.pval
    )
}

# Add in human friendly variable names
cvd_weights[, variable_name := fcase(
  variable_col == "CAD_NMR_score", "CHD NMR score", 
  variable_col == "Stroke_NMR_score", "Stroke NMR score", 
  variable_col == "CAD_metaGRS", "CHD PRS (PGS000018)", 
  variable_col == "Stroke_metaGRS", "Stroke PRS (PGS000039)",
  variable_col == "alb", "Albumin",
  variable_col == "alt", "ALT",
  variable_col == "alp", "ALP",
  variable_col == "apoa1", "ApoA1",
  variable_col == "apob", "ApoB",
  variable_col == "asp", "AST",
  variable_col == "dbili", "Bilirubin (direct)",
  variable_col == "tbili", "Bilirubin (total)",
  variable_col == "calcium", "Calcium",
  variable_col == "creat", "Creatinine",
  variable_col == "crp", "CRP",
  variable_col == "cyst", "Cystatin-C",
  variable_col == "ggt", "GGT",
  variable_col == "glucose", "Glucose",
  variable_col == "hba1c", "HbA1c",
  variable_col == "igf1", "IGF-1",
  variable_col == "lpa", "Lp(a)",
  variable_col == "ldl", "LDL cholesterol",
  variable_col == "oest", "Oestradiol",
  variable_col == "phos", "Phosphate",
  variable_col == "rheuf", "RF",
  variable_col == "shbg", "SHBG",
  variable_col == "testos", "Testosterone",
  variable_col == "protein", "Total protein",
  variable_col == "trig", "Triglycerides",
  variable_col == "uric", "Urate",
  variable_col == "urea", "Urea",
  variable_col == "vitd25", "Vitamin D"
)]

# Write out
fwrite(cvd_weights, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/multivariable_model_weights.txt")

