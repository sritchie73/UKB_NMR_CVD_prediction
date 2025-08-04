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

# Load and add predicted assay scores in the discovery data
assay_scores <- fread("analyses/assay_score_training/aggregate_test_assay_scores.txt")
dat <- dat[assay_scores, on = .(eid)]

# Fit multivariable models
models <- c("NMR scores", "Assay scores", "PRS", "NMR scores + PRS", "Assay scores + PRS", "NMR scores + Assay scores", "NMR scores + Assay scores + PRS")
cvd_weights <- foreach(this_sex = c("Male", "Female"), .combine=rbind) %:%
  foreach(this_model = models, .combine=rbind) %do% {
    # Build model formula
    mf <- "Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2_excl_UKB)"
    if (this_model %like% "NMR scores")
      mf <- paste(mf, "+ scale(CAD_NMR_score) + scale(Stroke_NMR_score)")
    if (this_model %like% "Assay scores") 
      mf <- paste(mf, "+ scale(CAD_assay_score) + scale(Stroke_assay_score)")
    if (this_model %like% "PRS")
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
  variable_col == "CAD_assay_score", "CHD assay score", 
  variable_col == "Stroke_assay_score", "Stroke assay score", 
  variable_col == "CAD_metaGRS", "CHD PRS (PGS000018)", 
  variable_col == "Stroke_metaGRS", "Stroke PRS (PGS000039)"
)]

# Write out
fwrite(cvd_weights, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/multivariable_model_weights_with_assay_scores.txt")

