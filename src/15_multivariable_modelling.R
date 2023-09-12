library(data.table)
library(foreach)
library(survival)
library(ggplot2)
library(cowplot)

# Make output directory
system("mkdir -p analyses/test", wait=TRUE)

# Load required data
dat <- fread("data/cleaned/analysis_cohort.txt", select=c("eid", "sex", "age", "incident_cvd_followup", "incident_cvd",  "smoking", "sbp", "hdl", "tchol", "SCORE2_excl_UKB", "CAD_metaGRS", "Stroke_metaGRS"))
setnames(dat, "SCORE2_excl_UKB", "SCORE2")

# Add in test scores
test_scores <- fread("analyses/nmr_score_training/aggregate_test_non_derived_NMR_scores.txt")
setnames(test_scores, gsub("NMR", "NMR", names(test_scores)))
dat <- dat[test_scores, on = .(eid)]

# Code smoking status factor
dat[, smoking := fcase(
  is.na(smoking), "other",
  smoking == FALSE, "other",
  smoking == TRUE, "current"
)]
dat[, smoking := factor(smoking, levels=c("other", "current"))]

# Fit multivariable model with SCORE2 as an offset
offset_hrs <- foreach(this_sex = c("Sex-stratified", "Males", "Females"), .combine=rbind) %do% {
  if (this_sex == "Sex-stratified") {
    cx <- coxph(data=dat,
      Surv(incident_cvd_followup, incident_cvd) ~ strata(sex) + offset(SCORE2) +
        scale(CAD_NMR_score) + scale(Stroke_NMR_score) + scale(CAD_metaGRS) + scale(Stroke_metaGRS)
    )
  } else {
    if (this_sex == "Males") {
      this_dat <- dat[sex == "Male"]
    } else {
      this_dat <- dat[sex == "Female"]
    }
    cx <- coxph(data=this_dat,
      Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) +
        scale(CAD_NMR_score) + scale(Stroke_NMR_score) + scale(CAD_metaGRS) + scale(Stroke_metaGRS)
    )
  }

  this_hrs <- as.data.table(coef(summary(cx)), keep.rownames=TRUE)
  setnames(this_hrs, c("variable", "logHR", "HR", "SE", "Z", "pval"))
  this_hrs[, variable := c("CAD NMR score", "Stroke NMR score", "CAD PRS", "Stroke PRS")]
  this_hrs[, c("L95", "U95") :=  as.data.table(confint(cx))]
  this_hrs[, .(sex=this_sex, variable, HR, L95, U95, pval)]
}
fwrite(offset_hrs, sep="\t", quote=FALSE, file="analyses/test/score_HRs_independent_of_score2.txt")

# Fit multivariable models with SCORE2 component risk factors
multi_hrs <- foreach(this_sex = c("Sex-stratified", "Males", "Females"), .combine=rbind) %do% {
  if (this_sex == "Sex-stratified") {
    cx <- coxph(data=dat,
      Surv(incident_cvd_followup, incident_cvd) ~ strata(sex) + 
        scale(age)*smoking + scale(age)*scale(sbp) + scale(age)*scale(tchol) + scale(age)*scale(hdl) +
        scale(CAD_NMR_score) + scale(Stroke_NMR_score) + scale(CAD_metaGRS) + scale(Stroke_metaGRS)
    )
  } else {
    if (this_sex == "Males") {
      this_dat <- dat[sex == "Male"]
    } else {
      this_dat <- dat[sex == "Female"]
    }
    cx <- coxph(data=this_dat,
      Surv(incident_cvd_followup, incident_cvd) ~
        scale(age)*smoking + scale(age)*scale(sbp) + scale(age)*scale(tchol) + scale(age)*scale(hdl) +
        scale(CAD_NMR_score) + scale(Stroke_NMR_score) + scale(CAD_metaGRS) + scale(Stroke_metaGRS)
    )
  }

  this_hrs <- as.data.table(coef(summary(cx)), keep.rownames=TRUE)
  setnames(this_hrs, c("variable", "logHR", "HR", "SE", "Z", "pval"))
  this_hrs[, variable := c(
    "Age", "Smoker", "SBP", "Total cholesterol", "HDL cholesterol",
    "CAD NMR score", "Stroke NMR score", "CAD PRS", "Stroke PRS",
    "Age x smoker", "Age x SBP", "Age x total cholesterol", "Age x HDL cholesterol"
  )]
  this_hrs[, c("L95", "U95") :=  as.data.table(confint(cx))]
  this_hrs[, .(sex=this_sex, variable, HR, L95, U95, pval)]
}
fwrite(multi_hrs, sep="\t", quote=FALSE, file="analyses/test/score_HRs_independent_of_score2_risk_factors.txt")

