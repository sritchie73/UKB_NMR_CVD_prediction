library(data.table)
library(survival)
library(foreach)

# First load in and subset the full UKB cohort to those where PRS can be modelled
dat <- fread('data/cleaned/full_UKB_analysis_cohort.txt')
dat <- dat[!(no_genetics) & !(prs_training_samples)]

# Standardise PRS
dat[, CAD_metaGRS := scale(CAD_metaGRS)]
dat[, Stroke_metaGRS := scale(Stroke_metaGRS)]

# Assess association in sex-specific analysis
prs_sex_assocs <- foreach(this_sex = c("Male", "Female"), .combine=rbind) %do% {
  this_dat <- dat[sex == this_sex]

  mf <- Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2_excl_UKB) + CAD_metaGRS + Stroke_metaGRS # formula
  cx <- coxph(mf, data=this_dat) # Fit cox proportional hazards model
  ci <- confint(cx) # get 95% confidence intervals

  # collate information about all model coefficients
  cxdt <- as.data.table(coef(summary(cx)), keep.rownames="coefficient")
  cidt <- as.data.table(ci, keep.rownames="coefficient")
  dt <- cxdt[cidt, on = .(coefficient)]
  
  dt[, .(Sex=this_sex, Samples=cx$n, Cases=cx$nevent, coefficient, HR=`exp(coef)`, L95=exp(`2.5 %`), U95=exp(`97.5 %`), Pvalue=`Pr(>|z|)`)]
}

fwrite(prs_sex_assocs, sep="\t", quote=FALSE, file="analyses/test/full_UKB_sex_specific_prs_assocs.txt")

# Now look at sex interaction terms

# Code sex factor with females as reference
dat[, sex := factor(sex, levels=c("Female", "Male"))]

cx <- coxph(Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2_excl_UKB) + sex*CAD_metaGRS + sex*Stroke_metaGRS, data=dat)
ci <- confint(cx) # get 95% confidence intervals

# collate information about all model coefficients
cxdt <- as.data.table(coef(summary(cx)), keep.rownames="coefficient")
cidt <- as.data.table(ci, keep.rownames="coefficient")
prs_sex_interactions <- cxdt[cidt, on = .(coefficient)]
prs_sex_interactions <- prs_sex_interactions[, .(Samples=cx$n, Cases=cx$nevent, coefficient, HR=`exp(coef)`, L95=exp(`2.5 %`), U95=exp(`97.5 %`), Pvalue=`Pr(>|z|)`)]
fwrite(prs_sex_interactions, sep="\t", quote=FALSE, file="analyses/test/full_UKB_prs_sex_interactions.txt")

