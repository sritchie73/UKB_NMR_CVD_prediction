library(data.table)
library(foreach)
library(survival)
source("src/utils/aki_absrisk.R")
source("src/utils/calibration.R")
source("src/utils/risk_recalibration.R")
source("src/utils/factor_by_size.R")

# Make output directory
system("mkdir -p analyses/public_health_modelling/risk_recalibration")

# Load test dataset
test <- fread("data/processed/test/processed_test_data.txt")

# Code factors, using non-risk/lower-risk group as reference
test[, sex := factor(sex, levels=c("Female", "Male"))]
test[, diabetes := factor(diabetes, levels=c("FALSE", "TRUE"))]
test[, smoking := factor(smoking, levels=c("FALSE", "TRUE"))]
test[, family_history_cvd := factor(family_history_cvd, levels=c("FALSE", "TRUE"))]
test[, assessment_centre := factor_by_size(assessment_centre)]
test[, earliest_hospital_nation := factor_by_size(earliest_hospital_nation)]
test[, latest_hospital_nation := factor_by_size(latest_hospital_nation)]

# Load model information
model_info <- fread("analyses/test/model_fit_information.txt")

# Build information about sample exclusions required for analysis
sample_info <- data.table(step="Test dataset", samples=test[,.N], exited=0)

# Risk recalibration works on five year age groups:
#
#    age_group    N
# 1:        40 3078
# 2:        45 3773
# 3:        50 4249
# 4:        55 4625
# 5:        60 5352
# 6:        65 3110
# 7:        70   71
#
# We'll drop people aged >= 70 to prevent this small group skewing the 
# model fit for recalibration
test[, age_group := age %/% 5 * 5]

test <- test[age >= 40]
sample_info <- rbind(sample_info, data.table(step=">= 40 years of age", samples=test[,.N], exited=sample_info[.N, samples] - test[,.N]))

test <- test[age < 70]
sample_info <- rbind(sample_info, data.table(step="< 70 years of age", samples=test[,.N], exited=sample_info[.N, samples] - test[,.N]))

# We also need to drop people censored before 10 years of follow-up for
# non-CVD reasons
test <- test[(incident_cvd) | incident_followup == 10]
sample_info <- rbind(sample_info, data.table(step="Right censored before 10 years for non-CVD reasons", samples=test[,.N], exited=sample_info[.N, samples] - test[,.N]))

# Write out sample information
fwrite(sample_info, sep="\t", quote=FALSE, file="analyses/public_health_modelling/sample_flowchart.txt")

# Get predicted risk, recalibrated risk, and CPRD-matched incidence rates for each person and model
risk <- foreach(midx = model_info[,.I], .combine=rbind) %do% {
  this_model <- model_info[midx]

  # Fit cox proportional hazards model
  cph <- coxph(as.formula(this_model$formula), data=test, x=TRUE)

  # Extract absolute risk
  pred_risk <- Coxar(cph, 10)

  # Extract individuals for which 10-year risk could be predicted
  dt <- test[as.integer(rownames(cph$x)), .(eid, sex, age, age_group, incident_cvd, incident_followup, predicted_risk=pred_risk)]

  # Recalibrate risk to sex-specific 5-year age group incidence rates from CPRD
  dt[, recalibrated_risk := recalibrate_risk(pred_risk, eid, floor(age), sex, male="Male")$recalibrated_risk]

  # Also get incidence rates for each person's age group from CPRD
  dt[, CPRD_incidence := recalibrate_risk(pred_risk, eid, floor(age), sex, male="Male")$CPRD_risk]

  # Add in model information
  this_info <- this_model[,.(name, lambda, PGS, long_name)]
  cbind(this_info, dt)
}

# Write out predicted and recalibrated risks
fwrite(risk, sep="\t", quote=FALSE, file="analyses/public_health_modelling/risk_recalibration/absolute_risks.txt")

# Summarise risk for each age group so we can generate diagnostic plots:
risk_summary <- risk[, .(
  samples=.N, cases=sum(incident_cvd),
  predicted=mean(predicted_risk), predicted.L95=quantile(predicted_risk, 0.025), predicted.U95=quantile(predicted_risk, 0.975),
  recalibrated=mean(recalibrated_risk), recalibrated.L95=quantile(recalibrated_risk, 0.025), recalibrated.U95=quantile(recalibrated_risk, 0.975),
  CPRD=mean(CPRD_incidence), CPRD.L95=quantile(CPRD_incidence, 0.025), CPRD.U95=quantile(CPRD_incidence, 0.975)
), by = .(name, lambda, PGS, long_name, sex, age_group)]
risk_summary[, age_group := sprintf("%s-%s", age_group, age_group+4)]

# Get observed risk for each five-year age group
obs_risk <- foreach(midx = model_info[,.I], .combine=rbind) %do% {
  this_model <- model_info[midx]

  # Get the predicted risk for this model
  this_risk <- risk[this_model, on = .(name, lambda, PGS, long_name)]

  # Extract subset of test dataset with no missing data for all model coefficients 
  # (i.e. participants in the 'risk' dataset for this model)
  this_dat <- test[eid %in% this_risk$eid]
  this_dat[, age_group := NULL]

  # Run model calibration function to extract observed risk per five year age group
  dt <- calibration.fit(this_model$formula, this_dat, 10, this_risk$recalibrated_risk, byage=TRUE)

  # Extract components of interest
  dt <- dt[,.(sex, age_group, observed=obs.risk, observed.L95=obs.L95, observed.U95=obs.U95)]
  
  # Add in model information
  this_info <- this_model[,.(name, lambda, PGS, long_name)]
  cbind(this_info, dt)
}

# Combine information
risk_summary <- obs_risk[risk_summary, on = .(name, lambda, PGS, long_name, sex, age_group)]

# Write out
fwrite(risk_summary, sep="\t", quote=FALSE, file="analyses/public_health_modelling/risk_recalibration/risk_by_five_yr_age.txt")

# Show diagnostic plots of risk recalibration for each model
for (midx in model_info[,.I]) {
  this_model <- model_info[midx]

  ggdt <- risk_summary[this_model, on = .(name, lambda, PGS, long_name)]

  ggdt <- rbind(idcol="metric",
    "Observed risk" = ggdt[,.(age_group, sex, risk=observed, L95=observed.L95, U95=observed.U95)],
    "Predicted risk" = ggdt[,.(age_group, sex, risk=predicted, L95=predicted.L95, U95=predicted.U95)],
    "Recalibrated risk" = ggdt[,.(age_group, sex, risk=recalibrated, L95=recalibrated.L95, U95=recalibrated.U95)],
    "CPRD incidence" = ggdt[,.(age_group, sex, risk=CPRD, L95=CPRD.L95, U95=CPRD.U95)])
  ggdt[, metric := factor(metric, levels=unique(metric))]

  g <- ggplot(ggdt) +
    aes(x=factor(age_group), y=risk, ymin=L95, ymax=U95, color=metric, group=metric) +
    geom_errorbar(width=0, alpha=0.5, position=position_dodge(width=0.4)) +
    geom_point(shape = 19, position=position_dodge(width=0.4)) +
    scale_color_manual(name="", values=c("Observed risk"="#fdbf6f", "Predicted risk"="#fb9a99", "Recalibrated risk"="#e31a1c", "CPRD incidence"="black")) +
    xlab("5 year age group") +
    ylab("10-year risk (95% CI)") +
    facet_wrap(~ sex, scales="free") +
    theme_bw() + theme(legend.position="bottom")

  ggsave(g, width=8, height=4, file=sprintf("analyses/public_health_modelling/risk_recalibration/%s%s%s.pdf",
    gsub(" \\+? ?", "_", this_model$name), ifelse(this_model$PGS, "_with_PGS", ""),
    ifelse(this_model$lambda == "", "", paste0("_", this_model$lambda))))
}

