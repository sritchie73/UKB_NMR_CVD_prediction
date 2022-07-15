library(data.table)
library(foreach)
source("src/utils/aki_absrisk.R")
source("src/utils/factor_by_size.R")
source("src/utils/risk_recalibration.R")

# Create output directory
out_dir <- "analyses/public_health_modelling/UK_population_generalised/targeted_screening"
system(sprintf("mkdir -p %s", out_dir))

# Load hypothetical population 
ons_pop <- fread("analyses/public_health_modelling/UK_population_generalised/ONS_hypothetical_100k_pop_by_age_sex.txt")
ons_pop_summary <- fread("analyses/public_health_modelling/UK_population_generalised/ONS_hypothetical_100k_pop.txt")

# Load test dataset - we need to refit models and recalibrate risk in subset of 
# participants with complete data for each model. In this instance, this means
# dropping people who are missing conventional risk factors (read: HDL and total
# cholesterol which are dropped for NMR biomarker models) and LDL cholesterol
# which is used to determine statin initation in intermediate risk individuals
test <- fread("data/processed/test/processed_test_data.txt")
test <- test[!is.na(ldl) & !is.na(hdl) & !is.na(tchol)]

# Code factors, using non-risk/lower-risk group as reference
test[, sex := factor(sex, levels=c("Female", "Male"))]
test[, diabetes := factor(diabetes, levels=c("FALSE", "TRUE"))]
test[, smoking := factor(smoking, levels=c("FALSE", "TRUE"))]
test[, family_history_cvd := factor(family_history_cvd, levels=c("FALSE", "TRUE"))]
test[, assessment_centre := factor_by_size(assessment_centre)]
test[, earliest_hospital_nation := factor_by_size(earliest_hospital_nation)]
test[, latest_hospital_nation := factor_by_size(latest_hospital_nation)]

# Reapply filtering used for risk recalibration
test[, age_group := age %/% 5 * 5]
test <- test[age >= 40 & age < 70]
test <- test[(incident_cvd) | incident_followup == 10]

# Load model information
model_info <- fread("analyses/test/model_fit_information.txt")

# Get recalibrated risk for each person and model (and for conventional risk
# factors alone in that subset)
pred <- foreach(midx = model_info[,.I], .combine=rbind) %do% {
  this_model <- model_info[midx]

  # Skip if conventional risk factors alone, doesn't make sense to do targeted re-assessment using
  # same model
  if (this_model$name == "Conventional RF" && !(this_model$PGS)) return(NULL)

  # Fit cox proportional hazards model
  cph <- coxph(as.formula(this_model$formula), data=test, x=TRUE)

  # Extract absolute risk
  pred_risk <- Coxar(cph, 10)

  # Extract individuals for which 10-year risk could be predicted
  dt <- test[as.integer(rownames(cph$x)), .(eid, sex, age, age_group, ldl, diabetes, incident_cvd, pred_risk)]

  # Recalibrate risk to sex-specific 5-year age group incidence rates from CPRD
  dt[, recalibrated_risk := recalibrate_risk(pred_risk, eid, age, sex, male="Male")$recalibrated_risk]

  # Also get the recalibrated risk for conventional risk factors alone in this cohort
  cph <- coxph(as.formula(model_info[name == "Conventional RF" & !(PGS), formula]), data=test[as.integer(rownames(cph$x))], x=TRUE)
  dt[, conv_rf_risk := recalibrate_risk(Coxar(cph, 10), eid, age, sex, male="Male")$recalibrated_risk]

  # No longer need predicted risk (or age)
  dt[, pred_risk := NULL]
  dt[, age := NULL]

  # Add in model information
  this_info <- this_model[,.(name, lambda, PGS, long_name)]
  cbind(this_info, dt)
}

# Stratify risk according to NICE and ACC/AHA guidelines
pred[, ACC.AHA.2019 := fcase(
  recalibrated_risk <= 0.05, "low",
  recalibrated_risk > 0.075, "high", 
  default="medium"
)]

pred[, NICE.2014 := fcase(
  recalibrated_risk <= 0.05, "low",
  recalibrated_risk > 0.10, "high",
  default="medium"
)]

# Also do the same for conventional risk factors
pred[, conv.rf.ACC.AHA.2019 := fcase(
  conv_rf_risk <= 0.05, "low",
  conv_rf_risk > 0.075, "high", 
  default="medium"
)]

pred[, conv.rf.NICE.2014 := fcase(
  conv_rf_risk <= 0.05, "low",
  conv_rf_risk > 0.10, "high",
  default="medium"
)]

# Melt
pred1 <- melt(pred, measure.vars=c("ACC.AHA.2019", "NICE.2014"), variable.name="guidelines", value.name="risk_group")
pred1[, c("conv.rf.ACC.AHA.2019", "conv.rf.NICE.2014") := NULL]
pred1 <- unique(pred1)

pred2 <- melt(pred, measure.vars=c("conv.rf.ACC.AHA.2019", "conv.rf.NICE.2014"), variable.name="guidelines", value.name="risk_group")
pred2[, c("ACC.AHA.2019", "NICE.2014") := NULL]
pred2 <- unique(pred2)
pred2[, guidelines := gsub("conv.rf.", "", guidelines)]

pred1[pred2, on = .(name, lambda, PGS, eid, guidelines), conv_rf_risk_group := i.risk_group]
pred <- pred1

# Get total number of cases and controls in each age and sex group - this is the denominator when computing
# % of cases and controls allocated to risk strata for a given model and thresholds
grp_totals <- pred[guidelines == "NICE.2014", # group totals agnostic to thresholds chosen
  .(.N, cases=sum(incident_cvd), controls=sum(!(incident_cvd))), 
  by=.(sex, age_group, name, lambda, PGS)]
grp_totals <- grp_totals[order(age_group)][order(sex)]

# Build empty table of all possible groups - this allows us to fill in 0s for groups with no participants when computing % allocated
empty <- expand.grid(sex=unique(pred$sex), age_group=unique(pred$age_group), guidelines=unique(pred$guidelines), risk_group=unique(pred$risk_group), stringsAsFactors=FALSE)
setDT(empty)
empty <- rbind(empty[risk_group == "low"], empty[risk_group == "medium"], empty[risk_group == "high"])
empty <- empty[order(age_group)][order(sex)][order(guidelines)]
empty[, c("pct_cases", "pct_controls") := 0]

# Determine % of cases and controls allocated to each risk strata using the blanket approach 
# (conventional risk factors alone, in the subset of individuals with complete data on both
# conventional risk factors, LDL, and the targeted assessment model)
models <- unique(pred[,.(name, lambda, PGS, long_name)])
conv_rf_strata_alloc <- foreach(mIdx = models[,.I], .combine=rbind) %do% {
  cbind(models[mIdx], empty)
}
to_fill <- pred[, .(cases=sum(incident_cvd), controls=sum(!(incident_cvd))),
  by=.(name, lambda, PGS, long_name, sex, age_group, guidelines, risk_group=conv_rf_risk_group)]
to_fill[grp_totals, on = .(sex, age_group, name, lambda, PGS), c("pct_cases", "pct_controls") := .(cases/i.cases, controls/i.controls)]
conv_rf_strata_alloc[to_fill, on = .(name, lambda, PGS, long_name, sex, age_group, guidelines, risk_group),
  c("pct_cases", "pct_controls") := .(i.pct_cases, i.pct_controls)]

# Generalise to hypothetical ONS population
ons_conv_rf <- conv_rf_strata_alloc[ons_pop, on = .(sex, age_group),
  .(name, lambda, PGS, long_name, sex, age_group, guidelines, risk_group,
    N = cases * pct_cases + controls * pct_controls,
    cases = cases * pct_cases, controls = controls * pct_controls)]
ons_conv_rf <- ons_conv_rf[order(guidelines)]

# Summarise to population totals for each risk group
ons_conv_rf_summary <- ons_conv_rf[,
  .(N=sum(cases) + sum(controls), cases=sum(cases), controls=sum(controls)),
  by=.(name, lambda, PGS, long_name, guidelines, risk_group)]

# Write out
fwrite(ons_conv_rf, sep="\t", quote=FALSE, file=sprintf("%s/conv_rf_risk_stratified_by_age_sex.txt", out_dir))
fwrite(ons_conv_rf_summary, sep="\t", quote=FALSE, file=sprintf("%s/conv_rf_risk_stratified.txt", out_dir))

# For people at intermediate risk, determine % of cases and controls that would be allocated treatment due to
# having history of diabetes or elevated LDL-C (>= 5.0 mmol/L)
intermed_treat <- foreach(mIdx = models[,.I], .combine=rbind) %do% {
  this_intermed_treat <- empty[risk_group != "low"]
  this_pred <- pred[models[mIdx], on = .(name, lambda, PGS, long_name)]
  this_grp_totals <- grp_totals[models[mIdx], on = .(name, lambda, PGS)]
  to_fill <- this_pred[conv_rf_risk_group == "medium",
    .(cases=sum(incident_cvd), controls=sum(!(incident_cvd))),
    by=.(name, lambda, PGS, sex, age_group, guidelines,
         risk_group = ifelse(diabetes == FALSE & ldl < 5, "medium", "high"))]
  to_fill[this_grp_totals, on = .(sex, age_group), c("pct_cases", "pct_controls") := .(cases/i.cases, controls/i.controls)]
  this_intermed_treat[to_fill,  on = .(sex, age_group, guidelines, risk_group), c("pct_cases", "pct_controls") := .(i.pct_cases, i.pct_controls)]
  cbind(models[mIdx], this_intermed_treat)
}

# Generalise to hypothetical ONS population
ons_intermed_treat <- intermed_treat[ons_pop, on = .(sex, age_group),
  .(name, lambda, PGS, long_name, sex, age_group, guidelines, risk_group,
    N = cases * pct_cases + controls * pct_controls,
    cases = cases * pct_cases, controls = controls * pct_controls)]
ons_intermed_treat <- ons_intermed_treat[order(guidelines)]

# Summarise to population totals for each risk group
ons_intermed_treat_summary <- ons_intermed_treat[,
  .(N=sum(cases) + sum(controls), cases=sum(cases), controls=sum(controls)),
  by=.(name, lambda, PGS, long_name, guidelines, risk_group)]

# Write out
fwrite(ons_intermed_treat, sep="\t", quote=FALSE, file=sprintf("%s/intermediate_risk_treat_alloc_by_age_sex.txt", out_dir))
fwrite(ons_intermed_treat_summary, sep="\t", quote=FALSE, file=sprintf("%s/intermediate_risk_treat_alloc.txt", out_dir))

# Determine % of cases and controls allocated to each risk strata using the targeted assessment
# subsequent to the above
pred_strata_alloc <- foreach(mIdx = models[,.I], .combine=rbind) %do% {
  cbind(models[mIdx], empty)
}
to_fill <- pred[conv_rf_risk_group == "medium" & diabetes == FALSE & ldl < 5,
  .(cases=sum(incident_cvd), controls=sum(!(incident_cvd))),
  by=.(name, lambda, PGS, long_name, sex, age_group, guidelines, risk_group)]
to_fill[grp_totals, on = .(sex, age_group, name, lambda, PGS), c("pct_cases", "pct_controls") := .(cases/i.cases, controls/i.controls)]
pred_strata_alloc[to_fill, on = .(name, lambda, PGS, long_name, sex, age_group, guidelines, risk_group),
  c("pct_cases", "pct_controls") := .(i.pct_cases, i.pct_controls)]

# Generalise to hypothetical ONS population
ons_pred <- pred_strata_alloc[ons_pop, on = .(sex, age_group),
  .(name, lambda, PGS, long_name, sex, age_group, guidelines, risk_group,
    N = cases * pct_cases + controls * pct_controls,
    cases = cases * pct_cases, controls = controls * pct_controls)]
ons_pred <- ons_pred[order(guidelines)]

# Summarise to population totals for each risk group
ons_pred_summary <- ons_pred[,
  .(N=sum(cases) + sum(controls), cases=sum(cases), controls=sum(controls)),
  by=.(name, lambda, PGS, long_name, guidelines, risk_group)]

# Write out
fwrite(ons_pred, sep="\t", quote=FALSE, file=sprintf("%s/targeted_biomarker_prs_stratified_by_age_sex.txt", out_dir))
fwrite(ons_pred_summary, sep="\t", quote=FALSE, file=sprintf("%s/targeted_biomarker_prs_stratified.txt", out_dir))

# Compute number of cases treated and prevented using the flowchart for each set of models
case_treat_prevent <- copy(models)
case_treat_prevent <- rbind(idcol="guidelines", "ACC.AHA.2019"=case_treat_prevent, "NICE.2014"=case_treat_prevent)
case_treat_prevent[ons_conv_rf_summary[risk_group == "high"], on = .(name, lambda, PGS, guidelines), conv_rf_treat := i.cases]
case_treat_prevent[ons_intermed_treat_summary[risk_group == "high"], on = .(name, lambda, PGS, guidelines), intermed_treat := i.cases]
case_treat_prevent[ons_pred_summary[risk_group == "high"], on = .(name, lambda, PGS, guidelines), pred_treat := i.cases]
case_treat_prevent[is.na(pred_treat), pred_treat := 0]
case_treat_prevent <- melt(case_treat_prevent, measure.vars=c("conv_rf_treat", "intermed_treat", "pred_treat"),
                           variable.name="treatment_reason", value.name="cases_treated")
case_treat_prevent[, treatment_reason := fcase(
  treatment_reason == "conv_rf_treat", "High risk, blanket screening with conventional risk factors", 
  treatment_reason == "intermed_treat", "Medium risk, diabetic or LDL-C >= 5 mmol/L",
  treatment_reason == "pred_treat", "Medium risk reclassified as high risk with targeted assessment of biomarkers/PRS"
)]
case_treat_prevent[, pct_cases_treated := cases_treated / ons_pop_summary$cases]

# Same for controls
control_treat <- copy(models)
control_treat <- rbind(idcol="guidelines", "ACC.AHA.2019"=control_treat, "NICE.2014"=control_treat)
control_treat[ons_conv_rf_summary[risk_group == "high"], on = .(name, lambda, PGS, guidelines), conv_rf_treat := i.controls]
control_treat[ons_intermed_treat_summary[risk_group == "high"], on = .(name, lambda, PGS, guidelines), intermed_treat := i.controls]
control_treat[ons_pred_summary[risk_group == "high"], on = .(name, lambda, PGS, guidelines), pred_treat := i.controls]
control_treat[is.na(pred_treat), pred_treat := 0]
control_treat <- melt(control_treat, measure.vars=c("conv_rf_treat", "intermed_treat", "pred_treat"),
                           variable.name="treatment_reason", value.name="controls_treated")
control_treat[, treatment_reason := fcase(
  treatment_reason == "conv_rf_treat", "High risk, blanket screening with conventional risk factors", 
  treatment_reason == "intermed_treat", "Medium risk, diabetic or LDL-C >= 5 mmol/L",
  treatment_reason == "pred_treat", "Medium risk reclassified as high risk with targeted assessment of biomarkers/PRS"
)]
control_treat[, pct_controls_treated := controls_treated / ons_pop_summary$controls]

# Combine and compute other public health statistics for each treatment reason, model, and threshold guidelines
phs <- case_treat_prevent[control_treat, on = .(name, lambda, PGS, long_name, guidelines, treatment_reason)]
phs[, events_prevented := cases_treated * 0.2] # assuming 20% reduction in risk due to statins
phs[, pct_events_prevented := pct_cases_treated * 0.2]

# Summarise to population level
phs_summary <- phs[, .(
  total_treated=sum(cases_treated) + sum(controls_treated),
  pct_total_treated = (sum(cases_treated) + sum(controls_treated)) / 100000,
  cases_treated=sum(cases_treated), pct_cases_treated=sum(pct_cases_treated),
  controls_treated=sum(controls_treated), pct_controls_treated=sum(pct_controls_treated)),
  by=.(name, lambda, PGS, long_name, guidelines)]
phs_summary[, events_prevented := cases_treated * 0.2]
phs_summary[, NNT := total_treated / events_prevented]
phs_summary[, NNS := 100000 / events_prevented]

# Write out
fwrite(phs, sep="\t", quote=FALSE, file=sprintf("%s/public_health_statistics_by_treatment_reason.txt", out_dir))
fwrite(phs_summary, sep="\t", quote=FALSE, file=sprintf("%s/public_health_statistics.txt", out_dir))

