library(openxlsx)
library(data.table)
library(foreach)
library(ggplot2)

# Create output directory
out_dir <- "analyses/public_health_modelling/UK_population_generalised/targeted_screening"
system(sprintf("mkdir -p %s", out_dir))

# Load hypothetical population
ons_pop <- fread("analyses/public_health_modelling/UK_population_generalised/ONS_hypothetical_100k_pop_by_age_sex.txt")
ons_pop_summary <- fread("analyses/public_health_modelling/UK_population_generalised/ONS_hypothetical_100k_pop.txt")

# Load in predicted risk levels for all models
pred <- fread("analyses/public_health_modelling/risk_recalibration/absolute_risks.txt")

# Drop columns we won't be using
pred[, age := NULL] # five year age group already available as 'age_group'
pred[, incident_followup := NULL]
pred[, predicted_risk := NULL]
pred[, CPRD_incidence := NULL]

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

# Melt
pred <- melt(pred, measure.vars=c("ACC.AHA.2019", "NICE.2014"), variable.name="guidelines", value.name="risk_group")

# Add in LDL-C and diabetes status
test <- fread("data/cleaned/test_data.txt")
pred[test, on = .(eid), diabetes := i.diabetes]
pred[test, on = .(eid), ldl := i.ldl] # from clinical chemistry

# Extract recalibrated predicted risk for conventional risk factors alone vs. other models
conv_rf <- pred[name == "Conventional RF" & !(PGS)]
pred <- pred[name != "Conventional RF" | (PGS)]
pred[conv_rf, on = .(eid, guidelines), conv_rf_risk_group := i.risk_group]

# Get total number of cases and controls in each age and sex group - this is the denominator when computing
# % of cases and controls allocated to risk strata for a given model and thresholds
grp_totals <- conv_rf[guidelines == "NICE.2014", # group totals identical regardless of thresholds
  .(.N, cases=sum(incident_cvd), controls=sum(!(incident_cvd))), 
  by=.(sex, age_group)]
grp_totals <- grp_totals[order(age_group)][order(sex)]

# Build empty table of all possible groups - this allows us to fill in 0s for groups with no participants when computing % allocated
empty <- expand.grid(sex=unique(pred$sex), age_group=unique(pred$age_group), guidelines=unique(pred$guidelines), risk_group=unique(pred$risk_group), stringsAsFactors=FALSE)
setDT(empty)
empty <- rbind(empty[risk_group == "low"], empty[risk_group == "medium"], empty[risk_group == "high"])
empty <- empty[order(age_group)][order(sex)][order(guidelines)]
empty[, c("pct_cases", "pct_controls") := 0]

# In each sex and age-group, determine the % of cases and controls allocated to each risk strata when using conventional risk
# factors alone
conv_rf_strata_alloc <- copy(empty)
to_fill <- conv_rf[, .(cases=sum(incident_cvd), controls=sum(!(incident_cvd))), 
 by=.(sex, age_group, guidelines, risk_group)]
to_fill[grp_totals, on = .(sex, age_group),  c("pct_cases", "pct_controls") := .(cases/i.cases, controls/i.controls)]
conv_rf_strata_alloc[to_fill, on = .(sex, age_group, guidelines, risk_group), c("pct_cases", "pct_controls") := .(i.pct_cases, i.pct_controls)]

# Now apply to hypothetical ONS population 
ons_conv_rf <- conv_rf_strata_alloc[ons_pop, on = .(sex, age_group),
  .(sex, age_group, guidelines, risk_group, 
    N = cases * pct_cases + controls * pct_controls,
    cases = cases * pct_cases, controls = controls * pct_controls)]
ons_conv_rf <- ons_conv_rf[order(guidelines)]

# Summarise to population totals for each risk group
ons_conv_rf_summary <- ons_conv_rf[, 
  .(N=sum(cases) + sum(controls), cases=sum(cases), controls=sum(controls)), 
  by=.(guidelines, risk_group)]

# Write out
fwrite(ons_conv_rf, sep="\t", quote=FALSE, file=sprintf("%s/conventional_risk_stratified_by_age_sex.txt", out_dir))
fwrite(ons_conv_rf_summary, sep="\t", quote=FALSE, file=sprintf("%s/conventional_risk_stratified.txt", out_dir))

# For people at intermediate risk, determine % of cases and controls that would be allocated treatment due to 
# having history of diabetes or elevated LDL-C (>= 5.0 mmol/L)
intermed_treat <- empty[risk_group != "low"]
to_fill <- conv_rf[risk_group == "medium", 
  .(cases=sum(incident_cvd), controls=sum(!(incident_cvd))),
  by=.(sex, age_group, guidelines, 
       risk_group = ifelse(!(diabetes) & (is.na(ldl) | ldl < 5), "medium", "high"))]
to_fill[grp_totals, on = .(sex, age_group), c("pct_cases", "pct_controls") := .(cases/i.cases, controls/i.controls)]
intermed_treat[to_fill,  on = .(sex, age_group, guidelines, risk_group), c("pct_cases", "pct_controls") := .(i.pct_cases, i.pct_controls)]

# Generalise to hypothetical ONS population
ons_intermed_treat <- intermed_treat[ons_pop, on = .(sex, age_group),
  .(sex, age_group, guidelines, risk_group,
    N = cases * pct_cases + controls * pct_controls,
    cases = cases * pct_cases, controls = controls * pct_controls)]
ons_intermed_treat <- ons_intermed_treat[order(guidelines)]

# Summarise to population totals for each risk group
ons_intermed_treat_summary <- ons_intermed_treat[,
  .(N=sum(cases) + sum(controls), cases=sum(cases), controls=sum(controls)), 
  by=.(guidelines, risk_group)]

# Write out
fwrite(ons_intermed_treat, sep="\t", quote=FALSE, file=sprintf("%s/intermediate_risk_treat_alloc_by_age_sex.txt", out_dir))
fwrite(ons_intermed_treat_summary, sep="\t", quote=FALSE, file=sprintf("%s/intermediate_risk_treat_alloc.txt", out_dir))

# For those at medium risk with conventional risk factors alone and not allocated treatment due to history of diabetes 
# or elevated LDL-C, determine % of cases and controls that would be reclassified to each risk strata by adding PRS
# and/or biomarkers 
models <- unique(pred[,.(name, lambda, PGS, long_name)])
pred_strata_alloc <- foreach(mIdx = models[,.I], .combine=rbind) %do% {
  cbind(models[mIdx], empty)
}
to_fill <- pred[conv_rf_risk_group == "medium" & !(diabetes) & (ldl < 5 | is.na(ldl)),
  .(cases=sum(incident_cvd), controls=sum(!(incident_cvd))),
  by=.(name, lambda, PGS, long_name, sex, age_group, guidelines, risk_group)]

# People with missing data in the relevant biomarkers are kept as the medium risk
missing <- foreach(mIdx = models[,.I], .combine=rbind) %do% {
  this_pred <- pred[models[mIdx], on = .(name, lambda, PGS, long_name)]
  this_pred_eid <- unique(this_pred$eid)
  this_pred_missing <- conv_rf[!(eid %in% this_pred_eid)]
  this_pred_missing <- this_pred_missing[risk_group == "medium" & !(diabetes) & (ldl < 5 | is.na(ldl))]
	this_pred_missing <- this_pred_missing[, .(cases=sum(incident_cvd), controls=sum(!(incident_cvd))), 
		by=.(sex, age_group, guidelines)]
  if (nrow(this_pred_missing) > 0) {
    return(cbind(models[mIdx], this_pred_missing))
  }
}
missing[, risk_group := "medium"]
to_fill[missing, on = .(name, lambda, PGS, long_name, sex, age_group, guidelines, risk_group), 
  c("cases", "controls") := .(cases + i.cases, controls + i.controls)]

to_fill[grp_totals, on = .(sex, age_group), c("pct_cases", "pct_controls") := .(cases/i.cases, controls/i.controls)]
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
fwrite(ons_pred, sep="\t", quote=FALSE, file=sprintf("%s/biomarker_prs_risk_stratified_by_age_sex.txt", out_dir))
fwrite(ons_pred_summary, sep="\t", quote=FALSE, file=sprintf("%s/biomarker_prs_risk_stratified_by_age_sex.txt", out_dir))

# Compute number of cases treated and prevented using the flowchart for each set of models
case_treat_prevent <- rbind(models[1], models)
case_treat_prevent[1, c("PGS", "long_name") := .(FALSE, "Conventional RF")]
case_treat_prevent <- rbind(idcol="guidelines", "ACC.AHA.2019"=case_treat_prevent, "NICE.2014"=case_treat_prevent)
case_treat_prevent[ons_conv_rf_summary[risk_group == "high"], on = .(guidelines), conv_rf := i.cases]
case_treat_prevent[ons_intermed_treat_summary[risk_group == "high"], on = .(guidelines), intermed_treat := i.cases]
case_treat_prevent[ons_pred_summary[risk_group == "high"], on = .(name, lambda, PGS, guidelines), pred_treat := i.cases]
case_treat_prevent[is.na(pred_treat), pred_treat := 0] 
case_treat_prevent <- melt(case_treat_prevent, measure.vars=c("conv_rf", "intermed_treat", "pred_treat"),
                           variable.name="treatment_reason", value.name="cases_treated")
case_treat_prevent[, treatment_reason := fcase(
  treatment_reason == "conv_rf", "High risk due to conventional risk factors",
  treatment_reason == "intermed_treat", "Medium risk, diabetic or LDL-C >= 5 mmol/L",
  treatment_reason == "pred_treat", "Reclassified as high risk from biomarkers/PRS"
)]
case_treat_prevent[, pct_cases_treated := cases_treated / ons_pop_summary$cases]

# Same for controls
control_treat <- rbind(models[1], models)
control_treat[1, c("PGS", "long_name") := .(FALSE, "Conventional RF")]
control_treat <- rbind(idcol="guidelines", "ACC.AHA.2019"=control_treat, "NICE.2014"=control_treat)
control_treat[ons_conv_rf_summary[risk_group == "high"], on = .(guidelines), conv_rf := i.controls]
control_treat[ons_intermed_treat_summary[risk_group == "high"], on = .(guidelines), intermed_treat := i.controls]
control_treat[ons_pred_summary[risk_group == "high"], on = .(name, lambda, PGS, guidelines), pred_treat := i.controls]
control_treat[is.na(pred_treat), pred_treat := 0]
control_treat <- melt(control_treat, measure.vars=c("conv_rf", "intermed_treat", "pred_treat"),
                           variable.name="treatment_reason", value.name="controls_treated")
control_treat[, treatment_reason := fcase(
  treatment_reason == "conv_rf", "High risk due to conventional risk factors",
  treatment_reason == "intermed_treat", "Medium risk, diabetic or LDL-C >= 5 mmol/L",
  treatment_reason == "pred_treat", "Reclassified as high risk from biomarkers/PRS"
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
