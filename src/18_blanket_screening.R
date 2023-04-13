library(data.table)
library(foreach)

# Create output directory
out_dir <- "analyses/public_health_modelling/UK_population_generalised/blanket_screening"
system(sprintf("mkdir -p %s", out_dir))

# Load hypothetical population 
ons_pop <- fread("analyses/public_health_modelling/UK_population_generalised/ONS_hypothetical_100k_pop_by_age_sex.txt")
ons_pop_summary <- fread("analyses/public_health_modelling/UK_population_generalised/ONS_hypothetical_100k_pop.txt")

# load precomputed absolute risks
abs_risk <- fread("analyses/public_health_modelling/risk_recalibration/absolute_risks.txt")

# Load test dataset and add in LDL cholesterol and history of diabetes
test <- fread("data/processed/test/processed_test_data.txt")
abs_risk[test, on = .(eid), ldl := i.ldl]
abs_risk[test, on = .(eid), diabetes := i.diabetes]

# Load model information
model_info <- fread("analyses/test/model_fit_information.txt")

# filter to models of interest
model_info <- model_info[name != "CRP" & lambda != "lambda.1se"]
abs_risk <- abs_risk[name != "CRP" & lambda != "lambda.1se"]

# Stratify risk according to NICE 2014 guidelines
abs_risk[, risk_group := fcase(
  recalibrated_risk <= 0.05, "low",
  recalibrated_risk > 0.10, "high",
  default="medium"
)]

# For people at intermediate risk, apply a second step, setting to high risk those with
# a history of diabetes or elevated LDL-C (> 5.0 mmol/L)
abs_risk[, risk_group2 := risk_group]
abs_risk[risk_group == "medium" & (diabetes | ldl > 5), risk_group2 := "high"]

# Get total number of cases and controls in each age and sex group - this is the denominator when computing
# % of cases and controls allocated to risk strata for a given model and thresholds
grp_totals <- abs_risk[long_name == "Conventional RF",
  .(.N, cases=sum(incident_cvd), controls=sum(!(incident_cvd))),
  by=.(age_group, sex)]
grp_totals <- grp_totals[order(age_group)][order(sex)]

# Write out to document sample sizes in UKB
fwrite(grp_totals, sep="\t", quote=FALSE, file=sprintf("%s/UKB_model_sample_size_by_age_sex.txt", out_dir))

# Build empty table of all possible groups - this allows us to fill in 0s for groups with no participants when computing % allocated
empty <- expand.grid(sex=unique(abs_risk$sex), age_group=unique(abs_risk$age_group), risk_group=unique(abs_risk$risk_group), stringsAsFactors=FALSE)
setDT(empty)
empty <- rbind(empty[risk_group == "low"], empty[risk_group == "medium"], empty[risk_group == "high"])
empty <- empty[order(age_group)][order(sex)]
empty[, pct_samples := 0]
empty[, pct_cases := 0]

# Determine % of samples and cases allocated to each risk strata by each model at the first screening step
models <- unique(abs_risk[,.(name, PGS, long_name)])
abs_risk_strata_alloc <- foreach(mIdx = models[,.I], .combine=rbind) %do% {
  cbind(models[mIdx], empty)
}
to_fill <- abs_risk[, .(samples=.N, cases=sum(incident_cvd)), by=.(name, PGS, long_name, sex, age_group, risk_group)]
to_fill[grp_totals, on = .(sex, age_group), pct_samples := samples/i.N]
to_fill[grp_totals, on = .(sex, age_group), pct_cases := cases/i.cases]
abs_risk_strata_alloc[to_fill, on = .(name, PGS, long_name, sex, age_group, risk_group), pct_samples := i.pct_samples]
abs_risk_strata_alloc[to_fill, on = .(name, PGS, long_name, sex, age_group, risk_group), pct_cases := i.pct_cases]

# Generalise to hypothetical ONS population. 
#
# Making sure the numbers add up is slightly tricky here, as simply multiplying out the percentages 
# gives fractional numbers (which don't make sense in the context of a hypothetical population). 
# 
# To handle this, we estimate the number of people allocated to high and intermediate risk groups,
# taking the floor of the number (so we don't end up with fractional samples), then allocate the 
# remainder to the low risk group.
# 
# Likewise, we then do the same for CVD cases, multiplying out the percentage and taking the floor,
# setting the remainder as controls. Cases not allocated to intermediate or high risk groups are allocated
# to the low risk group (to handle fractional cases)
blanket_screening <- abs_risk_strata_alloc[ons_pop, on = .(sex, age_group), 
  .(name, PGS, long_name, sex, age_group, risk_group)]
blanket_screening[ons_pop, on = .(sex, age_group), grp_total := i.N]
blanket_screening[ons_pop, on = .(sex, age_group), case_total := i.cases]

blanket_screening[abs_risk_strata_alloc[risk_group != "low"], on = .(name, PGS, long_name, sex, age_group, risk_group),
  allocated := floor(grp_total * pct_samples)]
alloc_risk <- blanket_screening[risk_group != "low", .(risk_group="low", allocated=sum(allocated)), by=.(name, PGS, long_name, sex, age_group)] 
blanket_screening[alloc_risk, on = .(name, PGS, long_name, sex, age_group, risk_group), allocated := grp_total - i.allocated]

blanket_screening[abs_risk_strata_alloc[risk_group != "low"], on = .(name, PGS, long_name, sex, age_group, risk_group), 
  cases := floor(pct_cases * case_total)]
alloc_cases <- blanket_screening[risk_group != "low", .(risk_group="low", allocated=sum(cases)), by=.(name, PGS, long_name, sex, age_group)]
blanket_screening[alloc_cases, on = .(name, PGS, long_name, sex, age_group, risk_group), cases := case_total - i.allocated]

blanket_screening[, controls := allocated - cases]
blanket_screening[, c("grp_total", "case_total") := NULL]

# Summarise to population totals for each risk group
blanket_screening_summary <- blanket_screening[,
  .(N=sum(cases) + sum(controls), cases=sum(cases), controls=sum(controls)),
  by=.(name, PGS, long_name, risk_group)]

# Write out
fwrite(blanket_screening, sep="\t", quote=FALSE, file=sprintf("%s/ONS_screening_step_1_stratified_by_age_sex.txt", out_dir))
fwrite(blanket_screening_summary, sep="\t", quote=FALSE, file=sprintf("%s/ONS_screening_step_1.txt", out_dir))

# Also document numbers in UKB
ukb_abs_risk <- foreach(mIdx = models[,.I], .combine=rbind) %do% {
  cbind(models[mIdx], empty)
}
setnames(ukb_abs_risk, c("pct_cases", "pct_samples"), c("cases", "samples"))
to_fill <- abs_risk[, .(cases=sum(incident_cvd), samples=.N),
  by=.(name, PGS, long_name, sex, age_group, risk_group)]
ukb_abs_risk[to_fill, on = .(name, PGS, long_name, sex, age_group, risk_group),
  c("cases", "samples") := .(i.cases, i.samples)]
ukb_abs_risk[, controls := samples - cases]

ukb_abs_risk_summary <- ukb_abs_risk[,
  .(cases=sum(cases), controls=sum(controls), samples=sum(samples)),
  by=.(name, PGS, long_name, risk_group)]

fwrite(ukb_abs_risk, sep="\t", quote=FALSE, file=sprintf("%s/UKB_screening_step_1_stratified_by_age_sex.txt", out_dir))
fwrite(ukb_abs_risk_summary, sep="\t", quote=FALSE, file=sprintf("%s/UKB_screening_step_1.txt", out_dir))

# -----------------------------------------------------------
# Now do the same thing for the second stage of screening
# -----------------------------------------------------------

# Get total number of cases and controls in each age and sex group that
# were input at the second stage of the screening
grp_totals2 <- abs_risk[risk_group == "medium",
  .(.N, cases=sum(incident_cvd), controls=sum(!(incident_cvd))),
  by=.(name, PGS, long_name, age_group, sex)]
grp_totals2 <- grp_totals2[order(age_group)][order(sex)][order(long_name)]

# Build empty table of all possible groups - this allows us to fill in 0s for groups with no participants when computing % allocated
empty2 <- empty[risk_group != "low"]
setnames(empty2, "risk_group", "risk_group2")

# Determine % of samples and cases allocated to each risk strata by each model at the first screening step
models <- unique(abs_risk[,.(name, PGS, long_name)])
abs_risk_strata_alloc2 <- foreach(mIdx = models[,.I], .combine=rbind) %do% {
  cbind(models[mIdx], empty2)
}
to_fill <- abs_risk[risk_group == "medium", .(samples=.N, cases=sum(incident_cvd)), by=.(name, PGS, long_name, sex, age_group, risk_group2)]
to_fill[grp_totals2, on = .(name, PGS, long_name, sex, age_group), pct_samples := samples/i.N]
to_fill[grp_totals2, on = .(name, PGS, long_name, sex, age_group), pct_cases := cases/i.cases]
abs_risk_strata_alloc2[to_fill, on = .(name, PGS, long_name, sex, age_group, risk_group2), pct_samples := i.pct_samples]
abs_risk_strata_alloc2[to_fill, on = .(name, PGS, long_name, sex, age_group, risk_group2), pct_cases := i.pct_cases]

# Generalise to hypothetical ONS population. 
intermed_treat <- abs_risk_strata_alloc2[ons_pop, on = .(sex, age_group), 
  .(name, PGS, long_name, sex, age_group, risk_group2)]
intermed_treat[blanket_screening[risk_group == "medium"], on = .(name, PGS, long_name, sex, age_group), grp_total := i.allocated]
intermed_treat[blanket_screening[risk_group == "medium"], on = .(name, PGS, long_name, sex, age_group), case_total := i.cases]

intermed_treat[abs_risk_strata_alloc2[risk_group2 == "high"], on = .(name, PGS, long_name, sex, age_group, risk_group2),
  allocated := floor(grp_total * pct_samples)]
alloc_risk <- intermed_treat[risk_group2 == "high", .(risk_group2="medium", allocated=sum(allocated)), by=.(name, PGS, long_name, sex, age_group)] 
intermed_treat[alloc_risk, on = .(name, PGS, long_name, sex, age_group, risk_group2), allocated := grp_total - i.allocated]

intermed_treat[abs_risk_strata_alloc2[risk_group2 == "high"], on = .(name, PGS, long_name, sex, age_group, risk_group2), 
  cases := floor(pct_cases * case_total)]
alloc_cases <- intermed_treat[risk_group2 == "high", .(risk_group2="medium", allocated=sum(cases)), by=.(name, PGS, long_name, sex, age_group)]
intermed_treat[alloc_cases, on = .(name, PGS, long_name, sex, age_group, risk_group2), cases := case_total - i.allocated]

intermed_treat[, controls := allocated - cases]
intermed_treat[, c("grp_total", "case_total") := NULL]

# Summarise to population totals for each risk group
intermed_treat_summary <- intermed_treat[,
  .(N=sum(cases) + sum(controls), cases=sum(cases), controls=sum(controls)),
  by=.(name, PGS, long_name, risk_group2)]

# Write out
fwrite(intermed_treat, sep="\t", quote=FALSE, file=sprintf("%s/ONS_screening_step_2_stratified_by_age_sex.txt", out_dir))
fwrite(intermed_treat_summary, sep="\t", quote=FALSE, file=sprintf("%s/ONS_screening_step_2_stratified.txt", out_dir))

# Also document numbers in UKB
ukb_abs_risk2 <- foreach(mIdx = models[,.I], .combine=rbind) %do% {
  cbind(models[mIdx], empty2)
}
setnames(ukb_abs_risk2, c("pct_cases", "pct_samples"), c("cases", "samples"))
to_fill <- abs_risk[risk_group == "medium", .(cases=sum(incident_cvd), samples=.N),
  by=.(name, PGS, long_name, sex, age_group, risk_group2)]
ukb_abs_risk2[to_fill, on = .(name, PGS, long_name, sex, age_group, risk_group2),
  c("cases", "samples") := .(i.cases, i.samples)]
ukb_abs_risk2[, controls := samples - cases]

ukb_abs_risk_summary2 <- ukb_abs_risk2[,
  .(cases=sum(cases), controls=sum(controls), samples=sum(samples)),
  by=.(name, PGS, long_name, risk_group2)]

fwrite(ukb_abs_risk2, sep="\t", quote=FALSE, file=sprintf("%s/UKB_screening_step_2_stratified_by_age_sex.txt", out_dir))
fwrite(ukb_abs_risk_summary2, sep="\t", quote=FALSE, file=sprintf("%s/UKB_screening_step_2.txt", out_dir))

# ----------------
# Compute number of cases treated and prevented using the flowchart for each set of models
case_treat_prevent <- copy(models)
case_treat_prevent[blanket_screening_summary[risk_group == "high"], on = .(name, PGS, long_name), step_1_treat := i.cases]
case_treat_prevent[intermed_treat_summary[risk_group2 == "high"], on = .(name, PGS, long_name), step_2_treat := i.cases]
case_treat_prevent <- melt(case_treat_prevent, id.vars=c("name", "PGS", "long_name"), variable.name="screening_step", value.name="cases_treated")
case_treat_prevent[, screening_step := as.integer(gsub("step_", "", gsub("_treat", "", screening_step)))]
case_treat_prevent[, pct_cases_treated := cases_treated / ons_pop_summary$cases]

# Same for controls
control_treat <- copy(models)
control_treat[blanket_screening_summary[risk_group == "high"], on = .(name, PGS, long_name), step_1_treat := i.controls]
control_treat[intermed_treat_summary[risk_group2 == "high"], on = .(name, PGS, long_name), step_2_treat := i.controls]
control_treat <- melt(control_treat, id.vars=c("name", "PGS", "long_name"), variable.name="screening_step", value.name="controls_treated")
control_treat[, screening_step := as.integer(gsub("step_", "", gsub("_treat", "", screening_step)))]
control_treat[, pct_controls_treated := controls_treated / ons_pop_summary$controls]

# Combine and compute other public health statistics for each treatment reason, model, and threshold guidelines
phs <- case_treat_prevent[control_treat, on = .(name, PGS, long_name, screening_step)]
phs[, events_prevented := floor(cases_treated * 0.2)] # assuming 20% reduction in risk due to statins
phs[, pct_events_prevented := pct_cases_treated * 0.2]

# Summarise to population level
phs_summary <- phs[, .(
  total_treated=sum(cases_treated) + sum(controls_treated),
  pct_total_treated = (sum(cases_treated) + sum(controls_treated)) / 100000,
  cases_treated=sum(cases_treated), pct_cases_treated=sum(pct_cases_treated),
  controls_treated=sum(controls_treated), pct_controls_treated=sum(pct_controls_treated),
  events_prevented=sum(events_prevented), pct_events_prevented=sum(pct_events_prevented)),
  by=.(name, PGS, long_name)]
phs_summary[, NNT := ceiling(total_treated / events_prevented)]
phs_summary[, NNS := ceiling(100000 / events_prevented)]

# Write out
fwrite(phs, sep="\t", quote=FALSE, file=sprintf("%s/ONS_public_health_statistics_by_screening_step.txt", out_dir))
fwrite(phs_summary, sep="\t", quote=FALSE, file=sprintf("%s/ONS_public_health_statistics.txt", out_dir))

