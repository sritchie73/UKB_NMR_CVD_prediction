library(data.table)
library(foreach)

# Create output directory
out_dir <- "analyses/public_health_modelling/UK_population_generalised/targeted_screening"
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

# filter to models of interest
abs_risk <- abs_risk[name != "CRP" & lambda != "lambda.1se"]

# Set up models of interest
strategies <- rbind(use.names=FALSE,
  data.table(blanket_screening="Conventional RF", targeted_screening="Conventional RF + NMR"),
  data.table("Conventional RF", "Conventional RF + PRS"),
  data.table("Conventional RF", "Conventional RF + NMR + PRS"),
  data.table("Conventional RF + PRS", "Conventional RF + NMR + PRS"),
  data.table("Conventional RF + NMR", "Conventional RF + NMR + PRS")
)

# map strategies to models
map <- rbind(
  data.table(strategy="Conventional RF", name="Conventional RF", PGS=FALSE),
  data.table(strategy="Conventional RF + NMR", name="NMR", PGS=FALSE),
  data.table(strategy="Conventional RF + PRS", name="Conventional RF", PGS=TRUE),
  data.table(strategy="Conventional RF + NMR + PRS", name="NMR", PGS=TRUE)
)

# Extract relavent absolute risk and stratify accordingly:
abs_risk <- foreach(sidx = strategies[,.I], .combine=rbind) %do% {
  # Extract absolute risks needed for blanket screening
  this_blanket <- strategies[sidx, blanket_screening]
  this_abs_risk <- abs_risk[map[strategy == this_blanket], on = .(name, PGS)]
  this_abs_risk <- this_abs_risk[,.(eid, sex, age_group, incident_cvd, blanket_screening_risk=recalibrated_risk, ldl, diabetes)]

  # Add in absolute risks needed for targeted screening
  this_targeted <- strategies[sidx, targeted_screening]
  updated_risk <- abs_risk[map[strategy == this_targeted], on = .(name, PGS)]
  this_abs_risk[updated_risk, on = .(eid), targeted_screening_risk := i.recalibrated_risk]

  # add in screening strategy information
  this_abs_risk <- cbind(strategies[sidx], this_abs_risk)

  # Apply the blanket screening stratification
  this_abs_risk[, risk_group := fcase(
    blanket_screening_risk <= 0.05, "low",
    blanket_screening_risk > 0.10, "high",
    default="medium"
  )]

  # For people at intermediate risk, apply the second step, setting to high risk those with
  # a history of diabetes or elevated LDL-C (> 5.0 mmol/L)
  this_abs_risk[, risk_group2 := risk_group]
  this_abs_risk[risk_group == "medium" & (diabetes | ldl > 5), risk_group2 := "high"]

  # Finally apply the targeted screening step
  this_abs_risk[, risk_group3 := risk_group2]
  this_abs_risk[risk_group2 == "medium", risk_group3 := fcase(
    targeted_screening_risk <= 0.05, "low",
    targeted_screening_risk > 0.10, "high",
    default="medium"
  )] 

  # Return
  return(this_abs_risk)
}

# ----------------------------------------
# generalise first step to ONS population
# ----------------------------------------

# Get total number of cases and controls in each age and sex group - this is the denominator when computing
# % of cases and controls allocated to risk strata for a given model and thresholds
grp_totals <- abs_risk[
  # Group totals identical for all strategies at this step
  blanket_screening == "Conventional RF" & targeted_screening == "Conventional RF + NMR",
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
abs_risk_strata_alloc <- foreach(sIdx = strategies[,.I], .combine=rbind) %do% {
  cbind(strategies[sIdx], empty)
}
to_fill <- abs_risk[, .(samples=.N, cases=sum(incident_cvd)), by=.(blanket_screening, targeted_screening, sex, age_group, risk_group)]
to_fill[grp_totals, on = .(sex, age_group), pct_samples := samples/i.N]
to_fill[grp_totals, on = .(sex, age_group), pct_cases := cases/i.cases]
abs_risk_strata_alloc[to_fill, on = .(blanket_screening, targeted_screening, sex, age_group, risk_group), pct_samples := i.pct_samples]
abs_risk_strata_alloc[to_fill, on = .(blanket_screening, targeted_screening, sex, age_group, risk_group), pct_cases := i.pct_cases]

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
  .(blanket_screening, targeted_screening, sex, age_group, risk_group)]
blanket_screening[ons_pop, on = .(sex, age_group), grp_total := i.N]
blanket_screening[ons_pop, on = .(sex, age_group), case_total := i.cases]

blanket_screening[abs_risk_strata_alloc[risk_group != "low"], on = .(blanket_screening, targeted_screening, sex, age_group, risk_group),
  allocated := floor(grp_total * pct_samples)]
alloc_risk <- blanket_screening[risk_group != "low", .(risk_group="low", allocated=sum(allocated)), by=.(blanket_screening, targeted_screening, sex, age_group)] 
blanket_screening[alloc_risk, on = .(blanket_screening, targeted_screening, sex, age_group, risk_group), allocated := grp_total - i.allocated]

blanket_screening[abs_risk_strata_alloc[risk_group != "low"], on = .(blanket_screening, targeted_screening, sex, age_group, risk_group), 
  cases := floor(pct_cases * case_total)]
alloc_cases <- blanket_screening[risk_group != "low", .(risk_group="low", allocated=sum(cases)), by=.(blanket_screening, targeted_screening, sex, age_group)]
blanket_screening[alloc_cases, on = .(blanket_screening, targeted_screening, sex, age_group, risk_group), cases := case_total - i.allocated]

blanket_screening[, controls := allocated - cases]
blanket_screening[, c("grp_total", "case_total") := NULL]

# Summarise to population totals for each risk group
blanket_screening_summary <- blanket_screening[,
  .(N=sum(cases) + sum(controls), cases=sum(cases), controls=sum(controls)),
  by=.(blanket_screening, targeted_screening, risk_group)]

# Write out
fwrite(blanket_screening, sep="\t", quote=FALSE, file=sprintf("%s/ONS_screening_step_1_stratified_by_age_sex.txt", out_dir))
fwrite(blanket_screening_summary, sep="\t", quote=FALSE, file=sprintf("%s/ONS_screening_step_1.txt", out_dir))

# Also document numbers in UKB
ukb_abs_risk <- foreach(sIdx = strategies[,.I], .combine=rbind) %do% {
  cbind(strategies[sIdx], empty)
}
setnames(ukb_abs_risk, c("pct_cases", "pct_samples"), c("cases", "samples"))
to_fill <- abs_risk[, .(cases=sum(incident_cvd), samples=.N),
  by=.(blanket_screening, targeted_screening, sex, age_group, risk_group)]
ukb_abs_risk[to_fill, on = .(blanket_screening, targeted_screening, sex, age_group, risk_group),
  c("cases", "samples") := .(i.cases, i.samples)]
ukb_abs_risk[, controls := samples - cases]

ukb_abs_risk_summary <- ukb_abs_risk[,
  .(cases=sum(cases), controls=sum(controls), samples=sum(samples)),
  by=.(blanket_screening, targeted_screening, risk_group)]

fwrite(ukb_abs_risk, sep="\t", quote=FALSE, file=sprintf("%s/UKB_screening_step_1_stratified_by_age_sex.txt", out_dir))
fwrite(ukb_abs_risk_summary, sep="\t", quote=FALSE, file=sprintf("%s/UKB_screening_step_1.txt", out_dir))

# -----------------------------------------------------------
# Now do the same thing for the second stage of screening
# -----------------------------------------------------------

# Get total number of cases and controls in each age and sex group that
# were input at the second stage of the screening
grp_totals2 <- abs_risk[risk_group == "medium",
  .(.N, cases=sum(incident_cvd), controls=sum(!(incident_cvd))),
  by=.(blanket_screening, targeted_screening, age_group, sex)]
grp_totals2 <- grp_totals2[order(age_group)][order(sex)][order(blanket_screening)][order(targeted_screening)]

# Build empty table of all possible groups - this allows us to fill in 0s for groups with no participants when computing % allocated
empty2 <- empty[risk_group != "low"]
setnames(empty2, "risk_group", "risk_group2")

# Determine % of samples and cases allocated to each risk strata by each model at the first screening step
abs_risk_strata_alloc2 <- foreach(sIdx = strategies[,.I], .combine=rbind) %do% {
  cbind(strategies[sIdx], empty2)
}
to_fill <- abs_risk[risk_group == "medium", .(samples=.N, cases=sum(incident_cvd)), by=.(blanket_screening, targeted_screening, sex, age_group, risk_group2)]
to_fill[grp_totals2, on = .(blanket_screening, targeted_screening, sex, age_group), pct_samples := samples/i.N]
to_fill[grp_totals2, on = .(blanket_screening, targeted_screening, sex, age_group), pct_cases := cases/i.cases]
abs_risk_strata_alloc2[to_fill, on = .(blanket_screening, targeted_screening, sex, age_group, risk_group2), pct_samples := i.pct_samples]
abs_risk_strata_alloc2[to_fill, on = .(blanket_screening, targeted_screening, sex, age_group, risk_group2), pct_cases := i.pct_cases]

# Generalise to hypothetical ONS population. 
intermed_treat <- abs_risk_strata_alloc2[ons_pop, on = .(sex, age_group), 
  .(blanket_screening, targeted_screening, sex, age_group, risk_group2)]
med_risk <- blanket_screening[risk_group == "medium"]
intermed_treat[med_risk, on = .(blanket_screening, targeted_screening, sex, age_group), grp_total := i.allocated]
intermed_treat[med_risk, on = .(blanket_screening, targeted_screening, sex, age_group), case_total := i.cases]

intermed_treat[abs_risk_strata_alloc2[risk_group2 == "high"], on = .(blanket_screening, targeted_screening, sex, age_group, risk_group2),
  allocated := floor(grp_total * pct_samples)]
alloc_risk <- intermed_treat[risk_group2 == "high", .(risk_group2="medium", allocated=sum(allocated)), by=.(blanket_screening, targeted_screening, sex, age_group)] 
intermed_treat[alloc_risk, on = .(blanket_screening, targeted_screening, sex, age_group, risk_group2), allocated := grp_total - i.allocated]

intermed_treat[abs_risk_strata_alloc2[risk_group2 == "high"], on = .(blanket_screening, targeted_screening, sex, age_group, risk_group2), 
  cases := floor(pct_cases * case_total)]
alloc_cases <- intermed_treat[risk_group2 == "high", .(risk_group2="medium", allocated=sum(cases)), by=.(blanket_screening, targeted_screening, sex, age_group)]
intermed_treat[alloc_cases, on = .(blanket_screening, targeted_screening, sex, age_group, risk_group2), cases := case_total - i.allocated]

intermed_treat[, controls := allocated - cases]
intermed_treat[, c("grp_total", "case_total") := NULL]

# Summarise to population totals for each risk group
intermed_treat_summary <- intermed_treat[,
  .(N=sum(cases) + sum(controls), cases=sum(cases), controls=sum(controls)),
  by=.(blanket_screening, targeted_screening, risk_group2)]

# Write out
fwrite(intermed_treat, sep="\t", quote=FALSE, file=sprintf("%s/ONS_screening_step_2_stratified_by_age_sex.txt", out_dir))
fwrite(intermed_treat_summary, sep="\t", quote=FALSE, file=sprintf("%s/ONS_screening_step_2_stratified.txt", out_dir))

# Also document numbers in UKB
ukb_abs_risk2 <- foreach(sIdx = strategies[,.I], .combine=rbind) %do% {
  cbind(strategies[sIdx], empty2)
}
setnames(ukb_abs_risk2, c("pct_cases", "pct_samples"), c("cases", "samples"))
to_fill <- abs_risk[risk_group == "medium", .(cases=sum(incident_cvd), samples=.N),
  by=.(blanket_screening, targeted_screening, sex, age_group, risk_group2)]
ukb_abs_risk2[to_fill, on = .(blanket_screening, targeted_screening, sex, age_group, risk_group2),
  c("cases", "samples") := .(i.cases, i.samples)]
ukb_abs_risk2[, controls := samples - cases]

ukb_abs_risk_summary2 <- ukb_abs_risk2[,
  .(cases=sum(cases), controls=sum(controls), samples=sum(samples)),
  by=.(blanket_screening, targeted_screening, risk_group2)]

fwrite(ukb_abs_risk2, sep="\t", quote=FALSE, file=sprintf("%s/UKB_screening_step_2_stratified_by_age_sex.txt", out_dir))
fwrite(ukb_abs_risk_summary2, sep="\t", quote=FALSE, file=sprintf("%s/UKB_screening_step_2.txt", out_dir))

# -----------------------------------------------------------
# Now do the same thing for the third stage of screening
# -----------------------------------------------------------

# Get total number of cases and controls in each age and sex group that
# were input at the third stage of the screening
grp_totals3 <- abs_risk[risk_group == "medium" & risk_group2 == "medium",
  .(.N, cases=sum(incident_cvd), controls=sum(!(incident_cvd))),
  by=.(blanket_screening, targeted_screening, age_group, sex)]
grp_totals3 <- grp_totals3[order(age_group)][order(sex)][order(blanket_screening)][order(targeted_screening)]

# Build empty table of all possible groups - this allows us to fill in 0s for groups with no participants when computing % allocated
empty3 <- copy(empty)
setnames(empty3, "risk_group", "risk_group3")

# Determine % of samples and cases allocated to each risk strata by each model at the first screening step
abs_risk_strata_alloc3 <- foreach(sIdx = strategies[,.I], .combine=rbind) %do% {
  cbind(strategies[sIdx], empty3)
}
to_fill <- abs_risk[risk_group == "medium" & risk_group2 == "medium", .(samples=.N, cases=sum(incident_cvd)), 
  by=.(blanket_screening, targeted_screening, sex, age_group, risk_group3)]
to_fill[grp_totals3, on = .(blanket_screening, targeted_screening, sex, age_group), pct_samples := samples/i.N]
to_fill[grp_totals3, on = .(blanket_screening, targeted_screening, sex, age_group), pct_cases := cases/i.cases]
abs_risk_strata_alloc3[to_fill, on = .(blanket_screening, targeted_screening, sex, age_group, risk_group3), pct_samples := i.pct_samples]
abs_risk_strata_alloc3[to_fill, on = .(blanket_screening, targeted_screening, sex, age_group, risk_group3), pct_cases := i.pct_cases]

# Generalise to hypothetical ONS population. 
targeted_screening <- abs_risk_strata_alloc3[ons_pop, on = .(sex, age_group), 
  .(blanket_screening, targeted_screening, sex, age_group, risk_group3)]
med_risk <- intermed_treat[risk_group2 == "medium"]
targeted_screening[med_risk, on = .(blanket_screening, targeted_screening, sex, age_group), grp_total := i.allocated]
targeted_screening[med_risk, on = .(blanket_screening, targeted_screening, sex, age_group), case_total := i.cases]

targeted_screening[abs_risk_strata_alloc3[risk_group3 != "low"], on = .(blanket_screening, targeted_screening, sex, age_group, risk_group3),
  allocated := floor(grp_total * pct_samples)]
alloc_risk <- targeted_screening[risk_group3 != "low", .(risk_group3="low", allocated=sum(allocated)), by=.(blanket_screening, targeted_screening, sex, age_group)] 
targeted_screening[alloc_risk, on = .(blanket_screening, targeted_screening, sex, age_group, risk_group3), allocated := grp_total - i.allocated]

targeted_screening[abs_risk_strata_alloc3[risk_group3 != "low"], on = .(blanket_screening, targeted_screening, sex, age_group, risk_group3), 
  cases := floor(pct_cases * case_total)]
alloc_cases <- targeted_screening[risk_group3 != "low", .(risk_group3="low", allocated=sum(cases)), by=.(blanket_screening, targeted_screening, sex, age_group)]
targeted_screening[alloc_cases, on = .(blanket_screening, targeted_screening, sex, age_group, risk_group3), cases := case_total - i.allocated]

targeted_screening[, controls := allocated - cases]
targeted_screening[, c("grp_total", "case_total") := NULL]

# Summarise to population totals for each risk group
targeted_screening_summary <- targeted_screening[,
  .(N=sum(cases) + sum(controls), cases=sum(cases), controls=sum(controls)),
  by=.(blanket_screening, targeted_screening, risk_group3)]

# Write out
fwrite(targeted_screening, sep="\t", quote=FALSE, file=sprintf("%s/ONS_screening_step_3_stratified_by_age_sex.txt", out_dir))
fwrite(targeted_screening_summary, sep="\t", quote=FALSE, file=sprintf("%s/ONS_screening_step_3_stratified.txt", out_dir))

# Also document numbers in UKB
ukb_abs_risk3 <- foreach(sIdx = strategies[,.I], .combine=rbind) %do% {
  cbind(strategies[sIdx], empty3)
}
setnames(ukb_abs_risk3, c("pct_cases", "pct_samples"), c("cases", "samples"))
to_fill <- abs_risk[risk_group == "medium" & risk_group2 == "medium", 
  .(cases=sum(incident_cvd), samples=.N),
  by=.(blanket_screening, targeted_screening, sex, age_group, risk_group3)]
ukb_abs_risk3[to_fill, on = .(blanket_screening, targeted_screening, sex, age_group, risk_group3),
  c("cases", "samples") := .(i.cases, i.samples)]
ukb_abs_risk3[, controls := samples - cases]

ukb_abs_risk_summary3 <- ukb_abs_risk3[,
  .(cases=sum(cases), controls=sum(controls), samples=sum(samples)),
  by=.(blanket_screening, targeted_screening, risk_group3)]

fwrite(ukb_abs_risk3, sep="\t", quote=FALSE, file=sprintf("%s/UKB_screening_step_3_stratified_by_age_sex.txt", out_dir))
fwrite(ukb_abs_risk_summary3, sep="\t", quote=FALSE, file=sprintf("%s/UKB_screening_step_3.txt", out_dir))

# ---------------------------------
# Compute public health statistics
# ---------------------------------

# Compute number of cases treated and prevented using the flowchart for each set of models
case_treat_prevent <- copy(strategies)
case_treat_prevent[blanket_screening_summary[risk_group == "high"], on = .(blanket_screening, targeted_screening), step_1_treat := i.cases]
case_treat_prevent[intermed_treat_summary[risk_group2 == "high"], on = .(blanket_screening, targeted_screening), step_2_treat := i.cases]
case_treat_prevent[targeted_screening_summary[risk_group3 == "high"], on = .(blanket_screening, targeted_screening), step_3_treat := i.cases]
case_treat_prevent <- melt(case_treat_prevent, id.vars=c("blanket_screening", "targeted_screening"), variable.name="screening_step", value.name="cases_treated")
case_treat_prevent[, screening_step := as.integer(gsub("step_", "", gsub("_treat", "", screening_step)))]
case_treat_prevent[, pct_cases_treated := cases_treated / ons_pop_summary$cases]

# Same for controls
control_treat <- copy(strategies)
control_treat[blanket_screening_summary[risk_group == "high"], on = .(blanket_screening, targeted_screening), step_1_treat := i.controls]
control_treat[intermed_treat_summary[risk_group2 == "high"], on = .(blanket_screening, targeted_screening), step_2_treat := i.controls]
control_treat[targeted_screening_summary[risk_group3 == "high"], on = .(blanket_screening, targeted_screening), step_3_treat := i.controls]
control_treat <- melt(control_treat, id.vars=c("blanket_screening", "targeted_screening"), variable.name="screening_step", value.name="controls_treated")
control_treat[, screening_step := as.integer(gsub("step_", "", gsub("_treat", "", screening_step)))]
control_treat[, pct_controls_treated := controls_treated / ons_pop_summary$controls]

# Combine and compute other public health statistics for each treatment reason, model, and threshold guidelines
phs <- case_treat_prevent[control_treat, on = .(blanket_screening, targeted_screening, screening_step)]
phs[, events_prevented := floor(cases_treated * 0.2)] # assuming 20% reduction in risk due to statins
phs[, pct_events_prevented := pct_cases_treated * 0.2]

# Summarise to population level
phs_summary <- phs[, .(
  total_treated=sum(cases_treated) + sum(controls_treated),
  pct_total_treated = (sum(cases_treated) + sum(controls_treated)) / 100000,
  cases_treated=sum(cases_treated), pct_cases_treated=sum(pct_cases_treated),
  controls_treated=sum(controls_treated), pct_controls_treated=sum(pct_controls_treated),
  events_prevented=sum(events_prevented), pct_events_prevented=sum(pct_events_prevented)),
  by=.(blanket_screening, targeted_screening)]
phs_summary[, NNT := ceiling(total_treated / events_prevented)]
phs_summary[, NNS := ceiling(100000 / events_prevented)]

# Write out
fwrite(phs, sep="\t", quote=FALSE, file=sprintf("%s/ONS_public_health_statistics_by_screening_step.txt", out_dir))
fwrite(phs_summary, sep="\t", quote=FALSE, file=sprintf("%s/ONS_public_health_statistics.txt", out_dir))

