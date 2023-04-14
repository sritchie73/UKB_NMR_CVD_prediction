library(data.table)
library(ggplot2)
library(ggforce)
library(scales)

# Load ONS pop so we know total cases
ons_pop <- fread("analyses/public_health_modelling/UK_population_generalised/ONS_hypothetical_100k_pop.txt")

# Load information about case identification and harmonize
blanket <- fread("analyses/public_health_modelling/UK_population_generalised/blanket_screening/ONS_public_health_statistics_by_screening_step.txt")
targeted <- fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening/ONS_public_health_statistics_by_screening_step.txt")

blanket[, blanket_screening := fcase(
  name == "Conventional RF" & !(PGS), "Conventional RF",
  name == "Conventional RF" & (PGS), "Conventional RF + PRS",
  name == "NMR" & !(PGS), "Conventional RF + NMR",
  name == "NMR" & (PGS), "Conventional RF + NMR + PRS"
)]
blanket[, targeted_screening := NA_character_]
blanket[, c("name", "PGS", "long_name") := NULL]

caseid <- rbind(targeted, blanket)

# Load information about number needed to screen/treat and harmonize
blanket <- fread("analyses/public_health_modelling/UK_population_generalised/blanket_screening/ONS_public_health_statistics.txt")
targeted <- fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening/ONS_public_health_statistics.txt")

blanket[, blanket_screening := fcase(
  name == "Conventional RF" & !(PGS), "Conventional RF",
  name == "Conventional RF" & (PGS), "Conventional RF + PRS",
  name == "NMR" & !(PGS), "Conventional RF + NMR",
  name == "NMR" & (PGS), "Conventional RF + NMR + PRS"
)]
blanket[, targeted_screening := NA_character_]
blanket[, c("name", "PGS", "long_name") := NULL]

phs <- rbind(targeted, blanket)

# Add plotting names for treatment reason
caseid[, treatment_reason := fcase(
  screening_step == 1, "High risk in blanket screening",
  screening_step == 2, "Medium risk but diabetic or LDL-C >= 5 mmol/L",
  screening_step == 3, "High risk in targeted screening"
)]
caseid[, treatment_reason := factor(treatment_reason, levels=rev(c("High risk in blanket screening", "Medium risk but diabetic or LDL-C >= 5 mmol/L", "High risk in targeted screening")))]

# Add in facet strata
caseid[, facet_group := blanket_screening]
caseid[is.na(targeted_screening), facet_group := "Blanket screening"]
caseid[, facet_group := factor(facet_group, levels=c("Blanket screening", "Conventional RF", "Conventional RF + PRS", "Conventional RF + NMR"))]

phs[, facet_group := blanket_screening]
phs[is.na(targeted_screening), facet_group := "Blanket screening"]
phs[, facet_group := factor(facet_group, levels=c("Blanket screening", "Conventional RF", "Conventional RF + PRS", "Conventional RF + NMR"))]

# Add in row names
caseid[, row_name := targeted_screening]
caseid[is.na(targeted_screening), row_name := blanket_screening]
caseid[, row_name := factor(row_name, levels=c("Conventional RF", "Conventional RF + PRS", "Conventional RF + NMR", "Conventional RF + NMR + PRS"))]

phs[, row_name := targeted_screening]
phs[is.na(targeted_screening), row_name := blanket_screening]
phs[, row_name := factor(row_name, levels=c("Conventional RF", "Conventional RF + PRS", "Conventional RF + NMR", "Conventional RF + NMR + PRS"))]

# Order rows
caseid <- caseid[order(-treatment_reason)][order(row_name)][order(facet_group)]
phs <- phs[order(row_name)][order(facet_group)]

# Select and order columns
caseid <- caseid[,.(facet_group, row_name, treatment_reason, cases_treated, pct_cases_treated)]
phs <- phs[,.(facet_group, row_name, NNT, NNS)]


# Plot number of cases identified
xlim <- c( 
  caseid[treatment_reason == "High risk in blanket screening", min(cases_treated)],
  caseid[, .(total_cases=sum(cases_treated)), by=.(facet_group, row_name)][, max(total_cases)]
)
xlim <- c(xlim[1] - diff(xlim)*0.1, xlim[2] + diff(xlim)*0.1)
xlim <- c(floor(xlim[1]/100)*100 - 30, ceiling(xlim[2]/100)*100 + 10)
g <- ggplot(caseid) +
  aes(x=cases_treated, y=row_name, fill=treatment_reason) +
  geom_bar(position="stack", stat="identity", color="black") +
  scale_fill_manual(name="", values=c(
    "High risk in blanket screening"="#bd0026",
    "Medium risk but diabetic or LDL-C >= 5 mmol/L"="#fc4e2a",
    "High risk in targeted screening"="#feb24c"
  )) + 
  facet_col(~ facet_group, space="free", scales="free_y") +
  scale_x_continuous(name="CVD cases caught", limits=xlim, oob=oob_keep, expand=c(0,0),
    sec.axis=sec_axis(~ . / ons_pop$cases * 100, name="% CVD cases caught")) + 
  ylab("") + 
  theme_bw() +
  theme(axis.text=element_text(size=7), axis.title=element_text(size=8), 
        panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
        legend.position="bottom", legend.text=element_text(size=8))
ggsave(g, width=10, height=6, file="analyses/public_health_modelling/UK_population_generalised/cases_identified.pdf")

# Plot NNS
xlim <- range(phs$NNS)
xlim <- c(xlim[1] - diff(xlim)*0.1, xlim[2] + diff(xlim)*0.1)
xlim <- c(floor(xlim[1]/10)*10 - 3, ceiling(xlim[2]/10)*10 + 1)
g <- ggplot(phs) +
  aes(x=NNS, y=row_name) +
  geom_col() + 
  facet_col(~ facet_group, space="free", scales="free_y") +
  scale_x_continuous(name="Number needed to screen", limits=xlim, oob=squish, expand=c(0,0)) +
  ylab("") + 
  theme_bw() +
  theme(axis.text=element_text(size=7), axis.title=element_text(size=8), 
        panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank())
ggsave(g, width=10, height=5, file="analyses/public_health_modelling/UK_population_generalised/NNS.pdf")

# Plot NNT
xlim <- range(phs$NNT)
xlim <- c(xlim[1] - diff(xlim)*0.05, xlim[2] + diff(xlim)*0.05)
xlim <- c(floor(xlim[1]/5)*5 - 3, ceiling(xlim[2]/5)*5 + 1)
g <- ggplot(phs) +
  aes(x=NNT, y=row_name) +
  geom_col() + 
  facet_col(~ facet_group, space="free", scales="free_y") +
  scale_x_continuous(name="Number needed to treat", limits=xlim, oob=squish, expand=c(0,0)) +
  ylab("") + 
  theme_bw() +
  theme(axis.text=element_text(size=7), axis.title=element_text(size=8), 
        panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank())
ggsave(g, width=10, height=5, file="analyses/public_health_modelling/UK_population_generalised/NNT.pdf")

stop()

# Now to build the wide table of numbers for the supp.
ukb_sample_size <- rbind(idcol="strategy",
  blanket = fread("analyses/public_health_modelling/UK_population_generalised/blanket_screening/UKB_model_sample_size.txt"),
  targeted = fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening/UKB_model_sample_size.txt"),
  targeted_above_PRS = fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening_above_PRS/UKB_model_sample_size.txt")
)

# Duplicate rows for the two sets of guidelines
ukb_sample_size <- rbind(idcol="guidelines", NICE.2014=ukb_sample_size, ACC.AHA.2019=ukb_sample_size)

# Enforce row ordering of models
ukb_sample_size[, name := factor(name, levels=c("Conventional RF", "NMR", "CRP", "NMR + Assays", "Assays"))]
ukb_sample_size[, PGS := factor(PGS, levels=c("FALSE", "TRUE"))]
ukb_sample_size <- ukb_sample_size[order(name)][order(PGS)]
ukb_sample_size <- ukb_sample_size[order(strategy)][order(-guidelines)]

# Add in long name
ukb_sample_size[phs, on = .(name, PGS, lambda), long_name := i.long_name]

# Filter to data we want to keep:
wide <- ukb_sample_size[lambda != "lambda.1se", .(guidelines, strategy, long_name, UKB_total=N, UKB_cases=cases, UKB_controls=controls)]

# Add in number of UKB participants allocated to low, medium, and high risk groups in blanket screening
ukb_blanket_alloc <- rbind(idcol="strategy",
  blanket = fread("analyses/public_health_modelling/UK_population_generalised/blanket_screening/UKB_biomarker_prs_risk_stratified.txt"),
  targeted = fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening/UKB_conv_rf_risk_stratified.txt"),
  targeted_above_PRS = fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening_above_PRS/UKB_conv_rf_prs_risk_stratified.txt")
)
wide[ukb_blanket_alloc[risk_group == "high"], on = .(long_name, guidelines, strategy), UKB_blanket_high_cases := i.cases]
wide[ukb_blanket_alloc[risk_group == "high"], on = .(long_name, guidelines, strategy), UKB_blanket_high_controls := i.controls]
wide[ukb_blanket_alloc[risk_group == "medium"], on = .(long_name, guidelines, strategy), UKB_blanket_medium_cases := i.cases]
wide[ukb_blanket_alloc[risk_group == "medium"], on = .(long_name, guidelines, strategy), UKB_blanket_medium_controls := i.controls]
wide[ukb_blanket_alloc[risk_group == "low"], on = .(long_name, guidelines, strategy), UKB_blanket_low_cases := i.cases]
wide[ukb_blanket_alloc[risk_group == "low"], on = .(long_name, guidelines, strategy), UKB_blanket_low_controls := i.controls]

# Add in number of UKB participants reclassified based on LDL and diabetes
ukb_intermed_treat <- rbind(idcol="strategy",
  blanket = fread("analyses/public_health_modelling/UK_population_generalised/blanket_screening/UKB_intermediate_risk_treat_alloc.txt"),
  targeted = fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening/UKB_intermediate_risk_treat_alloc.txt"),
  targeted_above_PRS = fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening_above_PRS/UKB_intermediate_risk_treat_alloc.txt")
)
wide[ukb_intermed_treat[risk_group == "high"], on = .(long_name, guidelines, strategy), UKB_intermed_high_cases := i.cases]
wide[ukb_intermed_treat[risk_group == "high"], on = .(long_name, guidelines, strategy), UKB_intermed_high_controls := i.controls]
wide[ukb_intermed_treat[risk_group == "medium"], on = .(long_name, guidelines, strategy), UKB_intermed_medium_cases := i.cases]
wide[ukb_intermed_treat[risk_group == "medium"], on = .(long_name, guidelines, strategy), UKB_intermed_medium_controls := i.controls]

# Add in number of UKB participants reclassified in targeted assessment
ukb_targeted <- rbind(idcol="strategy",
  targeted = fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening/UKB_targeted_biomarker_prs_risk_stratified.txt"),
  targeted_above_PRS = fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening_above_PRS/UKB_targeted_biomarker_risk_stratified.txt")
)
wide[ukb_targeted[risk_group == "high"], on = .(long_name, guidelines, strategy), UKB_targeted_high_cases := i.cases]
wide[ukb_targeted[risk_group == "high"], on = .(long_name, guidelines, strategy), UKB_targeted_high_controls := i.controls]
wide[ukb_targeted[risk_group == "medium"], on = .(long_name, guidelines, strategy), UKB_targeted_medium_cases := i.cases]
wide[ukb_targeted[risk_group == "medium"], on = .(long_name, guidelines, strategy), UKB_targeted_medium_controls := i.controls]
wide[ukb_targeted[risk_group == "low"], on = .(long_name, guidelines, strategy), UKB_targeted_low_cases := i.cases]
wide[ukb_targeted[risk_group == "low"], on = .(long_name, guidelines, strategy), UKB_targeted_low_controls := i.controls]

# Compute relevant public health statistics
wide[, UKB_high_risk_cases := UKB_blanket_high_cases + UKB_intermed_high_cases + ifelse(!is.na(UKB_targeted_high_cases), UKB_targeted_high_cases, 0)]
wide[, UKB_high_risk_cases_pct := UKB_high_risk_cases / UKB_cases]
wide[, UKB_high_risk_controls := UKB_blanket_high_controls + UKB_intermed_high_controls + ifelse(!is.na(UKB_targeted_high_controls), UKB_targeted_high_controls, 0)]
wide[, UKB_high_risk_controls_pct := UKB_high_risk_controls / UKB_controls]
wide[, UKB_cases_prevented := UKB_high_risk_cases * 0.2]
wide[, UKB_NNS := UKB_total / UKB_cases_prevented]
wide[, UKB_NNT := (UKB_high_risk_cases + UKB_high_risk_controls) / UKB_cases_prevented]

# Add in information about ONS hypothetical population
ons <- fread("analyses/public_health_modelling/UK_population_generalised/ONS_hypothetical_100k_pop.txt")
wide[, ONS_total := ons$N]
wide[, ONS_cases := ons$cases]
wide[, ONS_controls := ons$controls]

# Add in number of hypothetical 100K participants allocated to low, medium, and high risk groups in blanket screening
ons_blanket_alloc <- rbind(idcol="strategy",
  blanket = fread("analyses/public_health_modelling/UK_population_generalised/blanket_screening/ONS_biomarker_prs_risk_stratified.txt"),
  targeted = fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening/ONS_conv_rf_risk_stratified.txt"),
  targeted_above_PRS = fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening_above_PRS/ONS_conv_rf_prs_risk_stratified.txt")
)
wide[ons_blanket_alloc[risk_group == "high"], on = .(long_name, guidelines, strategy), ONS_blanket_high_cases := i.cases]
wide[ons_blanket_alloc[risk_group == "high"], on = .(long_name, guidelines, strategy), ONS_blanket_high_controls := i.controls]
wide[ons_blanket_alloc[risk_group == "medium"], on = .(long_name, guidelines, strategy), ONS_blanket_medium_cases := i.cases]
wide[ons_blanket_alloc[risk_group == "medium"], on = .(long_name, guidelines, strategy), ONS_blanket_medium_controls := i.controls]
wide[ons_blanket_alloc[risk_group == "low"], on = .(long_name, guidelines, strategy), ONS_blanket_low_cases := i.cases]
wide[ons_blanket_alloc[risk_group == "low"], on = .(long_name, guidelines, strategy), ONS_blanket_low_controls := i.controls]

# Add in number of hypothetical 100K participants reclassified based on LDL and diabetes
ons_intermed_treat <- rbind(idcol="strategy",
  blanket = fread("analyses/public_health_modelling/UK_population_generalised/blanket_screening/ONS_intermediate_risk_treat_alloc.txt"),
  targeted = fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening/ONS_intermediate_risk_treat_alloc.txt"),
  targeted_above_PRS = fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening_above_PRS/ONS_intermediate_risk_treat_alloc.txt")
)
wide[ons_intermed_treat[risk_group == "high"], on = .(long_name, guidelines, strategy), ONS_intermed_high_cases := i.cases]
wide[ons_intermed_treat[risk_group == "high"], on = .(long_name, guidelines, strategy), ONS_intermed_high_controls := i.controls]
wide[ons_intermed_treat[risk_group == "medium"], on = .(long_name, guidelines, strategy), ONS_intermed_medium_cases := i.cases]
wide[ons_intermed_treat[risk_group == "medium"], on = .(long_name, guidelines, strategy), ONS_intermed_medium_controls := i.controls]

# Add in number of hypothetical 100K participants reclassified in targeted assessment
ons_targeted <- rbind(idcol="strategy",
  targeted = fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening/ONS_targeted_biomarker_prs_risk_stratified.txt"),
  targeted_above_PRS = fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening_above_PRS/ONS_targeted_biomarker_risk_stratified.txt")
)
wide[ons_targeted[risk_group == "high"], on = .(long_name, guidelines, strategy), ONS_targeted_high_cases := i.cases]
wide[ons_targeted[risk_group == "high"], on = .(long_name, guidelines, strategy), ONS_targeted_high_controls := i.controls]
wide[ons_targeted[risk_group == "medium"], on = .(long_name, guidelines, strategy), ONS_targeted_medium_cases := i.cases]
wide[ons_targeted[risk_group == "medium"], on = .(long_name, guidelines, strategy), ONS_targeted_medium_controls := i.controls]
wide[ons_targeted[risk_group == "low"], on = .(long_name, guidelines, strategy), ONS_targeted_low_cases := i.cases]
wide[ons_targeted[risk_group == "low"], on = .(long_name, guidelines, strategy), ONS_targeted_low_controls := i.controls]

# Compute relevant public health statistics
wide[, ONS_high_risk_cases := ONS_blanket_high_cases + ONS_intermed_high_cases + ifelse(!is.na(ONS_targeted_high_cases), ONS_targeted_high_cases, 0)]
wide[, ONS_high_risk_cases_pct := ONS_high_risk_cases / ONS_cases]
wide[, ONS_high_risk_controls := ONS_blanket_high_controls + ONS_intermed_high_controls + ifelse(!is.na(ONS_targeted_high_controls), ONS_targeted_high_controls, 0)]
wide[, ONS_high_risk_controls_pct := ONS_high_risk_controls / ONS_controls]
wide[, ONS_cases_prevented := ONS_high_risk_cases * 0.2]
wide[, ONS_NNS := ONS_total / ONS_cases_prevented]
wide[, ONS_NNT := (ONS_high_risk_cases + ONS_high_risk_controls) / ONS_cases_prevented]

# Write out
fwrite(wide, sep="\t", quote=FALSE, file="analyses/public_health_modelling/UK_population_generalised/collated_public_health_statistics.txt")
