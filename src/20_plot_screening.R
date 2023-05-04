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

# -------------------------------

# Create table for supp
supp <- rbind(
  phs[facet_group == "Blanket screening", .(screening_step1 = row_name, screening_step3 = "")],
  phs[facet_group != "Blanket screening", .(screening_step1 = facet_group, screening_step3 = row_name)]
)

# How many UKB participants were allocated to each risk group
ukb_step1 <- fread("analyses/public_health_modelling/UK_population_generalised/blanket_screening/UKB_screening_step_1.txt")
ukb_step1[, short_name := fcase(
  name == "Conventional RF" & !(PGS), "Conventional RF",
  name == "Conventional RF" & (PGS), "Conventional RF + PRS",
  name == "NMR" & !(PGS), "Conventional RF + NMR",
  name == "NMR" & (PGS), "Conventional RF + NMR + PRS"
)]

supp[ukb_step1[risk_group == "high"], on = .(screening_step1 = short_name), UKB_step1_high_cases := i.cases]
supp[ukb_step1[risk_group == "high"], on = .(screening_step1 = short_name), UKB_step1_high_controls := i.controls]
supp[ukb_step1[risk_group == "medium"], on = .(screening_step1 = short_name), UKB_step1_medium_cases := i.cases]
supp[ukb_step1[risk_group == "medium"], on = .(screening_step1 = short_name), UKB_step1_medium_controls := i.controls]
supp[ukb_step1[risk_group == "low"], on = .(screening_step1 = short_name), UKB_step1_low_cases := i.cases]
supp[ukb_step1[risk_group == "low"], on = .(screening_step1 = short_name), UKB_step1_low_controls := i.controls]

ukb_step2 <- fread("analyses/public_health_modelling/UK_population_generalised/blanket_screening/UKB_screening_step_2.txt")
ukb_step2[, short_name := fcase(
  name == "Conventional RF" & !(PGS), "Conventional RF",
  name == "Conventional RF" & (PGS), "Conventional RF + PRS",
  name == "NMR" & !(PGS), "Conventional RF + NMR",
  name == "NMR" & (PGS), "Conventional RF + NMR + PRS"
)]

supp[ukb_step2[risk_group2 == "high"], on = .(screening_step1 = short_name), UKB_step2_high_cases := i.cases]
supp[ukb_step2[risk_group2 == "high"], on = .(screening_step1 = short_name), UKB_step2_high_controls := i.controls]
supp[ukb_step2[risk_group2 == "medium"], on = .(screening_step1 = short_name), UKB_step2_medium_cases := i.cases]
supp[ukb_step2[risk_group2 == "medium"], on = .(screening_step1 = short_name), UKB_step2_medium_controls := i.controls]

ukb_step3 <- fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening/UKB_screening_step_3.txt")

supp[ukb_step3[risk_group3 == "high"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening), UKB_step3_high_cases := i.cases]
supp[ukb_step3[risk_group3 == "high"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening), UKB_step3_high_controls := i.controls]
supp[ukb_step3[risk_group3 == "medium"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening), UKB_step3_medium_cases := i.cases]
supp[ukb_step3[risk_group3 == "medium"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening), UKB_step3_medium_controls := i.controls]
supp[ukb_step3[risk_group3 == "low"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening), UKB_step3_low_cases := i.cases]
supp[ukb_step3[risk_group3 == "low"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening), UKB_step3_low_controls := i.controls]

ukb_totals <- fread("analyses/public_health_modelling/sample_flowchart.txt")

supp[screening_step3 == "", UKB_high_cases := UKB_step1_high_cases + UKB_step2_high_cases]
supp[screening_step3 != "", UKB_high_cases := UKB_step1_high_cases + UKB_step2_high_cases + UKB_step3_high_cases]
supp[, UKB_high_cases_pct := UKB_high_cases / ukb_totals[.N, cases]]

supp[screening_step3 == "", UKB_high_controls := UKB_step1_high_controls + UKB_step2_high_controls]
supp[screening_step3 != "", UKB_high_controls := UKB_step1_high_controls + UKB_step2_high_controls + UKB_step3_high_controls]
supp[, UKB_high_controls_pct := UKB_high_controls / (ukb_totals[.N, samples] - ukb_totals[.N, cases])]

# Compute public health statistics
supp[, UKB_events_prevented := floor(UKB_high_cases / 5)]
supp[, UKB_NNS := floor(ukb_totals[.N, samples] / UKB_events_prevented)]
supp[, UKB_NNT := floor((UKB_high_cases + UKB_high_controls) / UKB_events_prevented)]

# add in same information generalized to ONS
ons_step1 <- fread("analyses/public_health_modelling/UK_population_generalised/blanket_screening/ONS_screening_step_1.txt")
ons_step1[, short_name := fcase(
  name == "Conventional RF" & !(PGS), "Conventional RF",
  name == "Conventional RF" & (PGS), "Conventional RF + PRS",
  name == "NMR" & !(PGS), "Conventional RF + NMR",
  name == "NMR" & (PGS), "Conventional RF + NMR + PRS"
)]

supp[ons_step1[risk_group == "high"], on = .(screening_step1 = short_name), ONS_step1_high_cases := i.cases]
supp[ons_step1[risk_group == "high"], on = .(screening_step1 = short_name), ONS_step1_high_controls := i.controls]
supp[ons_step1[risk_group == "medium"], on = .(screening_step1 = short_name), ONS_step1_medium_cases := i.cases]
supp[ons_step1[risk_group == "medium"], on = .(screening_step1 = short_name), ONS_step1_medium_controls := i.controls]
supp[ons_step1[risk_group == "low"], on = .(screening_step1 = short_name), ONS_step1_low_cases := i.cases]
supp[ons_step1[risk_group == "low"], on = .(screening_step1 = short_name), ONS_step1_low_controls := i.controls]

ons_step2 <- fread("analyses/public_health_modelling/UK_population_generalised/blanket_screening/ONS_screening_step_2_stratified.txt")
ons_step2[, short_name := fcase(
  name == "Conventional RF" & !(PGS), "Conventional RF",
  name == "Conventional RF" & (PGS), "Conventional RF + PRS",
  name == "NMR" & !(PGS), "Conventional RF + NMR",
  name == "NMR" & (PGS), "Conventional RF + NMR + PRS"
)]

supp[ons_step2[risk_group2 == "high"], on = .(screening_step1 = short_name), ONS_step2_high_cases := i.cases]
supp[ons_step2[risk_group2 == "high"], on = .(screening_step1 = short_name), ONS_step2_high_controls := i.controls]
supp[ons_step2[risk_group2 == "medium"], on = .(screening_step1 = short_name), ONS_step2_medium_cases := i.cases]
supp[ons_step2[risk_group2 == "medium"], on = .(screening_step1 = short_name), ONS_step2_medium_controls := i.controls]

ons_step3 <- fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening/ONS_screening_step_3_stratified.txt")

supp[ons_step3[risk_group3 == "high"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening), ONS_step3_high_cases := i.cases]
supp[ons_step3[risk_group3 == "high"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening), ONS_step3_high_controls := i.controls]
supp[ons_step3[risk_group3 == "medium"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening), ONS_step3_medium_cases := i.cases]
supp[ons_step3[risk_group3 == "medium"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening), ONS_step3_medium_controls := i.controls]
supp[ons_step3[risk_group3 == "low"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening), ONS_step3_low_cases := i.cases]
supp[ons_step3[risk_group3 == "low"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening), ONS_step3_low_controls := i.controls]

ons_totals <- fread("analyses/public_health_modelling/UK_population_generalised/ONS_hypothetical_100k_pop.txt")

supp[screening_step3 == "", ONS_high_cases := ONS_step1_high_cases + ONS_step2_high_cases]
supp[screening_step3 != "", ONS_high_cases := ONS_step1_high_cases + ONS_step2_high_cases + ONS_step3_high_cases]
supp[, ONS_high_cases_pct := ONS_high_cases / ons_totals[, cases]]

supp[screening_step3 == "", ONS_high_controls := ONS_step1_high_controls + ONS_step2_high_controls]
supp[screening_step3 != "", ONS_high_controls := ONS_step1_high_controls + ONS_step2_high_controls + ONS_step3_high_controls]
supp[, ONS_high_controls_pct := ONS_high_controls / (ons_totals[, N] - ons_totals[, cases])]

supp[, ONS_events_prevented := floor(ONS_high_cases / 5)]
supp[, ONS_NNS := floor(ons_totals[, N] / ONS_events_prevented)]
supp[, ONS_NNT := floor((ONS_high_cases + ONS_high_controls) / ONS_events_prevented)]

# Write out
fwrite(supp, sep="\t", quote=FALSE, file="analyses/public_health_modelling/UK_population_generalised/collated_public_health_statistics.txt")

# -------------------------------

# Stratify by age and sex
supp <- foreach(this_sex = c("Male", "Female"), .combine=rbind) %do% {
  foreach(this_age_group = sprintf("%s-%s", seq(40, 65, by=5), seq(45, 70, by=5)), .combine=rbind) %do% {
    cbind(sex=this_sex, age_group=this_age_group, 
      rbind(
        phs[facet_group == "Blanket screening", .(screening_step1 = row_name, screening_step3 = "")],
        phs[facet_group != "Blanket screening", .(screening_step1 = facet_group, screening_step3 = row_name)]
      ))
  }
}

# How many UKB participants were allocated to each risk group
ukb_step1 <- fread("analyses/public_health_modelling/UK_population_generalised/blanket_screening/UKB_screening_step_1_stratified_by_age_sex.txt")
ukb_step1[, short_name := fcase(
  name == "Conventional RF" & !(PGS), "Conventional RF",
  name == "Conventional RF" & (PGS), "Conventional RF + PRS",
  name == "NMR" & !(PGS), "Conventional RF + NMR",
  name == "NMR" & (PGS), "Conventional RF + NMR + PRS"
)]
ukb_step1[, age_group := sprintf("%s-%s", age_group, age_group + 5)]

supp[ukb_step1[risk_group == "high"], on = .(screening_step1 = short_name, age_group, sex), UKB_step1_high_cases := i.cases]
supp[ukb_step1[risk_group == "high"], on = .(screening_step1 = short_name, age_group, sex), UKB_step1_high_controls := i.controls]
supp[ukb_step1[risk_group == "medium"], on = .(screening_step1 = short_name, age_group, sex), UKB_step1_medium_cases := i.cases]
supp[ukb_step1[risk_group == "medium"], on = .(screening_step1 = short_name, age_group, sex), UKB_step1_medium_controls := i.controls]
supp[ukb_step1[risk_group == "low"], on = .(screening_step1 = short_name, age_group, sex), UKB_step1_low_cases := i.cases]
supp[ukb_step1[risk_group == "low"], on = .(screening_step1 = short_name, age_group, sex), UKB_step1_low_controls := i.controls]

ukb_step2 <- fread("analyses/public_health_modelling/UK_population_generalised/blanket_screening/UKB_screening_step_2_stratified_by_age_sex.txt")
ukb_step2[, short_name := fcase(
  name == "Conventional RF" & !(PGS), "Conventional RF",
  name == "Conventional RF" & (PGS), "Conventional RF + PRS",
  name == "NMR" & !(PGS), "Conventional RF + NMR",
  name == "NMR" & (PGS), "Conventional RF + NMR + PRS"
)]
ukb_step2[, age_group := sprintf("%s-%s", age_group, age_group + 5)]

supp[ukb_step2[risk_group2 == "high"], on = .(screening_step1 = short_name, age_group, sex), UKB_step2_high_cases := i.cases]
supp[ukb_step2[risk_group2 == "high"], on = .(screening_step1 = short_name, age_group, sex), UKB_step2_high_controls := i.controls]
supp[ukb_step2[risk_group2 == "medium"], on = .(screening_step1 = short_name, age_group, sex), UKB_step2_medium_cases := i.cases]
supp[ukb_step2[risk_group2 == "medium"], on = .(screening_step1 = short_name, age_group, sex), UKB_step2_medium_controls := i.controls]

ukb_step3 <- fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening/UKB_screening_step_3_stratified_by_age_sex.txt")
ukb_step3[, age_group := sprintf("%s-%s", age_group, age_group + 5)]

supp[ukb_step3[risk_group3 == "high"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening, age_group, sex), UKB_step3_high_cases := i.cases]
supp[ukb_step3[risk_group3 == "high"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening, age_group, sex), UKB_step3_high_controls := i.controls]
supp[ukb_step3[risk_group3 == "medium"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening, age_group, sex), UKB_step3_medium_cases := i.cases]
supp[ukb_step3[risk_group3 == "medium"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening, age_group, sex), UKB_step3_medium_controls := i.controls]
supp[ukb_step3[risk_group3 == "low"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening, age_group, sex), UKB_step3_low_cases := i.cases]
supp[ukb_step3[risk_group3 == "low"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening, age_group, sex), UKB_step3_low_controls := i.controls]

# add in same information generalized to ONS
ons_step1 <- fread("analyses/public_health_modelling/UK_population_generalised/blanket_screening/ONS_screening_step_1_stratified_by_age_sex.txt")
ons_step1[, short_name := fcase(
  name == "Conventional RF" & !(PGS), "Conventional RF",
  name == "Conventional RF" & (PGS), "Conventional RF + PRS",
  name == "NMR" & !(PGS), "Conventional RF + NMR",
  name == "NMR" & (PGS), "Conventional RF + NMR + PRS"
)]
ons_step1[, age_group := sprintf("%s-%s", age_group, age_group + 5)]

supp[ons_step1[risk_group == "high"], on = .(screening_step1 = short_name, age_group, sex), ONS_step1_high_cases := i.cases]
supp[ons_step1[risk_group == "high"], on = .(screening_step1 = short_name, age_group, sex), ONS_step1_high_controls := i.controls]
supp[ons_step1[risk_group == "medium"], on = .(screening_step1 = short_name, age_group, sex), ONS_step1_medium_cases := i.cases]
supp[ons_step1[risk_group == "medium"], on = .(screening_step1 = short_name, age_group, sex), ONS_step1_medium_controls := i.controls]
supp[ons_step1[risk_group == "low"], on = .(screening_step1 = short_name, age_group, sex), ONS_step1_low_cases := i.cases]
supp[ons_step1[risk_group == "low"], on = .(screening_step1 = short_name, age_group, sex), ONS_step1_low_controls := i.controls]

ons_step2 <- fread("analyses/public_health_modelling/UK_population_generalised/blanket_screening/ONS_screening_step_2_stratified_by_age_sex.txt")
ons_step2[, short_name := fcase(
  name == "Conventional RF" & !(PGS), "Conventional RF",
  name == "Conventional RF" & (PGS), "Conventional RF + PRS",
  name == "NMR" & !(PGS), "Conventional RF + NMR",
  name == "NMR" & (PGS), "Conventional RF + NMR + PRS"
)]
ons_step2[, age_group := sprintf("%s-%s", age_group, age_group + 5)]

supp[ons_step2[risk_group2 == "high"], on = .(screening_step1 = short_name, age_group, sex), ONS_step2_high_cases := i.cases]
supp[ons_step2[risk_group2 == "high"], on = .(screening_step1 = short_name, age_group, sex), ONS_step2_high_controls := i.controls]
supp[ons_step2[risk_group2 == "medium"], on = .(screening_step1 = short_name, age_group, sex), ONS_step2_medium_cases := i.cases]
supp[ons_step2[risk_group2 == "medium"], on = .(screening_step1 = short_name, age_group, sex), ONS_step2_medium_controls := i.controls]

ons_step3 <- fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening/ONS_screening_step_3_stratified_by_age_sex.txt")
ons_step3[, age_group := sprintf("%s-%s", age_group, age_group + 5)]

supp[ons_step3[risk_group3 == "high"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening, age_group, sex), ONS_step3_high_cases := i.cases]
supp[ons_step3[risk_group3 == "high"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening, age_group, sex), ONS_step3_high_controls := i.controls]
supp[ons_step3[risk_group3 == "medium"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening, age_group, sex), ONS_step3_medium_cases := i.cases]
supp[ons_step3[risk_group3 == "medium"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening, age_group, sex), ONS_step3_medium_controls := i.controls]
supp[ons_step3[risk_group3 == "low"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening, age_group, sex), ONS_step3_low_cases := i.cases]
supp[ons_step3[risk_group3 == "low"], on = .(screening_step1 = blanket_screening, screening_step3 = targeted_screening, age_group, sex), ONS_step3_low_controls := i.controls]

# order rows
ro <- rbind(
  phs[facet_group == "Blanket screening", .(screening_step1 = row_name, screening_step3 = "")],
  phs[facet_group != "Blanket screening", .(screening_step1 = facet_group, screening_step3 = row_name)]
)
supp <- supp[ro, on = .(screening_step1, screening_step3)]

# Write out
fwrite(supp, sep="\t", quote=FALSE, file="analyses/public_health_modelling/UK_population_generalised/collated_risk_strata_by_age_sex.txt")


