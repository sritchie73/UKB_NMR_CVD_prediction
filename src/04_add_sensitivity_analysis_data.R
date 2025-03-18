library(data.table)
source("src/utils/QRISK3.R")

# Load analysis cohort for filtering
dat <- fread("data/cleaned/analysis_cohort.txt", select="eid")

# Load in relevant anthropometric data for QRISK3
anthro <- fread("data/ukb/anthropometrics/output/anthropometrics.txt")
anthro <- anthro[visit_index == 0, .(eid, assessment_date, age, sex, height, weight, townsend, ethnicity_group, ethnicity_subgroup)] 
dat <- dat[anthro, on = .(eid), nomatch=0]

# Code ethnicity following QRISK3 algorithm
dat[, QRisk_ethnicity := fcase(
  ethnicity_group %in% c("White", "Prefer not to answer", ""), "White or not stated",
  ethnicity_subgroup == "Indian", "Indian",
  ethnicity_subgroup == "Pakistani", "Pakistani", 
  ethnicity_subgroup == "Bangladeshi", "Bangladeshi", 
  ethnicity_subgroup %in% c("Asian or Asian British", "Any other Asian background"), "Other Asian",
  ethnicity_subgroup == "Caribbean", "Black Caribbean",
  ethnicity_subgroup == "African", "Black African",
  ethnicity_group == "Chinese", "Chinese", 
  default="Other ethnic group")
]

# Add in SBP and SD of SBP
sbp <- fread("data/ukb/blood_pressure/output/blood_pressure.txt")
sbp <- sbp[visit_index == 0, .(eid, sbp, sd_sbp)]
dat <- merge(dat, sbp, by="eid", all.x=TRUE)

# Add in relevant biochemistry biomarkers
bio <- fread("data/ukb/biomarkers/output/biomarkers.txt")
bio <- bio[visit_index == 0L, .(eid, hdl, tchol)]
dat <- merge(dat, bio, by="eid", all.x=TRUE)

# Add in ratio of tchol:hdl
dat[, tchol_hdl_ratio := tchol / hdl]
dat[tchol_hdl_ratio < 1, tchol_hdl_ratio := 1]
dat[tchol_hdl_ratio > 11, tchol_hdl_ratio := 11]

# Add in smoking status
smoking <- fread("data/ukb/smoking/output/smoking_status.txt")
smoking <- smoking[visit_index == 0, .(eid, smoking_status, daily_cigarettes)]
dat <- merge(dat, smoking, by="eid", all.x=TRUE)

# Code smoking variable following QRISK3 algorithm
dat[, smoke := fcase(
  smoking_status == "Previous", "ex-smoker",
  smoking_status == "Current" & daily_cigarettes < 10, "light smoker",
  smoking_status == "Current" & daily_cigarettes >= 10 & daily_cigarettes < 20, "moderate smoker", 
  smoking_status == "Current" & daily_cigarettes >= 20, "heavy smoker", 
  smoking_status == "Current" & is.na(daily_cigarettes), "moderate smoker",
  default="non-smoker"
)]

# Add in medications
meds_touchscreen <- fread("data/ukb/medication/output/medications_simple.txt")
meds_interview <- fread("data/ukb/medication/output/detailed_medications_summarised.txt")
meds <- merge(meds_touchscreen, meds_interview, by=c("eid", "visit_index"), suffixes=c("_touchscreen", "_interview"), all=TRUE)
meds <- meds[visit_index == 0]
meds[, blood_pressure_treatment := blood_pressure_medication | hypertension_medication] # Touchsreen OR curated list of antihypertensive from nurse interview

dat[meds, on = .(eid), blood_pressure_treatment := i.blood_pressure_treatment]
dat[is.na(blood_pressure_treatment), blood_pressure_treatment := FALSE]
 
dat[meds, on = .(eid), atypical_antipsychotics := i.atypical_antipsychotics]
dat[is.na(atypical_antipsychotics), atypical_antipsychotics := FALSE]

dat[meds, on = .(eid), systematic_corticosteroids := i.systematic_corticosteroids]
dat[is.na(systematic_corticosteroids), systematic_corticosteroids := FALSE]

# Add in diabetes status (Eastwood 2016 algorithm)
diabetes <- fread("data/ukb/Eastwood_diabetes/output/prevalent_diabetes.txt")
dat[diabetes[visit_index == 0], on = .(eid), prevalent_t1d := fcase(
  adjudicated_diabetes %in% c("Probable type 1 diabetes",  "Possible type 1 diabetes"), TRUE, 
  default=FALSE)]
dat[diabetes[visit_index == 0], on = .(eid), prevalent_t2d := fcase(
  adjudicated_diabetes %in% c("Probable type 2 diabetes",  "Possible type 2 diabetes"), TRUE, 
  default=FALSE)]

# Add in erectile dysfunction (diagnosis or treatment)
ed <- fread("data/ukb/endpoints/endpoints/prevalent_erectile_dysfunction/events_and_followup.txt")
ed <- ed[visit_index == 0, .(eid, erectile_dysfunction=prevalent_event)]
ed[meds, on = .(eid), erectile_dysfunction := erectile_dysfunction | i.erectile_dysfunction] # Diagnosis or treatment
ed[is.na(erectile_dysfunction), erectile_dysfunction := FALSE]
dat <- merge(dat, ed, by="eid", all.x=TRUE)

# Other prevalent diseases
afib <- fread("data/ukb/endpoints/endpoints/prevalent_afib/events_and_followup.txt")
afib <- afib[visit_index == 0, .(eid, atrial_fibrillation=prevalent_event)]
afib[is.na(atrial_fibrillation), atrial_fibrillation := FALSE]
dat <- merge(dat, afib, by="eid", all.x=TRUE)

ckd <- fread("data/ukb/endpoints/endpoints/prevalent_CKD/events_and_followup.txt")
ckd <- ckd[visit_index == 0, .(eid, chronic_kidney_disease=prevalent_event)]
ckd[is.na(chronic_kidney_disease), chronic_kidney_disease := FALSE]
dat <- merge(dat, ckd, by="eid", all.x=TRUE)

migr <- fread("data/ukb/endpoints/endpoints/prevalent_migraine/events_and_followup.txt")
migr <- migr[visit_index == 0, .(eid, migraine=prevalent_event)]
migr[is.na(migraine), migraine := FALSE]
dat <- merge(dat, migr, by="eid", all.x=TRUE)

ra <- fread("data/ukb/endpoints/endpoints/prevalent_rheumatoid_arthritis/events_and_followup.txt")
ra <- ra[visit_index == 0, .(eid, rheumatoid_arthritis=prevalent_event)]
ra[is.na(rheumatoid_arthritis), rheumatoid_arthritis := FALSE]
dat <- merge(dat, ra, by="eid", all.x=TRUE)

sle <- fread("data/ukb/endpoints/endpoints/prevalent_SLE/events_and_followup.txt")
sle <- sle[visit_index == 0, .(eid, systemic_lupus_erythematosis=prevalent_event)]
sle[is.na(systemic_lupus_erythematosis), systemic_lupus_erythematosis := FALSE]
dat <- merge(dat, sle, by="eid", all.x=TRUE)

smi <- fread("data/ukb/endpoints/endpoints/prevalent_mental_illness/events_and_followup.txt")
smi <- smi[visit_index == 0, .(eid, severe_mental_illness=prevalent_event)]
smi[is.na(severe_mental_illness), severe_mental_illness := FALSE]
dat <- merge(dat, smi, by="eid", all.x=TRUE)

# Compute QRISK3 linear predictor
dat[, QRISK3 := qrisk3(
  sex, age, QRisk_ethnicity, townsend, smoke, sbp, sd_sbp, weight, height, tchol_hdl_ratio, FALSE,
  atrial_fibrillation, chronic_kidney_disease, erectile_dysfunction, severe_mental_illness, migraine, rheumatoid_arthritis, 
  systemic_lupus_erythematosis, prevalent_t1d, prevalent_t2d, atypical_antipsychotics, blood_pressure_treatment, systematic_corticosteroids,
  type="linear predictor"
)]

# Add in alternate CVD endpoint (CEU definition)
cvd <- fread("data/ukb/endpoints/endpoints/CEU_CVD_10yr/events_and_followup.txt")
cvd <- cvd[visit_index == 0]
dat[cvd, on = .(eid), incident_cvd2 := i.incident_event]
dat[cvd, on = .(eid), incident_cvd2_followup := i.incident_event_followup]
dat[cvd, on = .(eid), incident_cvd2_followup_date := i.incident_event_followup_date]
dat[cvd, on = .(eid), incident_cvd2_is_fatal := fcase(
  i.incident_event_type == "death", TRUE,
  i.incident_event_type == "hospitalisation", FALSE,
  i.incident_event_type == "" & !is.na(i.incident_event), FALSE,
  is.na(i.incident_event), NA
)]
dat[cvd, on = . (eid), cvd2_is_primary_cause := fcase(
  i.incident_cause_type == "primary", TRUE,
  i.incident_cause_type == "secondary", FALSE,
  i.incident_cause_type == "" & !is.na(i.incident_cause_type), FALSE,
  is.na(i.incident_cause_type), NA
)]

# Write out extended additional information
fwrite(dat, sep="\t", quote=FALSE, file="data/cleaned/sensitivity_analysis_extended_data.txt")

# subset to specific columns of interest and add to main dataset
sens <- dat[,.(eid, QRISK3, incident_cvd2, incident_cvd2_followup)]
dat <- fread("data/cleaned/analysis_cohort.txt")
dat <- merge(dat, sens, by="eid", all.x=TRUE)
fwrite(dat, sep="\t", quote=FALSE, file="data/cleaned/analysis_cohort.txt")

