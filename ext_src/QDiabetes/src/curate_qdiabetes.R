library(data.table)
library(QDiabetes)
library(lubridate)

# Load curated anthropometrics data and filter to relevant fields
anthro <- fread("data/curated/anthropometrics/anthropometrics.txt")

# Extract self-reported sex
sex <- anthro[, .(eid, visit_index, sex)]

# Extract age
age <- anthro[, .(eid, visit_index, age, age_decimal)]

# Extract BMI
bmi <- anthro[, .(eid, visit_index, bmi)]

# Extract height
height <- anthro[, .(eid, visit_index, height=height/100)] # QDiabetes package expects meters, not centimeters

# Extract weight
weight <- anthro[, .(eid, visit_index, weight)]

# Extract Townsend deprivation index
townsend <- anthro[, .(eid, visit_index, townsend)]

# Extract ethnicity
ethnicity <- anthro[,.(eid, visit_index, ethnicity=ethnicity_subgroup)]

# Recode ethnicity to groupings expected by QDiabetes package
ethnicity[, ethnicity := fcase(
  ethnicity == "White", "WhiteNA",
  ethnicity == "British", "WhiteNA",
  ethnicity == "Irish", "WhiteNA",
  ethnicity == "Any other white background", "WhiteNA",
  ethnicity == "Do not know", "WhiteNA",
  ethnicity == "Prefer not to answer", "WhiteNA",
  ethnicity == "", "WhiteNA",
  ethnicity == "Indian", "Indian",
  ethnicity == "Pakistani", "Pakistani",
  ethnicity == "Bangladeshi", "Bangladeshi",
  ethnicity == "Chinese", "Chinese",
  ethnicity == "Asian or Asian British", "OtherAsian",
  ethnicity == "Any other Asian background", "OtherAsian", 
  ethnicity == "African", "BlackAfrican",
  ethnicity == "Caribbean", "BlackCaribbean",
  ethnicity == "Mixed", "Other",
  ethnicity == "White and Black Caribbean", "Other",
  ethnicity == "White and Black African", "Other",
  ethnicity == "White and Asian", "Other",
  ethnicity == "Any other mixed background", "Other",
  ethnicity == "Black or Black British", "Other",
  ethnicity == "Any other Black background", "Other",
  ethnicity == "Other ethnic group", "Other",
  default = "Coding Error"
)]

# Load curated smoking status data
smoking <- fread("data/curated/smoking/smoking_status.txt")

# Curate categories as per QDiabetes package
smoking <- smoking[, .(eid, visit_index, smoking_status = fcase(
  smoking_status == "Never", "Non",
  smoking_status == "Previous", "Ex",
  smoking_status == "Current" & daily_cigarettes < 10, "Light",
  smoking_status == "Current" & daily_cigarettes >= 10 & daily_cigarettes < 20, "Moderate",
  smoking_status == "Current" & daily_cigarettes >= 20, "Heavy",
  smoking_status == "Current" & is.na(daily_cigarettes), "Moderate",
  smoking_status == "Prefer not to answer", NA_character_,
  default = "Coding Error"
))]

# Load family history of diabetes
family_history <- fread("data/curated/family_history/illness_of_first_degree_relatives.txt")
family_history_diabetes <- family_history[, .(eid, visit_index, family_history_diabetes=diabetes)]

# Load history of cardiovascular disease
history_cvd <- fread("data/adjudicated_medical_history/CVD/events_and_followup.txt")
history_cvd <- history_cvd[,.(eid, visit_index, history_cvd=prevalent_event)]

# Load history of learning difficulties - this is purely based on incidental codings in
# hospital records (i.e. no history of such was purposely collected by UKB). 
history_learning_difficulties <- fread("data/adjudicated_medical_history/learning_difficulties/events_and_followup.txt")
history_learning_difficulties <- history_learning_difficulties[,.(eid, visit_index, history_learning_difficulties=prevalent_event)]

# Load history of PCOS
history_pcos <- fread("data/adjudicated_medical_history/PCOS/events_and_followup.txt")
history_pcos[sex == "Male" & is.na(prevalent_event), prevalent_event := FALSE] # No ovaries; no PCOS
history_pcos <- history_pcos[,.(eid, visit_index, history_pcos=prevalent_event)]
 
# Load history of bipolar or schizophrenia disorders
history_bipolar_schizophrenia <- fread("data/adjudicated_medical_history/bipolar_schizophrenia/events_and_followup.txt")
history_bipolar_schizophrenia <- history_bipolar_schizophrenia[,.(eid, visit_index, history_bipolar_schizophrenia=prevalent_event)]

# Next, curate medication history
medications <- fread("data/curated/medication/detailed_medications_field_20003.txt")
medications_summarised <- fread("data/curated/medication/detailed_medications_summarised.txt")
med_survey <- fread("data/curated/medication/medications_simple.txt")

# History of statins (or other lipid lowering medications): 
# Either answering "Yes" to the touchscreen survey question on cholesterol medication, or
# having one of the lipid lowering medications in the list of medications coded by trained 
# nurse at biobank assessment (see README.txt for full list)
lipid_lowering_medication <- med_survey[medications_summarised, on = .(eid, visit_index), 
  .(eid, visit_index, lipid_lowering_medication = lipid_lowering_medication | cholesterol_medication)]

# History of treated hypertension
# Either answering "Yes" to touchscreen survey question on blood pressure medication, or
# having one of the hypertension medicaitons in the list of medications coded by trained
# nurse at biobank assessment (see README.txt for full list)
hypertension_medication <- med_survey[medications_summarised, on = .(eid, visit_index),
  .(eid, visit_index, hypertension_medication = hypertension_medication | blood_pressure_medication)]

# Curate data on second generation atypical antipsychotics
# List of medications obtained from QDiabetes 2018 paper
atypical_antipsychotics <- medications[,.(atypical_antipsychotics = ifelse(
  any(medication_code == 1141153490) |  # amisulpride
  any(medication_code == 1141195974) |  # aripiprazole
  any(medication_code == 1140867420) |  # clozapine
  any(medication_code == 1140928916) |  # olanzapine
  any(medication_code == 1141152848) |  # quetiapine
  any(medication_code == 1140867444) |  # risperidone
  any(medication_code == 1140927956) |  # sertindole
  any(medication_code == 1141169714),   # zotepine
  TRUE, FALSE)), by=.(eid, visit_index)]


# Curate data on systematic corticosteroid use (oral or injected)
#
# This list is my best-effort to curate the list that most closely matches what is listed 
# in the QDiabetes paper (i.e. those listed in the British National Formulary chapter 
# 6.3.2 or otherwise directly mentioned). 
#
# Here, history of corticosteroid usage is I believe is intended to capture long-term use
# of (strong) corticosteroids (oral or injected) not transient (e.g. over the counter topical
# treatments) or preventive inhalers for asthma. 
systematic_corticosteroids <- medications[,.(systematic_corticosteroids = ifelse(
  any(medication_code == 1140874790) |  # betamethasone
  any(medication_code == 1141145782) |  # deflazacort
  any(medication_code == 1140874816) |  # dexamethasone
  any(medication_code == 1140874896) |  # hydrocortisone
  any(medication_code == 1140874976) |  # methylprednisolone
  any(medication_code == 1140874930) |  # prednisolone
  any(medication_code == 1141157402) |  # prednisolone product
  any(medication_code == 1140868364) |  # prednisone
  any(medication_code == 1140868426),   # triamcinolone
  TRUE, FALSE)), by=.(eid, visit_index)]

# Load biomarkers
biomarkers <- fread("data/curated/biomarkers/biomarkers.txt")
biomarkers <- biomarkers[,.(eid, visit_index, hba1c, fasting_glucose, non_fasting_glucose=glucose)]

# Also load fasting time information
fasting <- fread("data/curated/biomarkers/fasting_time_hours.txt")

# Load history of diabetes adjudicted from self-report data using the Eastwood et al. algorithms
diabetes_self_report <- fread("data/adjudicated_medical_history/Eastwood_diabetes/EastwoodDiabetes_PrevalentDiabetes.csv")
diabetes_self_report <- diabetes_self_report[,.(eid=V1, visit_index=0, adjudicated_diabetes=diabetes_EastwoodAdjudicated)]
diabetes_self_report <- diabetes_self_report[,.(eid, visit_index, adjudicated_diabetes)]

# Load history of different diabetes types from retrospective hospital records
t1d_hes <- fread("data/adjudicated_medical_history/type_1_diabetes/events_and_followup.txt")
t2d_hes <- fread("data/adjudicated_medical_history/type_2_diabetes/events_and_followup.txt")
gd_hes <- fread("data/adjudicated_medical_history/gestational_diabetes/events_and_followup.txt")
ud_hes <- fread("data/adjudicated_medical_history/unspecified_diabetes/events_and_followup.txt")
udp_hes <- fread("data/adjudicated_medical_history/unspecified_diabetes_in_pregnancy/events_and_followup.txt")
dc_hes <- fread("data/adjudicated_medical_history/diabetes_complications/events_and_followup.txt")
mn_hes <- fread("data/adjudicated_medical_history/malnutrition_diabetes/events_and_followup.txt")

gd_hes[sex == "Male", prevalent_event := FALSE] # Male, no gestational diabetes
udp_hes[sex == "Male", prevalent_event := FALSE] # Male, no diabetes in pregnancy

t1d_hes <- t1d_hes[,.(eid, visit_index, type_1_diabetes=prevalent_event)]
t2d_hes <- t2d_hes[,.(eid, visit_index, type_2_diabetes=prevalent_event)]
gd_hes <- gd_hes[,.(eid, visit_index, gestational_diabetes=prevalent_event)]
ud_hes <- ud_hes[,.(eid, visit_index, unspecified_diabetes=prevalent_event)]
udp_hes <- udp_hes[,.(eid, visit_index, unspecified_diabetes_in_pregnancy=prevalent_event)]
dc_hes <- dc_hes[,.(eid, visit_index, diabetes_complications=prevalent_event)]
mn_hes <- mn_hes[,.(eid, visit_index, malnutrition_diabetes=prevalent_event)]

diabetes_hes <- t1d_hes[t2d_hes, on = .(eid, visit_index)]
diabetes_hes <- diabetes_hes[gd_hes, on = .(eid, visit_index)]
diabetes_hes <- diabetes_hes[ud_hes, on = .(eid, visit_index)]
diabetes_hes <- diabetes_hes[udp_hes, on = .(eid, visit_index)]
diabetes_hes <- diabetes_hes[dc_hes, on = .(eid, visit_index)]
diabetes_hes <- diabetes_hes[mn_hes, on = .(eid, visit_index)]

# Adjudicate diabetes type - follow flowchart guidance in Eastwood et al.
diabetes_hes[, no_diabetes := !(type_1_diabetes) & !(type_2_diabetes) & !(gestational_diabetes) &
  !(unspecified_diabetes) & !(unspecified_diabetes_in_pregnancy) & !(diabetes_complications) & !(malnutrition_diabetes)]

diabetes_hes[(unspecified_diabetes_in_pregnancy) & !(gestational_diabetes), gestational_diabetes := NA]

diabetes_hes[, adjudication := fcase(
  is.na(no_diabetes), "Withdrawn consent for hospital record linkage",
  (no_diabetes), "Diabetes unlikely",
  # Gestational diabetes if and only if no record also of type 1, type 2, or otherwise unspecified diabetes
  (gestational_diabetes) & !(type_1_diabetes) & !(type_2_diabetes) & !(unspecified_diabetes), "Possible gestational diabetes",
  # Type 2 diabetes only if no record of type 1 or otherwise unspecified diabetes
  (type_2_diabetes) & !(type_1_diabetes) & !(unspecified_diabetes), "Possible Type 2 diabetes",
  # Type 1 diabetes only if no record of type 2 or otherwise unspecified diabetes
  (type_1_diabetes) & !(type_2_diabetes) & !(unspecified_diabetes), "Possible Type 1 diabetes",
  # Otherwise not possible to determine specific diabetes 
  default = "Uncertain diabetes status"
)]

# Combine self-report and hospital record diabetes status adjudication
prevalent_diabetes <- merge(diabetes_self_report, diabetes_hes[,.(eid, visit_index, adjudication)], by=c("eid", "visit_index"), all=TRUE)
prevalent_diabetes <- prevalent_diabetes[visit_index == 0] # Self report adjudication only performed for baseline assessment currently
setnames(prevalent_diabetes, c("adjudicated_diabetes", "adjudication"), c("adjudication_self_report", "adjudication_hes"))

prevalent_diabetes[, adjudication_reason := fcase(
  is.na(adjudication_hes), "Withdrawn consent", # Missing data here because Eastwood adjudication run on older UKB dataset with fewer sample withdrawals
  adjudication_self_report == "Diabetes unlikely" & adjudication_hes == "Diabetes unlikely", "Diabetes unlikely",
  adjudication_self_report == "Diabetes unlikely" & adjudication_hes == "Possible Type 2 diabetes", "Possible Type 2 diabetes (hospital records)",
  adjudication_self_report == "Diabetes unlikely" & adjudication_hes == "Possible Type 1 diabetes", "Possible Type 1 diabetes (hospital records)",
  adjudication_self_report == "Diabetes unlikely" & adjudication_hes == "Possible gestational diabetes", "Possible gestational diabetes (hospital records)",
  adjudication_self_report == "Diabetes unlikely" & adjudication_hes == "Uncertain diabetes status", "Uncertain diabetes status (hospital records)",
  adjudication_self_report == "Diabetes unlikely" & adjudication_hes == "Withdrawn consent for hospital record linkage", "Withdrawn consent for hospital record linkage",
  adjudication_self_report == "Uncertain diabetes status" & adjudication_hes == "Diabetes unlikely", "Uncertain diabetes status (self-report)",
  adjudication_self_report == "Uncertain diabetes status" & adjudication_hes == "Possible Type 2 diabetes", "Possible Type 2 diabetes (hospital records)", # all missing data in self-report
  adjudication_self_report == "Possible gestational diabetes" & adjudication_hes == "Possible Type 2 diabetes", "Possible Type 2 diabetes (hospital records)", # T2D in HES at later date than gestation diabetes
  adjudication_self_report == "Possible gestational diabetes" & adjudication_hes == "Possible Type 1 diabetes", "Possible Type 1 diabetes (hospital records)", # T1D independent of pregnancy
  adjudication_self_report == "Possible gestational diabetes" & adjudication_hes == "Uncertain diabetes status", "Uncertain diabetes status (hospital records)", # Unspecified diabetes subsequent to or independent of pregnancy codings
  adjudication_self_report == "Possible gestational diabetes", "Possible gestational diabetes (self-report)", # Hospital records in agreement, or no diabetes codings in hospital records
  adjudication_self_report == "Possible Type 1 diabetes", "Possible Type 1 diabetes (self-report)",
  adjudication_self_report == "Probable Type 1 diabetes", "Probable Type 1 diabetes (self-report)",
  adjudication_self_report == "Possible Type 2 diabetes", "Possible Type 2 diabetes (self-report)",
  adjudication_self_report == "Probable Type 2 diabetes", "Probable Type 2 diabetes (self-report)",
  default = "Coding Error"
)]

# Create specific diabetes columns
prevalent_diabetes <- prevalent_diabetes[, .(eid, visit_index, diabetes_history_reason=adjudication_reason,
  prevalent_diabetes = fcase(
    adjudication_reason %like% "Type 1 diabetes", TRUE,
    adjudication_reason %like% "Type 2 diabetes", TRUE,
    adjudication_reason %like% "Uncertain diabetes", TRUE,
    adjudication_reason %like% "Withdrawn consent", NA,
    default = FALSE
  ),
  no_history_any_diabetes = fcase(
    adjudication_reason %like% "Withdrawn consent", NA,
    adjudication_reason == "Diabetes unlikely", TRUE,
    default = FALSE
  ), 
  type_1_diabetes = fcase(
    adjudication_reason %like% "Type 1 diabetes", TRUE,
    adjudication_reason %like% "Uncertain diabetes", NA, 
    adjudication_reason %like% "Withdrawn consent", NA,
    default = FALSE
  ), 
  type_2_diabetes = fcase(
    adjudication_reason %like% "Type 2 diabetes", TRUE,
    adjudication_reason %like% "Uncertain diabetes", NA, 
    adjudication_reason %like% "Withdrawn consent", NA,
    default = FALSE
  ), 
  uncertain_diabetes = fcase(
    adjudication_reason %like% "Uncertain diabetes", TRUE,
    adjudication_reason %like% "Withdrawn consent", NA,
    default = FALSE
  ),
  history_gestational_diabetes = fcase(
    adjudication_reason %like% "gestational diabetes", TRUE,
    adjudication_reason == "Diabetes unlikely", FALSE,
    default = NA
  )
)]

# Get follow-up information for prevalent disease
followup_history <- fread("data/curated/followup/followup.txt")
followup_history <- followup_history[,.(eid, visit_index, 
  lost_to_followup, lost_to_followup_reason, lost_to_followup_date, # Needed for incident disease as well - curate this later
  earliest_hospital_date, earliest_hospital_nation)]

assessment_information <- anthro[, .(eid, visit_index, assessment_date, assessment_centre, assessment_nation)]

####################################################################
# Combine baseline and prevalent disease data into a single dataset
####################################################################
dat <- copy(ethnicity)
dat <- merge(dat, sex, by=c("eid", "visit_index"), all.x=TRUE)
dat <- merge(dat, age, by=c("eid", "visit_index"), all.x=TRUE)
dat <- merge(dat, bmi, by=c("eid", "visit_index"), all.x=TRUE)
dat <- merge(dat, height, by=c("eid", "visit_index"), all.x=TRUE)
dat <- merge(dat, weight, by=c("eid", "visit_index"), all.x=TRUE)
dat <- merge(dat, townsend, by=c("eid", "visit_index"), all.x=TRUE)
dat <- merge(dat, smoking, by=c("eid", "visit_index"), all.x=TRUE)
dat <- merge(dat, family_history_diabetes, by=c("eid", "visit_index"), all.x=TRUE)
dat <- merge(dat, history_cvd, by=c("eid", "visit_index"), all.x=TRUE)
dat <- merge(dat, history_pcos, by=c("eid", "visit_index"), all.x=TRUE)
dat <- merge(dat, history_learning_difficulties, by=c("eid", "visit_index"), all.x=TRUE)
dat <- merge(dat, history_bipolar_schizophrenia, by=c("eid", "visit_index"), all.x=TRUE)
dat <- merge(dat, hypertension_medication, by=c("eid", "visit_index"), all.x=TRUE)
dat <- merge(dat, lipid_lowering_medication, by=c("eid", "visit_index"), all.x=TRUE)
dat <- merge(dat, atypical_antipsychotics, by=c("eid", "visit_index"), all.x=TRUE)
dat <- merge(dat, systematic_corticosteroids, by=c("eid", "visit_index"), all.x=TRUE)
dat <- merge(dat, biomarkers, by=c("eid", "visit_index"), all.x=TRUE)
dat <- merge(dat, fasting, by=c("eid", "visit_index"), all.x=TRUE)
dat <- merge(dat, prevalent_diabetes, by=c("eid", "visit_index"), all.x=TRUE)
dat <- merge(dat, followup_history, by=c("eid", "visit_index"), all.x=TRUE)
dat <- merge(dat, assessment_information, by=c("eid", "visit_index"), all.x=TRUE)

dat <- dat[visit_index == 0] # Eastwood et al. algorithms only run for baseline assessment

# For gestational diabetes, need to make sure self-reported Males with withdrawn consent for
# hospital linkage are set to FALSE
dat[sex == "Male", history_gestational_diabetes := FALSE]

# Modify lost to follow-up information - we only care about participants being lost to follow-up 
# where we can perform incident disease analysis (except where the reason lost to follow up is withdrawn consent
# for hospital record linkage)
dat[(prevalent_diabetes) & lost_to_followup_reason != "Participant has withdrawn consent for future linkage", lost_to_followup := NA]
dat[(prevalent_diabetes) & lost_to_followup_reason != "Participant has withdrawn consent for future linkage", lost_to_followup_date := NA]
dat[(prevalent_diabetes) & lost_to_followup_reason != "Participant has withdrawn consent for future linkage", lost_to_followup_reason := NA]

#############################################################################################
# Flag cases where its possible someone has undiagnosed diabetes based on HBA1C and glucose
#############################################################################################
# https://cks.nice.org.uk/topics/diabetes-type-2/diagnosis/diagnosis-in-adults/

# Using non-fasting glucose
dat[, non_fasting_undiagnosed_diabetes_reason := fcase(
  (type_2_diabetes), "Adjudicated T2D case",
  hba1c >= 48 & non_fasting_glucose >= 11.1, "Possible undiagnosed diabetes (tier 1) | HbA1c >= 48 mmol/mol and non-fasting glucose >= 11.1 mmol/L",
  hba1c >= 48 & is.na(non_fasting_glucose), "Possible undiagnosed diabetes (tier 2) | HbA1c >= 48 mmol/mol, missing non-fasting glucose measurement",
  non_fasting_glucose >= 11.1 & is.na(hba1c), "Possible undiagnosed diabetes (tier 2) | Non-fasting glucose >= 11.1 mmol/L, missing HbA1c measurement",
  hba1c >= 48 & non_fasting_glucose < 11.1, "Possible undiagnosed diabetes (tier 3) | HbA1c >= 48 mmol/mol but non-fasting glucose < 11.1 mmol/L",
  non_fasting_glucose >= 11.1 & hba1c < 48, "Possible undiagnosed diabetes (tier 3) | Non-fasting glucose >= 11.1 mmol/L but HbA1c < 48 mmol/mol",
  hba1c < 48 & non_fasting_glucose < 11.1, "No evidence of diabetes (tier 1) | HbA1c < 48 mmol/mol and non-fasting glucose < 11.1 mmol/L",
  is.na(hba1c) & non_fasting_glucose < 11.1, "No evidence of diabetes (tier 2) | Non-fasting glucose < 11.1 mmol/L but missing HbA1c measurement",
  is.na(non_fasting_glucose) & hba1c < 48, "No evidence of diabetes (tier 2) | HbA1c < 48 mmol/mol but missing non-fasting glucose measurement",
  is.na(hba1c) & is.na(non_fasting_glucose), "No evidence of diabetes (tier 3) | Missing HbA1c and non-fasting glucose measurements",
  default = "Coding Error"
)]

dat[, non_fasting_undiagnosed_diabetes := non_fasting_undiagnosed_diabetes_reason %like% "Possible undiagnosed diabetes"]

# Using fasting glucose
dat[, fasting_undiagnosed_diabetes_reason := fcase(
  (type_2_diabetes), "Adjudicated T2D case",
  hba1c >= 48 & fasting_glucose >= 7, "Possible undiagnosed diabetes (tier 1) | HbA1c >= 48 mmol/mol and fasting glucose >= 7 mmol/L",
  hba1c >= 48 & is.na(fasting_glucose), "Possible undiagnosed diabetes (tier 2) | HbA1c >= 48 mmol/mol, missing fasting glucose measurement",
  fasting_glucose >= 7 & is.na(hba1c), "Possible undiagnosed diabetes (tier 2) | Fasting glucose >= 7 mmol/L, missing HbA1c measurement",
  hba1c >= 48 & fasting_glucose < 7, "Possible undiagnosed diabetes (tier 3) | HbA1c >= 48 mmol/mol but fasting glucose < 7 mmol/L",
  fasting_glucose >= 7 & hba1c < 48, "Possible undiagnosed diabetes (tier 3) | Fasting glucose >= 7 mmol/L but HbA1c < 48 mmol/mol",
  hba1c < 48 & fasting_glucose < 7, "No evidence of diabetes (tier 1) | HbA1c < 48 mmol/mol and fasting glucose < 7 mmol/L",
  is.na(hba1c) & fasting_glucose < 7, "No evidence of diabetes (tier 2) | Fasting glucose < 7 mmol/L but missing HbA1c measurement",
  is.na(fasting_glucose) & hba1c < 48, "No evidence of diabetes (tier 2) | HbA1c < 48 mmol/mol but missing fasting glucose measurement",
  is.na(hba1c) & is.na(fasting_glucose), "No evidence of diabetes (tier 3) | Missing HbA1c and fasting glucose measurements",
  default = "Coding Error"
)]

dat[, fasting_undiagnosed_diabetes := fasting_undiagnosed_diabetes_reason %like% "Possible undiagnosed diabetes"]

#############################################################################################
# Flag cases of potential pre-diabetes
#############################################################################################

dat[, non_fasting_prediabetes_reason := fcase(
  (type_2_diabetes), "Adjudicated T2D case",
  hba1c >= 48 & non_fasting_glucose >= 11.1, "Possible undiagnosed diabetes (tier 1) | HbA1c >= 48 mmol/mol and non-fasting glucose >= 11.1 mmol/L",
  hba1c >= 48 & is.na(non_fasting_glucose), "Possible undiagnosed diabetes (tier 2) | HbA1c >= 48 mmol/mol, missing glucose measurement",
  non_fasting_glucose >= 11.1 & is.na(hba1c), "Possible undiagnosed diabetes (tier 2) | Non-fasting glucose >= 11.1 mmol/L, missing HbA1c measurement",
  hba1c >= 48 & non_fasting_glucose < 11.1, "Possible undiagnosed diabetes (tier 3)| HbA1c >= 48 mmol/mol but non-fasting glucose < 11.1 mmol/L",
  non_fasting_glucose >= 11.1 & hba1c < 48, "Possible undiagnosed diabetes (tier 3)| Non-fasting glucose >= 11.1 mmol/L but HbA1c < 48 mmol/mol",
  hba1c >= 42, "Possible pre-diabetes (tier 1) | HbA1c >= 42 mmol/mol",
  hba1c < 42, "No evidence of pre-diabetes (tier 1)| HbA1c < 42 mmol/mol",
  is.na(hba1c), "No evidence of pre-diabetes (tier 3) | Missing HbA1c measurement",
  default = "Coding Error"
)]

dat[, non_fasting_prediabetes := non_fasting_prediabetes_reason %like% "Possible pre-diabetes"] 

dat[, fasting_prediabetes_reason := fcase(
  (type_2_diabetes), "Adjudicated T2D case",
  hba1c >= 48 & fasting_glucose >= 7, "Possible undiagnosed diabetes (tier 1) | HbA1c >= 48 mmol/mol and fasting glucose >= 7 mmol/L",
  hba1c >= 48 & is.na(fasting_glucose), "Possible undiagnosed diabetes (tier 2) | HbA1c >= 48 mmol/mol, missing glucose measurement",
  fasting_glucose >= 7 & is.na(hba1c), "Possible undiagnosed diabetes (tier 2) | Fasting glucose >= 7 mmol/L, missing HbA1c measurement",
  hba1c >= 48 & fasting_glucose < 7, "Possible undiagnosed diabetes (tier 3)| HbA1c >= 48 mmol/mol but fasting glucose < 7 mmol/L",
  fasting_glucose >= 7 & hba1c < 48, "Possible undiagnosed diabetes (tier 3)| Fasting glucose >= 7 mmol/L but HbA1c < 48 mmol/mol",
  hba1c >= 42 & fasting_glucose >= 5.6, "Possible pre-diabetes (tier 1) | HbA1c >= 42 mmol/mol and fasting-glucose >= 5.6 mmol/L",
  hba1c >= 42 & is.na(fasting_glucose), "Possible pre-diabetes (tier 2) | HbA1c >= 42 mmol/mol, missing glucose measurement",
  is.na(hba1c) & fasting_glucose >= 5.6, "Possible pre-diabetes (tier 2) | Fasting glucose >= 5.6 mmol/L, missing HbA1c measurement",
  hba1c >= 42 & fasting_glucose < 5.6, "Possible pre-diabetes (tier 3) | HbA1c >= 42 mmol/mol but fasting glucose < 5.6 mmol/L",
  hba1c < 42 & fasting_glucose >= 5.6, "Possible pre-diabetes (tier 3) | Fasting glucose >= 5.6 mmol/L but HbA1c < 42 mmol/mol",
  hba1c < 42 & fasting_glucose < 5.6, "No evidence of pre-diabetes (tier 1)| HbA1c < 42 mmol/mol and fasting glucose < 5.6 mmol/L",
  is.na(hba1c) & fasting_glucose < 5.6, "No evidence of pre-diabetes (tier 2) | Fasting glucose < 5.6 mmol/L but missing HbA1c measurement",
  is.na(fasting_glucose) & hba1c < 42, "No evidence of pre-diabetes (tier 2) | HbA1c < 42 mmol/mol but missing glucose measurement",
  is.na(hba1c) & is.na(fasting_glucose), "No evidence of pre-diabetes (tier 3) | Missing HbA1c and glucose measurements",
  default = "Coding Error"
)]

dat[, fasting_prediabetes := fasting_prediabetes_reason %like% "Possible pre-diabetes"] 

###################################
# Curate incident diabetes data
###################################
source("src/functions/calendar_year_math.R")

incident_diabetes <- fread("data/adjudicated_medical_history/Eastwood_diabetes/EastwoodDiabetes_IncidentDiabetes.csv")
incident_diabetes <- incident_diabetes[,.(eid=V1, visit_index=0,
  incident_censor_years=T2D.Incident_Censor_Years,
  incident_censor_date=as.IDate(T2D.Incident_Censor_Date),
  incident_censor_reason=T2D.Incident_Censor_Reason)] # all same for different diabetes types

# Drop people who have prevalent t2d, t1d, or uncertain diabetes status now that we've included the hospital records
drop_inci <- prevalent_diabetes[(prevalent_diabetes) | is.na(prevalent_diabetes)]
incident_diabetes <- incident_diabetes[!drop_inci, on = .(eid, visit_index)]

# Add in additional information we need to curate incident disease follow-up
followup <- fread("data/curated/followup/followup.txt")
incident_diabetes[followup, on = .(eid, visit_index), max_death_date := latest_mortality_date]
incident_diabetes[followup, on = .(eid, visit_index), censor_hospital_nation := latest_hospital_nation]
incident_diabetes[followup, on = .(eid, visit_index), max_censor_date_by_nation := latest_hospital_date]
incident_diabetes[followup, on = .(eid, visit_index), assessment_date := i.assessment_date]
incident_diabetes[followup, on = .(eid, visit_index), age_decimal := i.age_decimal]
incident_diabetes[followup, on = .(eid, visit_index), lost_to_followup_date := i.lost_to_followup_date]
incident_diabetes[followup, on = .(eid, visit_index), lost_to_followup_reason := i.lost_to_followup_reason]

# Flag cases where lack of diabetes is participant being lost to followup
incident_diabetes[lost_to_followup_date < incident_censor_date, 
  c("incident_censor_years", "incident_censor_date", "incident_censor_reason") := 
  .(years_between(assessment_date, lost_to_followup_date), lost_to_followup_date, 
    ifelse(incident_censor_reason %like% "No in-patient data",
      "No in-patient data (HES) | Lost to follow-up (see 'lost_to_followup_reason')",
      "No evidence of diabetes (HES) | Lost to follow-up (see 'lost_to_followup_reason')"))]

# Flag cases where follow-up is truncated to prevent confounding from pandemic
# (i.e. due to change in behaviour or environment due to lockdowns, or SARS-CoV2 exposure)
incident_diabetes[incident_censor_date == "2020-02-01" & 
  !(incident_censor_reason %like% "Incident") & 
  !(incident_censor_reason %like% "Death"),
    incident_censor_reason := ifelse(incident_censor_reason %like% "No in-patient data",
      "No in-patient data (HES) | Max follow-date (2020-02-01; pandemic cut-off)",
      "No evidence of diabetes (HES) | Max follow-date (2020-02-01; pandemic cut-off)")]

# Where the participant is diabetes free at the max censor date, make sure to update this with the 
# inferred nation of residence at follow-up (rather than nation at baseline) and update with latest
# follow-up cut-off dates available (cut-off dates are hard-coded in incident diabetes notebook, not
# detected from the data)

# Some of these people are those in Wales at baseline who have moved to England/Scotland but don't 
# have any diabetes records in the HES data or death records
incident_diabetes[incident_censor_reason %like% "CENSOR DATE" & max_censor_date_by_nation > "2020-02-01",
  c("incident_censor_years", "incident_censor_date", "incident_censor_reason") :=
  .(years_between(assessment_date, as.IDate("2020-02-01")), as.IDate("2020-02-01"), 
    ifelse(incident_censor_reason %like% "No in-patient data",
      "No in-patient data (HES) | Max follow-date (2020-02-01; pandemic cut-off)",
      "No evidence of diabetes (HES) | Max follow-date (2020-02-01; pandemic cut-off)"))]

# Check if max_censor_date_by_nation matches max_death_date - in this case, these are people
# who have died, but have no diabetes in their HES/DEATH records
incident_diabetes[incident_censor_reason %like% "CENSOR DATE" & max_censor_date_by_nation == max_death_date,
  c("incident_censor_years", "incident_censor_date", "incident_censor_reason") :=
  .(years_between(assessment_date, max_censor_date_by_nation), max_censor_date_by_nation,
    gsub("CENSOR DATE (Wales)", "Death Record (no diabetes)", incident_censor_reason, fixed=TRUE))]

# Remaining are people in Wales, but make sure to update with appropriate censor date in the data
incident_diabetes[incident_censor_reason %like% "CENSOR DATE",
  c("incident_censor_years", "incident_censor_date", "incident_censor_reason") :=
  .(years_between(assessment_date, max_censor_date_by_nation), max_censor_date_by_nation,
    gsub("CENSOR DATE (Wales)", "Max follow-up available (Hospitals in Wales)", incident_censor_reason, fixed=TRUE))]

# Truncate follow-up at 10 years to match QDiabetes time horizon
incident_diabetes[, incident_censor_years := years_between(assessment_date, incident_censor_date)]
incident_diabetes[incident_censor_years > 10,
  c("incident_censor_years", "incident_censor_date", "incident_censor_reason") :=
  .(10, add_years(assessment_date, 10), ifelse(incident_censor_reason %like% "No in-patient data",
    "No in-patient data (HES) | Max follow time (10 years)",
    "No evidence of diabetes (HES) | Max follow time (10 years)"))]

# Relabel Diabetes_unspecified
incident_diabetes[, incident_censor_reason := gsub("Diabetes_unspecified", "Uncertain Diabetes", incident_censor_reason)]

# Compute age at censor
incident_diabetes[age, on = .(eid, visit_index), age_decimal := i.age_decimal]
incident_diabetes[, incident_censor_age := floor(age_decimal + incident_censor_years)]

# Collate TRUE/FALSE columns for incident diabetes case status
incident_diabetes[, incident_type_2_diabetes := fcase(
  incident_censor_reason %like% "Incident T2D", TRUE,
  incident_censor_reason %like% "Incident Uncertain Diabetes", NA,
  default = FALSE)]
incident_diabetes[, incident_type_1_diabetes := fcase(
  incident_censor_reason %like% "Incident T1D", TRUE,
  incident_censor_reason %like% "Incident Uncertain Diabetes", NA,
  default = FALSE)]
incident_diabetes[, incident_uncertain_diabetes := ifelse(incident_censor_reason %like% "Incident Uncertain Diabetes", TRUE, FALSE)]
incident_diabetes[, incident_any_diabetes := ifelse(incident_censor_reason %like% "Incident", TRUE, FALSE)]

# Determine all-cause mortality
incident_diabetes[, death_at_censor_date := incident_censor_date == max_death_date | lost_to_followup_reason %like% "Death"]

# Filter to columns we want to add to main 'dat'
incident_diabetes <- incident_diabetes[,.(eid, visit_index, 
  incident_censor_years, incident_censor_date, incident_censor_reason, incident_censor_age,
  incident_type_2_diabetes, incident_type_1_diabetes, incident_uncertain_diabetes, incident_any_diabetes,
  censor_hospital_nation, death_at_censor_date)]

# Add to 'dat'
dat <- merge(dat, incident_diabetes, by=c("eid", "visit_index"), all.x=TRUE)

# Drop additional sample withdrawal dropped in incident diabetes file, but not yet dropped in others
# (Note analysts are expected to drop latest withdrawals downstream, as these change between data releases)
dat <- dat[!(is.na(incident_any_diabetes) & diabetes_history_reason == "Diabetes unlikely")]

#####################################
# Compute QDiabetes score(s)
#####################################

dat[!(prevalent_diabetes) & townsend < 11 &
    !is.na(smoking_status) & !is.na(bmi) & !is.na(townsend) & !is.na(family_history_diabetes) & 
    !is.na(family_history_diabetes) & !is.na(hypertension_medication) & !is.na(history_cvd) &
    !is.na(systematic_corticosteroids) & height >= 1.4 & height <= 2.1 & weight >= 40 & weight <= 180,
			QDiabetes2013 := QDR2013(surv = 10, # 10 year risk to match 2018 score time horizons
				sex = sex, age = age, bmi = bmi, ethn = ethnicity, smoke = smoking_status, tds = townsend,
				fhdm = family_history_diabetes, htn = hypertension_medication, cvd = history_cvd, ster = systematic_corticosteroids)]

dat[!(prevalent_diabetes) &
    !is.na(smoking_status) & !is.na(bmi) & !is.na(townsend) & !is.na(family_history_diabetes) & 
    !is.na(hypertension_medication) & !is.na(history_cvd) & !is.na(history_gestational_diabetes) &
    !is.na(history_pcos) & !is.na(history_learning_difficulties) & !is.na(history_bipolar_schizophrenia) &
    !is.na(systematic_corticosteroids) & !is.na(lipid_lowering_medication) & !is.na(atypical_antipsychotics) &
    height >= 1.4 & height <= 2.1 & weight >= 40 & weight <= 180,
      QDiabetes2018A := QDR2018A(
        sex = sex, age = age, bmi = bmi, ethn = ethnicity, smoke = smoking_status, tds = townsend,
        fhdm = family_history_diabetes, htn = hypertension_medication, cvd = history_cvd, 
        gdm = history_gestational_diabetes, pcos = history_pcos, learn = history_learning_difficulties,
        psy = history_bipolar_schizophrenia, ster = systematic_corticosteroids, stat = lipid_lowering_medication,
        apsy = atypical_antipsychotics)]   

dat[!(prevalent_diabetes) &
    !is.na(smoking_status) & !is.na(bmi) & !is.na(townsend) & !is.na(family_history_diabetes) & 
    !is.na(hypertension_medication) & !is.na(history_cvd) & !is.na(history_gestational_diabetes) &
    !is.na(history_pcos) & !is.na(history_learning_difficulties) & !is.na(history_bipolar_schizophrenia) &
    !is.na(systematic_corticosteroids) & !is.na(lipid_lowering_medication) & !is.na(atypical_antipsychotics) &
    height >= 1.4 & height <= 2.1 & weight >= 40 & weight <= 180 & fasting_glucose >= 2 & fasting_glucose < 7, 
      QDiabetes2018B_fasting := QDR2018B(
        sex = sex, age = age, bmi = bmi, ethn = ethnicity, smoke = smoking_status, tds = townsend,
        fhdm = family_history_diabetes, htn = hypertension_medication, cvd = history_cvd, 
        gdm = history_gestational_diabetes, pcos = history_pcos, learn = history_learning_difficulties,
        psy = history_bipolar_schizophrenia, ster = systematic_corticosteroids, stat = lipid_lowering_medication,
        apsy = atypical_antipsychotics, fpg = fasting_glucose)]   

dat[!(prevalent_diabetes) &
    !is.na(smoking_status) & !is.na(bmi) & !is.na(townsend) & !is.na(family_history_diabetes) & 
    !is.na(hypertension_medication) & !is.na(history_cvd) & !is.na(history_gestational_diabetes) &
    !is.na(history_pcos) & !is.na(history_learning_difficulties) & !is.na(history_bipolar_schizophrenia) &
    !is.na(systematic_corticosteroids) & !is.na(lipid_lowering_medication) & !is.na(atypical_antipsychotics) &
    height >= 1.4 & height <= 2.1 & weight >= 40 & weight <= 180 & non_fasting_glucose >= 2 & non_fasting_glucose < 7, 
      QDiabetes2018B_non_fasting := QDR2018B(
        sex = sex, age = age, bmi = bmi, ethn = ethnicity, smoke = smoking_status, tds = townsend,
        fhdm = family_history_diabetes, htn = hypertension_medication, cvd = history_cvd, 
        gdm = history_gestational_diabetes, pcos = history_pcos, learn = history_learning_difficulties,
        psy = history_bipolar_schizophrenia, ster = systematic_corticosteroids, stat = lipid_lowering_medication,
        apsy = atypical_antipsychotics, fpg = non_fasting_glucose)]   

dat[!(prevalent_diabetes) &
    !is.na(smoking_status) & !is.na(bmi) & !is.na(townsend) & !is.na(family_history_diabetes) & 
    !is.na(hypertension_medication) & !is.na(history_cvd) & !is.na(history_gestational_diabetes) &
    !is.na(history_pcos) & !is.na(history_learning_difficulties) & !is.na(history_bipolar_schizophrenia) &
    !is.na(systematic_corticosteroids) & !is.na(lipid_lowering_medication) & !is.na(atypical_antipsychotics) &
    height >= 1.4 & height <= 2.1 & weight >= 40 & weight <= 180 & hba1c >= 15 & hba1c < 48,
      QDiabetes2018C := QDR2018C(
        sex = sex, age = age, bmi = bmi, ethn = ethnicity, smoke = smoking_status, tds = townsend,
        fhdm = family_history_diabetes, htn = hypertension_medication, cvd = history_cvd, 
        gdm = history_gestational_diabetes, pcos = history_pcos, learn = history_learning_difficulties,
        psy = history_bipolar_schizophrenia, ster = systematic_corticosteroids, stat = lipid_lowering_medication,
        apsy = atypical_antipsychotics, hba1c = hba1c)]

#################################################
# Impose organisation on rows and columns
#################################################

dat <- dat[anthro[,.(eid, visit_index)], on = .(eid, visit_index), nomatch=0] # restore row order

dat <- dat[, .(
  # Participant and visit indentifiers
  eid, visit_index,
  # QDiabetes scores
  QDiabetes2018A, QDiabetes2018B_fasting, QDiabetes2018B_non_fasting, QDiabetes2018C, QDiabetes2013,
  # QDiabetes risk factors
  ethnicity, sex, age, bmi, height, weight, smoking_status, townsend, family_history_diabetes, 
  history_cvd, history_gestational_diabetes, history_pcos, history_learning_difficulties,
  history_bipolar_schizophrenia, hypertension_medication, lipid_lowering_medication, 
  systematic_corticosteroids, atypical_antipsychotics, hba1c, fasting_glucose,
  # Other glucose information
  fasting_time, non_fasting_glucose,
  # Prevalent diabetes
  prevalent_diabetes, no_history_any_diabetes, 
  type_2_diabetes, type_1_diabetes, uncertain_diabetes,
  diabetes_history_reason, 
  # Prevalent follow-up in hospital records
  earliest_hospital_date, earliest_hospital_nation,
  # Possible undiagnosed diabetes from HbA1c or Glucose
  fasting_undiagnosed_diabetes, fasting_undiagnosed_diabetes_reason,
  non_fasting_undiagnosed_diabetes, non_fasting_undiagnosed_diabetes_reason,
  # Possible pre-diabetes from HbA1c or Glucose
  fasting_prediabetes, fasting_prediabetes_reason,
  non_fasting_prediabetes, non_fasting_prediabetes_reason,
  # Incident diabetes
  incident_type_2_diabetes, incident_type_1_diabetes, incident_uncertain_diabetes,
  incident_any_diabetes, incident_censor_date, incident_censor_years, incident_censor_age,
  incident_censor_reason,
  # Follow-up information 
  assessment_date, assessment_centre, assessment_nation,
  death_at_censor_date, censor_hospital_nation, 
  lost_to_followup, lost_to_followup_date, lost_to_followup_reason
)]
  
# Write out 
system("mkdir -p output", wait=TRUE)
fwrite(dat, sep="\t", quote=FALSE, file="output/qdiabetes.txt")

# Curate column information
info <- rbind(use.names=FALSE,
  data.table(column_name="eid", description="Participant ID in project 7439"),
  data.table("visit_index", "Assessment visit (e.g. '0' for baseline assessment)"),
  data.table("QDiabetes2018A", "QDiabetes 2018 model A; does not use glucose or HbA1c measures"),
  data.table("QDiabetes2018B_fasting", "QDiabetes 2018 model B; uses glucose adjusted for fasting time (see 'fasting_glucose')"),
  data.table("QDiabetes2018B_non_fasting", "QDiabetes 2018 model B; uses glucose as-is (note: model expects fasting glucose)"),
  data.table("QDiabetes2018C", "QDiabetes 2018 model C; includes HbA1c"),
  data.table("QDiabetes2013", "QDiabetes 2013 score with 10 year time horizon"),
  data.table("ethnicity", "Ethnicity, mapped to QDiabetes categories. Admixed are grouped into 'Other' and missing/not reported into 'WhiteNA'"),
  data.table("sex", "Self-reported sex at baseline assessment"),
  data.table("age", "Age in years (integer) at baseline assessment"),
  data.table("bmi", "Body mass index"),
  data.table("height", "Standing height (in meters)"),
  data.table("weight", "Weight (in kilograms)"),
  data.table("smoking_status", "Smoking status mapped to QDiabetes categories. Light: < 10 cigarettes per day, Moderate: < 20 cigarettes per day (or missing data), Heavy: >= 20 cigarettes per day."),
  data.table("townsend", "Townsend deprivation index"),
  data.table("family_history_diabetes", "History of diabetes in first degree relatives (father/mother/siblings)."),
  data.table("history_cvd", "History of CVD (ischaemic heart disease, stroke, or transient ischaemic attack)"),
  data.table("history_gestational_diabetes", "History of gestational diabetes; adjudicated from self-report and hospital records (see README.txt)"),
  data.table("history_pcos", "History of polycystic ovary syndrome; self-reported or in hospital records"),
  data.table("history_learning_difficulties", "History of learning difficulties incidentally in hospital records"),
  data.table("history_bipolar_schizophrenia", "History of bipolar of schizophrenia disorders; self-reported or in hospital records"),
  data.table("hypertension_medication", "History of blood pressure medication, self-reported; used as proxy for history of treated hypertension"),
  data.table("lipid_lowering_medication", "History of lipid loewring medication, self-reported; used as proxy/umbrella term for statin usage"),
  data.table("systematic_corticosteroids", "History of systematic corticosteroid usage (oral or injected)"),
  data.table("atypical_antipsychotics", "History of second generation atypical antipsychotics"),
  data.table("hba1c", "Glycated haemoglobin (HbA1c) (mmol/mol)"),
  data.table("fasting_glucose", "Glucose (mmol/L) corrected for sex-specific differences in medians at 0 hours fasting, 1 hours, and 2 hours, vs. 3 or more hours based on findings of Moebus et al. 2011"),
  data.table("fasting_time", "Hours since last meal"),
  data.table("non_fasting_glucose", "Glucose (mmol/L) without any adjustments for fasting time"), 
  data.table("prevalent_diabetes", "TRUE where participant has any history of type 1, type 2, or uncertain diabetes; adjudicated from self-report and hospital records (see README.txt)"),
  data.table("no_history_any_diabetes", "TRUE where no history of any type of diabetes, including gestational; adjudicated from self-report and hospital records (see README.txt)"),
  data.table("type_2_diabetes", "Type 2 diabetes cases (TRUE) and controls (FALSE) at assessment; adjudicated from self-report and hospital records (see README.txt)"),
  data.table("type_1_diabetes", "Type 1 diabetes cases (TRUE) and controls (FALSE) at assessment; adjudicated from self-report and hospital records (see README.txt)"),
  data.table("uncertain_diabetes", "TRUE where participant has history of diabetes that could not be classified as type 1, type 2, or gestational; adjudicated from self-report and hospital records (see README.txt)"),
  data.table("type_2_diabetes", "Type 2 diabetes cases (TRUE) and controls (FALSE) at assessment; adjudicated from self-report and hospital records (see README.txt)"),
  data.table("diabetes_history_reason", "Details on adjudication of diabetes history"),
  data.table("earliest_hospital_date", "Earliest date at which hospital records could hypothetically be found for this participant based on their 'earliest_hospital_nation'"),
  data.table("earliest_hospital_nation", "Inferred nation of historical residence that determines the earliest date at which hospital records could be found based on the earliest retrospective linkage available for each nation"),
  data.table("fasting_undiagnosed_diabetes", "Possible undiagnosed diabetes based on HbA1c and/or fasting glucose concentrations"),
  data.table("fasting_undiagnosed_diabetes_reason", "Details on adjudication of fasting_undiagnosed_diabetes"),
  data.table("non_fasting_undiagnosed_diabetes", "Possible undiagnosed diabetes based on HbA1c and/or non-fasting glucose concentrations"),
  data.table("non_fasting_undiagnosed_diabetes_reason", "Details on adjudication of non_fasting_undiagnosed_diabetes"),
  data.table("fasting_prediabetes", "Possible prediabetes based on HbA1c and/or fasting glucose concentrations"),
  data.table("fasting_prediabetes_reason", "Details on adjudication of fasting_prediabetes"),
  data.table("non_fasting_prediabetes", "Possible prediabetes based on HbA1c alone"),
  data.table("non_fasting_prediabetes_reason", "Details on adjudication of non_fasting_prediabetes"),
  data.table("incident_type_2_diabetes", "TRUE where participant develops type 2 diabetes within 10 years of UK Biobank assessment; adjudicated from hospital (and death) records using Eastwood et al. 2016 algorithm"),
  data.table("incident_type_1_diabetes", "TRUE where participant develops type 1 diabetes within 10 years of UK Biobank assessment; adjudicated from hospital (and death) records using Eastwood et al. 2016 algorithm"),
  data.table("incident_uncertain_diabetes", "TRUE where participant develops uncertain diabetes within 10 years of UK Biobank assessment; adjudicated from hospital (and death) records using Eastwood et al. 2016 algorithm"),
  data.table("incident_any_diabetes", "TRUE where participant develops type 1, type 2, or uncertain diabetes within 10 years of UK Biobank assessment"),
  data.table("incident_censor_date", "Date type 1, type 2, or uncertain diabetes case status is adjudicated to begin, or maximum date of follow-up"),
  data.table("incident_censor_years", "Years after UK Biobank assessment diabetes the censor date occurs"),
  data.table("incident_censor_age", "Age (integer) of participant at censor date, accounting for month (and approximate) date of birth"),
  data.table("incident_censor_reason", "Details on censor reason at the censor date (e.g. end of follow-up, diabetes case occurs, etc.)"),
  data.table("assessment_date", "Date of assessment at UK Biobank"),
  data.table("assessment_centre", "Location of assessment centre attended by participant"),
  data.table("assessment_nation", "UK nation assessment centre is located in"),
  data.table("death_at_censor_date", "TRUE where the reason follow-up ends is due to participant death"),
  data.table("censor_hospital_nation", "Nation of hospital that records may most recently be found for each participant; note hospitals in Wales have less follow-up than England or Scotland"),
  data.table("lost_to_followup", "TRUE where participant is recorded by UK Biobank as being lost to follow-up"),
  data.table("lost_to_followup_date", "Latest date at which follow-up is available for the participant"),
  data.table("lost_to_followup_reason", "Reason recorded by UK Biobank participant was lost to follow-up")
)
fwrite(info, sep="\t", quote=FALSE, file="output/column_headers.txt")
