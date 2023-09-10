library(data.table)
source('src/utils/SCORE2.R')

#####################################################################
# Define analysis cohort without filtering to samples with NMR data
#####################################################################

# Load basic cohort information
dat <- fread("data/ukb_7439/anthropometrics/output/anthropometrics.txt")

# Filter to baseline assessment
dat <- dat[visit_index == 0]

# Extract relevant columns:
dat <- dat[, .(eid, assessment_date, assessment_centre, age, age_decimal, sex, bmi)]

# Drop latest set of sample withdrawals
withdrawals <- fread("data/ukb_7439/latest_withdrawals/output/latest_withdrawals.txt")
dat <- dat[!withdrawals, on = .(eid)]

# Add in systolic blood pressure
bp <- fread("data/ukb_7439/blood_pressure/output/blood_pressure.txt")
dat[bp[visit_index == 0], on = .(eid), sbp := i.sbp]

# Add in medication usage
meds_touchscreen <- fread("data/ukb_7439/medication/output/medications_simple.txt")
meds_interview <- fread("data/ukb_7439/medication/output/detailed_medications_summarised.txt")
meds <- merge(meds_touchscreen, meds_interview, by=c("eid", "visit_index"), all=TRUE)

dat[meds[visit_index == 0], on = .(eid), cholesterol_medication := i.cholesterol_medication | i.lipid_lowering_medication]

# Add in smoking status
smoking <- fread("data/ukb_7439/smoking/output/smoking_status.txt")
dat[smoking[visit_index == 0], on = .(eid), smoking := i.current_smoker]

# Add in diabetes status (Eastwood 2016 algorithm)
diab_sr <- fread("data/ukb_7439/Eastwood_diabetes/output/prevalent_diabetes.txt")
diab_hes <- fread("data/ukb_7439/Eastwood_diabetes/output/incident_diabetes.txt")
diabetes <- merge(diab_sr, diab_hes, by=c("eid", "visit_index"), all=TRUE, suffixes=c("_self_report", "_hes"))
dat[diabetes[visit_index == 0], on = .(eid), prevalent_diabetes_mellitus := fcase(
  adjudicated_diabetes_self_report == "Probable type 2 diabetes", TRUE,
  adjudicated_diabetes_self_report == "Probable type 1 diabetes", TRUE,
  adjudicated_diabetes_self_report == "Possible type 2 diabetes", TRUE,
  adjudicated_diabetes_self_report == "Possible type 1 diabetes", TRUE,
  adjudicated_diabetes_hes == "Prevalent diabetes", TRUE,
  default = FALSE
)]
# Add in prevalent CKD (UKB algorithmically defined outcome)
ado <- fread("data/ukb_7439/algorithmically_defined_outcomes/output/algorithmically_defined_outcomes.txt")
dat[ado, on = .(eid), prevalent_CKD := esrd_date < assessment_date]
dat[is.na(prevalent_CKD), prevalent_CKD := FALSE]

# Add in earliest onset of prevalent vascular disease and CAD (as defined by the Dutch Lipic Clinical Network) 
# and determine where the event onset is premature, to be used in conjunction with LDL cholesterol 
# to predict familiar hypercholesterolemia
DLCN_vasc_disease <- fread("data/ukb_7439/endpoints/endpoints/DLCN_premature_vascular_disease/events_and_followup.txt")
dat[DLCN_vasc_disease[visit_index == 0], on = .(eid), premature_vascular_disease := fcase(
  prevalent_event & sex == "Male" & prevalent_event_age < 55, TRUE,
  prevalent_event & sex == "Female" & prevalent_event_age < 60, TRUE,
  prevalent_event & sex == "Male" & is.na(prevalent_event_age) & x.age < 55, TRUE,
  prevalent_event & sex == "Female" & is.na(prevalent_event_age) & x.age < 60, TRUE,
  default = FALSE
)]

DLCN_cad <- fread("data/ukb_7439/endpoints/endpoints/DLCN_premature_CAD/events_and_followup.txt")
dat[DLCN_cad[visit_index == 0], on = .(eid), premature_cad := fcase(
  prevalent_event & sex == "Male" & prevalent_event_age < 55, TRUE,
  prevalent_event & sex == "Female" & prevalent_event_age < 60, TRUE,
  prevalent_event & sex == "Male" & is.na(prevalent_event_age) & x.age < 55, TRUE,
  prevalent_event & sex == "Female" & is.na(prevalent_event_age) & x.age < 60, TRUE,
  default = FALSE
)]

# Add in established atherosclerotic cardiovascular disease: the ESC 2021 guidelines
# (Table 4) pretty closely match CEU's prevalent vascular disease definition
ascvd <- fread("data/ukb_7439/endpoints/endpoints/CEU_prevalent_vascular_disease/events_and_followup.txt")
dat[ascvd[visit_index == 0], on = .(eid), ASCVD := i.prevalent_event]

# Add in CVD endpoint (SCORE2 definition)
cvd <- fread("data/ukb_7439/endpoints/endpoints/SCORE2_CVD/events_and_followup.txt")
cvd <- cvd[visit_index == 0]
dat[cvd, on = .(eid), incident_cvd := i.incident_event]
dat[cvd, on = .(eid), incident_cvd_followup := i.incident_event_followup]
dat[cvd, on = .(eid), incident_cvd_followup_date := i.incident_event_followup_date]
dat[cvd, on = .(eid), incident_cvd_is_fatal := fcase(
  i.incident_event_type == "death", TRUE,
  i.incident_event_type == "hospitalisation", FALSE,
  i.incident_event_type == "" & !is.na(i.incident_event), FALSE,
  is.na(i.incident_event), NA
)]
dat[cvd, on = . (eid), cvd_is_primary_cause := fcase(
  i.incident_cause_type == "primary", TRUE,
  i.incident_cause_type == "secondary", FALSE,
  i.incident_cause_type == "" & !is.na(i.incident_cause_type), FALSE,
  is.na(i.incident_cause_type), NA
)]

# Get information on where participants resided at baseline assessment, 
# maximum, and minimum follow-up available in hospital records (different
# hospital systems have different follow-up time available depending on
# nation of hospital)
dat[cvd, on = .(eid), latest_hospital_nation := i.latest_hospital_nation]
dat[ascvd, on = .(eid), earliest_hospital_nation := i.earliest_hospital_nation]

# Get additional information on follow-up
dat[cvd, on = .(eid), mortality_at_cvd_followup := i.mortality_at_followup_date]
dat[cvd, on = .(eid), lost_at_cvd_followup := i.lost_to_followup_reason != ""]
dat[cvd, on = .(eid), lost_at_cvd_followup_reason := i.lost_to_followup_reason]
dat[, cvd_follow_lt10_Wales := !(incident_cvd) & !(lost_at_cvd_followup) & !(mortality_at_cvd_followup) & incident_cvd_followup < 10]
stopifnot(all(dat[(cvd_follow_lt10_Wales), latest_hospital_nation == "Wales"]))

# Add in biochemistry biomarker data
bio <- fread("data/ukb_7439/biomarkers/output/biomarkers.txt")
bio <- bio[visit_index == 0L]
bio[, visit_index := NULL]
dat <- merge(dat, bio, by="eid", all.x=TRUE)

# Flag people missing blood biochemistry data
bio_sinfo <- fread("data/ukb_7439/biomarkers/output/samples_not_measured.txt")
bio_sinfo <- bio_sinfo[visit_index == 0L]
bio_sinfo[, visit_index := NULL]
dat[, no_blood_sample := FALSE]
dat[, no_urine_sample := FALSE]
dat[bio_sinfo, on = .(eid), no_blood_sample := i.no_blood_sample]
dat[bio_sinfo, on = .(eid), no_urine_sample := i.no_urine_sample]

blood_bio <- setdiff(names(bio), c("eid", "uriacc", "urianac", "uriamac", "uriakc"))
urine_bio <- setdiff(names(bio), c("eid", blood_bio))
dat[, no_blood_biomarkers := apply(as.matrix(dat[,blood_bio,with=FALSE]), 1, function(rr) { all(is.na(rr)) })]
dat[, no_urine_biomarkers := apply(as.matrix(dat[,urine_bio,with=FALSE]), 1, function(rr) { all(is.na(rr)) })]

# Score FH status according to Dutch Lipid Clinic Network diagnosit criteria:
dat[, DLCN_FH_score := 0]
dat[(premature_cad), DLCN_FH_score := DLCN_FH_score + 2]
dat[(premature_vascular_disease), DLCN_FH_score := DLCN_FH_score + 1]
dat[ldl >= 8.5, DLCN_FH_score := DLCN_FH_score + 8]
dat[ldl >= 6.5 & ldl < 8.5, DLCN_FH_score := DLCN_FH_score + 5]
dat[ldl >= 5.0 & ldl < 6.5, DLCN_FH_score := DLCN_FH_score + 3]
dat[ldl >= 4.0 & ldl < 5.0, DLCN_FH_score := DLCN_FH_score + 1]

# Add PRSs
PRSs <- rbind(idcol="PRS", fill=TRUE,
  CAD_metaGRS = fread("data/ukb_7439/PRS/CAD_metaGRS/CAD_metaGRS_PGS000018_b097e681_P7439_from_dosage.sscore.gz"),
  Stroke_metaGRS = fread("data/ukb_7439/PRS/Stroke_metaGRS/Stroke_metaGRS_PGS000039_6a7832a2_P7439_from_dosage.sscore.gz")
)
PRSs <- dcast(PRSs, IID ~ PRS, value.var="score_sum")
dat[PRSs, on = .(eid = IID), CAD_metaGRS := i.CAD_metaGRS]
dat[PRSs, on = .(eid = IID), Stroke_metaGRS := i.Stroke_metaGRS]
dat[, no_genetics := ifelse(is.na(CAD_metaGRS), TRUE, FALSE)]

# Set PRS to missing for people used to train the PRSs
prs_training <- unique(rbind(
  fread("data/ukb_7439/PRS/sample_splits/CAD_metaGRS_training_samples.txt"),
  fread("data/ukb_7439/PRS/sample_splits/Stroke_metaGRS_training_samples.txt")
))
setnames(prs_training, "eid")
dat[, prs_training_samples := FALSE]
dat[prs_training, on = .(eid), prs_training_samples := TRUE]

# -------------------------------------------------------------------
# Now do sample exclusions to derive analysis cohort
# -------------------------------------------------------------------
# Build second sample flowchart for full cohort
sample_info <- data.table(step="Baseline (excl. withdrawals)",
  samples=dat[,.N], CVD=NA_real_, exited=NA_real_, exited_cvd=NA_real_)

# Function to update sample information
update_sample_info <- function(step_name, dataset, last_dataset) {
  if (missing(dataset)) {
    dataset <- dat
  }
  current_samples <- dataset[,.N]
  current_cvd <- sum(dataset$incident_cvd, na.rm=TRUE)
  if (missing(last_dataset)) {
    last_samples <- sample_info[.N, samples]
    last_cvd <- sample_info[.N, CVD]
  } else {
    last_samples <- last_dataset[,.N]
    last_cvd <- last_dataset[, sum(incident_cvd, na.rm=TRUE)]
  }
  new_row <- data.table(
    step=step_name,
    samples=current_samples,
    CVD=current_cvd,
    exited=last_samples - current_samples,
    exited_cvd=ifelse(is.na(last_cvd), NA, last_cvd - current_cvd)
  )
  sample_info <<- rbind(sample_info, new_row)
}

# Flag withdrawal from EHR linkage
dat <- dat[!cvd[(ehr_linkage_withdrawn)], on = .(eid)] # no EHR linkage
update_sample_info("With EHR linkage")

# Assess eligibility for SCORE2 screening
dat_cpy <- copy(dat)
dat <- dat[age >= 40]
update_sample_info("40 years or older")

dat <- dat[age <= 69]
update_sample_info("69 years or younger")
update_sample_info("Eligible age for SCORE2 risk prediction", dat, dat_cpy)

dat <- dat[!(ASCVD) | is.na(ASCVD)]
update_sample_info("Without established atherosclerotic CVD")

dat_cpy2 <- copy(dat)
dat <- dat[!(prevalent_diabetes_mellitus) | is.na(prevalent_diabetes_mellitus)]
update_sample_info("Without type 1 or type 1 diabetes", dat, dat_cpy2)

dat_cpy2 <- copy(dat)
dat <- dat[!(prevalent_CKD) | is.na(prevalent_CKD)]
update_sample_info("Without CKD", dat, dat_cpy2)

dat_cpy2 <- copy(dat)
update_sample_info("Definite FH (LDL ≥ 8.5 mmol/L and premature CAD or vascular disease):", dat[DLCN_FH_score > 8], dat[DLCN_FH_score > 8])
update_sample_info("Probable FH (LDL ≥ 8.5 mmol/L)", dat[DLCN_FH_score >= 6 & !(premature_vascular_disease)], dat[DLCN_FH_score >= 6 & !(premature_vascular_disease)])
update_sample_info("Probable FH (LDL ≥ 6.5 mmol/L and premature CAD or vascular disease)", dat[DLCN_FH_score >= 6 & DLCN_FH_score <= 8 & premature_vascular_disease], dat[DLCN_FH_score >= 6 & DLCN_FH_score <= 8 & premature_vascular_disease])
update_sample_info("Possible FH (LDL ≥ 5 mmol/L)", dat[DLCN_FH_score >= 3 & DLCN_FH_score <= 5], dat[DLCN_FH_score >= 3 & DLCN_FH_score <= 5])
update_sample_info("No FH (LDL < 5 mmol/L or missing)", dat[DLCN_FH_score < 3], dat[DLCN_FH_score < 3])
dat <- dat[DLCN_FH_score < 6]
update_sample_info("Without probable FH", dat, dat_cpy2)

update_sample_info("Eligible for screening with SCORE2 according to ESC 2021 guidelines", dat, dat_cpy)

# Exclude people on statins
dat <- dat[!(cholesterol_medication) | is.na(cholesterol_medication)]
update_sample_info("Already treated with statins")

# Drop people missing quantitative conventional risk factors
dat_cpy <- copy(dat)

dat <- dat[!is.na(sbp)]
update_sample_info("With known SBP", dat, dat_cpy)

dat <- dat[!(no_blood_sample)]
update_sample_info("With blood sample taken")

dat <- dat[!(no_blood_biomarkers)]
update_sample_info("With data on any clinical biochemistry assays")

# Flag issues with clinical biochemistry assays
update_sample_info("Missing HDL cholesterol:", dat[is.na(hdl)], dat[is.na(hdl)])
update_sample_info("Missing Total cholesterol:", dat[is.na(tchol)], dat[is.na(tchol)])
update_sample_info("Missing HDL and total cholesterol:", dat[is.na(hdl) & is.na(tchol) & is.na(ldl)], dat[is.na(hdl) & is.na(tchol) & is.na(ldl)])

dat_cpy2 <- copy(dat)
dat <- dat[!is.na(hdl) & !is.na(tchol) & !is.na(ldl)]
update_sample_info("With non-missing data for HDL and total cholesterol", dat, dat_cpy2)

update_sample_info("With non-missing quantitative SCORE2 risk factors", dat, dat_cpy)


# Format flowchart and add percentages
sample_info[, exited_cvd := ifelse(
  is.na(exited_cvd), NA_character_,
  ifelse(exited == 0L, "0 (-%)",
  sprintf("%s (%s%%)", format(exited_cvd, big.mark=","), round(exited_cvd/(exited_cvd + CVD)*100, digits=2))
))]

sample_info[, exited := ifelse(
  is.na(exited), NA_character_,
  ifelse(exited == 0L, "0 (-%)",
  sprintf("%s (%s%%)", format(exited, big.mark=","), round(exited/(exited + samples)*100, digits=2))
))]

sample_info[, CVD := ifelse(
  is.na(CVD), NA_character_,
  ifelse(samples == 0L, "0 (-%)",
  sprintf("%s (%s%%)", format(CVD, big.mark=","), round(CVD/samples*100, digits=2))
))]

sample_info[, samples := ifelse(is.na(samples), NA_character_, format(samples, big.mark=","))]

# Write out sample information
fwrite(sample_info, sep="\t", quote=FALSE, file="analyses/full_UKB_sample_flowchart.txt")

# Add SCORE2
dat[, SCORE2 := score2(sex, age, smoking, sbp, tchol, hdl, type="linear predictor")]
dat[, SCORE2_excl_UKB := score2(sex, age, smoking, sbp, tchol, hdl, type="linear predictor", weights="excluding UK Biobank")]

# Write out analysis cohort
fwrite(dat, sep="\t", quote=FALSE, file="data/cleaned/full_UKB_analysis_cohort.txt")

