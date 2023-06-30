library(data.table)
library(ukbnmr)
library(ggplot2)
library(caret)
source('src/utils/SCORE2.R')

# Make output directories
system("mkdir -p data/cleaned/")
system("mkdir -p analyses")

# Load basic cohort information
dat <- fread("data/ukb/anthropometrics/output/anthropometrics.txt")

# Filter to baseline assessment
dat <- dat[visit_index == 0]

# Extract relevant columns:
dat <- dat[, .(eid, assessment_date, assessment_centre, age, age_decimal, sex, bmi)]

# Drop latest set of sample withdrawals
withdrawals <- fread("data/ukb/latest_withdrawals/output/latest_withdrawals.txt")
dat <- dat[!withdrawals, on = .(eid)]

# Add in systolic blood pressure
bp <- fread("data/ukb/blood_pressure/output/blood_pressure.txt")
dat[bp[visit_index == 0], on = .(eid), sbp := i.sbp]

# Add in medication usage
meds_touchscreen <- fread("data/ukb/medication/output/medications_simple.txt")
meds_interview <- fread("data/ukb/medication/output/detailed_medications_summarised.txt")
meds <- merge(meds_touchscreen, meds_interview, by=c("eid", "visit_index"), all=TRUE)

dat[meds[visit_index == 0], on = .(eid), cholesterol_medication := i.cholesterol_medication | i.lipid_lowering_medication]

# Add in smoking status
smoking <- fread("data/ukb/smoking/output/smoking_status.txt")
dat[smoking[visit_index == 0], on = .(eid), smoking := i.current_smoker]

# Add in diabetes status (Eastwood 2016 algorithm)
diab_sr <- fread("data/ukb/Eastwood_diabetes/output/prevalent_diabetes.txt")
diab_hes <- fread("data/ukb/Eastwood_diabetes/output/incident_diabetes.txt")
diabetes <- merge(diab_sr, diab_hes, by=c("eid", "visit_index"), all=TRUE, suffixes=c("_self_report", "_hes"))
dat[diabetes[visit_index == 0], on = .(eid), prevalent_diabetes_mellitus := fcase(
  adjudicated_diabetes_self_report == "Probable type 2 diabetes", TRUE,
  adjudicated_diabetes_self_report == "Probable type 1 diabetes", TRUE,
  adjudicated_diabetes_self_report == "Possible type 2 diabetes", TRUE,
  adjudicated_diabetes_self_report == "Possible type 1 diabetes", TRUE,
  adjudicated_diabetes_hes == "Prevalent diabetes", TRUE,
  default = FALSE
)]
dat[diabetes[visit_index == 0], on = .(eid), history_gestational_diabetes := fcase(
  adjudicated_diabetes_self_report == "Possible gestational diabetes", TRUE,
  default = FALSE
)]

# Add in prevalent CKD (UKB algorithmically defined outcome)
ado <- fread("data/ukb/algorithmically_defined_outcomes/output/algorithmically_defined_outcomes.txt")
dat[ado, on = .(eid), prevalent_CKD := esrd_date < assessment_date]
dat[is.na(prevalent_CKD), prevalent_CKD := FALSE]

# Add in earliest onset of prevalent vascular disease and CAD (as defined by the Dutch Lipic Clinical Network) 
# and determine where the event onset is premature, to be used in conjunction with LDL cholesterol 
# to predict familiar hypercholesterolemia
DLCN_vasc_disease <- fread("data/ukb/endpoints/endpoints/DLCN_premature_vascular_disease/events_and_followup.txt")
dat[DLCN_vasc_disease[visit_index == 0], on = .(eid), premature_vascular_disease := fcase(
  prevalent_event & sex == "Male" & prevalent_event_age < 55, TRUE,
  prevalent_event & sex == "Female" & prevalent_event_age < 60, TRUE,
  prevalent_event & sex == "Male" & is.na(prevalent_event_age) & x.age < 55, TRUE,
  prevalent_event & sex == "Female" & is.na(prevalent_event_age) & x.age < 60, TRUE,
  default = FALSE
)]

DLCN_cad <- fread("data/ukb/endpoints/endpoints/DLCN_premature_CAD/events_and_followup.txt")
dat[DLCN_cad[visit_index == 0], on = .(eid), premature_cad := fcase(
  prevalent_event & sex == "Male" & prevalent_event_age < 55, TRUE,
  prevalent_event & sex == "Female" & prevalent_event_age < 60, TRUE,
  prevalent_event & sex == "Male" & is.na(prevalent_event_age) & x.age < 55, TRUE,
  prevalent_event & sex == "Female" & is.na(prevalent_event_age) & x.age < 60, TRUE,
  default = FALSE
)]

# Add in established atherosclerotic cardiovascular disease: the ESC 2021 guidelines
# (Table 4) pretty closely match CEU's prevalent vascular disease definition
ascvd <- fread("data/ukb/endpoints/endpoints/CEU_prevalent_vascular_disease/events_and_followup.txt")
dat[ascvd[visit_index == 0], on = .(eid), ASCVD := i.prevalent_event]

# Add in CVD endpoint (SCORE2 definition)
cvd <- fread("data/ukb/endpoints/endpoints/SCORE2_CVD/events_and_followup.txt")
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

# Add in Ischaemic Stroke and CAD endpoints, to be used for model training. Note these have
# separate follow-up columns as we don't treat CVD as a competing risk, only all-cause mortality.

# For Stroke we use UK Biobank's algorithmically defined outcome to match the 2019 Stroke metaGRS
# paper (here with additional curation for follow-up time in non-cases)
stroke <- fread("data/ukb/endpoints/endpoints/ADO_Stroke_IS_10year/events_and_followup.txt")
stroke <- stroke[visit_index == 0]
dat[stroke, on = .(eid), incident_stroke := i.incident_event]
dat[stroke, on = .(eid), incident_stroke_followup := i.incident_event_followup]
dat[stroke, on = .(eid), incident_stroke_followup_date := i.incident_event_followup_date]
dat[stroke, on = .(eid), incident_stroke_is_fatal := fcase(
  i.incident_event_type == "death", TRUE,
  i.incident_event_type == "hospitalisation", FALSE,
  i.incident_event_type == "" & !is.na(i.incident_event), FALSE,
  is.na(i.incident_event), NA
)]
dat[stroke, on = . (eid), stroke_is_primary_cause := fcase(
  i.incident_cause_type == "primary", TRUE,
  i.incident_cause_type == "secondary", FALSE,
  i.incident_cause_type == "" & !is.na(i.incident_cause_type), FALSE,
  is.na(i.incident_cause_type), NA
)]
dat[stroke, on = .(eid), prevalent_stroke := i.prevalent_event]

# For CAD we the same definition as the 2018 CAD metaGRS paper
cad <- fread("data/ukb/endpoints/endpoints/MetaGRS_CAD_10yr/events_and_followup.txt")
cad <- cad[visit_index == 0]
dat[cad, on = .(eid), incident_cad := i.incident_event]
dat[cad, on = .(eid), incident_cad_followup := i.incident_event_followup]
dat[cad, on = .(eid), incident_cad_followup_date := i.incident_event_followup_date]
dat[cad, on = .(eid), incident_cad_is_fatal := fcase(
  i.incident_event_type == "death", TRUE,
  i.incident_event_type == "hospitalisation", FALSE,
  i.incident_event_type == "operation", FALSE,
  i.incident_event_type == "" & !is.na(i.incident_event), FALSE,
  is.na(i.incident_event), NA
)]
dat[cad, on = . (eid), cad_is_primary_cause := fcase(
  i.incident_cause_type == "primary", TRUE,
  i.incident_cause_type == "secondary", FALSE,
  i.incident_cause_type == "" & !is.na(i.incident_cause_type), FALSE,
  is.na(i.incident_cause_type), NA
)]
dat[cad, on = .(eid), prevalent_cad := i.prevalent_event]
 
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

dat[stroke, on = .(eid), mortality_at_stroke_followup := i.mortality_at_followup_date]
dat[stroke, on = .(eid), lost_at_stroke_followup := i.lost_to_followup_reason != ""]
dat[stroke, on = .(eid), lost_at_stroke_followup_reason := i.lost_to_followup_reason]
dat[stroke, on = .(eid), stroke_follow_lt10_Wales := !(incident_event) & incident_event_followup < 10 & incident_event_followup_date == latest_hospital_date & latest_hospital_nation == "Wales"]
dat[, stroke_follow_lt10_Wales := !(incident_stroke) & !(lost_at_stroke_followup) & !(mortality_at_stroke_followup) & incident_stroke_followup < 10]
stopifnot(all(dat[(stroke_follow_lt10_Wales), latest_hospital_nation == "Wales"]))

dat[cad, on = .(eid), mortality_at_cad_followup := i.mortality_at_followup_date]
dat[cad, on = .(eid), lost_at_cad_followup := i.lost_to_followup_reason != ""]
dat[cad, on = .(eid), lost_at_cad_followup_reason := i.lost_to_followup_reason]
dat[cad, on = .(eid), cad_follow_lt10_Wales := !(incident_event) & incident_event_followup < 10 & incident_event_followup_date == latest_hospital_date & latest_hospital_nation == "Wales"]
dat[, cad_follow_lt10_Wales := !(incident_cad) & !(lost_at_cad_followup) & !(mortality_at_cad_followup) & incident_cad_followup < 10]
stopifnot(all(dat[(cad_follow_lt10_Wales), latest_hospital_nation == "Wales"]))

# Add in biochemistry biomarker data
bio <- fread("data/ukb/biomarkers/output/biomarkers.txt")
bio <- bio[visit_index == 0L]
bio[, visit_index := NULL]
dat <- merge(dat, bio, by="eid", all.x=TRUE)

# Flag people missing blood biochemistry data
bio_sinfo <- fread("data/ukb/biomarkers/output/samples_not_measured.txt")
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

# Add in NMR data
nmr <- fread("data/ukb/NMR_metabolomics/biomarker_measurements.txt")
nmr <- nmr[visit_index == 0] # baseline assessment only
nmr[, visit_index := NULL]
dat <- merge(dat, nmr, by="eid", all.x=TRUE)
dat[, no_nmr_data := FALSE]
dat[!nmr, on = .(eid), no_nmr_data := TRUE]

# Add PRSs
PRSs <- rbind(idcol="PRS", fill=TRUE,
  CAD_metaGRS = fread("data/ukb/PRS/CAD_metaGRS/CAD_metaGRS_PGS000018_7457de2d_UKBv3.sscore.gz"),
  Stroke_metaGRS = fread("data/ukb/PRS/Stroke_metaGRS/Stroke_metaGRS_PGS000039_519864bc_UKBv3.sscore.gz")
)
PRSs <- dcast(PRSs, IID ~ PRS, value.var="score_sum")
dat[PRSs, on = .(eid = IID), CAD_metaGRS := i.CAD_metaGRS]
dat[PRSs, on = .(eid = IID), Stroke_metaGRS := i.Stroke_metaGRS]
dat[, no_genetics := ifelse(is.na(CAD_metaGRS), TRUE, FALSE)]

# Set PRS to missing for people used to train the PRSs
prs_training <- unique(rbind(
  fread("data/ukb/PRS/sample_splits/CAD_metaGRS_training_samples.txt"),
  fread("data/ukb/PRS/sample_splits/Stroke_metaGRS_training_samples.txt")
))
setnames(prs_training, "eid")
dat[, prs_training_samples := FALSE]
dat[prs_training, on = .(eid), prs_training_samples := TRUE]

# -----------------------------------------------------
# Now do sample exclusions to derive analysis cohort
# -----------------------------------------------------

# Start building table of sample flowchart
sample_info <- data.table(step="Baseline (excl. withdrawals)", 
  samples=dat[,.N], CVD=NA_real_, CAD=NA_real_, Stroke=NA_real_,
  exited=NA_real_, exited_cvd=NA_real_, exited_cad=NA_real_, exited_stroke=NA_real_)

# Function to update sample information
update_sample_info <- function(step_name, dataset, last_dataset) {
  if (missing(dataset)) {
    dataset <- dat
  }
  current_samples <- dataset[,.N]
  current_cvd <- sum(dataset$incident_cvd, na.rm=TRUE)
  current_cad <- sum(dataset$incident_cad, na.rm=TRUE)
  current_stroke <- sum(dataset$incident_stroke, na.rm=TRUE)
  if (missing(last_dataset)) {
		last_samples <- sample_info[.N, samples]
		last_cvd <- sample_info[.N, CVD]
		last_cad <- sample_info[.N, CAD]
		last_stroke <- sample_info[.N, Stroke]
  } else {
    last_samples <- last_dataset[,.N]
    last_cvd <- last_dataset[, sum(incident_cvd, na.rm=TRUE)]
    last_cad <- last_dataset[, sum(incident_cad, na.rm=TRUE)]
    last_stroke <- last_dataset[, sum(incident_stroke, na.rm=TRUE)]
  }
  new_row <- data.table(
    step=step_name, 
    samples=current_samples, 
    CVD=current_cvd, 
    CAD=current_cad,
    Stroke=current_stroke,
    exited=last_samples - current_samples,
    exited_cvd=ifelse(is.na(last_cvd), NA, last_cvd - current_cvd),
    exited_cad=ifelse(is.na(last_cad), NA, last_cad - current_cad),
    exited_stroke=ifelse(is.na(last_stroke), NA, last_stroke - current_stroke)
  )
  sample_info <<- rbind(sample_info, new_row) 
}

# Filter to people with requisite data
dat_cpy <- copy(dat)

# Drop people with no linkage in electronic health records and update flow-chart
dat <- dat[!cvd[(ehr_linkage_withdrawn)], on = .(eid)]
update_sample_info("With EHR linkage")

# Filter to people with blood samples
dat <- dat[!(no_blood_sample)]
update_sample_info("With blood samples")

# Filter to people with NMR data
dat <- dat[!(no_nmr_data)]
update_sample_info("With NMR data")

# Drop people without genetic data
dat <- dat[!(no_genetics)]
update_sample_info("With linked genotype information")

# Add entry summarising steps
update_sample_info("With EHR, NMR, and genotype data", dat, dat_cpy)

# Flag data quality issues or other ineligibility criteria
dat_cpy <- copy(dat)

dat <- dat[!(prs_training_samples)]
update_sample_info("Not used for PRS training")

# Get information on NMR data missingness in non-derived biomarkers
non_derived <- ukbnmr::nmr_info[Type == "Non-derived", Biomarker]
nmr_miss <- apply(dat[,.SD,.SDcols=non_derived], 1, function(rr) { sum(is.na(rr)) / length(non_derived) })
miss <- dat[, .(eid, nmr_missingness=nmr_miss)]

ggdt <- rbind(
  data.table(type="NMR metabolomics", pct=miss$nmr_missingness)
)

g <- ggplot(ggdt, aes(x=pct, color=type)) +
  geom_density(trim=TRUE) +
  xlab("% sample measurements missing") +
  facet_wrap( ~ type, scales="free") +
  theme_bw() +
  theme(legend.position="bottom")
ggsave(g, width=6, height=4, units="in", file="data/cleaned/sample_missingness.png")

miss_thresh <- data.table(
  biomarkers_missing = 0:length(non_derived),
  pct_biomarkers_missing = 0:length(non_derived) / length(non_derived) * 100,
  samples_included = sapply(0:length(non_derived), function(ii) {  ggdt[pct <= ii/length(non_derived), .N] }),
  pct_samples_included = sapply(0:length(non_derived), function(ii) {  ggdt[pct <= ii/length(non_derived), .N] / ggdt[,.N] * 100 })
)
fwrite(miss_thresh, sep="\t", quote=FALSE, file="data/cleaned/missingness_threshold_consequences.txt")

# Get biomarker missingness reason
nmr_qc_flags <- fread("data/ukb/NMR_metabolomics/measurement_qc_flags.txt")
nmr_miss2 <- melt(dat, id.vars="eid", measure.vars=non_derived, variable.name="biomarker")[is.na(value)]
nmr_miss2[, value := NULL]
nmr_miss2 <- nmr_qc_flags[nmr_miss2, on = .(eid, biomarker)]
nmr_miss2 <- nmr_miss2[,.(QC_tag=strsplit(QC_tag, "; ")[[1]]), by=.(eid, biomarker)] # some NAs have multiple tags
nmr_miss_tags <- nmr_miss2[, .(n_missing=.N, pct_missing=round(.N/nmr_miss2[,.N,by=.(eid,biomarker)][,.N]*100, digits=2)), by=QC_tag]
nmr_miss_tags <- nmr_miss_tags[order(-pct_missing)] # note % missing is relative to number of NAs, not number of tags
fwrite(nmr_miss_tags, sep="\t", quote=FALSE, file="data/cleaned/missingness_reasons.txt")

sample_miss <- nmr_miss2[,.(n_missing=length(unique(biomarker)), pct_missing=round(length(unique(biomarker))/length(non_derived)*100, digits=2)), by=eid]
sample_miss_tags <- nmr_miss2[sample_miss, on = .(eid)][, .(.N, pct=round(.N/n_missing*100, digits=2)), by=.(eid, n_missing, pct_missing, QC_tag)] 
sample_miss_tags <- sample_miss_tags[order(-pct)][order(eid)][order(-pct_missing)]
fwrite(sample_miss_tags, sep="\t", quote=FALSE, file="data/cleaned/missingness_reason_by_sample.txt")

# Set missingness threshold as 5% (6 or more biomarkers with missing data): 
# This is the point where the majority of missing data arise due to outlier 
# plates of non-biological origin (93% of excluded samples) or technical 
# error (6% of excluded samples). Setting the threshold as 10% (as we had originally)
# only includes 231 additional samples 
dat_cpy2 <- copy(dat)
update_sample_info("Subset with complete NMR data:", dat[eid %in% miss[nmr_missingness == 0, eid]], dat[eid %in% miss[nmr_missingness == 0, eid]])
dat <- dat[!miss[nmr_missingness > 0.05], on = .(eid)]
update_sample_info("With <5% missing NMR data", dat, dat_cpy2)

# Flag issues with clinical biochemistry assays
dat <- dat[!(no_blood_biomarkers)]
update_sample_info("Clinical chemistry failed")

# Add new row flagging requisite data for analysis
update_sample_info("Passing data analysis QC", dat, dat_cpy)

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
update_sample_info("History of gestational diabetes included:", dat[(history_gestational_diabetes)], dat[(history_gestational_diabetes)]) 

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

sample_info[, exited_cad := ifelse(
  is.na(exited_cad), NA_character_,
  ifelse(exited == 0L, "0 (-%)",
  sprintf("%s (%s%%)", format(exited_cad, big.mark=","), round(exited_cad/(exited_cad + CVD)*100, digits=2))
))]

sample_info[, exited_stroke := ifelse(
  is.na(exited_stroke), NA_character_,
  ifelse(exited == 0L, "0 (-%)",
  sprintf("%s (%s%%)", format(exited_stroke, big.mark=","), round(exited_stroke/(exited_stroke + CVD)*100, digits=2))
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

sample_info[, CAD := ifelse(
  is.na(CAD), NA_character_,
  ifelse(samples == 0L, "0 (-%)",
  sprintf("%s (%s%%)", format(CAD, big.mark=","), round(CAD/samples*100, digits=2))
))]

sample_info[, Stroke := ifelse(
  is.na(Stroke), NA_character_,
  ifelse(samples == 0L, "0 (-%)",
  sprintf("%s (%s%%)", format(Stroke, big.mark=","), round(Stroke/samples*100, digits=2))
))]

sample_info[, samples := ifelse(is.na(samples), NA_character_, format(samples, big.mark=","))]

# Write out sample information
fwrite(sample_info, sep="\t", quote=FALSE, file="analyses/sample_flowchart.txt")

# Add SCORE2
dat[, SCORE2 := score2(sex, age, smoking, sbp, tchol, hdl, type="linear predictor")]

# We will perform nested cross-validation. Here, split the data into 5-folds balancing CVD cases status and sex.
# We will later split each 4/5ths of the data into 10-folds for elasticnet cross-validation, balancing by CAD or 
# Stroke case status.
dat[, foldgrp := paste(incident_cvd, sex)]
dat[, prediction_cv_foldid := createFolds(foldgrp, k=5, list=FALSE)]
dat[, foldgrp := NULL]

# Write out analysis cohort
fwrite(dat, sep="\t", quote=FALSE, file="data/cleaned/analysis_cohort.txt")

# Tabulate cohort information by CVD case status and sex
dat[, CAD_metaGRS := scale(CAD_metaGRS)]
dat[, Stroke_metaGRS := scale(Stroke_metaGRS)]
cvd_cohort_info <- dat[,.(
  samples = sprintf("%s (%s%%)", format(.N, big.mark=","), round(.N/dat[,.N]*100, digits=1)),
  assumed_ascvd_free = sprintf("%s (%s%%)", format(sum(is.na(ASCVD)), big.mark=","), round(sum(is.na(ASCVD))/.N*100, digits=2)),
  assumed_medication_free = sprintf("%s (%s%%)", format(sum(is.na(cholesterol_medication)), big.mark=","), round(sum(is.na(cholesterol_medication))/.N*100, digits=2)),
  men = sprintf("%s (%s%%)", format(sum(sex == "Male"), big.mark=","), round(sum(sex == "Male")/.N*100, digits=2)),
  age = sprintf("%s (%s)", round(median(age), digits=1), round(sd(age), digits=2)),
  BMI = sprintf("%s (%s)", round(median(na.omit(bmi)), digits=1), round(sd(na.omit(bmi)), digits=2)),
  SBP = sprintf("%s (%s)", round(median(sbp), digits=1), round(sd(sbp), digits=2)),
  smokers = sprintf("%s (%s%%)", format(sum(na.omit(smoking)), big.mark=","), round(sum(na.omit(smoking))/.N*100, digits=2)),
  assumed_non_smoking = sprintf("%s (%s%%)", format(sum(is.na(smoking)), big.mark=","), round(sum(is.na(smoking))/.N*100, digits=2)),
  history_gestational_diabetes = sprintf("%s (%s%%)", format(sum(na.omit(history_gestational_diabetes)), big.mark=","), round(sum(na.omit(history_gestational_diabetes))/dat[sex == "Female", .N]*100, digits=2)),
  hdl_cholesterol = sprintf("%s (%s)", round(median(hdl, na.rm=TRUE), digits=2), round(sd(hdl, na.rm=TRUE), digits=2)),
  total_cholesterol = sprintf("%s (%s)", round(median(tchol, na.rm=TRUE), digits=2), round(sd(tchol, na.rm=TRUE), digits=2)),
  ldl_cholesterol = sprintf("%s (%s)", round(median(ldl, na.rm=TRUE), digits=2), round(sd(ldl, na.rm=TRUE), digits=2)),
  SCORE2 = sprintf("%s (%s)", round(median(SCORE2), digits=2), round(sd(SCORE2), digits=2)),
  CAD_metaGRS = sprintf("%s (%s)", round(median(CAD_metaGRS, na.rm=TRUE), digits=2), round(sd(CAD_metaGRS, na.rm=TRUE), digits=2)),
  Stroke_metaGRS = sprintf("%s (%s)", round(median(Stroke_metaGRS, na.rm=TRUE), digits=2), round(sd(Stroke_metaGRS, na.rm=TRUE), digits=2)),
  Median_followup = sprintf("%s (%s)", round(median(incident_cvd * incident_cvd_followup), digits=1), round(sd(incident_cvd * incident_cvd_followup), digits=2)),
  primary_cause = sprintf("%s (%s%%)", format(sum(cvd_is_primary_cause), big.mark=","), round(sum(cvd_is_primary_cause)/sum(incident_cvd)*100, digits=2)),
  fatal = sprintf("%s (%s%%)", format(sum(incident_cvd_is_fatal), big.mark=","), round(sum(incident_cvd_is_fatal)/sum(incident_cvd)*100, digits=2)),
  Censored_lt_10yr = sprintf("%s (%s%%)", format(sum(!incident_cvd & incident_cvd_followup < 10), big.mark=","), round(sum(!incident_cvd & incident_cvd_followup < 10)/.N*100, digits=2)),
  Censored_fatal = sprintf("%s (%s%%)", format(sum(!incident_cvd & incident_cvd_followup < 10 & mortality_at_cvd_followup), big.mark=","), 
                                          round(sum(!incident_cvd & incident_cvd_followup < 10 & mortality_at_cvd_followup)/sum(!incident_cvd & incident_cvd_followup < 10)*100, digits=2)),
  Censored_max_Wales = sprintf("%s (%s%%)", format(sum(!incident_cvd & cvd_follow_lt10_Wales), big.mark=","), round(sum(!incident_cvd & cvd_follow_lt10_Wales)/sum(!incident_cvd & incident_cvd_followup < 10)*100, digits=2)),
  Censored_lost = sprintf("%s (%s%%)", format(sum(!incident_cvd & lost_at_cvd_followup  & !mortality_at_cvd_followup), big.mark=","), 
                                         round(sum(!incident_cvd & lost_at_cvd_followup  & !mortality_at_cvd_followup)/sum(!incident_cvd & incident_cvd_followup < 10)*100, digits=2))
), by=.(case_status=fcase(incident_cvd, "cvd", default="cvd non-case"))]

cad_cohort_info <- dat[!(prevalent_cad) | is.na(prevalent_cad),.(
  samples = sprintf("%s (%s%%)", format(.N, big.mark=","), round(.N/dat[,.N]*100, digits=1)),
  assumed_ascvd_free = sprintf("%s (%s%%)", format(sum(is.na(ASCVD)), big.mark=","), round(sum(is.na(ASCVD))/.N*100, digits=2)),
  assumed_medication_free = sprintf("%s (%s%%)", format(sum(is.na(cholesterol_medication)), big.mark=","), round(sum(is.na(cholesterol_medication))/.N*100, digits=2)),
  men = sprintf("%s (%s%%)", format(sum(sex == "Male"), big.mark=","), round(sum(sex == "Male")/.N*100, digits=2)),
  age = sprintf("%s (%s)", round(median(age), digits=1), round(sd(age), digits=2)),
  BMI = sprintf("%s (%s)", round(median(na.omit(bmi)), digits=1), round(sd(na.omit(bmi)), digits=2)),
  SBP = sprintf("%s (%s)", round(median(sbp), digits=1), round(sd(sbp), digits=2)),
  smokers = sprintf("%s (%s%%)", format(sum(na.omit(smoking)), big.mark=","), round(sum(na.omit(smoking))/.N*100, digits=2)),
  assumed_non_smoking = sprintf("%s (%s%%)", format(sum(is.na(smoking)), big.mark=","), round(sum(is.na(smoking))/.N*100, digits=2)),
  history_gestational_diabetes = sprintf("%s (%s%%)", format(sum(na.omit(history_gestational_diabetes)), big.mark=","), round(sum(na.omit(history_gestational_diabetes))/dat[sex == "Female", .N]*100, digits=2)),
  hdl_cholesterol = sprintf("%s (%s)", round(median(hdl, na.rm=TRUE), digits=2), round(sd(hdl, na.rm=TRUE), digits=2)),
  total_cholesterol = sprintf("%s (%s)", round(median(tchol, na.rm=TRUE), digits=2), round(sd(tchol, na.rm=TRUE), digits=2)),
  ldl_cholesterol = sprintf("%s (%s)", round(median(ldl, na.rm=TRUE), digits=2), round(sd(ldl, na.rm=TRUE), digits=2)),
  SCORE2 = sprintf("%s (%s)", round(median(SCORE2), digits=2), round(sd(SCORE2), digits=2)),
  CAD_metaGRS = sprintf("%s (%s)", round(median(CAD_metaGRS, na.rm=TRUE), digits=2), round(sd(CAD_metaGRS, na.rm=TRUE), digits=2)),
  Stroke_metaGRS = sprintf("%s (%s)", round(median(Stroke_metaGRS, na.rm=TRUE), digits=2), round(sd(Stroke_metaGRS, na.rm=TRUE), digits=2)),
  Median_followup = sprintf("%s (%s)", round(median(incident_cad * incident_cad_followup), digits=1), round(sd(incident_cad * incident_cad_followup), digits=2)),
  primary_cause = sprintf("%s (%s%%)", format(sum(cad_is_primary_cause), big.mark=","), round(sum(cad_is_primary_cause)/sum(incident_cad)*100, digits=2)),
  fatal = sprintf("%s (%s%%)", format(sum(incident_cad_is_fatal), big.mark=","), round(sum(incident_cad_is_fatal)/sum(incident_cad)*100, digits=2)),
  Censored_lt_10yr = sprintf("%s (%s%%)", format(sum(!incident_cad & incident_cad_followup < 10), big.mark=","), round(sum(!incident_cad & incident_cad_followup < 10)/.N*100, digits=2)),
  Censored_fatal = sprintf("%s (%s%%)", format(sum(!incident_cad & incident_cad_followup < 10 & mortality_at_cad_followup), big.mark=","), 
                                          round(sum(!incident_cad & incident_cad_followup < 10 & mortality_at_cad_followup)/sum(!incident_cad & incident_cad_followup < 10)*100, digits=2)),
  Censored_max_Wales = sprintf("%s (%s%%)", format(sum(!incident_cad & cad_follow_lt10_Wales), big.mark=","), round(sum(!incident_cad & cad_follow_lt10_Wales)/sum(!incident_cad & incident_cad_followup < 10)*100, digits=2)),
  Censored_lost = sprintf("%s (%s%%)", format(sum(!incident_cad & lost_at_cad_followup  & !mortality_at_cad_followup), big.mark=","), 
                                         round(sum(!incident_cad & lost_at_cad_followup  & !mortality_at_cad_followup)/sum(!incident_cad & incident_cad_followup < 10)*100, digits=2))
), by=.(case_status=fcase(incident_cad, "cad", default="cad non-case"))]

stroke_cohort_info <- dat[!(prevalent_stroke) | is.na(prevalent_stroke),.(
  samples = sprintf("%s (%s%%)", format(.N, big.mark=","), round(.N/dat[,.N]*100, digits=1)),
  assumed_ascvd_free = sprintf("%s (%s%%)", format(sum(is.na(ASCVD)), big.mark=","), round(sum(is.na(ASCVD))/.N*100, digits=2)),
  assumed_medication_free = sprintf("%s (%s%%)", format(sum(is.na(cholesterol_medication)), big.mark=","), round(sum(is.na(cholesterol_medication))/.N*100, digits=2)),
  men = sprintf("%s (%s%%)", format(sum(sex == "Male"), big.mark=","), round(sum(sex == "Male")/.N*100, digits=2)),
  age = sprintf("%s (%s)", round(median(age), digits=1), round(sd(age), digits=2)),
  BMI = sprintf("%s (%s)", round(median(na.omit(bmi)), digits=1), round(sd(na.omit(bmi)), digits=2)),
  SBP = sprintf("%s (%s)", round(median(sbp), digits=1), round(sd(sbp), digits=2)),
  smokers = sprintf("%s (%s%%)", format(sum(na.omit(smoking)), big.mark=","), round(sum(na.omit(smoking))/.N*100, digits=2)),
  assumed_non_smoking = sprintf("%s (%s%%)", format(sum(is.na(smoking)), big.mark=","), round(sum(is.na(smoking))/.N*100, digits=2)),
  history_gestational_diabetes = sprintf("%s (%s%%)", format(sum(na.omit(history_gestational_diabetes)), big.mark=","), round(sum(na.omit(history_gestational_diabetes))/dat[sex == "Female", .N]*100, digits=2)),
  hdl_cholesterol = sprintf("%s (%s)", round(median(hdl, na.rm=TRUE), digits=2), round(sd(hdl, na.rm=TRUE), digits=2)),
  total_cholesterol = sprintf("%s (%s)", round(median(tchol, na.rm=TRUE), digits=2), round(sd(tchol, na.rm=TRUE), digits=2)),
  ldl_cholesterol = sprintf("%s (%s)", round(median(ldl, na.rm=TRUE), digits=2), round(sd(ldl, na.rm=TRUE), digits=2)),
  SCORE2 = sprintf("%s (%s)", round(median(SCORE2), digits=2), round(sd(SCORE2), digits=2)),
  CAD_metaGRS = sprintf("%s (%s)", round(median(CAD_metaGRS, na.rm=TRUE), digits=2), round(sd(CAD_metaGRS, na.rm=TRUE), digits=2)),
  Stroke_metaGRS = sprintf("%s (%s)", round(median(Stroke_metaGRS, na.rm=TRUE), digits=2), round(sd(Stroke_metaGRS, na.rm=TRUE), digits=2)),
  Median_followup = sprintf("%s (%s)", round(median(incident_stroke * incident_stroke_followup), digits=1), round(sd(incident_stroke * incident_stroke_followup), digits=2)),
  primary_cause = sprintf("%s (%s%%)", format(sum(stroke_is_primary_cause), big.mark=","), round(sum(stroke_is_primary_cause)/sum(incident_stroke)*100, digits=2)),
  fatal = sprintf("%s (%s%%)", format(sum(incident_stroke_is_fatal), big.mark=","), round(sum(incident_stroke_is_fatal)/sum(incident_stroke)*100, digits=2)),
  Censored_lt_10yr = sprintf("%s (%s%%)", format(sum(!incident_stroke & incident_stroke_followup < 10), big.mark=","), round(sum(!incident_stroke & incident_stroke_followup < 10)/.N*100, digits=2)),
  Censored_fatal = sprintf("%s (%s%%)", format(sum(!incident_stroke & incident_stroke_followup < 10 & mortality_at_stroke_followup), big.mark=","), 
                                          round(sum(!incident_stroke & incident_stroke_followup < 10 & mortality_at_stroke_followup)/sum(!incident_stroke & incident_stroke_followup < 10)*100, digits=2)),
  Censored_max_Wales = sprintf("%s (%s%%)", format(sum(!incident_stroke & stroke_follow_lt10_Wales), big.mark=","), round(sum(!incident_stroke & stroke_follow_lt10_Wales)/sum(!incident_stroke & incident_stroke_followup < 10)*100, digits=2)),
  Censored_lost = sprintf("%s (%s%%)", format(sum(!incident_stroke & lost_at_stroke_followup  & !mortality_at_stroke_followup), big.mark=","), 
                                         round(sum(!incident_stroke & lost_at_stroke_followup  & !mortality_at_stroke_followup)/sum(!incident_stroke & incident_stroke_followup < 10)*100, digits=2))
), by=.(case_status=fcase(incident_stroke, "stroke", default="stroke non-case"))]

cohort_info <- rbind(cvd_cohort_info, cad_cohort_info, stroke_cohort_info)
cohort_info <- as.data.table(t(cohort_info), keep.rownames=TRUE)
cohort_info <- cohort_info[, c(1,2,3,5,4,7,6), with=FALSE]

fwrite(cohort_info, sep="\t", quote=FALSE, col.names=FALSE, file="analyses/cohort_information.txt")

# Note disease subtype analysis exclusions
excl_info <- rbind(
  data.table("Prevalent events excluded from CAD analyses", sprintf("%s (%s%%)", format(dat[(prevalent_cad), .N], big.mark=","), round(dat[(prevalent_cad), .N]/dat[,.N]*100, digits=2))),
  data.table("Assumed CAD free", sprintf("%s (%s%%)", format(dat[is.na(prevalent_cad), .N], big.mark=","), round(dat[is.na(prevalent_cad), .N]/dat[,.N]*100, digits=2))),
  data.table("Prevalent events excluded from stroke analyses", sprintf("%s (%s%%)", format(dat[(prevalent_stroke), .N], big.mark=","), round(dat[(prevalent_stroke), .N]/dat[,.N]*100, digits=2))),
  data.table("Assumed stroke free", sprintf("%s (%s%%)", format(dat[is.na(prevalent_stroke), .N], big.mark=","), round(dat[is.na(prevalent_stroke), .N]/dat[,.N]*100, digits=2)))
)
fwrite(excl_info, sep="\t", quote=FALSE, col.names=FALSE, file="analyses/sub_cohort_information.txt")

