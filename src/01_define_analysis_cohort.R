library(data.table)
library(ukbnmr)
library(ggplot2)
library(caret)

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

# Start building table of sample flowchart
sample_info <- data.table(step="Baseline (excl. withdrawals)", samples=dat[,.N], CVD=NA_real_, exited=NA_real_, exited_cases=NA_real_)

# Add in systolic blood pressure
bp <- fread("data/ukb/blood_pressure/output/blood_pressure.txt")
dat[bp[visit_index == 0], on = .(eid), sbp := i.sbp]

# Add in medication usage
meds_touchscreen <- fread("data/ukb/medication/output/medications_simple.txt")
meds_interview <- fread("data/ukb/medication/output/detailed_medications_summarised.txt")
meds <- merge(meds_touchscreen, meds_interview, by=c("eid", "visit_index"), all=TRUE)

dat[meds[visit_index == 0], on = .(eid), cholesterol_medication := i.cholesterol_medication | i.lipid_lowering_medication]
dat[meds[visit_index == 0], on = .(eid), blood_pressure_medication := i.blood_pressure_medication | i.hypertension_medication]

# Add in smoking status
smoking <- fread("data/ukb/smoking/output/smoking_status.txt")
dat[smoking[visit_index == 0], on = .(eid), smoking := i.current_smoker]

# Add in diabetes status
diabetes <- fread("data/ukb/QDiabetes/output/qdiabetes.txt")
diabetes <- diabetes[visit_index == 0, .(eid, prevalent_t2d=type_2_diabetes, prevalent_t1d=type_1_diabetes, uncertain_diabetes, history_gestational_diabetes)]
dat <- merge(dat, diabetes, by="eid", all.x=TRUE)

# Add in prevalent vascular disease
prev_cvd <- fread('data/ukb/endpoints/endpoints/CEU_prevalent_vascular_disease/events_and_followup.txt')
dat[prev_cvd[visit_index == 0], on = .(eid), prevalent_vascular_disease := i.prevalent_event]

# Add in incident CVD
cvd <- fread("data/ukb/endpoints/endpoints/CEU_CVD_10year/events_and_followup.txt")
cvd <- cvd[visit_index == 0]
dat[cvd, on = .(eid), incident_cvd := i.incident_event]
dat[cvd, on = .(eid), incident_followup := i.incident_event_followup]
dat[cvd, on = .(eid), incident_followup_date := i.incident_event_followup_date]
dat[cvd, on = .(eid), cvd_is_fatal := fcase(
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

# Determine whether CVD was (primarily) due to CHD or Stroke
dat[cvd, on = .(eid), cvd_primarily_chd := fcase(
  substr(i.incident_code, 0, 2) == "I2", TRUE,
  substr(i.incident_code, 0, 2) %in% c("I6", "F0"), FALSE,
  substr(i.incident_code, 0, 2) == "" & !is.na(i.incident_code), FALSE,
  is.na(i.incident_code), NA
)]

dat[cvd, on = .(eid), cvd_primarily_stroke := fcase(
  substr(i.incident_code, 0, 2) %in% c("I6", "F0"), TRUE,
  substr(i.incident_code, 0, 2) == "I2", FALSE,
  substr(i.incident_code, 0, 2) == "" & !is.na(i.incident_code), FALSE,
  is.na(i.incident_code), NA
)]

# Get information on where participants resided at baseline assessment, 
# maximum, and minimum follow-up available in hospital records (different
# hospital systems have different follow-up time available depending on
# nation of hospital)
dat[cvd, on = .(eid), latest_hospital_nation := i.latest_hospital_nation]
dat[prev_cvd, on = .(eid), earliest_hospital_nation := i.earliest_hospital_nation]

# Get additional information
dat[cvd, on = .(eid), all_cause_mortality := i.mortality_at_followup_date]
dat[cvd, on = .(eid), lost_to_followup := i.lost_to_followup_reason]

# Function to update sample information
update_sample_info <- function(step_name, dataset, last_dataset) {
  if (missing(dataset)) {
    dataset <- dat
  }
  current_samples <- dataset[,.N]
  current_cases <- sum(dataset$incident_cvd, na.rm=TRUE)
  if (missing(last_dataset)) {
		last_samples <- sample_info[.N, samples]
		last_cases <- sample_info[.N, CVD]
  } else {
    last_samples <- last_dataset[,.N]
    last_cases <- last_dataset[, sum(incident_cvd, na.rm=TRUE)]
  }
  new_row <- data.table(
    step=step_name, samples=current_samples, CVD=current_cases, 
    exited=last_samples - current_samples,
    exited_cases=ifelse(is.na(last_cases), NA, last_cases - current_cases)
  )
  sample_info <<- rbind(sample_info, new_row) 
}

# Drop people with no linkage in electronic health records and update flow-chart
dat <- dat[!cvd[(ehr_linkage_withdrawn)], on = .(eid)]
update_sample_info("With EHR linkage")

# Add in NMR data and filter to participants with that data
nmr <- fread("data/ukb/NMR_metabolomics/biomarker_measurements.txt")
nmr <- nmr[visit_index == 0] # baseline assessment only
nmr[, visit_index := NULL]
dat <- dat[nmr, on = .(eid), nomatch=0]
update_sample_info("With NMR data")

# Only include people in the age range eligible for SCORE2
dat_cpy <- copy(dat) # for tracking sample and case exits aggregating multiple steps
dat <- dat[age >= 40]
update_sample_info("40 years or older")

dat <- dat[age <= 69]
update_sample_info("69 years or younger")
update_sample_info("Eligible age for SCORE2 risk prediction", dat, dat_cpy) 

# Drop people otherwise ineligible for SCORE2 risk prediction due to disease or medication
dat_cpy <- copy(dat) # for tracking sample and case exits aggregating multiple steps
dat <- dat[!(prevalent_vascular_disease) | is.na(prevalent_vascular_disease)]
update_sample_info("Without prevalent vascular disease")
update_sample_info("Missing/unclassifiable treated as non-prevalent", dat[is.na(prevalent_vascular_disease)], dat[is.na(prevalent_vascular_disease)])

dat_cpy2 <- copy(dat)
dat <- dat[!(cholesterol_medication) | is.na(cholesterol_medication)]
update_sample_info("Not on lipid lowering medication", dat, dat_cpy2)
update_sample_info("Missing data treated as medication-free", dat[is.na(cholesterol_medication)], dat[is.na(cholesterol_medication)])

dat_cpy2 <- copy(dat)
dat <- dat[!(prevalent_t2d) | is.na(prevalent_t2d)]
update_sample_info("Without type 2 diabetes", dat, dat_cpy2)
update_sample_info("Type 1 diabetics included:", dat[(prevalent_t1d)], dat[(prevalent_t1d)])
update_sample_info("Uncertain diabetes included:", dat[(uncertain_diabetes)], dat[(uncertain_diabetes)])
update_sample_info("History of gestational diabetes included:", dat[(history_gestational_diabetes)], dat[(history_gestational_diabetes)]) 

update_sample_info("Eligible disease and medication history for SCORE2 risk prediction", dat, dat_cpy)

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

# Drop people missing quantitative conventional risk factors
dat_cpy <- copy(dat)

dat <- dat[!is.na(sbp)]
update_sample_info("With known SBP", dat, dat_cpy)

dat <- dat[!(no_blood_sample)]
update_sample_info("With blood sample for biochemistry assays")

update_sample_info("Missing HDL cholesterol:", dat[is.na(hdl)], dat[is.na(hdl)])
update_sample_info("Missing Total cholesterol:", dat[is.na(tchol)], dat[is.na(tchol)])
update_sample_info("Missing LDL cholesterol:", dat[is.na(ldl)], dat[is.na(ldl)])
update_sample_info("Missing HDL, LDL, and total cholesterol:", dat[is.na(hdl) & is.na(tchol) & is.na(ldl)], dat[is.na(hdl) & is.na(tchol) & is.na(ldl)])
update_sample_info("Missing all blood biochemistry biomarkers:", dat[(no_blood_biomarkers)], dat[(no_blood_biomarkers)])
update_sample_info("Missing 2/3 cholesterols:", 
  dat[(is.na(hdl) & is.na(ldl) & !is.na(tchol)) | (!is.na(hdl) & is.na(ldl) & is.na(tchol)) | (is.na(hdl) & !is.na(ldl) & is.na(tchol))],
  dat[(is.na(hdl) & is.na(ldl) & !is.na(tchol)) | (!is.na(hdl) & is.na(ldl) & is.na(tchol)) | (is.na(hdl) & !is.na(ldl) & is.na(tchol))])
update_sample_info("Missing HDL and total cholesterol:", dat[is.na(hdl) & is.na(tchol) & !is.na(ldl)], dat[is.na(hdl) & is.na(tchol) & !is.na(ldl)])
update_sample_info("Missing HDL and LDL cholesterol:", dat[is.na(hdl) & !is.na(tchol) & is.na(ldl)], dat[is.na(hdl) & !is.na(tchol) & is.na(ldl)])
update_sample_info("Missing LDL and total cholesterol:", dat[!is.na(hdl) & is.na(tchol) & is.na(ldl)], dat[!is.na(hdl) & is.na(tchol) & is.na(ldl)])
update_sample_info("Missing 1/3 cholesterols:", 
  dat[(is.na(hdl) & !is.na(ldl) & !is.na(tchol)) | (!is.na(hdl) & is.na(ldl) & !is.na(tchol)) | (!is.na(hdl) & !is.na(ldl) & is.na(tchol))],
  dat[(is.na(hdl) & !is.na(ldl) & !is.na(tchol)) | (!is.na(hdl) & is.na(ldl) & !is.na(tchol)) | (!is.na(hdl) & !is.na(ldl) & is.na(tchol))])
update_sample_info("Missing only HDL cholesterol:", dat[is.na(hdl) & !is.na(tchol) & !is.na(ldl)], dat[is.na(hdl) & !is.na(tchol) & !is.na(ldl)])
update_sample_info("Missing only total cholesterol:", dat[!is.na(hdl) & is.na(tchol) & !is.na(ldl)], dat[!is.na(hdl) & is.na(tchol) & !is.na(ldl)])
update_sample_info("Missing only LDL cholesterol:", dat[!is.na(hdl) & !is.na(tchol) & is.na(ldl)], dat[!is.na(hdl) & !is.na(tchol) & is.na(ldl)])

dat_cpy2 <- copy(dat)
dat <- dat[!is.na(hdl) & !is.na(tchol) & !is.na(ldl)]
update_sample_info("With non-missing data for HDL, total, and LDL cholesterol", dat, dat_cpy2)

update_sample_info("With non-missing quantitative clinical risk factors", dat, dat_cpy)

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

# Remove samples with >10% missing (i.e. 10 or more) NMR biomarkers (~2% of samples - 90% have complete data)
update_sample_info("Subset with complete NMR data:", dat[eid %in% miss[nmr_missingness == 0, eid]], dat)
dat_cpy <- copy(dat)
dat <- dat[!miss[nmr_missingness > 0.1], on = .(eid)]
update_sample_info("With <10% missing NMR data", dat, dat_cpy)

# Load and add PCs and add
pcs <- fread("data/ukb/genetics/reference_files/ukb_sqc_v2.txt")
pcs <- pcs[,.(eid, chip=genotyping.array, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)]
dat <- merge(dat, pcs, by="eid", all.x=TRUE)

# Add PRSs
PRSs <- rbind(idcol="PRS", fill=TRUE,
  CAD_metaGRS = fread("data/ukb/PRS/CAD_metaGRS/CAD_metaGRS_PGS000018_7457de2d_UKBv3.sscore.gz"),
  Stroke_metaGRS = fread("data/ukb/PRS/Stroke_metaGRS/Stroke_metaGRS_PGS000039_519864bc_UKBv3.sscore.gz")
)
PRSs <- dcast(PRSs, IID ~ PRS, value.var="score_sum")
dat[PRSs, on = .(eid = IID), CAD_metaGRS := i.CAD_metaGRS]
dat[PRSs, on = .(eid = IID), Stroke_metaGRS := i.Stroke_metaGRS]

# Set PRS to missing for people used to train the PRSs
prs_training <- unique(rbind(
  fread("data/ukb/PRS/sample_splits/CAD_metaGRS_training_samples.txt"),
  fread("data/ukb/PRS/sample_splits/Stroke_metaGRS_training_samples.txt")
))
setnames(prs_training, "eid")
dat[prs_training, on = .(eid), CAD_metaGRS := NA]
dat[prs_training, on = .(eid), Stroke_metaGRS := NA]

# Drop people who can't be jointly analysed with PRS
dat_cpy <- copy(dat)
dat <- dat[!is.na(PC1)]
update_sample_info("With linked genotype information")

dat <- dat[!is.na(CAD_metaGRS)]
update_sample_info("Not used for PRS training")
update_sample_info("Can be jointly analyzed with PRS", dat, dat_cpy)

# Add percentages to sample flowchart
sample_info[, exited_cases := ifelse(
  is.na(exited_cases), NA_character_,
  ifelse(exited == 0L, "0 (0%)",
  sprintf("%s (%s%%)", format(exited_cases, big.mark=","), round(exited_cases/(exited_cases + CVD)*100, digits=2))
))]

sample_info[, exited := ifelse(
  is.na(exited), NA_character_,
  ifelse(exited == 0L, "0 (0%)",
  sprintf("%s (%s%%)", format(exited, big.mark=","), round(exited/(exited + samples)*100, digits=2))
))]

sample_info[, CVD := ifelse(
  is.na(CVD), NA_character_,
  sprintf("%s (%s%%)", format(CVD, big.mark=","), round(CVD/samples*100, digits=2))
)]

sample_info[, samples := ifelse(is.na(samples), NA_character_, format(samples, big.mark=","))]

# Write out sample information
fwrite(sample_info, sep="\t", quote=FALSE, file="analyses/sample_flowchart.txt")

# Split data for nested cross-validation, balancing case control status and sex
dat[, foldgrp := paste(incident_cvd, cvd_primarily_stroke, sex)]
dat[, foldid := createFolds(foldgrp, k=5*10, list=FALSE) - 1]
dat[, elasticnet_cv_foldid := foldid %% 10]
dat[, prediction_cv_foldid := floor(foldid / 10)]

# Write out analysis cohort
fwrite(dat, file="data/cleaned/analysis_cohort.txt")

# Tabulate cohort information by CVD case status and sex
dat[, CAD_metaGRS := scale(CAD_metaGRS)]
dat[, Stroke_metaGRS := scale(Stroke_metaGRS)]
cohort_info <- dat[,.(
  samples = sprintf("%s (%s%%)", format(.N, big.mark=","), round(.N/dat[,.N]*100, digits=1)),
  men = sprintf("%s (%s%%)", format(sum(sex == "Male"), big.mark=","), round(sum(sex == "Male")/.N*100, digits=2)),
  age = sprintf("%s (%s)", round(median(age), digits=1), round(sd(age), digits=2)),
  BMI = sprintf("%s (%s)", round(median(na.omit(bmi)), digits=1), round(sd(na.omit(bmi)), digits=2)),
  SBP = sprintf("%s (%s)", round(median(sbp), digits=1), round(sd(sbp), digits=2)),
  smokers = sprintf("%s (%s%%)", format(sum(na.omit(smoking)), big.mark=","), round(sum(na.omit(smoking))/.N*100, digits=2)),
  assumed_medication_free = sprintf("%s (%s%%)", format(sum(is.na(cholesterol_medication)), big.mark=","), round(sum(is.na(cholesterol_medication))/.N*100, digits=2)),
  assumed_cvd_free = sprintf("%s (%s%%)", format(sum(is.na(prevalent_vascular_disease)), big.mark=","), round(sum(is.na(prevalent_vascular_disease))/.N*100, digits=2)),
  assumed_non_smoking = sprintf("%s (%s%%)", format(sum(is.na(smoking)), big.mark=","), round(sum(is.na(smoking))/.N*100, digits=2)),
  uncertain_diabetes = sprintf("%s (%s%%)", format(sum(uncertain_diabetes), big.mark=","), round(sum(uncertain_diabetes)/.N*100, digits=2)),
  prevalent_t1d = sprintf("%s (%s%%)", format(sum(na.omit(prevalent_t1d)), big.mark=","), round(sum(na.omit(prevalent_t1d))/.N*100, digits=2)),
  history_gestational_diabetes = sprintf("%s (%s%%)", format(sum(na.omit(history_gestational_diabetes)), big.mark=","), round(sum(na.omit(history_gestational_diabetes))/dat[sex == "Female", .N]*100, digits=2)),
  hdl_cholesterol = sprintf("%s (%s)", round(median(hdl, na.rm=TRUE), digits=2), round(sd(hdl, na.rm=TRUE), digits=2)),
  total_cholesterol = sprintf("%s (%s)", round(median(tchol, na.rm=TRUE), digits=2), round(sd(tchol, na.rm=TRUE), digits=2)),
  ldl_cholesterol = sprintf("%s (%s)", round(median(ldl, na.rm=TRUE), digits=2), round(sd(ldl, na.rm=TRUE), digits=2)),
  CAD_metaGRS = sprintf("%s (%s)", round(median(CAD_metaGRS, na.rm=TRUE), digits=2), round(sd(CAD_metaGRS, na.rm=TRUE), digits=2)),
  Stroke_metaGRS = sprintf("%s (%s)", round(median(Stroke_metaGRS, na.rm=TRUE), digits=2), round(sd(Stroke_metaGRS, na.rm=TRUE), digits=2)),
  fatal = sprintf("%s (%s%%)", format(sum(cvd_is_fatal), big.mark=","), round(sum(cvd_is_fatal)/sum(incident_cvd)*100, digits=2)),
  primary_cause = sprintf("%s (%s%%)", format(sum(cvd_is_primary_cause), big.mark=","), round(sum(cvd_is_primary_cause)/sum(incident_cvd)*100, digits=2)),
  Median_followup = sprintf("%s (%s)", round(median(incident_followup), digits=1), round(sd(incident_followup), digits=2)),
  Censored_lt_10yr = sprintf("%s (%s%%)", format(sum(incident_followup < 10), big.mark=","), round(sum(incident_followup < 10)/.N*100, digits=2)),
  Censored_fatal = sprintf("%s (%s%%)", format(sum(all_cause_mortality & incident_followup < 10), big.mark=","), round(sum(all_cause_mortality & incident_followup < 10)/.N*100, digits=2)),
  Censored_lost = sprintf("%s (%s%%)", format(sum(lost_to_followup != "" & incident_followup < 10 & !all_cause_mortality), big.mark=","), 
                                       round(sum(lost_to_followup != "" & incident_followup < 10 & !all_cause_mortality)/.N*100, digits=2)),
  Censored_max_Wales = sprintf("%s (%s%%)", format(sum(latest_hospital_nation == "Wales" & incident_followup < 10 & !all_cause_mortality & lost_to_followup == ""), big.mark=","), 
                                            round(sum(latest_hospital_nation == "Wales" & incident_followup < 10 & !all_cause_mortality & lost_to_followup == "")/.N*100, digits=2))
), by=.(case_status=fcase(cvd_primarily_chd, "chd", cvd_primarily_stroke, "stroke", default="non-case"))]
cohort_info <- cohort_info[order(factor(case_status, levels=c("chd", "stroke", "non-case")))]
cohort_info <- as.data.table(t(cohort_info), keep.rownames=TRUE)

fwrite(cohort_info, sep="\t", quote=FALSE, col.names=FALSE, file="analyses/cohort_information.txt")
