library(data.table)
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
dat <- dat[, .(eid, assessment_date, assessment_centre, age, age_decimal, sex, bmi, genetic_white_british)]

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

# Add in family history of CVD
famhist <- fread("data/ukb/family_history/output/illness_of_first_degree_relatives.txt")
dat[famhist[visit_index == 0], on = .(eid), family_history_cvd := i.cardiovascular_disease]

# Add in diabetes status
diabetes <- fread("data/ukb/QDiabetes/output/qdiabetes.txt")
dat[diabetes, on = .(eid), diabetes := i.prevalent_diabetes] # T2D, T1D, or uncertain

# Add in prevalent vascular disease and incident CVD
cvd <- fread("data/ukb/endpoints/endpoints/CEU_CVD_10year/events_and_followup.txt")
cvd <- cvd[visit_index == 0]
dat[cvd, on = .(eid), prevalent_vascular_disease := i.prevalent_event]
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
dat[cvd, on = .(eid), earliest_hospital_nation := i.earliest_hospital_nation]

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

# Add in NMR release phase
nmr_sinfo <- fread("data/ukb/NMR_metabolomics/sample_processing_information.txt")
nmr_sinfo <- nmr_sinfo[visit_index == 0]
nmr_sinfo[, phase := ifelse(Shipment.Batch %like% "Phase_2", 2, 1)]
dat[nmr_sinfo, on = .(eid), NMR_release := fcase(
  i.phase == 1, "Phase 1 public",
  i.phase == 2, "Phase 2 pre-release")]

# Get information on NMR data missingness
nmr_miss <- apply(dat[,.SD,.SDcols=names(nmr)[-1]], 1, function(rr) { sum(is.na(rr)) / (ncol(nmr)-1) })
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

# Flag samples with excess missingness
dat[miss, on = .(eid), nmr_excess_miss := ifelse(nmr_missingness < 0.1, FALSE, TRUE)]
dat <- dat[!(nmr_excess_miss)]
update_sample_info("With <10% missing NMR data")

# Drop people already at high risk for cardiovascular disease
dat_cpy <- copy(dat) # for tracking sample and case exits aggregating multiple steps
dat <- dat[!(prevalent_vascular_disease) | is.na(prevalent_vascular_disease)]
update_sample_info("Without prevalent vascular disease")

dat <- dat[!(cholesterol_medication) | is.na(cholesterol_medication)]
update_sample_info("Not on lipid lowering medication")

dat <- dat[!(blood_pressure_medication) | is.na(blood_pressure_medication)]
update_sample_info("Not on blood pressure lowering medication")

update_sample_info("Not had CVD or evaluated as having high risk of CVD", dat, dat_cpy)

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

# Adjust PRS for PCs
dat[!is.na(CAD_metaGRS), CAD_metaGRS := scale(lm(CAD_metaGRS ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)$residuals)]
dat[!is.na(Stroke_metaGRS), Stroke_metaGRS := scale(lm(Stroke_metaGRS ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)$residuals)]

# Drop people who can't be jointly analysed with PRS
dat_cpy <- copy(dat)
dat <- dat[!is.na(PC1)]
update_sample_info("With linked genotype information")

dat <- dat[!is.na(CAD_metaGRS)]
update_sample_info("Not used for PRS training")
update_sample_info("Can be jointly analyzed with PRS", dat, dat_cpy)

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
update_sample_info("With known SBP")

dat <- dat[!(no_blood_sample)]
update_sample_info("With blood sample for biochemistry assays")

dat <- dat[!(no_blood_biomarkers)]
update_sample_info("With non-missing data for any blood biochemistry biomarkers")

dat <- dat[!is.na(hdl)]
update_sample_info("With known HDL cholesterol")

dat <- dat[!is.na(tchol)]
update_sample_info("With known total cholesterol")

dat <- dat[!is.na(crp)]
update_sample_info("With known CRP")

update_sample_info("With non-missing quantitative clinical risk factors", dat, dat_cpy)

# For survey-based clinical risk factors, set to FALSE those with missing data
dat[is.na(smoking), smoking := FALSE]
dat[is.na(family_history_cvd), family_history_cvd := FALSE]

# With-hold phase 2 data for testing
dat_p2 <- dat[NMR_release == "Phase 2 pre-release"]
dat_p1 <- dat[NMR_release == "Phase 1 public"]

update_sample_info("Phase 1 public release", dat_p1)
sample_info[.N, c("exited", "exited_cases") := .(NA, NA)]
update_sample_info("Phase 2 pre-release", dat_p2)
sample_info[.N, c("exited", "exited_cases") := .(NA, NA)]

# Split phase 1 data into training and test data
# Balance split by case/control status, type, and sex. Any more factors
# and the groups get too small to be meaningfully useful
dat_p1[, foldgrp := paste(incident_cvd, cvd_is_primary_cause, cvd_primarily_stroke, cvd_is_fatal, sex)]
dat_p1[, foldid := createFolds(foldgrp, k=2, list=FALSE)]

dat_p1[, partition := ifelse(foldid == 1, "train", "test")]
dat_p1[, c("foldgrp", "foldid") := NULL]
dat_p2[, partition := "test"]

train <- dat_p1[partition == "train"]
test <- rbind(dat_p1[partition == "test"], dat_p2)

dat[dat_p1, on = .(eid), partition := i.partition]
dat[dat_p2, on = .(eid), partition := i.partition]

update_sample_info("Training data", train)
sample_info[.N, c("exited", "exited_cases") := .(NA, NA)]
update_sample_info("Test data", test)
sample_info[.N, c("exited", "exited_cases") := .(NA, NA)]

# Write out cleaned data
fwrite(train, file="data/cleaned/training_data.txt")
fwrite(test, file="data/cleaned/test_data.txt")

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

# Tabulate cohort information by CVD case status and phase 1 / phase 2 data release
cohort_info <- dat[,.(
  samples = sprintf("%s (%s%%)", format(.N, big.mark=","), round(.N/dat[NMR_release == cohort, .N]*100, digits=1)),
  men = sprintf("%s (%s%%)", format(sum(sex == "Male"), big.mark=","), round(sum(sex == "Male")/.N*100, digits=1)),
  age = sprintf("%s (%s)", round(median(age), digits=1), round(sd(age), digits=2)),
  BMI = sprintf("%s (%s)", round(median(na.omit(bmi)), digits=1), round(sd(na.omit(bmi)), digits=2)),
  SBP = sprintf("%s (%s)", round(median(sbp), digits=1), round(sd(sbp), digits=2)),
  diabetics = sprintf("%s (%s%%)", format(sum(diabetes), big.mark=","), round(sum(diabetes)/.N*100, digits=2)),
  smokers = sprintf("%s (%s%%)", format(sum(smoking), big.mark=","), round(sum(smoking)/.N*100, digits=2)),
  total_cholesterol = sprintf("%s (%s)", round(median(tchol, na.rm=TRUE), digits=2), round(sd(tchol, na.rm=TRUE), digits=2)),
  hdl_cholesterol = sprintf("%s (%s)", round(median(hdl, na.rm=TRUE), digits=2), round(sd(hdl, na.rm=TRUE), digits=2)),
  ldl_cholesterol = sprintf("%s (%s)", round(median(ldl, na.rm=TRUE), digits=2), round(sd(ldl, na.rm=TRUE), digits=2)),
  CRP = sprintf("%s (%s)", round(median(crp, na.rm=TRUE), digits=2), round(sd(crp, na.rm=TRUE), digits=2)),
  CAD_metaGRS = sprintf("%s (%s)", round(median(CAD_metaGRS, na.rm=TRUE), digits=2), round(sd(CAD_metaGRS, na.rm=TRUE), digits=2)),
  Stroke_metaGRS = sprintf("%s (%s)", round(median(Stroke_metaGRS, na.rm=TRUE), digits=2), round(sd(Stroke_metaGRS, na.rm=TRUE), digits=2)),
  CVD_fatal = sprintf("%s (%s%%)", format(sum(cvd_is_fatal), big.mark=","), round(sum(cvd_is_fatal)/sum(incident_cvd)*100, digits=2)),
  primary_cause = sprintf("%s (%s%%)", format(sum(cvd_is_primary_cause), big.mark=","), round(sum(cvd_is_primary_cause)/sum(incident_cvd)*100, digits=2)),
  CHD = sprintf("%s (%s%%)", format(sum(cvd_primarily_chd), big.mark=","), round(sum(cvd_primarily_chd)/sum(incident_cvd)*100, digits=2)),
  CHD_fatal = sprintf("%s (%s%%)", format(sum(cvd_primarily_chd & cvd_is_fatal), big.mark=","), round(sum(cvd_primarily_chd & cvd_is_fatal)/sum(cvd_primarily_chd)*100, digits=2)),
  Stroke = sprintf("%s (%s%%)", format(sum(cvd_primarily_stroke), big.mark=","), round(sum(cvd_primarily_stroke)/sum(incident_cvd)*100, digits=2)),
  Stroke_fatal = sprintf("%s (%s%%)", format(sum(cvd_primarily_stroke & cvd_is_fatal), big.mark=","), round(sum(cvd_primarily_stroke & cvd_is_fatal)/sum(cvd_primarily_stroke)*100, digits=2)),
  Median_followup = sprintf("%s (%s)", round(median(incident_followup), digits=1), round(sd(incident_followup), digits=2)),
  Censored_lt_10yr = sprintf("%s (%s%%)", format(sum(incident_followup < 10), big.mark=","), round(sum(incident_followup < 10)/.N*100, digits=2)),
  Censored_fatal = sprintf("%s (%s%%)", format(sum(all_cause_mortality & incident_followup < 10), big.mark=","), round(sum(all_cause_mortality & incident_followup < 10)/.N*100, digits=2)),
  Censored_lost = sprintf("%s (%s%%)", format(sum(lost_to_followup != "" & incident_followup < 10 & !all_cause_mortality), big.mark=","), 
                                       round(sum(lost_to_followup != "" & incident_followup < 10 & !all_cause_mortality)/.N*100, digits=2)),
  Censored_max_Wales = sprintf("%s (%s%%)", format(sum(latest_hospital_nation == "Wales" & incident_followup < 10 & !all_cause_mortality & lost_to_followup == ""), big.mark=","), 
                                            round(sum(latest_hospital_nation == "Wales" & incident_followup < 10 & !all_cause_mortality & lost_to_followup == "")/.N*100, digits=2))
), by=.(cohort=NMR_release, case_status=ifelse(incident_cvd, "case", "control"))]
cohort_info <- cohort_info[order(case_status)][order(-cohort)]

fwrite(as.data.table(t(cohort_info), keep.rownames=TRUE), sep="\t", quote=FALSE, col.names=FALSE, file="analyses/cohort_information_by_case_status_and_NMR_release.txt")

# Tabulate cohort information by CVD case status and dataset split
cohort_info <- dat[,.(
  samples = sprintf("%s (%s%%)", format(.N, big.mark=","), round(.N/dat[partition == cohort, .N]*100, digits=1)),
  men = sprintf("%s (%s%%)", format(sum(sex == "Male"), big.mark=","), round(sum(sex == "Male")/.N*100, digits=1)),
  age = sprintf("%s (%s)", round(median(age), digits=1), round(sd(age), digits=2)),
  BMI = sprintf("%s (%s)", round(median(na.omit(bmi)), digits=1), round(sd(na.omit(bmi)), digits=2)),
  SBP = sprintf("%s (%s)", round(median(sbp), digits=1), round(sd(sbp), digits=2)),
  diabetics = sprintf("%s (%s%%)", format(sum(diabetes), big.mark=","), round(sum(diabetes)/.N*100, digits=2)),
  smokers = sprintf("%s (%s%%)", format(sum(smoking), big.mark=","), round(sum(smoking)/.N*100, digits=2)),
  total_cholesterol = sprintf("%s (%s)", round(median(tchol, na.rm=TRUE), digits=2), round(sd(tchol, na.rm=TRUE), digits=2)),
  hdl_cholesterol = sprintf("%s (%s)", round(median(hdl, na.rm=TRUE), digits=2), round(sd(hdl, na.rm=TRUE), digits=2)),
  ldl_cholesterol = sprintf("%s (%s)", round(median(ldl, na.rm=TRUE), digits=2), round(sd(ldl, na.rm=TRUE), digits=2)),
  CRP = sprintf("%s (%s)", round(median(crp, na.rm=TRUE), digits=2), round(sd(crp, na.rm=TRUE), digits=2)),
  CAD_metaGRS = sprintf("%s (%s)", round(median(CAD_metaGRS, na.rm=TRUE), digits=2), round(sd(CAD_metaGRS, na.rm=TRUE), digits=2)),
  Stroke_metaGRS = sprintf("%s (%s)", round(median(Stroke_metaGRS, na.rm=TRUE), digits=2), round(sd(Stroke_metaGRS, na.rm=TRUE), digits=2)),
  CVD_fatal = sprintf("%s (%s%%)", format(sum(cvd_is_fatal), big.mark=","), round(sum(cvd_is_fatal)/sum(incident_cvd)*100, digits=2)),
  primary_cause = sprintf("%s (%s%%)", format(sum(cvd_is_primary_cause), big.mark=","), round(sum(cvd_is_primary_cause)/sum(incident_cvd)*100, digits=2)),
  CHD = sprintf("%s (%s%%)", format(sum(cvd_primarily_chd), big.mark=","), round(sum(cvd_primarily_chd)/sum(incident_cvd)*100, digits=2)),
  CHD_fatal = sprintf("%s (%s%%)", format(sum(cvd_primarily_chd & cvd_is_fatal), big.mark=","), round(sum(cvd_primarily_chd & cvd_is_fatal)/sum(cvd_primarily_chd)*100, digits=2)),
  Stroke = sprintf("%s (%s%%)", format(sum(cvd_primarily_stroke), big.mark=","), round(sum(cvd_primarily_stroke)/sum(incident_cvd)*100, digits=2)),
  Stroke_fatal = sprintf("%s (%s%%)", format(sum(cvd_primarily_stroke & cvd_is_fatal), big.mark=","), round(sum(cvd_primarily_stroke & cvd_is_fatal)/sum(cvd_primarily_stroke)*100, digits=2)),
  Median_followup = sprintf("%s (%s)", round(median(incident_followup), digits=1), round(sd(incident_followup), digits=2)),
  Censored_lt_10yr = sprintf("%s (%s%%)", format(sum(incident_followup < 10), big.mark=","), round(sum(incident_followup < 10)/.N*100, digits=2)),
  Censored_fatal = sprintf("%s (%s%%)", format(sum(all_cause_mortality & incident_followup < 10), big.mark=","), round(sum(all_cause_mortality & incident_followup < 10)/.N*100, digits=2)),
  Censored_lost = sprintf("%s (%s%%)", format(sum(lost_to_followup != "" & incident_followup < 10 & !all_cause_mortality), big.mark=","),
                                       round(sum(lost_to_followup != "" & incident_followup < 10 & !all_cause_mortality)/.N*100, digits=2)),
  Censored_max_Wales = sprintf("%s (%s%%)", format(sum(latest_hospital_nation == "Wales" & incident_followup < 10 & !all_cause_mortality & lost_to_followup == ""), big.mark=","), 
                                            round(sum(latest_hospital_nation == "Wales" & incident_followup < 10 & !all_cause_mortality & lost_to_followup == "")/.N*100, digits=2))
), by=.(cohort=partition, case_status=ifelse(incident_cvd, "case", "control"))]
cohort_info <- cohort_info[order(case_status)][order(-cohort)]

fwrite(as.data.table(t(cohort_info), keep.rownames=TRUE), sep="\t", quote=FALSE, col.names=FALSE, file="analyses/cohort_information_by_case_status_and_dataset_split.txt")

# Tabulate cohort information by CVD case status dateset split, and phase 1 / phase 2 data release
cohort_info <- dat[,.(
  samples = sprintf("%s (%s%%)", format(.N, big.mark=","), round(.N/dat[NMR_release == cohort2 & partition == cohort1, .N]*100, digits=1)),
  men = sprintf("%s (%s%%)", format(sum(sex == "Male"), big.mark=","), round(sum(sex == "Male")/.N*100, digits=1)),
  age = sprintf("%s (%s)", round(median(age), digits=1), round(sd(age), digits=2)),
  BMI = sprintf("%s (%s)", round(median(na.omit(bmi)), digits=1), round(sd(na.omit(bmi)), digits=2)),
  SBP = sprintf("%s (%s)", round(median(sbp), digits=1), round(sd(sbp), digits=2)),
  diabetics = sprintf("%s (%s%%)", format(sum(diabetes), big.mark=","), round(sum(diabetes)/.N*100, digits=2)),
  smokers = sprintf("%s (%s%%)", format(sum(smoking), big.mark=","), round(sum(smoking)/.N*100, digits=2)),
  total_cholesterol = sprintf("%s (%s)", round(median(tchol, na.rm=TRUE), digits=2), round(sd(tchol, na.rm=TRUE), digits=2)),
  hdl_cholesterol = sprintf("%s (%s)", round(median(hdl, na.rm=TRUE), digits=2), round(sd(hdl, na.rm=TRUE), digits=2)),
  ldl_cholesterol = sprintf("%s (%s)", round(median(ldl, na.rm=TRUE), digits=2), round(sd(ldl, na.rm=TRUE), digits=2)),
  CRP = sprintf("%s (%s)", round(median(crp, na.rm=TRUE), digits=2), round(sd(crp, na.rm=TRUE), digits=2)),
  CAD_metaGRS = sprintf("%s (%s)", round(median(CAD_metaGRS, na.rm=TRUE), digits=2), round(sd(CAD_metaGRS, na.rm=TRUE), digits=2)),
  Stroke_metaGRS = sprintf("%s (%s)", round(median(Stroke_metaGRS, na.rm=TRUE), digits=2), round(sd(Stroke_metaGRS, na.rm=TRUE), digits=2)),
  CVD_fatal = sprintf("%s (%s%%)", format(sum(cvd_is_fatal), big.mark=","), round(sum(cvd_is_fatal)/sum(incident_cvd)*100, digits=2)),
  primary_cause = sprintf("%s (%s%%)", format(sum(cvd_is_primary_cause), big.mark=","), round(sum(cvd_is_primary_cause)/sum(incident_cvd)*100, digits=2)),
  CHD = sprintf("%s (%s%%)", format(sum(cvd_primarily_chd), big.mark=","), round(sum(cvd_primarily_chd)/sum(incident_cvd)*100, digits=2)),
  CHD_fatal = sprintf("%s (%s%%)", format(sum(cvd_primarily_chd & cvd_is_fatal), big.mark=","), round(sum(cvd_primarily_chd & cvd_is_fatal)/sum(cvd_primarily_chd)*100, digits=2)),
  Stroke = sprintf("%s (%s%%)", format(sum(cvd_primarily_stroke), big.mark=","), round(sum(cvd_primarily_stroke)/sum(incident_cvd)*100, digits=2)),
  Stroke_fatal = sprintf("%s (%s%%)", format(sum(cvd_primarily_stroke & cvd_is_fatal), big.mark=","), round(sum(cvd_primarily_stroke & cvd_is_fatal)/sum(cvd_primarily_stroke)*100, digits=2)),
  Median_followup = sprintf("%s (%s)", round(median(incident_followup), digits=1), round(sd(incident_followup), digits=2)),
  Censored_lt_10yr = sprintf("%s (%s%%)", format(sum(incident_followup < 10), big.mark=","), round(sum(incident_followup < 10)/.N*100, digits=2)),
  Censored_fatal = sprintf("%s (%s%%)", format(sum(all_cause_mortality & incident_followup < 10), big.mark=","), round(sum(all_cause_mortality & incident_followup < 10)/.N*100, digits=2)),
  Censored_lost = sprintf("%s (%s%%)", format(sum(lost_to_followup != "" & incident_followup < 10 & !all_cause_mortality), big.mark=","), 
                                       round(sum(lost_to_followup != "" & incident_followup < 10 & !all_cause_mortality)/.N*100, digits=2)),
  Censored_max_Wales = sprintf("%s (%s%%)", format(sum(latest_hospital_nation == "Wales" & incident_followup < 10 & !all_cause_mortality & lost_to_followup == ""), big.mark=","), 
                                            round(sum(latest_hospital_nation == "Wales" & incident_followup < 10 & !all_cause_mortality & lost_to_followup == "")/.N*100, digits=2))
), by=.(cohort1=partition, cohort2=NMR_release, case_status=ifelse(incident_cvd, "case", "control"))]
cohort_info <- cohort_info[order(case_status)][order(-cohort1)][order(cohort2)]

fwrite(as.data.table(t(cohort_info), keep.rownames=TRUE), sep="\t", quote=FALSE, col.names=FALSE, file="analyses/cohort_information_by_case_status_dataset_split_and_NMR_release.txt")
