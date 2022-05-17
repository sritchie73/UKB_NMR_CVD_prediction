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

# Start building table of sample flowchart
sample_info <- data.table(step="Baseline (excl. withdrawals)", samples=dat[,.N], CVD=NA_real_, exited=NA_real_, exited_cases=NA_real_)

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
update_sample_info <- function(step_name, dataset) {
  if (missing(dataset)) {
    dataset <- dat
  }
  current_samples <- dataset[,.N]
  current_cases <- sum(dataset$incident_cvd, na.rm=TRUE)
  last_samples <- sample_info[.N, samples]
  last_cases <- sample_info[.N, CVD]
  new_row <- data.table(
    step=step_name, samples=current_samples, CVD=current_cases, 
    exited=last_samples - current_samples,
    exited_cases=fcase(
      is.na(last_cases), NA_real_,
      default = last_cases - current_cases
    )
  )
  sample_info <<- rbind(sample_info, new_row) 
}

# Drop people with no linkage in electronic health records and update flow-chart
dat <- dat[!cvd[(ehr_linkage_withdrawn)], on = .(eid)]
update_sample_info("With EHR linkage")

# Add in NMR data and filter to participants with that data
nmr <- fread("data/ukb/NMR_metabolomics/output/nmr_techadj.txt")
nmr <- nmr[visit_index == 0L] # baseline assessment only
nmr[, visit_index := NULL]
dat <- dat[nmr, on = .(eid), nomatch=0]
update_sample_info("With NMR data")

# Add in biomarker data
bio <- fread("data/ukb/biomarkers/output/biomarkers.txt")
bio <- bio[visit_index == 0L]
bio[, visit_index := NULL]
bio[, hba1c_pct := NULL] # just keep hba1c
bio[, fasting_glucose := NULL] # correction of glucose for fasting time not in scope for this project
bio[, c("uriacc", "urianac", "uriamac", "uriakc") := NULL] # drop urine biomarkers
bio <- bio[apply(bio, 1, function(rr) { sum(is.na(rr)) != (ncol(bio)-1) })] # remove N=7,651 people with no blood biomarkers
dat <- dat[bio, on = .(eid), nomatch=0]
update_sample_info("With biomarker data")

# Generate updated plots on sample missingness
bio_miss <- apply(dat[,.SD,.SDcols=names(bio)[-1]], 1, function(rr) { sum(is.na(rr)) / (ncol(bio)-1) })
nmr_miss <- apply(dat[,.SD,.SDcols=names(nmr)[-1]], 1, function(rr) { sum(is.na(rr)) / (ncol(nmr)-1) })
miss <- dat[, .(eid, nmr_missingness=nmr_miss, biomarker_missingness=bio_miss)]

ggdt <- rbind(
  data.table(type="Clinical biochemistry", pct=bio_miss),
  data.table(type="NMR metabolomics", pct=nmr_miss)
)

g <- ggplot(ggdt, aes(x=pct, color=type)) +
  geom_density(trim=TRUE) +
  xlab("% sample measurements missing") +
  facet_wrap( ~ type, scales="free") +
  theme_bw() +
  theme(legend.position="bottom")
ggsave(g, width=13, height=4, units="in", file="data/cleaned/sample_missingness.png")

# Drop samples with excess missingness in NMR data (i.e. due to outlier plates of
# non-biological origin)
dat <- dat[miss[nmr_missingness < 0.1, .(eid)], on = .(eid), nomatch=0]
update_sample_info("< 10% missingness in NMR data")

# Drop samples with excess missingness in clinical biochemistry data
dat <- dat[miss[biomarker_missingness < 0.1, .(eid)], on = .(eid), nomatch=0]
update_sample_info("< 10% missingness in biomarker data")

# Drop people already at high risk for cardiovascular disease
dat <- dat[!is.na(prevalent_vascular_disease)]
update_sample_info("With known history for prevalent vascular disease")

dat <- dat[!(prevalent_vascular_disease)]
update_sample_info("Without prevalent vascular disease")

dat <- dat[!is.na(cholesterol_medication)]
update_sample_info("With non-missing lipid lowering medication status")

dat <- dat[!(cholesterol_medication)]
update_sample_info("Not on lipid lowering medication")

dat <- dat[!is.na(blood_pressure_medication)]
update_sample_info("With non-missing blood pressure medication status")

dat <- dat[!(blood_pressure_medication)]
update_sample_info("Not on blood pressure lowering medication")

# Drop people missing conventional risk factors
dat <- dat[!is.na(age)]
update_sample_info("With known age")

dat <- dat[!is.na(sex)]
update_sample_info("With known sex")

dat <- dat[!is.na(tchol)]
update_sample_info("With non-missing total cholesterol")

dat <- dat[!is.na(hdl)]
update_sample_info("With non-missing HDL cholesterol")

dat <- dat[!is.na(sbp)]
update_sample_info("With known SBP")
dat <- dat[!is.na(diabetes)]
update_sample_info("With known diabetes status")

dat <- dat[!is.na(smoking)]
update_sample_info("With known smoking status")

dat <- dat[!is.na(family_history_cvd)]
update_sample_info("With known family history")

dat <- dat[!is.na(bmi)]
update_sample_info("With known BMI")

update_sample_info("With non-missing conventional risk factors")

# Load and add PCs (and thereby filter to people with genetic data)
pcs <- fread("data/ukb/genetics/reference_files/ukb_sqc_v2.txt")
pcs <- pcs[,.(eid, chip=genotyping.array, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)]
dat <- dat[pcs, on = .(eid), nomatch=0]
update_sample_info("With genetics")

# Filter to white british ancestry subset
dat <- dat[(genetic_white_british)]
update_sample_info("In White British ancestry subset")

# Add PRSs
PRSs <- rbind(idcol="PRS", fill=TRUE,
  CAD_metaGRS = fread("data/ukb/PRS/CAD_metaGRS/CAD_metaGRS_PGS000018_b097e681_P7439_from_dosage.sscore.gz"),
  Stroke_metaGRS = fread("data/ukb/PRS/Stroke_metaGRS/Stroke_metaGRS_PGS000039_6a7832a2_P7439_from_dosage.sscore.gz")
)
PRSs <- dcast(PRSs, IID ~ PRS, value.var="score_sum")
dat <- dat[PRSs, on = .(eid = IID), nomatch=0]
update_sample_info("With PGS")

# Drop people used to train the PRSs
prs_training <- unique(rbind(
  fread("data/ukb/PRS/sample_splits/CAD_metaGRS_training_samples.txt"),
  fread("data/ukb/PRS/sample_splits/Stroke_metaGRS_training_samples.txt")
))
setnames(prs_training, "eid")
dat <- dat[!prs_training, on = .(eid)]
update_sample_info("Not used for PGS training")

# Prune out first-degree relatives. When choosing from pairs of related samples, prioritise keeping CVD cases
# over non-cases.
kinship <- fread("data/ukb/genetics/reference_files/kinship_relatedness.txt")
kinship <- kinship[Kinship > 0.0884] # cutoff from KING manual http://people.virginia.edu/~wc9c/KING/manual.html
kinship <- kinship[ID1 %in% dat$eid & ID2 %in% dat$eid]
kinship[dat, on = .(ID1=eid), ID1_CVD := i.incident_cvd]
kinship[dat, on = .(ID2=eid), ID2_CVD := i.incident_cvd]
kinship[, to_drop := ifelse(ID1_CVD & !(ID2_CVD), "ID1", "ID2")]
kinship[to_drop == "ID1", to_drop_eid := ID1]
kinship[to_drop == "ID2", to_drop_eid := ID2]
dat <- dat[!kinship, on = .(eid=to_drop_eid)]
update_sample_info("Pruned for first-degree relatives")

# Adjust PRS for PCs
dat[, CAD_metaGRS := scale(lm(CAD_metaGRS ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)$residuals)]
dat[, Stroke_metaGRS := scale(lm(Stroke_metaGRS ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)$residuals)]

# Split into training and test data
# Balance split by case/control status, sex, prevalent diabetes, smoking status, recruitment centre, and nation
dat[, foldgrp := paste(incident_cvd, cvd_is_primary_cause, cvd_is_fatal, cvd_primarily_chd, cvd_primarily_stroke, 
                       earliest_hospital_nation, latest_hospital_nation, assessment_centre, sex, diabetes, smoking, family_history_cvd)]
dat[, foldid := createFolds(foldgrp, k=2, list=FALSE)]

dat[, partition := ifelse(foldid == 1, "train", "test")]

test <- dat[partition == "test"]
train <- dat[partition == "train"]

update_sample_info("Training data", train)
sample_info[.N, c("exited", "exited_cases") := .(NA, NA)]
update_sample_info("Test data", test)
sample_info[.N, c("exited", "exited_cases") := .(NA, NA)]

# Remove foldgrp columns
train[, c("foldgrp", "foldid") := NULL]
test[, c("foldgrp", "foldid") := NULL]
 
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

# Tabulate cohort information by sex
cohort_info <- dat[,.(
  samples = sprintf("%s (%s%%)", format(.N, big.mark=","), round(.N/dat[partition == cohort, .N]*100, digits=1)),
  age = sprintf("%s (%s)", round(median(age), digits=1), round(sd(age), digits=2)),
  BMI = sprintf("%s (%s)", round(median(na.omit(bmi)), digits=1), round(sd(na.omit(bmi)), digits=2)),
  SBP = sprintf("%s (%s)", round(median(sbp), digits=1), round(sd(sbp), digits=2)),
  diabetics = sprintf("%s (%s%%)", format(sum(diabetes), big.mark=","), round(sum(diabetes)/.N*100, digits=2)),
  smokers = sprintf("%s (%s%%)", format(sum(smoking), big.mark=","), round(sum(smoking)/.N*100, digits=2)),
  total_cholesterol = sprintf("%s (%s)", round(median(tchol), digits=2), round(sd(tchol), digits=2)),
  hdl_cholesterol = sprintf("%s (%s)", round(median(hdl), digits=2), round(sd(hdl), digits=2)),
  CVD = sprintf("%s (%s%%)", format(sum(incident_cvd), big.mark=","), round(sum(incident_cvd)/.N*100, digits=2)),
  CVD_fatal = sprintf("%s (%s%%)", format(sum(cvd_is_fatal), big.mark=","), round(sum(cvd_is_fatal)/sum(incident_cvd)*100, digits=2)),
  primary_cause = sprintf("%s (%s%%)", format(sum(cvd_is_primary_cause), big.mark=","), round(sum(cvd_is_primary_cause)/sum(incident_cvd)*100, digits=2)),
  CHD = sprintf("%s (%s%%)", format(sum(cvd_primarily_chd), big.mark=","), round(sum(cvd_primarily_chd)/sum(incident_cvd)*100, digits=2)),
  CHD_fatal = sprintf("%s (%s%%)", format(sum(cvd_primarily_chd & cvd_is_fatal), big.mark=","), round(sum(cvd_primarily_chd & cvd_is_fatal)/sum(cvd_primarily_chd)*100, digits=2)),
  Stroke = sprintf("%s (%s%%)", format(sum(cvd_primarily_stroke), big.mark=","), round(sum(cvd_primarily_stroke)/sum(incident_cvd)*100, digits=2)),
  Stroke_fatal = sprintf("%s (%s%%)", format(sum(cvd_primarily_stroke & cvd_is_fatal), big.mark=","), round(sum(cvd_primarily_stroke & cvd_is_fatal)/sum(cvd_primarily_stroke)*100, digits=2)),
  Median_followup = sprintf("%s (%s)", round(median(incident_followup), digits=1), round(sd(incident_followup), digits=2)),
  Censored_lt_10yr = sprintf("%s (%s%%)", format(sum(incident_followup < 10), big.mark=","), round(sum(incident_followup < 10)/.N*100, digits=2)),
  Censored_fatal = sprintf("%s (%s%%)", format(sum(all_cause_mortality & incident_followup < 10), big.mark=","), round(sum(all_cause_mortality & incident_followup < 10)/.N*100, digits=2)),
  Censored_pandemic = sprintf("%s (%s%%)", format(sum(incident_followup_date == "2020-03-01" & incident_followup < 10 & !all_cause_mortality), big.mark=","), 
                                           round(sum(incident_followup_date == "2020-03-01" & incident_followup < 10 & !all_cause_mortality)/.N*100, digits=2)),
  Censored_lost = sprintf("%s (%s%%)", format(sum(lost_to_followup != "" & incident_followup < 10 & !all_cause_mortality), big.mark=","), 
                                       round(sum(lost_to_followup != "" & incident_followup < 10 & !all_cause_mortality)/.N*100, digits=2)),
  Censored_max_Wales = sprintf("%s (%s%%)", format(sum(latest_hospital_nation == "Wales" & incident_followup < 10 & !all_cause_mortality & lost_to_followup == ""), big.mark=","), 
                                            round(sum(latest_hospital_nation == "Wales" & incident_followup < 10 & !all_cause_mortality & lost_to_followup == "")/.N*100, digits=2))
), by=.(cohort=partition, sex)]

cohort_info <- cohort_info[order(sex)][order(-cohort)]

fwrite(as.data.table(t(cohort_info), keep.rownames=TRUE), sep="\t", quote=FALSE, col.names=FALSE, file="analyses/cohort_information_by_sex.txt")

# Also partition by CVD cases status
cohort_info <- dat[,.(
  samples = sprintf("%s (%s%%)", format(.N, big.mark=","), round(.N/dat[partition == cohort, .N]*100, digits=1)),
  men = sprintf("%s (%s%%)", format(sum(sex == "Male"), big.mark=","), round(sum(sex == "Male")/.N*100, digits=1)),
  age = sprintf("%s (%s)", round(median(age), digits=1), round(sd(age), digits=2)),
  BMI = sprintf("%s (%s)", round(median(na.omit(bmi)), digits=1), round(sd(na.omit(bmi)), digits=2)),
  SBP = sprintf("%s (%s)", round(median(sbp), digits=1), round(sd(sbp), digits=2)),
  diabetics = sprintf("%s (%s%%)", format(sum(diabetes), big.mark=","), round(sum(diabetes)/.N*100, digits=2)),
  smokers = sprintf("%s (%s%%)", format(sum(smoking), big.mark=","), round(sum(smoking)/.N*100, digits=2)),
  total_cholesterol = sprintf("%s (%s)", round(median(tchol), digits=2), round(sd(tchol), digits=2)),
  hdl_cholesterol = sprintf("%s (%s)", round(median(hdl), digits=2), round(sd(hdl), digits=2)),
  CVD_fatal = sprintf("%s (%s%%)", format(sum(cvd_is_fatal), big.mark=","), round(sum(cvd_is_fatal)/sum(incident_cvd)*100, digits=2)),
  primary_cause = sprintf("%s (%s%%)", format(sum(cvd_is_primary_cause), big.mark=","), round(sum(cvd_is_primary_cause)/sum(incident_cvd)*100, digits=2)),
  CHD = sprintf("%s (%s%%)", format(sum(cvd_primarily_chd), big.mark=","), round(sum(cvd_primarily_chd)/sum(incident_cvd)*100, digits=2)),
  CHD_fatal = sprintf("%s (%s%%)", format(sum(cvd_primarily_chd & cvd_is_fatal), big.mark=","), round(sum(cvd_primarily_chd & cvd_is_fatal)/sum(cvd_primarily_chd)*100, digits=2)),
  Stroke = sprintf("%s (%s%%)", format(sum(cvd_primarily_stroke), big.mark=","), round(sum(cvd_primarily_stroke)/sum(incident_cvd)*100, digits=2)),
  Stroke_fatal = sprintf("%s (%s%%)", format(sum(cvd_primarily_stroke & cvd_is_fatal), big.mark=","), round(sum(cvd_primarily_stroke & cvd_is_fatal)/sum(cvd_primarily_stroke)*100, digits=2)),
  Median_followup = sprintf("%s (%s)", round(median(incident_followup), digits=1), round(sd(incident_followup), digits=2)),
  Censored_lt_10yr = sprintf("%s (%s%%)", format(sum(incident_followup < 10), big.mark=","), round(sum(incident_followup < 10)/.N*100, digits=2)),
  Censored_fatal = sprintf("%s (%s%%)", format(sum(all_cause_mortality & incident_followup < 10), big.mark=","), round(sum(all_cause_mortality & incident_followup < 10)/.N*100, digits=2)),
  Censored_pandemic = sprintf("%s (%s%%)", format(sum(incident_followup_date == "2020-03-01" & incident_followup < 10 & !all_cause_mortality), big.mark=","), 
                                           round(sum(incident_followup_date == "2020-03-01" & incident_followup < 10 & !all_cause_mortality)/.N*100, digits=2)),
  Censored_lost = sprintf("%s (%s%%)", format(sum(lost_to_followup != "" & incident_followup < 10 & !all_cause_mortality), big.mark=","), 
                                       round(sum(lost_to_followup != "" & incident_followup < 10 & !all_cause_mortality)/.N*100, digits=2)),
  Censored_max_Wales = sprintf("%s (%s%%)", format(sum(latest_hospital_nation == "Wales" & incident_followup < 10 & !all_cause_mortality & lost_to_followup == ""), big.mark=","), 
                                            round(sum(latest_hospital_nation == "Wales" & incident_followup < 10 & !all_cause_mortality & lost_to_followup == "")/.N*100, digits=2))
), by=.(cohort=partition, case_status=ifelse(incident_cvd, "case", "control"))]

fwrite(as.data.table(t(cohort_info), keep.rownames=TRUE), sep="\t", quote=FALSE, col.names=FALSE, file="analyses/cohort_information_by_case_status.txt")


