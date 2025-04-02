library(data.table)
library(foreach)
source('src/utils/SCORE2.R')

# Read in analysis cohort
p12 <- fread("data/cleaned/analysis_cohort.txt")
p3 <- fread("data/cleaned/phase3_analysis_cohort.txt")
dat <- rbind(p12, p3, fill=TRUE)

# Standardise PRSs 
dat[, CAD_metaGRS := scale(CAD_metaGRS)]
dat[, Stroke_metaGRS := scale(Stroke_metaGRS)]

# Compute absolute SCORE2 risk
dat[, SCORE2_uncalibrated_risk := score2_absrisk(sex, SCORE2)]
dat[, SCORE2_uk_risk := score2_recalibration(sex, SCORE2_uncalibrated_risk, "low")]
dat[, SCORE2_excl_UKB_uncalibrated_risk := score2_absrisk(sex, SCORE2_excl_UKB)]
dat[, SCORE2_excl_UKB_uk_risk := score2_recalibration(sex, SCORE2_excl_UKB_uncalibrated_risk, "low")]

# Duplicate with a "total" so we can easily tabulate males, females, and overall cohort
dat <- rbind(idcol="duplicate", dat, dat)
dat[duplicate == 2, sex := "Total"]

# Code sex as factor for output
dat[sex != "Total", sex := paste0(sex, "s")]
dat[, sex := factor(sex, levels=c("Males", "Females", "Total"))]

# Tabulate cohort information by sex
cvd_cohort_info <- dat[order(sex),.(
  samples = sprintf("%s (%s%%)", format(.N, big.mark=","), round(.N/dat[,.N]*100, digits=1)),
  assumed_ascvd_free = sprintf("%s (%s%%)", format(sum(is.na(ASCVD)), big.mark=","), round(sum(is.na(ASCVD))/.N*100, digits=2)),
  assumed_medication_free = sprintf("%s (%s%%)", format(sum(is.na(cholesterol_medication)), big.mark=","), round(sum(is.na(cholesterol_medication))/.N*100, digits=2)),
  age = sprintf("%s (%s)", round(mean(age)), round(sd(age), digits=1)),
  BMI = sprintf("%s (%s)", round(mean(na.omit(bmi)), digits=1), round(sd(na.omit(bmi)), digits=1)),
  SBP = sprintf("%s (%s)", round(mean(sbp)), round(sd(sbp))),
  smokers = sprintf("%s (%s%%)", format(sum(na.omit(smoking)), big.mark=","), round(sum(na.omit(smoking))/.N*100, digits=1)),
  assumed_non_smoking = sprintf("%s (%s%%)", format(sum(is.na(smoking)), big.mark=","), round(sum(is.na(smoking))/.N*100, digits=2)),
  history_gestational_diabetes = sprintf("%s (%s%%)", format(sum(na.omit(history_gestational_diabetes)), big.mark=","), round(sum(na.omit(history_gestational_diabetes))/dat[sex == "Female", .N]*100, digits=2)),
  hdl_cholesterol = sprintf("%s (%s)", round(mean(hdl, na.rm=TRUE), digits=1), round(sd(hdl, na.rm=TRUE), digits=1)),
  total_cholesterol = sprintf("%s (%s)", round(mean(tchol, na.rm=TRUE), digits=1), round(sd(tchol, na.rm=TRUE), digits=1)),
  ldl_cholesterol = sprintf("%s (%s)", round(mean(ldl, na.rm=TRUE), digits=1), round(sd(ldl, na.rm=TRUE), digits=1)),
  SCORE2 = sprintf("%s%% (%s%%)", round(mean(SCORE2_uk_risk*100), digits=1), round(sd(SCORE2_uk_risk * 100), digits=1)),
  SCORE2_excl_UKB = sprintf("%s%% (%s%%)", round(mean(SCORE2_excl_UKB_uk_risk * 100), digits=1), round(sd(SCORE2_excl_UKB_uk_risk * 100), digits=1)),
  CAD_metaGRS = sprintf("%s (%s)", round(mean(CAD_metaGRS, na.rm=TRUE), digits=1), round(sd(CAD_metaGRS, na.rm=TRUE), digits=1)),
  Stroke_metaGRS = sprintf("%s (%s)", round(mean(Stroke_metaGRS, na.rm=TRUE), digits=1), round(sd(Stroke_metaGRS, na.rm=TRUE), digits=1)),
  CVD = sprintf("%s (%s%%)", format(sum(incident_cvd), big.mark=","), round(sum(incident_cvd)/.N*100, digits=1)),
  onset = sprintf("%s (%s-%s)", round(median(incident_cvd_followup[incident_cvd]), digits=1),
                              round(quantile(incident_cvd_followup[incident_cvd], 0.25), digits=1), 
                              round(quantile(incident_cvd_followup[incident_cvd], 0.75), digits=1)),
  primary_cause = sprintf("%s (%s%%)", format(sum(cvd_is_primary_cause), big.mark=","), round(sum(cvd_is_primary_cause)/sum(incident_cvd)*100, digits=1)),
  fatal = sprintf("%s (%s%%)", format(sum(incident_cvd_is_fatal), big.mark=","), round(sum(incident_cvd_is_fatal)/sum(incident_cvd)*100, digits=1)),
  competing_fatal = sprintf("%s (%s%%)", format(sum(!incident_cvd & incident_cvd_followup < 10 & mortality_at_cvd_followup), big.mark=","), 
                                          round(sum(!incident_cvd & incident_cvd_followup < 10 & mortality_at_cvd_followup)/.N*100, digits=1)),
  fatal_onset = sprintf("%s (%s-%s)", round(median(incident_cvd_followup[!incident_cvd & incident_cvd_followup < 10 & mortality_at_cvd_followup]), digits=1),
                                    round(quantile(incident_cvd_followup[!incident_cvd & incident_cvd_followup < 10 & mortality_at_cvd_followup], 0.25), digits=1), 
                                    round(quantile(incident_cvd_followup[!incident_cvd & incident_cvd_followup < 10 & mortality_at_cvd_followup], 0.75), digits=1)),
  Censored_lt_10yr = sprintf("%s (%s%%)", format(sum(!incident_cvd & incident_cvd_followup < 10 & !mortality_at_cvd_followup), big.mark=","), 
                                           round(sum(!incident_cvd & incident_cvd_followup < 10 & !mortality_at_cvd_followup)/.N*100, digits=1)),
  Censored_max_Wales = sprintf("%s (%s%%)", format(sum(!incident_cvd & cvd_follow_lt10_Wales & !mortality_at_cvd_followup), big.mark=","), 
                                             round(sum(!incident_cvd & cvd_follow_lt10_Wales & !mortality_at_cvd_followup)/sum(!incident_cvd & incident_cvd_followup < 10 & !mortality_at_cvd_followup)*100, digits=1)),
  Wales_follow = sprintf("%s (%s-%s)", round(median(incident_cvd_followup[!incident_cvd & cvd_follow_lt10_Wales & !mortality_at_cvd_followup]), digits=1),
                                     round(quantile(incident_cvd_followup[!incident_cvd & cvd_follow_lt10_Wales & !mortality_at_cvd_followup], 0.25), digits=1), 
                                     round(quantile(incident_cvd_followup[!incident_cvd & cvd_follow_lt10_Wales & !mortality_at_cvd_followup], 0.75), digits=1)),
  Censored_lost = sprintf("%s (%s%%)", format(sum(!incident_cvd & lost_at_cvd_followup & !mortality_at_cvd_followup), big.mark=","), 
                                         round(sum(!incident_cvd & lost_at_cvd_followup & !mortality_at_cvd_followup)/sum(!incident_cvd & incident_cvd_followup < 10 & !mortality_at_cvd_followup)*100, digits=1)),
  lost_follow = sprintf("%s (%s-%s)", round(median(incident_cvd_followup[!incident_cvd & lost_at_cvd_followup & !mortality_at_cvd_followup]), digits=1),
                                    round(quantile(incident_cvd_followup[!incident_cvd & lost_at_cvd_followup & !mortality_at_cvd_followup], 0.25), digits=1), 
                                    round(quantile(incident_cvd_followup[!incident_cvd & lost_at_cvd_followup & !mortality_at_cvd_followup], 0.75), digits=1))
), by=.(sex)]

cvd_cohort_info <- as.data.table(t(cvd_cohort_info), keep.rownames=TRUE)
fwrite(cvd_cohort_info, col.names=FALSE, sep="\t", quote=FALSE, file="analyses/pooled_cohort_characteristics.txt")

# Tabulate missing biomarker data
bio_info <- fread("data/ukb/biomarkers/output/biomarker_info.txt")
bio_info <- bio_info[sample_type != "Urine" & !is.na(UKB.Field.ID)]
nmr_info <- fread("data/ukb/NMR_metabolomics/biomarker_information.txt")
nmr_info <- nmr_info[!is.na(UKB.Field.ID)]
biomarkers <- rbind(idcol="type", "assays"=bio_info[,.(biomarker=var)], "NMR"=nmr_info[,.(biomarker=Biomarker)])

dat <- rbind(idcol="cohort", "discovery"=p12, "replication"=p3, fill=TRUE)
stop()

miss <- foreach(this_cohort = c("discovery", "replication", "pooled"), .combine=rbind) %:% 
  foreach(bioIdx = biomarkers[,.I], .combine=rbind) %do% {
   
}






