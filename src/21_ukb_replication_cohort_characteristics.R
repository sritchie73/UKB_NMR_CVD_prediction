library(data.table)
source('src/utils/SCORE2.R')

# Read in analysis cohort
dat <- fread("data/cleaned/phase3_analysis_cohort.txt")

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
fwrite(cvd_cohort_info, col.names=FALSE, sep="\t", quote=FALSE, file="analyses/phase3_cohort_characteristics.txt")

# Get characteristics of CAD and Stroke followup
cad_cohort_info <- dat[order(sex),.(
  events = sprintf("%s (%s%%)", format(sum(incident_cad), big.mark=","), round(sum(incident_cad)/.N*100, digits=1)),
  onset = sprintf("%s (%s-%s)", round(median(incident_cad_followup[incident_cad]), digits=1),
                              round(quantile(incident_cad_followup[incident_cad], 0.25), digits=1), 
                              round(quantile(incident_cad_followup[incident_cad], 0.75), digits=1)),
  primary_cause = sprintf("%s (%s%%)", format(sum(cad_is_primary_cause), big.mark=","), round(sum(cad_is_primary_cause)/sum(incident_cad)*100, digits=1)),
  fatal = sprintf("%s (%s%%)", format(sum(incident_cad_is_fatal), big.mark=","), round(sum(incident_cad_is_fatal)/sum(incident_cad)*100, digits=1)),
  competing_fatal = sprintf("%s (%s%%)", format(sum(!incident_cad & incident_cad_followup < 10 & mortality_at_cad_followup), big.mark=","), 
                                          round(sum(!incident_cad & incident_cad_followup < 10 & mortality_at_cad_followup)/.N*100, digits=1)),
  fatal_onset = sprintf("%s (%s-%s)", round(median(incident_cad_followup[!incident_cad & incident_cad_followup < 10 & mortality_at_cad_followup]), digits=1),
                                    round(quantile(incident_cad_followup[!incident_cad & incident_cad_followup < 10 & mortality_at_cad_followup], 0.25), digits=1), 
                                    round(quantile(incident_cad_followup[!incident_cad & incident_cad_followup < 10 & mortality_at_cad_followup], 0.75), digits=1)),
  Censored_lt_10yr = sprintf("%s (%s%%)", format(sum(!incident_cad & incident_cad_followup < 10 & !mortality_at_cad_followup), big.mark=","), 
                                           round(sum(!incident_cad & incident_cad_followup < 10 & !mortality_at_cad_followup)/.N*100, digits=1)),
  Censored_max_Wales = sprintf("%s (%s%%)", format(sum(!incident_cad & cad_follow_lt10_Wales & !mortality_at_cad_followup), big.mark=","), 
                                             round(sum(!incident_cad & cad_follow_lt10_Wales & !mortality_at_cad_followup)/sum(!incident_cad & incident_cad_followup < 10 & !mortality_at_cad_followup)*100, digits=1)),
  Wales_follow = sprintf("%s (%s-%s)", round(median(incident_cad_followup[!incident_cad & cad_follow_lt10_Wales & !mortality_at_cad_followup]), digits=1),
                                     round(quantile(incident_cad_followup[!incident_cad & cad_follow_lt10_Wales & !mortality_at_cad_followup], 0.25), digits=1), 
                                     round(quantile(incident_cad_followup[!incident_cad & cad_follow_lt10_Wales & !mortality_at_cad_followup], 0.75), digits=1)),
  Censored_lost = sprintf("%s (%s%%)", format(sum(!incident_cad & lost_at_cad_followup & !mortality_at_cad_followup), big.mark=","), 
                                         round(sum(!incident_cad & lost_at_cad_followup & !mortality_at_cad_followup)/sum(!incident_cad & incident_cad_followup < 10 & !mortality_at_cvd_followup)*100, digits=1)),
  lost_follow = sprintf("%s (%s-%s)", round(median(incident_cad_followup[!incident_cad & lost_at_cad_followup & !mortality_at_cad_followup]), digits=1),
                                    round(quantile(incident_cad_followup[!incident_cad & lost_at_cad_followup & !mortality_at_cad_followup], 0.25), digits=1), 
                                    round(quantile(incident_cad_followup[!incident_cad & lost_at_cad_followup & !mortality_at_cad_followup], 0.75), digits=1))
), by=.(sex)]

cad_cohort_info <- as.data.table(t(cad_cohort_info), keep.rownames=TRUE)
fwrite(cad_cohort_info, col.names=FALSE, sep="\t", quote=FALSE, file="analyses/phase3_CAD_followup.txt")

stroke_cohort_info <- dat[order(sex),.(
  events = sprintf("%s (%s%%)", format(sum(incident_stroke), big.mark=","), round(sum(incident_stroke)/.N*100, digits=1)),
  onset = sprintf("%s (%s-%s)", round(median(incident_stroke_followup[incident_stroke]), digits=1),
                              round(quantile(incident_stroke_followup[incident_stroke], 0.25), digits=1), 
                              round(quantile(incident_stroke_followup[incident_stroke], 0.75), digits=1)),
  primary_cause = sprintf("%s (%s%%)", format(sum(stroke_is_primary_cause), big.mark=","), round(sum(stroke_is_primary_cause)/sum(incident_stroke)*100, digits=1)),
  fatal = sprintf("%s (%s%%)", format(sum(incident_stroke_is_fatal), big.mark=","), round(sum(incident_stroke_is_fatal)/sum(incident_stroke)*100, digits=1)),
  competing_fatal = sprintf("%s (%s%%)", format(sum(!incident_stroke & incident_stroke_followup < 10 & mortality_at_stroke_followup), big.mark=","), 
                                          round(sum(!incident_stroke & incident_stroke_followup < 10 & mortality_at_stroke_followup)/.N*100, digits=1)),
  fatal_onset = sprintf("%s (%s-%s)", round(median(incident_stroke_followup[!incident_stroke & incident_stroke_followup < 10 & mortality_at_stroke_followup]), digits=1),
                                    round(quantile(incident_stroke_followup[!incident_stroke & incident_stroke_followup < 10 & mortality_at_stroke_followup], 0.25), digits=1), 
                                    round(quantile(incident_stroke_followup[!incident_stroke & incident_stroke_followup < 10 & mortality_at_stroke_followup], 0.75), digits=1)),
  Censored_lt_10yr = sprintf("%s (%s%%)", format(sum(!incident_stroke & incident_stroke_followup < 10 & !mortality_at_stroke_followup), big.mark=","), 
                                           round(sum(!incident_stroke & incident_stroke_followup < 10 & !mortality_at_stroke_followup)/.N*100, digits=1)),
  Censored_max_Wales = sprintf("%s (%s%%)", format(sum(!incident_stroke & stroke_follow_lt10_Wales & !mortality_at_stroke_followup), big.mark=","), 
                                             round(sum(!incident_stroke & stroke_follow_lt10_Wales & !mortality_at_stroke_followup)/sum(!incident_stroke & incident_stroke_followup < 10 & !mortality_at_stroke_followup)*100, digits=1)),
  Wales_follow = sprintf("%s (%s-%s)", round(median(incident_stroke_followup[!incident_stroke & stroke_follow_lt10_Wales & !mortality_at_stroke_followup]), digits=1),
                                     round(quantile(incident_stroke_followup[!incident_stroke & stroke_follow_lt10_Wales & !mortality_at_stroke_followup], 0.25), digits=1), 
                                     round(quantile(incident_stroke_followup[!incident_stroke & stroke_follow_lt10_Wales & !mortality_at_stroke_followup], 0.75), digits=1)),
  Censored_lost = sprintf("%s (%s%%)", format(sum(!incident_stroke & lost_at_stroke_followup & !mortality_at_stroke_followup), big.mark=","), 
                                         round(sum(!incident_stroke & lost_at_stroke_followup & !mortality_at_stroke_followup)/sum(!incident_stroke & incident_stroke_followup < 10 & !mortality_at_cvd_followup)*100, digits=1)),
  lost_follow = sprintf("%s (%s-%s)", round(median(incident_stroke_followup[!incident_stroke & lost_at_stroke_followup & !mortality_at_stroke_followup]), digits=1),
                                    round(quantile(incident_stroke_followup[!incident_stroke & lost_at_stroke_followup & !mortality_at_stroke_followup], 0.25), digits=1), 
                                    round(quantile(incident_stroke_followup[!incident_stroke & lost_at_stroke_followup & !mortality_at_stroke_followup], 0.75), digits=1))
), by=.(sex)]

stroke_cohort_info <- as.data.table(t(stroke_cohort_info), keep.rownames=TRUE)
fwrite(stroke_cohort_info, col.names=FALSE, sep="\t", quote=FALSE, file="analyses/phase3_stroke_followup.txt")


