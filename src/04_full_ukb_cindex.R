library(data.table)
library(foreach)
library(survival)
library(ggplot2)
library(ggh4x)
library(forcats)
source('src/utils/SCORE2.R')
source('src/utils/score_cindex.R')
registerDoMC(10)

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

# Now do sample exclusions to derive analysis cohort
dat <- dat[!cvd[(ehr_linkage_withdrawn)], on = .(eid)] # no EHR linkage
dat <- dat[age >= 40 & age <= 69] # eligible age for SCORE2 risk prediction
dat <- dat[!(ASCVD) | is.na(ASCVD)] # without established atherosclerotic CVD
dat <- dat[!(prevalent_diabetes_mellitus) | is.na(prevalent_diabetes_mellitus)] # without prevalent diabetes
dat <- dat[!(prevalent_CKD) | is.na(prevalent_CKD)] # without CKD
dat <- dat[DLCN_FH_score < 6] # without probable familial hypercholesteremia
dat <- dat[!(cholesterol_medication) | is.na(cholesterol_medication)] # without lipid lowering medications
dat <- dat[!is.na(sbp)] # missing SBP for SCORE2 
dat <- dat[!is.na(hdl) & !is.na(tchol)] # missing HDL and total cholesterol for SCORE2

# Add SCORE2
dat[, SCORE2 := score2(sex, age, smoking, sbp, tchol, hdl, type="linear predictor")]
dat[, SCORE2_excl_UKB := score2(sex, age, smoking, sbp, tchol, hdl, type="linear predictor", weights="excluding UK Biobank")]

# Sanity check C-indices
score2_cind <- rbind(idcol="SCORE2_method",
  "Weights derived from all datasets"=dat[, score_cindex(Surv(incident_cvd_followup, incident_cvd) ~ SCORE2, data=.SD), by=.(sex=paste0(sex, "s"))],
  "Weights derived from all datasets"=dat[, .(sex="Sex-stratified", score_cindex(Surv(incident_cvd_followup, incident_cvd) ~ SCORE2 + strata(sex), data=.SD))],
  "Weights derived excluding UK Biobank"=dat[, score_cindex(Surv(incident_cvd_followup, incident_cvd) ~ SCORE2_excl_UKB, data=.SD), by=.(sex=paste0(sex, "s"))],
  "Weights derived excluding UK Biobank"=dat[, .(sex="Sex-stratified", score_cindex(Surv(incident_cvd_followup, incident_cvd) ~ SCORE2_excl_UKB + strata(sex), data=.SD))]
)
fwrite(score2_cind, sep="\t", quote=FALSE, file="analyses/test/cindex_by_SCORE2_method_full_UKB.txt")

score2_cind[, sex := factor(sex, levels=c("Males", "Females", "Sex-stratified"))]
score2_cind[, SCORE2_method := factor(SCORE2_method, levels=c("Weights derived excluding UK Biobank", "Weights derived from all datasets"))]
g <- ggplot(score2_cind) + 
  aes(x=C.index, xmin=L95, xmax=U95, y=SCORE2_method, color=SCORE2_method) +
  facet_wrap(~ sex, scales="free_x") +
  geom_errorbarh(height=0) +
  geom_point(shape=18) +
  scale_color_manual(values=c("Weights derived excluding UK Biobank"="#2166ac", "Weights derived from all datasets"="#b2182b")) +
  xlab("C-index (95% CI)") +
  guides(color=guide_legend(title="SCORE2 model", reverse=TRUE)) +
  theme_bw() +
  theme(
    axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    strip.background=element_blank(), strip.text=element_text(size=8, face="bold"), 
    legend.text=element_text(size=6), legend.title=element_text(size=8)
  )
ggsave(g, width=7.2, height=1.5, file="analyses/test/cindex_by_SCORE2_method_full_UKB.pdf")

# Write out analysis cohort
fwrite(dat, sep="\t", quote=FALSE, file="data/cleaned/full_UKB_analysis_cohort.txt")

#####################################################################
# Fit C-index for individual biochemistry biomarkers
#####################################################################

# Code risk factors as per SCORE2
dat[, age := (age - 60)/5]
dat[is.na(smoking), smoking := FALSE]
dat[, smoking := factor(smoking, levels=c(FALSE, TRUE))]
dat[, sbp := (sbp - 120)/20]
dat[, tchol := (tchol - 6)/1]
dat[, hdl := (hdl - 1.3)/0.5]

# Load biomarker information sheets
bio_info <- fread("data/ukb/biomarkers/output/biomarker_info.txt")

# Code factors, using non-risk/lower-risk group as reference
dat[, sex := factor(sex, levels=c("Female", "Male"))]

# Get list of biomarkers to test
test_assay <- bio_info[!is.na(UKB.Field.ID) & sample_type != "Urine" & var != "tchol" & var != "hdl", var]

# Build set of models to test
models <- rbind(
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ age", type="risk_factors", name="age"),
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ smoking", type="risk_factors", name="smoking"),
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ sbp", type="risk_factors", name="sbp"),
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ tchol", type="risk_factors", name="tchol"),
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ hdl", type="risk_factors", name="hdl"),
  data.table(formula=sprintf("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2_excl_UKB) + scale(%s)", test_assay), type="assays", name=test_assay)
)

# Compute C-indices
cinds <- foreach(this_sex=c("Males", "Females", "Sex-stratified"), .combine=rbind) %do% {

  # Extract subset of data needed
  if (this_sex == "Males") {
    this_dat <- dat[sex == "Male"]
  } else if (this_sex == "Females") {
    this_dat <- dat[sex == "Female"]
  } else {
    this_dat <- dat
  }

  # Fit individual models
  foreach(midx = models[,.I], .combine=rbind) %dopar% {
    # Extract model information
    this_model <- models[midx]

    # Set up formula
    mf <- this_model$formula
    if (this_sex == "Sex-stratified") {
      mf <- paste(mf, "+ strata(sex)")
    }

    # Fit Cox proportional hazards model
    cph <- coxph(as.formula(mf), data=this_dat)

    # Get C-index and its 95% CI 
    cindex <- summary(cph)$concordance[1]
    cindex.se <- summary(cph)$concordance[2]
    cindex.l95 <- cindex - qnorm(0.975)*cindex.se
    cindex.u95 <- cindex + qnorm(0.975)*cindex.se

    # Return results
    cbind(this_model, sex=this_sex, samples=cph$n, cases=cph$nevent,
      C.index=cindex, L95=cindex.l95, U95=cindex.u95)
  }
}

# Load in and add SCORE2 C-index
score2_cindex <- fread("analyses/test/cindex_by_SCORE2_method_full_UKB.txt")
score2_cindex <- score2_cindex[SCORE2_method == "Weights derived excluding UK Biobank"]
score2_cindex[, formula := "Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2_excl_UKB)"]
score2_cindex[, type := "SCORE2"]
score2_cindex[, name := "SCORE2"]
score2_cindex[sex == "Males", c("samples", "cases") := dat[sex == "Male", .(.N, sum(incident_cvd))]]
score2_cindex[sex == "Females", c("samples", "cases") := dat[sex == "Female", .(.N, sum(incident_cvd))]]
score2_cindex[sex == "Females", c("samples", "cases") := dat[sex == "Female", .(.N, sum(incident_cvd))]]
score2_cindex[sex == "Sex-stratified", c("samples", "cases") := dat[, .(.N, sum(incident_cvd))]]
score2_cindex[, c("SCORE2_method", "SE") := NULL]
cinds <- rbind(cinds, score2_cindex)

# Add in missing strata term to documented formula
cinds[sex == "Sex-stratified" & name != "sex", formula := paste(formula, "+ strata(sex)")]

# Add in delta change compared to SCORE2
score2_ref <- cinds[name == "SCORE2"]
score2_ref[, type := "assays"]
cinds[score2_ref, on = .(type, sex), c("deltaC", "deltaC.L95", "deltaC.U95") := .(C.index - i.C.index, L95 - i.C.index, U95 - i.C.index)]

# Add human friendly display name
cinds[, display_name := name]
cinds[name == "age", display_name := "Age"]
cinds[name == "sex", display_name := "Sex"]
cinds[name == "smoking", display_name := "Smoker"]
cinds[name == "sbp", display_name := "SBP"]
cinds[name == "hdl", display_name := "HDL cholesterol"]
cinds[name == "tchol", display_name := "Total cholesterol"]
cinds[name %in% test_assay, display_name := fcase(
  name == "alb", "Albumin",
  name == "alt", "ALT",
  name == "alp", "ALP",
  name == "apoa1", "ApoA1",
  name == "apob", "ApoB",
  name == "asp", "AST",
  name == "dbili", "Bilirubin (direct)",
  name == "tbili", "Bilirubin (total)",
  name == "calcium", "Calcium",
  name == "creat", "Creatinine",
  name == "crp", "CRP",
  name == "cyst", "Cystatin-C",
  name == "ggt", "GGT",
  name == "glucose", "Glucose",
  name == "hba1c", "HbA1c",
  name == "igf1", "IGF-1",
  name == "lpa", "Lp(a)",
  name == "ldl", "LDL cholesterol",
  name == "oest", "Oestradiol",
  name == "phos", "Phosphate",
  name == "rheuf", "RF",
  name == "shbg", "SHBG",
  name == "testos", "Testosterone",
  name == "protein", "Total protein",
  name == "trig", "Triglycerides",
  name == "uric", "Urate",
  name == "urea", "Urea",
  name == "vitd25", "Vitamin D"
)]

# Write out
fwrite(cinds, sep="\t", quote=FALSE, file="analyses/univariate/cindices_full_UKB.txt")

# Plot
ggdt <- rbind(
  cinds[type == "SCORE2"],
  cinds[type == "assays"][order(-C.index)]
)
ggdt[type != "SCORE2", display_name := paste("SCORE2 +", display_name)]
ggdt[, sex := factor(sex, c("Males", "Females", "Sex-stratified"))]
ggdt[, type := factor(type, levels=c("SCORE2", "assays"))]
ggdt <- ggdt[order(sex)]
ggdt[, rank := factor(.I)]

g <- ggplot(ggdt) +
  aes(y=fct_rev(rank), x=C.index, xmin=L95, xmax=U95, color=type) +
  facet_grid2(. ~ sex, scales="free", independent="all") +
  geom_errorbarh(height=0) +
  geom_point(shape=23, fill="white", size=1.2) +
  scale_color_manual(values=c("SCORE2"="#2166ac", "assays"="black")) +
  geom_vline(data=ggdt[type == "SCORE2"], aes(xintercept=C.index), linetype=2, color="#4393c3") +
  scale_y_discrete(labels=structure(ggdt$display_name, names=as.character(ggdt$rank))) +
  xlab("C-index (95% CI)") +
  theme_bw() +
  theme(
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    axis.text.y=element_text(size=6, color="black"), axis.title.y=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    strip.background=element_blank(), strip.text=element_text(size=8, face="bold"),
    legend.position="none"
  )
ggsave(g, width=7.2, height=4, file="analyses/univariate/assays.pdf")

# Format table for output
dt <- cinds[type %in% c("SCORE2", "assays")]
dt <- dcast(dt, display_name ~ sex, value.var=c("samples", "cases", "C.index", "L95", "U95", "deltaC", "deltaC.L95", "deltaC.U95"))
dt <- dt[,.SD, .SDcols=c("display_name",
  sapply(c("_Males", "_Females", "_Sex-stratified"), function(f) {
    paste0(c("samples", "cases", "C.index", "L95", "U95", "deltaC", "deltaC.L95", "deltaC.U95"), f)
  })
)]

dt <- dt[order(-C.index_Males)]
dt <- rbind(dt[display_name == "SCORE2"], dt[display_name != "SCORE2"])
dt[display_name != "SCORE2", display_name := paste("SCORE2 +", display_name)]
fwrite(dt, sep="\t", quote=FALSE, file="analyses/univariate/assays_cindex_supp.txt")

