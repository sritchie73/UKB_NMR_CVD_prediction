library(data.table)

# Following flowchart in Figure 2 of Eastwood et al. Algorithms for the Capture and 
# Adjudication of Prevalent and Incident Diabetes in UK Biobank. PLoS One (2016).

system("mkdir -p output")

########################################################
# Load extracted/curated variables to use in algorithm
########################################################

# Start with basic participant characteristics
anthro <- fread("data/curated/anthropometrics/output/anthropometrics.txt")
dat <- anthro[, .(eid, visit_index, sex)]

# Need to determine whether participant is South Asian or African Caribbean ethnicity
# as age of onset is younger in this group in the UK.
eth <- anthro[, .(eid, visit_index, ethnicity_SAorAC=fcase(
    ethnicity_subgroup == "Asian or Asian British", TRUE,
    ethnicity_subgroup == "Indian", TRUE,
    ethnicity_subgroup == "Pakistani", TRUE,
    ethnicity_subgroup == "Bangladeshi", TRUE,
    ethnicity_subgroup == "Black or Black British", TRUE,
    ethnicity_subgroup == "African", TRUE,
    ethnicity_subgroup == "Caribbean", TRUE,
    ethnicity_subgroup == "Any other Black background", TRUE,
    ethnicity_subgroup == "Prefer not to answer", NA,
    ethnicity_subgroup == "Do not know", NA,
    ethnicity_subgroup == "", NA,
    default = FALSE))]
dat <- dat[eth, on = .(eid, visit_index)]

# load self-report data
self_report <- fread("data/curated/self_report/self_report_codes.txt")

# Diabetes diagnosis (touchscreen)
diab_ts_dat <- self_report[field_id == 2443]
diab_ts_dat[, diab_ts := fcase(
  code == 1, TRUE, # Yes
  code == 0, FALSE,  # No
  code == -2, FALSE, # Not applicable (i.e. answered 'No' to any history of diabetes (field #2443) or is Male)
  code == -1, NA,   # Do not know
  code == -3, NA)]  # Prefer not to answer
dat[diab_ts_dat, on = .(eid, visit_index), diab_ts := i.diab_ts]

# Age at diabetes diagnosis (touchscreen, derived from field 2976)
dat[diab_ts_dat, on = .(eid, visit_index), age_diab_ts := i.age]

# Gestational diabetes (touchscreen)
gdm_ts_dat <- self_report[field_id == 4041]
gdm_ts_dat[, gdm_ts := fcase(
  code == 1, TRUE, # Yes
  code == 0, FALSE,  # No
  code == -2, FALSE, # Not applicable (i.e. answered 'No' to any history of diabetes (field #2443) or is Male)
  code == -1, NA,   # Do not know
  code == -3, NA)]  # Prefer not to answer
dat[gdm_ts_dat, on = .(eid, visit_index), gdm_ts := i.gdm_ts]
dat[sex == "Male", gdm_ts := FALSE]

# Load medications (touchscreen)
meds_ts <- fread("data/curated/medication/output/medications_simple.txt")

# insulin medication (touchscreen fields 6177 [males] and 6153 [females])
# codings curated and merged in respective directory
dat[meds_ts, on = .(eid, visit_index), insulin_ts := i.insulin_medication]

# insulin within 1 year of diabetes diagnosis (touchscreen)
insulin_1yr_diab_ts_dat <- self_report[field_id == 2986]
insulin_1yr_diab_ts_dat[, insulin_1yr_diab_ts := fcase(
  code == 1, TRUE, # Yes
  code == 0, FALSE,  # No
  code == -1, NA,   # Do not know
  code == -3, NA)]  # Prefer not to answer
dat[insulin_1yr_diab_ts_dat, on = .(eid, visit_index), insulin_1yr_diab_ts := i.insulin_1yr_diab_ts]
dat[!(diab_ts), insulin_1yr_diab_ts := FALSE]

# Get diabetes diagnoses and age of onset from interview with nurse
nurse_interview <- self_report[field_id == 20002]

dat[nurse_interview, on = .(eid, visit_index), c("diab_ns_ni", "age_diab_ns_ni") := .(FALSE, NA_real_)]
dat[nurse_interview[code == 1220], on = .(eid, visit_index), c("diab_ns_ni", "age_diab_ns_ni") := .(TRUE, i.age)]

dat[nurse_interview, on = .(eid, visit_index), c("gdm_ni", "age_gdm_ni") := .(FALSE, NA_real_)]
dat[nurse_interview[code == 1221], on = .(eid, visit_index), c("gdm_ni", "age_gdm_ni") := .(TRUE, i.age)]

dat[nurse_interview, on = .(eid, visit_index), c("t1d_ni", "age_t1d_ni") := .(FALSE, NA_real_)]
dat[nurse_interview[code == 1222], on = .(eid, visit_index), c("t1d_ni", "age_t1d_ni") := .(TRUE, i.age)]

dat[nurse_interview, on = .(eid, visit_index), c("t2d_ni", "age_t2d_ni") := .(FALSE, NA_real_)]
dat[nurse_interview[code == 1223], on = .(eid, visit_index), c("t2d_ni", "age_t2d_ni") := .(TRUE, i.age)]

# Fix coding error for gestational diabetes (set to nonspecific diabetes)
dat[sex == "Male" & gdm_ni, 
  c("gdm_ni", "age_gdm_ni", "diab_ns_ni", "age_diab_ns_ni") := 
  .(FALSE, NA, TRUE, age_gdm_ni)]

# People with t1d and t2d in nurse interview set to nonspecific (see Appendix A of Eastwood et al.)
dat[t2d_ni & t1d_ni,
  c("t2d_ni", "age_t2d_ni", "t1d_ni", "age_t1d_ni", "diab_ns_ni", "age_diab_ns_ni") := 
  .(FALSE, NA, FALSE, NA, TRUE, pmin(age_t1d_ni, age_t2d_ni))]

# Build collated field for any diabetes (excluding gestational)
dat[, diab_ni := diab_ns_ni | t1d_ni | t2d_ni]
dat[, age_diab_ni := pmin(age_diab_ns_ni, age_t1d_ni, age_t2d_ni, na.rm=TRUE)]

dat[, diab := diab_ni | diab_ts]
dat[, age_diab := age_diab_ts]
dat[!is.na(age_diab_ni), age_diab := age_diab_ni]

# Get medications from interview with nurse
meds_ni <- fread("data/curated/medication/output/detailed_medications_field_20003.txt")

dat[meds_ni, on = .(eid, visit_index), insulin_ni := FALSE]
dat[meds_ni[medication_code == 1140883066], on = .(eid, visit_index), insulin_ni := TRUE]

dat[meds_ni, on = .(eid, visit_index), metformin_ni := FALSE]
dat[meds_ni[medication_code %in% c(1140884600, 1140874686, 1141189090)], on = .(eid, visit_index), metformin_ni := TRUE]

dat[meds_ni, on = .(eid, visit_index), other_diab_med_ni := FALSE]
dat[meds_ni[medication_code %in% c(
    # Current non-metformin oral anti-diabetic receipt (Sulfonylureas)
    1140874718, 1140874744, 1140874746, 1141152590, 1141156984, 1140874646, 
    1141157284 , 1140874652, 1140874674, 1140874728,
    # Taking other oral anti-diabetic (acarbose, guar gum)
    1140868902, 1140868908, 1140857508,
    # Meglitinides
    1141173882, 1141173786, 1141168660,
    # Glitazones
    1141171646, 1141171652, 1141153254, 1141177600, 1141177606
  )], on = .(eid, visit_index), other_diab_med_ni := TRUE]

# Build collated fields
dat[, diab_med_ni := insulin_ni | metformin_ni | other_diab_med_ni]
dat[, diab_med_ts := insulin_ts | insulin_1yr_diab_ts]
dat[, diab_med := diab_med_ts | diab_med_ni]

########################################################
# Run algorithm 1 to identify possible diabetes
########################################################

cat("Starting algorithm 1: determining diabetes presence or absence, and sorting diabetes types...\n\n")

dat[, adjudicated_diabetes := NA_character_]
dat[, adjudication_finalized := NA_character_]
dat[, adjudication_note := NA_character_]

totals <- dat[,.N,by=visit_index]

# 1.1 Identify participants for whom diabetes is possible.
#
# Justification from Appendix A:
#
# Based on presence of self-reported diabetes diagnoses or medication. 
# The exception was diabetes (non-gestational) report by touchscreen alone, 
# which was not used to indicate possible diabetes, since only 3% reported 
# other markers of diabetes, compared with 73% of those with diabetes diagnoses 
# by touchscreen and nurse interview together, and 16% of those with diagnoses 
# in the nurse interview data alone.
#
# Presence of diabetes complications was not used to rule in/ out cases, due to 
# low prevalence (0.3%).
continue_1.1 <- dat[gdm_ts | gdm_ni | diab_ns_ni | t1d_ni | t2d_ni | diab_med]
exited_1.1 <- dat[!continue_1.1, on = .(eid, visit_index)]
dat <- continue_1.1

# Here, also additionally classify those with insufficient information to call
# diabetes unlikely\
exited_1.1[, adjudicated_diabetes := "Diabetes unlikely"]
exited_1.1[, adjudication_finalized := "Step 1.1"]
exited_1.1[is.na(gdm_ts) | is.na(gdm_ni) | is.na(diab_med) | is.na(diab), 
  adjudication_note := "1.1: Some missing data."]

cat("Exited at step 1.1 with adjudication reason \"Diabetes unlikely\":\n")
cat(" ", format(exited_1.1[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_1.1[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_1.1[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_1.1[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

# 1.2 Identify participants with possible gestational diabetes
#
# Justification from Appendix A:
#
# Based on self-report of gestational diabetes from touchscreen or nurse interview. 
# Participants classified as possible gestational diabetes at this step also must 
# not have any evidence of current diabetes (medication or self-report of other 
# types of diabetes), i.e. those who may have had pre -existing diabetes or 
# developed type 1 or 2 diabetes subsequent to gestational diabetes continue in 
# the algorithm. This included 32% women reporting gestational diabetes by touchscreen
# and 18% in nurse interview with evidence of current diabetes who were not currently 
# pregnant. In addition, women who reported age at diagnosis of >= 50 years also 
# continued in the algorithm.
exited_1.2 <- dat[gdm_ts | (gdm_ni & age_gdm_ni < 50)]
exited_1.2 <- exited_1.2[!t1d_ni & !t2d_ni] # drop concurrent t1d or t2d (no NAs)
exited_1.2 <- exited_1.2[!diab_med | is.na(diab_med)] # drop people with any medication in nurse interview or touchscreen (keeping also NAs)
dat <- dat[!exited_1.2, on = .(eid, visit_index)]

# But break down possible gestational diabetes status by medication missingingess
exited_1.2[, adjudicated_diabetes := "Possible gestational diabetes"]
exited_1.2[, adjudication_finalized := "Step 1.2"]
exited_1.2[is.na(diab_med_ts) & is.na(diab_med_ni), adjudication_note := "1.2: All medication data missing."]
exited_1.2[is.na(diab_med_ts) & !is.na(diab_med_ni), adjudication_note := "1.2: Touchscreen medication data missing."]
exited_1.2[!is.na(diab_med_ts) & is.na(diab_med_ni), adjudication_note := "1.2: Nurse interview medication data missing."]

cat("Exited at step 1.2 with adjudication reason \"Possible gestational diabetes\":\n")
cat(" ", format(exited_1.2[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_1.2[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_1.2[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_1.2[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

# 1.3 Possible T2D 
# 
# Justification from Appendix A:
#
# Use of this non-metformin and non-insulin diabaetes medication is likely 
# to be exclusive to type 2 diabetes.
dat[(other_diab_med_ni), adjudicated_diabetes := "Possible type 2 diabetes"]

exited_1.3 <- dat[!is.na(adjudicated_diabetes)]
dat <- dat[!exited_1.3, on = .(eid, visit_index)]
cat("Exited at step 1.3 with adjudication reason \"Possible type 2 diabetes\":\n")
cat(" ", format(exited_1.3[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_1.3[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_1.3[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_1.3[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

# 1.4 Determine T1D vs. T2D likelihood based on age of diagnosis
#
# Justification from Appendix A:
#
# Age of onset of diabetes is younger in the UK’s largest ethnic minority groups; 
# South Asians and African Caribbeans. A priori-determined ethnic-specific cut-points 
# are used to define older vs. younger age at diagnosis – indicating higher likelihood 
# of type 2 vs. type 1 diagnoses respectively. Participants with no age of diagnosis 
# information are classified into the T2D arm.
#
# Other notes:
# 
# Since age of diagnosis is numeric (e.g. after interpolation), we floor the age so that
# e.g. a 36.5 is treated as a 36 for the purposes of the algorithm (e.g. where the cutoff 
# for T1D is age of diagnosis <= 36)
#
# People with missing ethnicity data are treated as the majority group for the purposes of
# the cutoff
t1d_1.4 <- dat[(ethnicity_SAorAC & floor(age_diab) <= 30) | floor(age_diab) <= 36]
t2d_1.4 <- dat[!t1d_1.4, on = .(eid, visit_index)]

t1d_1.4[, adjudicated_diabetes := "Possible type 1 diabetes"]
t1d_1.4[is.na(ethnicity_SAorAC), adjudication_note := "1.4: Missing ethnicity for stratification of T1D/T2D by age of diagnosis."]

t2d_1.4[, adjudicated_diabetes := "Possible type 2 diabetes"]
t2d_1.4[is.na(age_diab) & is.na(adjudication_note), adjudication_note := "1.4: Missing age of diagnosis for stratification of T1D/T2D"]
t2d_1.4[!is.na(age_diab) & is.na(ethnicity_SAorAC), adjudication_note := "1.4: Missing ethnicity for stratification of T1D/T2D by age of diagnosis."]

exited_1.4 <- t2d_1.4
dat <- t1d_1.4

cat("Exited at step 1.4 with adjudication reason \"Possible type 2 diabetes\":\n")
cat(" ", format(exited_1.4[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_1.4[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_1.4[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_1.4[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

# 1.5 Further back up possibility of T1D diagnosis with medication data
#
# Justfication from Appendix A:
#
# Almost all participants with type 1 diabetes will be on insulin; coupled 
# with a lack of non-metformin oral anti-diabetic drugs and a younger age 
# at onset (see rules 1.3 and 1.4), specificity for type 1 diabetes is increased.
#
# Immediate commencement of insulin after diagnosis is obligatory for type 
# 1 diabetes. We included the criterion of commencement of insulin at less 
# than 12 months post-diagnosis in addition to current insulin use as 96 
# people entering rule 1.5 reported insulin use at less than 12 months but 
# no current insulin. These people are then carried forward to the next 
# flowchart for finalising type 1 diabetes status, and are labelled as 
# "possible type 1 diabetes". Otherwise they would be excluded, but it is 
# possible to have been on insulin initially and not currently, e.g. due to 
# pancreatic transplant/ incorrect diagnosis.
#
# Characteristics of participants self-reporting (nurse interview) type 1 
# diabetes at this stage in the algorithm appear to be congruent with the 
# diagnosis - their average age at diagnosis was 21 +/- 9 years, 98% (touchscreen) 
# or 95% (nurse interview) were on current insulin, only 10% were on metformin 
# and 92% started insulin within 1 year of diagnosis. However, due to the 
# low number of participants specifying this diagnosis at nurse interview 
# overall (0.09%), we sought other evidence of type 1 diabetes in addition. 
# NB: there was no option to specify diabetes type on the touchscreen question.
exited_1.5 <- dat[insulin_ts | insulin_ni | insulin_1yr_diab_ts | t1d_ni]
exited_1.5_t2d <- dat[!exited_1.5, on = .(eid, visit_index)]

exited_1.5_t2d[, adjudicated_diabetes := "Possible type 2 diabetes"]
exited_1.5_t2d[is.na(insulin_ts) | is.na(insulin_1yr_diab_ts), adjudication_note := paste(adjudication_note, 
  "1.5: Missing touchscreen medication data to support possibility of type 1 diabetes.")]
exited_1.5_t2d[, adjudication_note := gsub("^NA ", "", adjudication_note)]

cat("Exited at step 1.5 with adjudication reason \"Possible type 1 diabetes\":\n")
cat(" ", format(exited_1.5[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_1.5[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_1.5[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_1.5[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

cat("Exited at step 1.5 with adjudication reason \"Possible type 2 diabetes\":\n")
cat(" ", format(exited_1.5_t2d[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_1.5_t2d[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_1.5_t2d[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_1.5_t2d[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

# Collate diabetes status
diabetes_unlikely <- exited_1.1
possible_gdm <- exited_1.2
possible_t2d <- rbind(exited_1.3, exited_1.4, exited_1.5_t2d)
possible_t1d <- exited_1.5

# Print status at end of algorithm 1
cat("Algorithm 1 finished with following adjudication allocations:\n\n")

cat("  \"Diabetes unlikely\":\n")
cat("   ", format(diabetes_unlikely[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat("   ", format(diabetes_unlikely[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat("   ", format(diabetes_unlikely[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat("   ", format(diabetes_unlikely[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

cat("  \"Possible gestational diabetes\":\n")
cat("   ", format(possible_gdm[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat("   ", format(possible_gdm[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat("   ", format(possible_gdm[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat("   ", format(possible_gdm[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

cat("  \"Possible type 2 diabetes\":\n")
cat("   ", format(possible_t2d[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat("   ", format(possible_t2d[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat("   ", format(possible_t2d[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat("   ", format(possible_t2d[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

cat("  \"Possible type 1 diabetes\":\n")
cat("   ", format(possible_t1d[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat("   ", format(possible_t1d[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat("   ", format(possible_t1d[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat("   ", format(possible_t1d[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

########################################################
# Run algorithm 2 to finalize T1D diagnosis
########################################################

cat("Starting algorithm 2: determining certainty of type 1 diabetes diagnosis...\n\n")

possible_t1d[, adjudicated_diabetes := NA]

totals <- possible_t1d[,.N,by=visit_index]

# 2.1 Probable type 1 diabetes based on nurse interview
#
# Justification from Appendix A:
#
# After identifying possible type 1 diabetes cases based on age of diagnosis
# and supporting information from insulin medication usage/prescription, the
# self-report of type 1 diabetes in the nurse interview is deemed to be 
# likely correct upgrading the certainty from "possible" to "probable".
#
# NB: 10 people co-reporting type 1 and type 2 diabetes in nurse interview
# were reclassified as non-specific diabetes (see data loading steps above)
possible_t1d[(t1d_ni), adjudicated_diabetes := "Probable type 1 diabetes"]
exited_2.1 <- possible_t1d[!is.na(adjudicated_diabetes)]
exited_2.1[, adjudication_finalized := "Step 2.1"]
possible_t1d <- possible_t1d[!exited_2.1, on = .(eid, visit_index)]

cat("Exited at step 2.1 with adjudication reason \"Probable type 1 diabetes\":\n")
cat(" ", format(exited_2.1[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_2.1[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_2.1[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_2.1[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

# 2.2 Probable T1D from medication usage
#
# Justification from Appendix A:
#
# Amongst those entering this algorithm without nurse report of type 1 diabetes, 
# 78% received insulin within a year of diagnosis (touchscreen), 97% self-reported 
# current insulin use on touchscreen, and 95% reported current insulin use to the 
# nurse. 
# 
# Metformin use was reported in 52% of those with “possible” versus 11% of those 
# with “probable” type 1 diabetes. Similarly, proportions co-reporting type 2 
# diabetes were 9% in people with “possible” versus 0.5% in people with “probable” 
# type 1 diabetes. 
#
# Other notes:
#
# While the flowchart implies the conditional statement should be
#  (insulin_1yr_diab_ts & insulin_ts) | insulin_ni
# this leads to almost all remaining participants (97%) being classed as 
# "probable" whereas switching the AND and OR better matches the population
# split shown on the flowchart (74% probable) 
possible_t1d[insulin_1yr_diab_ts & (insulin_ts | insulin_ni), adjudicated_diabetes := "Probable type 1 diabetes"]
exited_2.2 <- possible_t1d[!is.na(adjudicated_diabetes)]
exited_2.2[, adjudication_finalized := "Step 2.2"]
possible_t1d <- possible_t1d[!exited_2.2, on = .(eid, visit_index)]
cat("Exited at step 2.2 with adjudication reason \"Probable type 1 diabetes\":\n")
cat(" ", format(exited_2.2[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_2.2[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_2.2[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_2.2[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

# 2.2 Possible T1D (no medication)
possible_t1d[, adjudicated_diabetes := "Possible type 1 diabetes"]
possible_t1d[, adjudication_finalized := "Step 2.2"]
possible_t1d[(insulin_ni | insulin_ts) & is.na(insulin_1yr_diab_ts), adjudication_note := paste(adjudication_note,
  "2.2: self-reported insulin, but missing data on whether prescribed within 1 year of diagnosis.")]
possible_t1d[, adjudication_note := gsub("^NA ", "", adjudication_note)]

cat("Exited at step 2.2 with adjudication reason \"Possible type 1 diabetes\":\n")
cat(" ", format(possible_t1d[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(possible_t1d[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(possible_t1d[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(possible_t1d[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

# Collate diabetes status
probable_t1d <- rbind(exited_2.1, exited_2.2)

# Print status at end of algorithm 1
cat("Algorithm 2 finished with following adjudication allocations:\n\n")

cat("  \"Probable type 1 diabetes\":\n")
cat("   ", format(probable_t1d[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat("   ", format(probable_t1d[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat("   ", format(probable_t1d[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat("   ", format(probable_t1d[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

cat("  \"Possible type 1 diabetes\":\n")
cat("   ", format(possible_t1d[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat("   ", format(possible_t1d[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat("   ", format(possible_t1d[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat("   ", format(possible_t1d[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

########################################################
# Run algorithm 3 to finalize T2D diagnosis
########################################################

cat("Starting algorithm 3: determining certainty of type 2 diabetes diagnosis...\n\n")

possible_t2d[, adjudicated_diabetes := NA]

totals <- possible_t2d[,.N,by=visit_index]

# 3.1 triage metformin users
#
# Justification from Appendix A:
#
# Metformin can be used to treat conditions other than diabetes, 
# e.g. polycystic ovarian syndrome, therefore treatment with 
# metformin only is less specific for diabetes than if used in 
# conjunction with other anti-diabetic drugs.
exited_3.1 <- possible_t2d[metformin_ni & !insulin_ni & !other_diab_med_ni]
possible_t2d <- possible_t2d[!exited_3.1, on = .(eid, visit_index)]

cat("Exited at step 3.1 to triage metformin usage cases:\n")
cat(" ", format(exited_3.1[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_3.1[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_3.1[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_3.1[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

# 3.2 identify people on metformin but no other indication of diabetes
#
# Justification from Appendix A:
#
# With the exception of the small number of participants (7%) taking metformin
# who reported diabetes on touchscreen, no other diabetes self-report data was 
# present (e.g. in nurse interview), so they were deemed unlikely to have 
# diabetes at this stage.
exited_3.2 <- exited_3.1[!(diab_ni)]
exited_3.2[, adjudicated_diabetes := "Diabetes unlikely"]
exited_3.2[, adjudication_finalized := "Step 3.2"]
exited_3.1 <- exited_3.1[!exited_3.2, on = .(eid, visit_index)]

cat("Exited at step 3.2 with adjudication reason \"Diabetes unlikely\":\n")
cat(" ", format(exited_3.2[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_3.2[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_3.2[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_3.2[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

cat("Individuals exited at step 3.1 returned to algorithm at step 3.2:\n")
cat(" ", format(exited_3.1[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_3.1[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_3.1[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_3.1[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

possible_t2d <- rbind(possible_t2d, exited_3.1)

# 3.3 probable t2d based on other non-metformin diabetes medication
#
# Justification from Appendix A:
#
# Participants assigned "probable type 2 diabetes" status at this stage were 
# more likely than those who continued in the algorithm to have developed 
# diabetes at an older age (55 +/- 9 vs. 52 +/- 9 years, p<0.001) and to be 
# taking metformin concurrently (79% vs. 51%, p<0.001), and less likely to 
# be on insulin (8% vs. 51%, p<0.001).
possible_t2d[(other_diab_med_ni), adjudicated_diabetes := "Probable type 2 diabetes"]
exited_3.3 <- possible_t2d[!is.na(adjudicated_diabetes)]
exited_3.3[, adjudication_finalized := "Step 3.3"]
possible_t2d <- possible_t2d[!exited_3.3, on = .(eid, visit_index)]

cat("Exited at step 3.3 with adjudication reason \"Probable type 2 diabetes\":\n")
cat(" ", format(exited_3.3[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_3.3[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_3.3[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_3.3[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

# 3.4 probable t2d based on absence of insulin
#
# Justification from Appendix A:
#
# Participants who were assigned "probable type 2 diabetes" status (i.e. were not 
# receiving insulin) at this stage were more likely than those who continued in 
# the algorithm to have developed diabetes at an older age (56 +/- 9 vs. 
# 46 +/- 9 years, p<0.001), though rates of metformin use were similar (51% vs. 52%,
# p=0.2).
possible_t2d[(!insulin_ts | is.na(insulin_ts)) & (!insulin_ni | is.na(insulin_ni)), 
  adjudicated_diabetes := "Probable type 2 diabetes"]
exited_3.4 <- possible_t2d[!is.na(adjudicated_diabetes)]
exited_3.4[, adjudication_finalized := "Step 3.4"]
possible_t2d <- possible_t2d[!exited_3.4, on = .(eid, visit_index)]

exited_3.4[is.na(insulin_ts) & is.na(insulin_ni), adjudication_note := paste(adjudication_note,
  "3.4: missing medication data to assess possibility of type 1 diabetes.")]
exited_3.4[is.na(insulin_ts) & !is.na(insulin_ni), adjudication_note := paste(adjudication_note,
  "3.4: missing touchscreen medication data to assess possibility of type 1 diabetes.")]
exited_3.4[!is.na(insulin_ts) & is.na(insulin_ni), adjudication_note := paste(adjudication_note,
  "3.4: missing nurse interview medication data to assess possibility of type 1 diabetes.")]
exited_3.4[, adjudication_note := gsub("^NA ", "", adjudication_note)]

cat("Exited at step 3.4 with adjudication reason \"Probable type 2 diabetes\":\n")
cat(" ", format(exited_3.4[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_3.4[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_3.4[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_3.4[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

# 3.5 possible t2d based on absence of report of t1d
#
# Justification from Appendix A:
#
# Participants who were assigned "probable type 1 diabetes" status at this stage 
# were more likely than those who were assigned "possible type 2 diabetes" status 
# to have developed diabetes at a younger age (48 +/- 8 vs. 50 +/- 7 years, p<0.001) 
# and to have commenced insulin at less than one year after diagnosis (72% vs. 39%, 
# p<0.001), and less likely to be taking metformin currently (27% vs. 53%, p<0.001).
#
# Participants were assigned "possible type 2 diabetes" rather than "probable type 
# 2 diabetes" status at this stage due to insulin use and relatively high rates of 
# commencement of insulin within a year of diagnosis (39%) casting doubt over a
# type 2 diagnosis.
possible_t2d[!(t1d_ni), adjudicated_diabetes := "Possible type 2 diabetes"]
exited_3.5 <- possible_t2d[!is.na(adjudicated_diabetes)]
exited_3.5[, adjudication_finalized := "Step 3.5"]
possible_t2d <- possible_t2d[!exited_3.5, on = .(eid, visit_index)]

cat("Exited at step 3.5 with adjudication reason \"Possible type 2 diabetes\":\n")
cat(" ", format(exited_3.5[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_3.5[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_3.5[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_3.5[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

# 3.5 probable t1d cases
exited_3.5_t1d <- possible_t2d
exited_3.5_t1d[, adjudicated_diabetes := "Probable type 1 diabetes"]
exited_3.5_t1d[, adjudication_finalized := "Step 3.5"]

cat("Exited at step 3.5 with adjudication reason \"Probable type 1 diabetes\":\n")
cat(" ", format(exited_3.5_t1d[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_3.5_t1d[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_3.5_t1d[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_3.5_t1d[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

# Collate diabetes status
probable_t1d_from_t2d <- exited_3.5_t1d
probable_t2d <- rbind(exited_3.3, exited_3.4)
possible_t2d <- exited_3.5
diabetes_unlikely_from_t2d <- exited_3.2

# Print status at end of algorithm 1
cat("Algorithm 3 finished with following adjudication allocations:\n\n")

cat("  \"Diabetes unlikely\":\n")
cat("   ", format(diabetes_unlikely_from_t2d[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat("   ", format(diabetes_unlikely_from_t2d[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat("   ", format(diabetes_unlikely_from_t2d[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat("   ", format(diabetes_unlikely_from_t2d[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

cat("  \"Probable type 2 diabetes\":\n")
cat("   ", format(probable_t2d[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat("   ", format(probable_t2d[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat("   ", format(probable_t2d[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat("   ", format(probable_t2d[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

cat("  \"Possible type 2 diabetes\":\n")
cat("   ", format(possible_t2d[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat("   ", format(possible_t2d[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat("   ", format(possible_t2d[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat("   ", format(possible_t2d[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

cat("  \"Probable type 1 diabetes\":\n")
cat("   ", format(probable_t1d_from_t2d[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat("   ", format(probable_t1d_from_t2d[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat("   ", format(probable_t1d_from_t2d[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat("   ", format(probable_t1d_from_t2d[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

########################################################
# Collate all results
########################################################

diabetes_unlikely <- rbind(diabetes_unlikely, diabetes_unlikely_from_t2d)
probable_t1d <- rbind(probable_t1d, probable_t1d_from_t2d)

dat <- rbind(diabetes_unlikely, possible_gdm, probable_t2d, possible_t2d, probable_t1d, possible_t1d)
totals <- dat[,.N,by=visit_index]

# Print status after all algorithms finished
cat("All algorithms finished with final adjudication allocations:\n\n")

cat("  \"Diabetes unlikely\":\n")
cat("   ", format(diabetes_unlikely[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat("   ", format(diabetes_unlikely[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat("   ", format(diabetes_unlikely[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat("   ", format(diabetes_unlikely[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

cat("  \"Possible gestational diabetes\":\n")
cat("   ", format(possible_gdm[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat("   ", format(possible_gdm[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat("   ", format(possible_gdm[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat("   ", format(possible_gdm[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

cat("  \"Probable type 2 diabetes\":\n")
cat("   ", format(probable_t2d[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat("   ", format(probable_t2d[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat("   ", format(probable_t2d[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat("   ", format(probable_t2d[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

cat("  \"Possible type 2 diabetes\":\n")
cat("   ", format(possible_t2d[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat("   ", format(possible_t2d[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat("   ", format(possible_t2d[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat("   ", format(possible_t2d[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

cat("  \"Probable type 1 diabetes\":\n")
cat("   ", format(probable_t1d[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat("   ", format(probable_t1d[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat("   ", format(probable_t1d[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat("   ", format(probable_t1d[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

cat("  \"Possible type 1 diabetes\":\n")
cat("   ", format(possible_t1d[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat("   ", format(possible_t1d[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat("   ", format(possible_t1d[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat("   ", format(possible_t1d[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

# extract columns we care about and reorder to match original data
final <- anthro[,.(eid, visit_index)]
dat <- dat[, .(eid, visit_index, adjudicated_diabetes, adjudication_finalized, adjudication_note)]
final <- dat[final, on = .(eid, visit_index)]

fwrite(final, sep="\t", quote=FALSE, file="output/prevalent_diabetes.txt")

