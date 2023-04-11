library(data.table)
library(lubridate)

# Following flowchart in Figure 3b of Eastwood et al. Algorithms for the Capture and
# Adjudication of Prevalent and Incident Diabetes in UK Biobank. PLoS One (2016).

#####################################
# Load data 
#####################################

cat("Loading electronic health records...\n\n")

# Load curated hospital inpatient records, death records, and follow-up information
hes <- fread("data/curated/hospital_records/hospital_records.txt")
deaths <- fread("data/curated/deaths/deaths.txt")
followup <- fread("data/curated/followup/followup.txt")

# Load assessment visit information
assessments <- fread("data/curated/anthropometrics/output/anthropometrics.txt")
assessments <- assessments[,.(eid, visit_index, assessment_date)]

###################################################
# Create sequence of events for later determining 
# date of diagnosis or maximum diabetes free date
###################################################

cat("Creating sequence of events...\n\n")

# Construct sequence of events for each participant, including hospitalisations
# (but not operations), UK Biobank assessments, deaths, and lost to follow-up 
# dates
events <- unique(assessments[,.(eid, event="UK Biobank assessment", 
  start_date=assessment_date, end_date=assessment_date)]) 

events <- rbind(events, unique(hes[event_type == "hospitalisation",
  .(eid, event="Hospital episode", start_date=event_date, end_date=episode_end)]))

# For deaths we use the follow-up table rather than deaths as the 
# deaths table misses mortality events with no recorded causes 
# available
events <- rbind(events, unique(followup[(all_cause_mortality) & !(lost_to_followup),
  .(eid, event="Death", start_date=latest_mortality_date, 
    end_date=latest_mortality_date)]))

events <- rbind(events, unique(followup[(lost_to_followup), 
  .(eid, event="Lost to followup", start_date=lost_to_followup_date,
    end_date=lost_to_followup_date)]))

events <- rbind(events, unique(followup[, .(eid, 
    event="Min hospital followup", start_date=earliest_hospital_date,
    end_date=earliest_hospital_date)]))

events <- rbind(events, unique(followup[!lost_to_followup & !all_cause_mortality,
  .(eid, event="Max hospital followup", start_date=latest_hospital_date,
    end_date=latest_hospital_date)]))

events <- rbind(events, unique(followup[!lost_to_followup & !all_cause_mortality,
  .(eid, event="Max deaths followup", start_date=latest_mortality_date,
    end_date=latest_mortality_date)]))

# For hospital episodes that co-occur with a mortality event, the episode end date
# is sometimes listed as the day after the hospitalisation + death. For ordering 
# purposes, here we set the episode end to the same date and redo the ordering.
events[, event_grp := .GRP, by=.(event)]
events <- events[order(end_date)][order(start_date)][order(eid)]
events[, event_idx := 1:.N, by=.(eid)]

death_idx <- events[event == "Death", .(eid, event_idx, start_date)]
hes_fix <- events[death_idx, on = .(eid, event_idx > event_idx, start_date), nomatch = 0]
events[hes_fix, on = .(eid, event, start_date), end_date := start_date]

# Also fix hospital records so we can match events later
hes_fix[, code_type := "ICD-10"]
hes_fix[, event_type := "hospitalisation"]
hes[hes_fix, on = .(eid, event_type, code_type, event_date=start_date), episode_end := event_date]

events[, event_idx := NULL]
events <- unique(events) # need to deduplicate some events now
hes <- unique(hes) 

events <- events[order(event_grp)][order(end_date)][order(start_date)][order(eid)]
events[, event_idx := 1:.N, by=.(eid)]

# Sanitise - remove any events after a death
# (currently 0 events)
death_idx <- events[event == "Death", .(eid, event_idx)]
after_death <- events[death_idx, on = .(eid, event_idx > event_idx), nomatch=0, .(eid, event_idx = x.event_idx)]
events <- events[!after_death, on = .(eid, event_idx)]

# Sanitise - remove any events after lost to followup date. First, fix any instances
# where a hospital episode occurs on the lost to followup date but has an episode end
# date after that
lost_idx <- events[event == "Lost to followup", .(eid, event_idx, start_date)]
hes_fix <- events[lost_idx, on = .(eid, event_idx > event_idx, start_date), nomatch = 0]
events[hes_fix, on = .(eid, event, start_date), end_date := start_date]

events <- events[order(event_grp)][order(end_date)][order(start_date)][order(eid)]
events[, event_idx := 1:.N, by=.(eid)]

# Now we can remove any events after lost to followup dates (usually someone returning to
# the UK for biobank imaging assessment)
lost_idx <- events[event == "Lost to followup", .(eid, event_idx)]
after_lost <- events[lost_idx, on = .(eid, event_idx > event_idx), nomatch=0, .(eid, event_idx = x.event_idx)] 
events <- events[!after_lost, on = .(eid, event_idx)]

###############################################
# Determine diabetes type at each event
###############################################

cat("Determining diabetes status at each event...\n\n")

# Function for matching ICD or OPCS codes
# Returns TRUE for each element of x if it starts
# with any pattern in code_list. E.g. returns TRUE for each
# hospital episode (x) if the icd10 code starts with any of the
# codes in the vec list.
#
# E.g. I20 matches I20, I20.1, I20.2 and so on.
#
`%starts_with%` <- function(x, code_list) {
  if (length(code_list) == 0) {
    return(rep(FALSE, length(x)))
  }
  code_list <- gsub("\\.", "", code_list)
  code_list <- paste0("^", code_list)
  foreach(i = seq_along(code_list), .combine=`|`) %do% {
    tolower(x) %like% tolower(code_list[i])
  }
}

# Now we can classify each event as having a diabetes encoding or not, which will help us
# determine date of diagnosis (or maximum diabetes-free date). To do so, we need to set 
# unknown events (deaths without cause record available, maximum follow-up date) to NA, so
# we can interpolate last known diabetes free date as halfway between that date and the last
# record with diagnoses.
events[, c("t1d", "t2d", "diab_ns") := NA]
events[event == "UK Biobank assessment", c("t1d", "t2d", "diab_ns") := FALSE] # handled by prevalence algorithm

hes_known_cause <- hes[!is.na(code) & event_type == "hospitalisation", .SD[1], 
  by=.(eid, start_date=event_date, end_date=episode_end)]
hes_known_cause[, event := "Hospital episode"]
events[hes_known_cause, on = .(eid, event), c("t1d", "t2d", "diab_ns") := FALSE]

hes_t1d <- hes[code_type == "ICD-10" & code %starts_with% "E10", 
  .(eid, event="Hospital episode", start_date=event_date, end_date=episode_end)]
events[hes_t1d, on = .(eid, event, start_date, end_date), t1d := TRUE]

hes_t2d <- hes[code_type == "ICD-10" & code %starts_with% "E11", 
  .(eid, event="Hospital episode", start_date=event_date, end_date=episode_end)]
events[hes_t2d, on = .(eid, event, start_date, end_date), t2d := TRUE]

hes_diab_ns <- hes[code_type == "ICD-10" & code %starts_with% c("E13", "E14"), 
  .(eid, event="Hospital episode", start_date=event_date, end_date=episode_end)]
events[hes_diab_ns, on = .(eid, event, start_date, end_date), diab_ns := TRUE]

deaths_known_cause <- deaths[,.SD[1],by=.(eid)][,.(eid, event="Death")]
events[deaths_known_cause, on = .(eid, event), c("t1d", "t2d", "diab_ns") := FALSE]

death_t1d <- deaths[cause_icd10 %starts_with% "E10", .(eid, event="Death", 
  start_date=date_of_death, end_date=date_of_death)]
events[death_t1d, on = .(eid, event, start_date, end_date), t1d := TRUE]

death_t2d <- deaths[cause_icd10 %starts_with% "E11", .(eid, event="Death", 
  start_date=date_of_death, end_date=date_of_death)]
events[death_t2d, on = .(eid, event, start_date, end_date), t2d := TRUE]

death_diab_ns <- deaths[cause_icd10 %starts_with% c("E13", "E14"), .(eid,
  event="Death", start_date=date_of_death, end_date=date_of_death)]
events[death_diab_ns, on = .(eid, event, start_date, end_date), diab_ns := TRUE]

events[, diab := t1d | t2d | diab_ns ]

# Recode instances where diabetes type is uncertain (e.g. where an
# event has both a T1D and T2D diagnosis, or an uncertain diabetes
# diagnosis co-occuring with a T2D or T1D diagnosis)
events[t1d + t2d + diab_ns > 1, c("t1d", "t2d", "diab_ns") := .(FALSE, FALSE, TRUE)]

#########################################################################
# Interpolate date of diabetes diagnosis or maximum diabetes free date
#########################################################################

cat("Interpolating date of diabetes diagnosis or maximum diabetes free data...\n\n")

diab_followup <- followup[,.(eid, visit_index, assessment_date, 
  interpolated_follow_date=as.IDate(NA_character_), interpolated_follow_years=NA_real_,
  max_follow_date=as.IDate(NA_character_), prev_record_date=as.IDate(NA_character_),
  any_diab=NA, t1d_first=NA, t2d_first=NA, diab_ns_first=NA,
  any_diab_death=NA, t1d_death=NA, t2d_death=NA, diab_ns_death=NA,
  t1d_records=NA_integer_, t2d_records=NA_integer_, diab_ns_records=NA_integer_, any_diab_records=NA_integer_)]

diab_followup[events, on = .(eid), c("any_diab", "t1d_first", "t2d_first", "diab_ns_first") := FALSE]
diab_followup[events, on = .(eid), c("t1d_records", "t2d_records", "diab_ns_records", "any_diab_records") := 0]

diab_followup[!followup[(all_cause_mortality)], on = .(eid, visit_index), diab_death_record := FALSE]
diab_followup[deaths, on = .(eid), c("any_diab_death", "t1d_death", "t2d_death", "diab_ns_death") := FALSE] # NAs are now deaths with no available cause information

# Add in first diagnosis of any diabetes
first_diab <- events[(diab), .SD[1], by=.(eid)]
diab_followup[first_diab, on = .(eid), c("t1d_first", "t2d_first", "diab_ns_first", "any_diab") := .(i.t1d, i.t2d, i.diab_ns, i.diab)]

# Add in information on whether there is diabetes in the death record
diab_followup[events[event == "Death" & (diab)], on = .(eid), any_diab_death:= TRUE]
diab_followup[events[event == "Death" & (t1d)], on = .(eid), t1d_death:= TRUE]
diab_followup[events[event == "Death" & (t2d)], on = .(eid), t2d_death:= TRUE]
diab_followup[events[event == "Death" & (diab_ns)], on = .(eid), diab_ns_death:= TRUE]

# Add in number of records for each type of diabetes
diab_followup[events[(diab), .N, by=eid], on = .(eid), any_diab_records := i.N]
diab_followup[events[(t1d), .N, by=eid], on = .(eid), t1d_records := i.N]
diab_followup[events[(t2d), .N, by=eid], on = .(eid), t2d_records := i.N]
diab_followup[events[(diab_ns), .N, by=eid], on = .(eid), diab_ns_records := i.N]

# For a diabetes diagnosis, the interpolated date of diagnosis is the midpoint between
# the first diabetes diagnosis and last diabetes-free event
first_diab[, prev_event_idx := event_idx - 1]
first_diab[events, on = .(eid, prev_event_idx=event_idx), last_diab_free_date := i.end_date]
first_diab[, interpolated_diagnosis_date := round(start_date - (start_date - last_diab_free_date)/2)]
diab_followup[first_diab, on = .(eid), 
  c("interpolated_follow_date", "max_follow_date", "prev_record_date") :=
  .(i.interpolated_diagnosis_date, i.last_diab_free_date, i.start_date)]

# For those without diabetes diagnosis who have died, the max follow-up date the 
# date of death.
no_diab <- events[!diab_followup[(any_diab)], on = .(eid)]
non_diab_death <- no_diab[event == "Death" & !is.na(diab)]
diab_followup[non_diab_death, on = .(eid), 
  c("interpolated_follow_date", "max_follow_date", "prev_record_date") :=
  .(i.start_date, i.start_date, i.start_date)]

# For those without any diabetes diagnosis who have died, but for whom there is
# no information available on cause of death, we can't rule out a diabetes diagnosis
# (although unlikely) at the stage, so we set these to NA
miss_death <- no_diab[event == "Death" & is.na(diab)]
diab_followup[miss_death, on = .(eid), c("t1d_first", "t2d_first", "diab_ns_first", "any_diab") := NA]

# For people still alive at maximum follow-up, or who have died but with no available
# information on cause of death, the maximum diabetes-free date is the midpoint between
# the last diabetes free timepoint and the next record with missing diagnosis data
# (one of maximum hospital record follow-up date, maximum death record follow-up date,
# lost to follow-up date, or date of death, whichever occurs earliest).
no_diab <- no_diab[!non_diab_death, on = .(eid)]
last_diab_free <- no_diab[!is.na(diab)][, .SD[.N], by=eid]
last_diab_free[, next_event_idx := event_idx + 1]
last_diab_free[events, on = .(eid, next_event_idx=event_idx), last_follow_date := i.start_date]
last_diab_free[, max_diabetes_free_date := round(end_date + (last_follow_date - end_date)/2)]
diab_followup[last_diab_free, on = .(eid), 
  c("interpolated_follow_date", "max_follow_date", "prev_record_date") :=
  .(i.max_diabetes_free_date, i.last_follow_date, i.end_date)]

# Compute follow-up time in years

# Function to compute years between two dates - note this preserves
# human notions of whole years, i.e:
#
#   2019-07-25 - 2009-07-25 = 10 years
#
# At the expense of (potentially) leading to small inaccuracies in
# rank ordering due to leap days, i.e. in pure terms of days,
#
#   2019-07-25 - 2009-07-25 = 9.9986 years
#
# As this makes it harder to accurately truncate follow-up (e.g.
# for testing model calibration)
years_between <- function(d1, d2) {
  year_diff <- as.period(interval(as.Date(d1), as.Date(d2))) / years(1)

  # When d1 or d2 is a feb 29th, we want to return a whole number for
  # number of years between d1 and d2 if the other date is (1) Feb 28th,
  # and (2) not also a leap year.
  leap_day_d1 <- which(month(d1) == 2 & day(d1) == 29)
  leap_day_d2 <- which(month(d2) == 2 & day(d1) == 29)
  feb_28_d1 <- which(month(d1) == 2 & day(d2) == 28 & !(leap_year(year(d1))))
  feb_28_d2 <- which(month(d2) == 2 & day(d2) == 28 & !(leap_year(year(d2))))

  to_correct <- union(intersect(leap_day_d1, feb_28_d2), intersect(feb_28_d1, leap_day_d2))
  year_diff[to_correct] <- year(d2[to_correct]) - year(d1[to_correct])
  return(year_diff)
}

diab_followup[, interpolated_follow_years := years_between(assessment_date, interpolated_follow_date)]

#####################################
# Curate evidence of diabetes
#####################################

cat("Curating evidence of incident diabetes using Eastwood algorithm...\n\n")

# Create empty dataset for adjudicating incident diabetes
dat <- followup[, .(eid, visit_index, assessment_date, 
  adjudicated_diabetes = NA_character_, adjudication_finalized = NA_character_, 
  adjudication_note = NA_character_)]

# compute total numbers per visit for diagnostic output
totals <- dat[,.N,by=visit_index]

# 1.1 Exclude people with prevalent diabetes
#
# Note: this only includes those withe "probable" type 1 or type 2
# diabetes as adjudicated by the Eastwood et al. algorithm (see 
# src/01_adjudicate_prevalent_diabetes.R)

# Load prevalent diabetes adjudication
prev_diab <- fread("output/prevalent_diabetes.txt")

prob_t1d <- prev_diab[adjudicated_diabetes == "Probable type 1 diabetes"]
prob_t2d <- prev_diab[adjudicated_diabetes == "Probable type 2 diabetes"]

dat[rbind(prob_t1d, prob_t2d), on = .(eid, visit_index), adjudicated_diabetes := "Prevalent diabetes"]
dat[rbind(prob_t1d, prob_t2d), on = .(eid, visit_index), adjudication_finalized := "Step 1.1"]
dat[prob_t1d, on = .(eid, visit_index), adjudication_note := "1.1: Prevalent type 1 diabetes (probable)."]
dat[prob_t2d, on = .(eid, visit_index), adjudication_note := "1.1: Prevalent type 2 diabetes (probable)."]
exited_1.1 <- dat[!is.na(adjudicated_diabetes)]
dat <- dat[!exited_1.1, on = .(eid, visit_index)]

cat("Exited at step 1.1 with adjudication reason \"Prevalent diabetes\" (self-report):\n")
cat(" ", format(exited_1.1[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_1.1[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_1.1[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_1.1[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

# Drop clear cut prevalent cases from follow-up table
diab_followup <- diab_followup[!exited_1.1, on = .(eid, visit_index)]

# Add note warning of possible diabetes
dat[prev_diab[adjudicated_diabetes == "Possible type 1 diabetes"], on = .(eid, visit_index),
  adjudication_note := "1.1: Prevalent type 1 diabetes (possible)."]
dat[prev_diab[adjudicated_diabetes == "Possible type 2 diabetes"], on = .(eid, visit_index),
  adjudication_note := "1.1: Prevalent type 2 diabetes (possible)."]
dat[prev_diab[adjudicated_diabetes == "Possible gestational diabetes"], on = .(eid, visit_index),
  adjudication_note := "1.1: Prevalent gestational diabetes (possible)."]

# 1.2 remove people without inpatient data
#
# Note: this step doesn't exclude as many people as the paper.
# Eastwood et al. exclude anyone without inpatient data at this step,
# however, a person may not have inpatient data due to good health 
# (i.e. no hospital visits) rather than lack of follow-up. Supporting
# this, even if we use the same criteria, about half as many samples 
# are excluded due to the increase in follow-up data available now
# (i.e. hospitalisations since ~2014/2015 when the Eastwood analysis
# would have been done). Further, UKB separately curates information
# on presence/absence of follow-up which we can use to identify people
# who are missing hospital linkage.

no_followup <- dat[!diab_followup, on = .(eid)] # N = 0
linkage_withdrawn <- dat[eid %in% followup[lost_to_followup_reason == "Participant has withdrawn consent for future linkage", eid]]

no_followup[, adjudicated_diabetes := "No evidence of diabetes"]
no_followup[, adjudication_finalized := "Step 1.2"]
no_followup[, adjudication_note := paste(adjudication_note, "1.2: No health record linkage.")]
no_followup[, adjudication_note := gsub("^NA ", "", adjudication_note)]

linkage_withdrawn[, adjudicated_diabetes := "No evidence of diabetes"]
linkage_withdrawn[, adjudication_finalized := "Step 1.2"]
linkage_withdrawn[, adjudication_note := paste(adjudication_note, "1.2: Withdrawn consent for health record linkage.")]
linkage_withdrawn[, adjudication_note := gsub("^NA ", "", adjudication_note)]

exited_1.2 <- rbind(no_followup, linkage_withdrawn)
dat <- dat[!exited_1.2, on = .(eid, visit_index)]

cat("Exited at step 1.2 with adjudication reason \"No evidence of diabetes (no health record linkage)\":\n")
cat(" ", format(exited_1.2[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_1.2[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_1.2[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_1.2[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

# Also drop from diabetes follow-up table (i.e. those with consent for linkage withdrawn)
diab_followup <- diab_followup[!exited_1.2, on = .(eid, visit_index)]

# 1.3 Remove people with no diabetes in linked health records
no_diab <- diab_followup[any_diab_records == 0]

exited_1.3 <- dat[no_diab[,.(eid, visit_index)], on = .(eid, visit_index)]
exited_1.3[, adjudicated_diabetes := "No evidence of diabetes"]
exited_1.3[, adjudication_finalized := "Step 1.3"]
exited_1.3[followup[lost_to_followup_date < assessment_date], on = .(eid, visit_index), adjudication_note := paste(adjudication_note,
  sprintf("1.3: Lost to follow-up prior to UK Biobank assessment visit with reason \"%s\" and no evidence of prevalent diabetes.", lost_to_followup_reason))] 
exited_1.3[followup[all_cause_mortality & !lost_to_followup], on = .(eid, visit_index), adjudication_note := paste(adjudication_note, 
  "1.3: Died from non-diabetes related causes, no diabetes in health records.")]
exited_1.3[followup[(lost_to_followup) & lost_to_followup_date >= assessment_date], on = .(eid, visit_index), adjudication_note := paste(adjudication_note,
  sprintf("1.3: No diabetes until lost to follow-up with reason \"%s\".", lost_to_followup_reason))]
exited_1.3[followup[!all_cause_mortality & !lost_to_followup], on = .(eid, visit_index), adjudication_note := paste(adjudication_note,
  sprintf("1.3: Alive at maximum follow-up date (%s for hospitals in %s) with no diabetes in health records.",
    latest_hospital_date, latest_hospital_nation))] 
exited_1.3[, adjudication_note := gsub("^NA ", "", adjudication_note)]
dat <- dat[!exited_1.3, on = .(eid, visit_index)]

cat("Exited at step 1.3 with adjudication reason \"No evidence of diabetes\":\n")
cat(" ", format(exited_1.3[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_1.3[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_1.3[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_1.3[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

# 1.4 There are no people without diagnosis date
no_date <- diab_followup[(any_diab) & is.na(interpolated_follow_date)]
exited_1.4 <- dat[no_date[,.(eid, visit_index)], on = .(eid, visit_index)]
exited_1.4[, adjudicated_diabetes := "Diabetes with no diagnosis date"]
exited_1.4[, adjudication_finalized := "Step 1.4"]
exited_1.4[no_date, on = .(eid, visit_index), adjudication_note := paste(adjudication_note, 
  sprintf("1.4: %s records of type 1 diabetes, %s records of type 2 diabetes, %s records of non-specific diabetes.",
    t1d_records, t2d_records, diab_ns_records))]
exited_1.4[, adjudication_note := gsub("^NA ", "", adjudication_note)]
dat <- dat[!exited_1.4, on = .(eid, visit_index)]

cat("Exited at step 1.4 with adjudication reason \"Diabetes with no diagnosis date\":\n")
cat(" ", format(exited_1.4[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_1.4[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_1.4[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_1.4[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

# 1.5 We've already interpolated the date of diagnosis above. 
# 
# Apart from that, step 1.5 doesn't appear to have an conditional statement?
# The flowchart in the paper suggests filtering on participants with >1 in-patient
# episodes containing diabetes codes, but does not list the number passing/failing
# this condition, and makes no other reference to this in the paper. 
#
# In subsequent steps we'll therefore make a compromise and triage cases into
# "possible" and "probable" based on the number of reports (1, or >1, respectively)
# indicating the confidence in the diagnosis code.
#
# We'll also treat death records as more reliable, so if there is a diabetes diagnosis
# in the death records we'll set the diagnosis as "probable" instead of "possible"
# if its the only record.

cat("No conditional statement in step 1.5, all remaining participants continue in algorithm...\n\n")

# 1.6 identify additional prevalent cases from the hospital inpatient records.
prev_diab_hes <- diab_followup[interpolated_follow_years < 0]
exited_1.6 <- dat[prev_diab_hes, on = .(eid, visit_index), nomatch=0]
exited_1.6[, adjudicated_diabetes := "Prevalent diabetes"]
exited_1.6[, adjudication_finalized := "Step 1.6"]
exited_1.6[, adjudication_note := paste(adjudication_note, sprintf(
  "1.6: %s first, %s type 1 diabetes records, %s type 2 diabetes records, %s non-specific diabetes records.",
   fcase(t1d_first, "type 1 diabetes", t2d_first, "type 2 diabetes", diab_ns_first, "non-specific diabetes"),
   t1d_records, t2d_records, diab_ns_records))]
exited_1.6[, adjudication_note := gsub("^NA ", "", adjudication_note)]
exited_1.6 <- exited_1.6[, .SD, .SDcols=names(dat)]
dat <- dat[!exited_1.6, on = .(eid, visit_index)]
diab_followup <- diab_followup[!exited_1.6, on = .(eid, visit_index)]
diab_followup <- diab_followup[!(interpolated_follow_years < 0 & any_diab_records == 0)] # Left UK, but returned for imaging assessment.

cat("Exited at step 1.6 with adjudication reason \"Prevalent diabetes\" (in hospital records):\n")
cat(" ", format(exited_1.6[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_1.6[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_1.6[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_1.6[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

# 1.7 Identify incident type 1 diabetes
exited_1.7 <- dat[diab_followup[(t1d_first)], on = .(eid, visit_index), nomatch=0]
exited_1.7[t1d_records > 1 | (t1d_death) | 
  (adjudication_note %like% "1.1 Prevalent type 1 diabetes (possible)"),
  adjudicated_diabetes := "Probable incident type 1 diabetes"]
exited_1.7[is.na(adjudicated_diabetes), adjudicated_diabetes := "Possible incident type 1 diabetes"]
exited_1.7[, adjudication_finalized := "Step 1.7"]
exited_1.7[, adjudication_note := paste(adjudication_note, "1.7:")]
exited_1.7[(t1d_death), adjudication_note := paste(adjudication_note, "type 1 diabetes death record,")]
exited_1.7[any_diab_death & !t1d_death, adjudication_note := paste(adjudication_note, "death record with other diabetes,")]
exited_1.7[events[event == "Death" & !(diab)], on = .(eid), adjudication_note := paste(adjudication_note, "death record with no diabetes,")]
exited_1.7[, adjudication_note := paste(adjudication_note, sprintf(
  "%s type 1 diabetes records, %s type 2 diabetes records, %s non-specific diabetes records.",
   t1d_records, t2d_records, diab_ns_records))]
exited_1.7[, adjudication_note := gsub("^NA ", "", adjudication_note)]
exited_1.7 <- exited_1.7[, .SD, .SDcols=names(dat)]
dat <- dat[!exited_1.7, on = .(eid, visit_index)]

exited_1.7_probable <- exited_1.7[adjudicated_diabetes == "Probable incident type 1 diabetes"]
exited_1.7_possible <- exited_1.7[adjudicated_diabetes == "Possible incident type 1 diabetes"]

cat("Exited at step 1.7 with adjudication reason \"Probable incident type 1 diabetes\":\n")
cat(" ", format(exited_1.7_probable[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_1.7_probable[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_1.7_probable[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_1.7_probable[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

cat("Exited at step 1.7 with adjudication reason \"Possible incident type 1 diabetes\":\n")
cat(" ", format(exited_1.7_possible[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_1.7_possible[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_1.7_possible[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_1.7_possible[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

# 1.8 Identify incident type 2 diabetes
exited_1.8 <- dat[diab_followup[(t2d_first)], on = .(eid, visit_index), nomatch=0]
exited_1.8[t2d_records > 1 | (t2d_death) | 
  (adjudication_note %like% "1.1 Prevalent type 2 diabetes (possible)"),
  adjudicated_diabetes := "Probable incident type 2 diabetes"]
exited_1.8[is.na(adjudicated_diabetes), adjudicated_diabetes := "Possible incident type 2 diabetes"]
exited_1.8[, adjudication_finalized := "Step 1.8"]
exited_1.8[, adjudication_note := paste(adjudication_note, "1.8:")]
exited_1.8[(t2d_death), adjudication_note := paste(adjudication_note, "type 2 diabetes death record,")]
exited_1.8[any_diab_death & !t2d_death, adjudication_note := paste(adjudication_note, "death record with other diabetes,")]
exited_1.8[events[event == "Death" & !(diab)], on = .(eid), adjudication_note := paste(adjudication_note, "death record with no diabetes,")]
exited_1.8[, adjudication_note := paste(adjudication_note, sprintf(
  "%s type 1 diabetes records, %s type 2 diabetes records, %s non-specific diabetes records.",
   t1d_records, t2d_records, diab_ns_records))]
exited_1.8[, adjudication_note := gsub("^NA ", "", adjudication_note)]
exited_1.8 <- exited_1.8[, .SD, .SDcols=names(dat)]
dat <- dat[!exited_1.8, on = .(eid, visit_index)]

exited_1.8_probable <- exited_1.8[adjudicated_diabetes == "Probable incident type 2 diabetes"]
exited_1.8_possible <- exited_1.8[adjudicated_diabetes == "Possible incident type 2 diabetes"]

cat("Exited at step 1.8 with adjudication reason \"Probable incident type 2 diabetes\":\n")
cat(" ", format(exited_1.8_probable[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_1.8_probable[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_1.8_probable[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_1.8_probable[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

cat("Exited at step 1.8 with adjudication reason \"Possible incident type 2 diabetes\":\n")
cat(" ", format(exited_1.8_possible[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_1.8_possible[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_1.8_possible[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_1.8_possible[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

# 1.9 Identify incident non-specific diabetes
exited_1.9 <- dat[diab_followup[(diab_ns_first)], on = .(eid, visit_index), nomatch=0]
exited_1.9[is.na(adjudicated_diabetes), adjudicated_diabetes := "Incident diabetes of uncertain type"]
exited_1.9[, adjudication_finalized := "Step 1.9"]
exited_1.9[, adjudication_note := paste(adjudication_note, "1.9:")]
exited_1.9[(t1d_death), adjudication_note := paste(adjudication_note, "type 1 diabetes death record,")]
exited_1.9[(t2d_death), adjudication_note := paste(adjudication_note, "type 2 diabetes death record,")]
exited_1.9[(diab_ns_death), adjudication_note := paste(adjudication_note, "death record with diabetes of uncertain type")]
exited_1.9[events[event == "Death" & !(diab)], on = .(eid), adjudication_note := paste(adjudication_note, "death record with no diabetes,")]
exited_1.9[, adjudication_note := paste(adjudication_note, sprintf(
  "%s type 1 diabetes records, %s type 2 diabetes records, %s non-specific diabetes records.",
   t1d_records, t2d_records, diab_ns_records))]
exited_1.9[, adjudication_note := gsub("^NA ", "", adjudication_note)]
exited_1.9 <- exited_1.9[, .SD, .SDcols=names(dat)]
dat <- dat[!exited_1.9, on = .(eid, visit_index)]

cat("Exited at step 1.9 with adjudication reason \"Incident diabetes of uncertain type\":\n")
cat(" ", format(exited_1.9[visit_index == 0, .N], big.mark=","), "of", format(totals[visit_index == 0, N], big.mark=","), "at basline assessment (2006-2010)\n")
cat(" ", format(exited_1.9[visit_index == 1, .N], big.mark=","), "of", format(totals[visit_index == 1, N], big.mark=","), "at first repeat assessment (2009-2014)\n")
cat(" ", format(exited_1.9[visit_index == 2, .N], big.mark=","), "of", format(totals[visit_index == 2, N], big.mark=","), "at first imaging assessment (2014-2020)\n")
cat(" ", format(exited_1.9[visit_index == 3, .N], big.mark=","), "of", format(totals[visit_index == 3, N], big.mark=","), "at repeat imaging assessment (2019-ongoing)\n")
cat("\n")

########################################################
# Collate all results
########################################################

final <- rbind(exited_1.1, exited_1.2, exited_1.3, exited_1.4, 
  exited_1.6, exited_1.7, exited_1.8, exited_1.9)

final <- merge(
  diab_followup[,.(eid, visit_index, assessment_date, 
    interpolated_follow_years, interpolated_follow_date,
    max_follow_date, prev_record_date)], 
  final, by=c("eid", "visit_index", "assessment_date"), 
  all.y=TRUE)

final <- final[assessments, on = .(eid, visit_index, assessment_date), nomatch=0]
diab_followup <- diab_followup[assessments, on = .(eid, visit_index, assessment_date), nomatch=0]

fwrite(final, sep="\t", quote=FALSE, file="output/incident_diabetes.txt")
fwrite(diab_followup, sep="\t", quote=FALSE, file="output/incident_diabetes_extended_info.txt")

