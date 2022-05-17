library(data.table)
library(lubridate)

system("mkdir -p data/curated/ukbiobank/followup/", wait=TRUE)

# Load previously curated anthropometric data
dt <- fread("data/curated/ukbiobank/anthropometrics/anthropometrics.txt")

# Filter to columns relevant for curating follow-up details
dt <- dt[, .(eid, visit_index, sex, age, age_decimal, assessment_date, assessment_centre, assessment_nation)]

# Load and information regarding lost to follow-up information
lost <- fread("data/raw/ukbiobank/extracted/followup/followup.csv")

lost <- lost[, .(eid, lost_to_followup=FALSE, 
  lost_to_followup_reason=`190-0.0`,
  lost_to_followup_date=`191-0.0`)]

lost[!is.na(lost_to_followup_reason) | !is.na(lost_to_followup_date), lost_to_followup := TRUE]

lost[, lost_to_followup_reason := fcase(
  lost_to_followup_reason == 1, "Death reported to UK Biobank by a relative",
  lost_to_followup_reason == 2, "NHS records indicate they are lost to follow-up",
  lost_to_followup_reason == 3, "NHS records indicate they have left the UK",
  lost_to_followup_reason == 4, "UK Biobank sources report they have left the UK",
  lost_to_followup_reason == 5, "Participant has withdrawn consent for future linkage"
)]

# Add to table
dt <- merge(dt, lost, by = "eid", all.x=TRUE)

# Load hospital records
hes <- fread("data/raw/ukbiobank/HES/hesin_20210727.txt")
hes[, epistart := as.IDate(epistart, format="%d/%m/%Y")]
hes[, epiend := as.IDate(epiend, format="%d/%m/%Y")]
hes[, elecdate := as.IDate(elecdate, format="%d/%m/%Y")]
hes[, admidate := as.IDate(admidate, format="%d/%m/%Y")]
hes[, disdate := as.IDate(disdate, format="%d/%m/%Y")]

# For the HES data, get the date of each event. This is the episode start date,
# or if missing the earlier of the date of decision to admit to hospital (elecdate)
# or actual hospital admission date (admidate)
hes[, event_date := epistart]
hes[is.na(event_date), event_date := pmin(admidate, elecdate, na.rm=TRUE)]

# Load death records
deaths <- fread("data/raw/ukbiobank/Death/death_20210727.txt")
deaths[, date_of_death := as.IDate(date_of_death, format="%d/%m/%Y")]

# Compute follow-up available per nation - they each have different cutoff and start
# dates
hes[, nation := fcase(
  dsource == "HES", "England",
  dsource == "PEDW", "Wales",
  dsource == "SMR", "Scotland"
)]
hes_follow <- hes[!is.na(event_date),.(min_follow=min(event_date), max_follow=max(event_date)), by=nation]

# Compute follow-up available per nation
death_follow <- deaths[,.(max_follow=max(date_of_death)), by=dsource]
death_follow <- rbind(idcol="nation",
  "England"=death_follow[dsource == "E/W"],
  "Wales"=death_follow[dsource == "E/W"],
  "Scotland"=death_follow[dsource == "SCOT"]
)

# Output information
sink("data/curated/ukbiobank/followup/max_followup_by_nation.txt")
cat("Death records (ONS):\n")
print(death_follow[])
cat("\nHospital records (HES):\n")
print(hes_follow[])
cat("\n")
sink()

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

# Add information on any hospital events
dt[, any_hospitalisations := FALSE]
dt[hes, on = .(eid), any_hospitalisations := TRUE]

# Use the hospital records to identify earliest and most recent hospital nation
# for each person so we know which follow-up date cutoffs to use
earliest <- hes[,.SD[which.min(event_date)], by=.(eid)][,.(eid, event_date, nation)]
latest <- hes[,.SD[which.max(event_date)], by=.(eid)][,.(eid, event_date, nation)]

# Add in retrospective hospital follow-up information
dt[earliest, on = .(eid, assessment_date > event_date), earliest_hospital_nation := i.nation]
dt[is.na(earliest_hospital_nation), earliest_hospital_nation := assessment_nation] # i.e. people with no hospitalisations prior to assessment
dt[dt[visit_index == 0], on = .(eid), earliest_hospital_nation := i.earliest_hospital_nation] # baseline overrules subsequent visits, i.e. to retrospective nation change between assessments

dt[hes_follow, on = .(earliest_hospital_nation=nation), earliest_hospital_date := min_follow] # note all min dates at least 10 years later than date of birth
dt[, earliest_hospital_followup := years_between(assessment_date, earliest_hospital_date)]

# Add prospective hospital follow-up information
dt[latest, on = .(eid, assessment_date < event_date), latest_hospital_nation := i.nation]

# Assessment visits that happened after the latest hospital visit overwrite nation -
# e.g. to prevent someone who moved England (baseline) -> Wales (latest hospital record) -> England (latest assessment)
# from being incorrectly classified
latest_assessment <- dt[,.SD[which.max(visit_index)], by=.(eid)][,.(eid, assessment_date, assessment_nation)]
latest_assessment <- latest_assessment[latest[, .(eid, event_date)], on = .(eid, assessment_date > event_date), nomatch=0]
dt[latest_assessment, on = .(eid), latest_hospital_nation := i.assessment_nation]

# For people with no prospective hospital records, use the most recent assessment visit nation
latest_assessment <- dt[,.SD[which.max(visit_index)], by=.(eid)][,.(eid, assessment_date, assessment_nation)]
no_hosp <- dt[is.na(latest_hospital_nation), .(eid, visit_index)]
no_hosp[latest_assessment, on = .(eid), latest_hospital_nation := assessment_nation]
dt[no_hosp, on = .(eid, visit_index), latest_hospital_nation := i.latest_hospital_nation]

# Add in date and follow-up
dt[hes_follow, on = .(latest_hospital_nation=nation), latest_hospital_date := max_follow] 
dt[, latest_hospital_followup := years_between(assessment_date, latest_hospital_date)]

# Add all cause mortality column
dt[, all_cause_mortality := FALSE]
dt[deaths, on = .(eid), all_cause_mortality := TRUE]

# Add death nation - mostly so we can cross-check hospital nation
dt[deaths, on = .(eid), death_source := fcase(
  i.dsource == "E/W", "NHS Digital records of death within England or Wales",
  i.dsource == "SCOT", "NHS Central Register of deaths within Scotland"
)]

# Add max follow-up for mortality based on nation
dt[death_follow, on = .(latest_hospital_nation=nation), latest_mortality_date := i.max_follow]

# For people who have died, truncate at date of death
dt[deaths, on = .(eid), latest_mortality_date := date_of_death]

# compute follow-up time
dt[, latest_mortality_followup := years_between(assessment_date, latest_mortality_date)]

# For people who have died, roll back max hospital record follow-up to that date
dt[(all_cause_mortality) & latest_mortality_date < latest_hospital_date,
   c("latest_hospital_date", "latest_hospital_followup") :=
   .(latest_mortality_date, latest_mortality_followup)]

# Cross-check hospital record nation against death records.
# In the case death_source == "E/W", we don't know if the death record came from England or Wales.
# In these cases we currently set the nation to England for several reasons:
#
#  (1) Balance of probabilities. England is much more populous than Wales, and ~90% of UK Biobank
#      participants are in England, so it's most likely these people are located in England at date
#      of death.
#
#  (2) None of the people satisfying this condition were located in Wales at assessment
#  
#  (3) Wales has shorter follow-up time in the hospital records than England. This leaves two possibilites: 
#       (1) they're in England and there's one or more hospital record after the maximum follow-up available in Wales.
#       (2) they're in Wales, but the hospital events are all prior to maximum follow-up in either nation, so being
#           less precise doesn't matter
#
#  (4) if we inspect the individual hospital records for these individuals, all events happen well before
#      the maximum follow-up available for hospitals in Wales (or England).
dt[death_source == "NHS Central Register of deaths within Scotland" & 
   latest_hospital_nation != "Scotland", 
   latest_hospital_nation := "Scotland"]
dt[death_source == "NHS Digital records of death within England or Wales" & 
   latest_hospital_nation == "Scotland", 
   latest_hospital_nation := "England"]

dt[hes_follow, on = .(latest_hospital_nation=nation, latest_hospital_date > max_follow),
   c("latest_hospital_date", "latest_hospital_followup") := 
   .(max_follow, years_between(assessment_date, max_follow))]

# Remove follow-up for people who have withdrawn consent for linkage
dt[lost_to_followup_reason == "Participant has withdrawn consent for future linkage",
   c("any_hospitalisations", "earliest_hospital_nation", "earliest_hospital_date", "earliest_hospital_followup",
     "latest_hospital_nation", "latest_hospital_date", "latest_hospital_followup", "all_cause_mortality", 
     "death_source", "latest_mortality_date", "latest_mortality_followup") := 
   .(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)]

# Add in death information reported by relative:
# 
# Note these are not necessarily reflected by ONS:
#
# - All but one have no death records in ONS, which could be explained e.g. by death abroad.
# - There is one instance where ONS reports the death many years later (2014 vs. 2021). This
#   could either be a linkage error (i.e. relative is correct), or incorrect reporting from
#   relative. In either case, we take the most conservative (earliest) date.
dt[lost_to_followup_reason == "Death reported to UK Biobank by a relative",
   c("latest_hospital_date", "latest_hospital_followup", "all_cause_mortality",
     "death_source", "latest_mortality_date", "latest_mortality_followup") := 
   .(lost_to_followup_date, years_between(assessment_date, lost_to_followup_date), TRUE,
     "Reported by relative to UK Biobank", lost_to_followup_date, years_between(assessment_date, lost_to_followup_date))]

# Others are reported as lost to follow-up either due to leaving the UK
# (reported by relative or NHS) or unspecified (reported by NHS). Some 
# of these people, however, have hospital or death records after these
# dates. Nevertheless, we truncate follow-up at these dates to take the
# most conservative approach.
dt[lost_to_followup_reason %in% c("NHS records indicate they have left the UK",
  "UK Biobank sources report they have left the UK", "NHS records indicate they are lost to follow-up"),
   c("latest_hospital_date", "latest_hospital_followup", "all_cause_mortality",
     "death_source", "latest_mortality_date", "latest_mortality_followup") := 
   .(lost_to_followup_date, years_between(assessment_date, lost_to_followup_date), FALSE,
     NA, lost_to_followup_date, years_between(assessment_date, lost_to_followup_date))]

# Update any_hospitalisations field for these people
hes <- hes[!dt[is.na(any_hospitalisations)], on = .(eid)] # drop withdrawn linkage 
hes <- unique(rbind(
  hes[unique(dt[, .(eid, earliest_hospital_date, latest_hospital_date)]), 
		  on = .(eid, event_date > earliest_hospital_date, event_date < latest_hospital_date),
		  nomatch=0, .(eid)],
  hes[is.na(event_date), .(eid)]
))
dt[!is.na(any_hospitalisations), any_hospitalisations := FALSE]
dt[hes, on = .(eid), any_hospitalisations := TRUE]

# Write out
fwrite(dt, sep="\t", quote=FALSE, file="data/curated/ukbiobank/followup/followup.txt")

# Curate information about each column:
info <- rbind(
  data.table(var="sex", description="Self-reported sex (field #31)"),
  data.table(var="age", description="Age at assessment (whole years) (field #21003)"),
  data.table(var="age_decimal", description=paste(
   "Age at assessment (decimal) computed from birth year (field #34) and birth month (field #52)",
   "assuming day of birth was at the mid-point (15th) of the month (specific date not made available by UK Biobank for privacy reasons)"
  )),
  data.table(var="assessment_date", description="Date of attending assessment centre (field #53)"),
  data.table(var="assessment_centre", description="Assessment centre (field #54)"),
  data.table(var="assessment_nation", description="Nation within the UK assessment centre is located in"),
  data.table(var="lost_to_followup", description="Participant has any records in fields #190 or #191 indicating lost to followup"),
  data.table(var="lost_to_followup_reason", description="Reason participant lost to follow-up (field #190)"),
  data.table(var="lost_to_followup_date", description="Date participant lost to follow-up (field #191)"),
  data.table(var="any_hospitalisation", description="TRUE where any hospital records exist for the participant, NA where linkage withdrawn ('lost_to_followup_reason')"),
  data.table(var="earliest_hospital_nation", description="Nation of earliest retrospective hospital record (or most likely nation)"),
  data.table(var="earliest_hospital_date", description="Earliest possible retrospective date hospital records might be available based on 'earliest_hospital_nation'"),
  data.table(var="earliest_hospital_followup", description="Time between 'assessment_date' and 'earliest_hospital_date' in years (decimal)"),
  data.table(var="latest_hospital_nation", description="Nation of most recent hospital record (or most likely nation)"),
  data.table(var="latest_hospital_date", description="Latest possible prospective date hospital records might be available based on 'latest_hospital_nation'"),
  data.table(var="latest_hospital_followup", description="Time between 'assessment_date' and 'latest_hospital_date' in years (decimal)"),
  data.table(var="all_cause_mortality", description="TRUE where any death records exist for the participant, NA where linkage withdrawn ('lost_to_followup_reason')"),
  data.table(var="death_source", description=paste(
    "One of 'NHS Digital records of death within England or Wales', 'NHS Central Register of deaths within Scotland',",
    "'Reported by relative to UK Biobank', or NA where linkage withdrawn"
  )),
  data.table(var="latest_mortality_date", description="Date of death, or latest possible date participant is not known to be dead"),
  data.table(var="latest_mortality_followup", description="Time between 'assessment_date' and 'mortality_date' in years (decimal)")
)

# Write out
fwrite(info, sep="\t", quote=FALSE, file="data/curated/ukbiobank/followup/column_info.txt")


