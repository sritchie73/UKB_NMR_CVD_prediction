library(data.table)
library(lubridate)

system("mkdir -p data/curated/ukbiobank/hospital_records/", wait=TRUE)

# UK Biobank documentation:
#
# Overview of hospital inpatient data:
# https://biobank.ndph.ox.ac.uk/showcase/ukb/docs/HospitalEpisodeStatistics.pdf
#
# Data dictionary for columns:
# https://biobank.ctsu.ox.ac.uk/crystal/ukb/docs/HESDataDic.xlsx
#

# Load hospital episode records
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

# Also get the event end. This is either the end of the episode, or if missing the
# discharge date. We use this to help identify episodes that contain fatal events.
hes[, episode_end := epiend]
hes[is.na(episode_end), episode_end := disdate]
hes[is.na(episode_end), episode_end := event_date]

# Obtain nation of event
hes[, hospital_nation := fcase(
  dsource == "HES", "England",
  dsource == "PEDW", "Wales",
  dsource == "SMR", "Scotland"
)]

# Filter to columns of interest
hes <- hes[,.(eid, ins_index, event_date, episode_end, hospital_nation)]

# Load diagnosis codes for each event
hes_diagnoses <- fread("data/raw/ukbiobank/HES/hesin_diag_20210727.txt") 

# Melt to long, keeping relevant columns
hes_diagnoses <- hes_diagnoses[, .(eid, ins_index, level, ICD9=diag_icd9, ICD10=diag_icd10)]
hes_diagnoses <- melt(hes_diagnoses, id.vars=c("eid", "ins_index", "level"), variable.name="code_type", value.name="code")
hes_diagnoses <- hes_diagnoses[!is.na(code) & code != ""]

# Add details on event dates, episode ends, and nation of hospital
hes_diagnoses <- merge(hes_diagnoses, hes, by=c("eid", "ins_index"), all.x=TRUE)

# Load operations performed in hospital
hes_operations <- fread("data/raw/ukbiobank/HES/hesin_oper_20210727.txt") 

# Melt to long, keeping relevant columns. Note event date here is the date of 
# the operation
hes_operations[, opdate := as.IDate(opdate, format="%d/%m/%Y")]
hes_operations <- hes_operations[, .(eid, ins_index, level, event_date=opdate, OPCS3=oper3, OPCS4=oper4)]
hes_operations <- melt(hes_operations, id.vars=c("eid", "ins_index", "level", "event_date"), variable.name="code_type", value.name="code")
hes_operations <- hes_operations[!is.na(code) & code != ""]

# For operations without operation dates, add in dates from the hospital episodes table
no_opdate <- hes_operations[is.na(event_date)]
no_opdate[hes, on = .(eid, ins_index), event_date := i.event_date]
hes_operations <- rbind(hes_operations[!is.na(event_date)], no_opdate)

# Add in episode end dates
hes_operations[hes, on = .(eid, ins_index), episode_end := i.episode_end]

# Add in nation of hospital
hes_operations[hes, on = .(eid, ins_index), hospital_nation := i.hospital_nation]

# Combine both tables
curated <- rbind(idcol="event_type", "hospitalisation"=hes_diagnoses, "operation"=hes_operations)

# Determine primary/secondary events/operations and external causes
curated[, level := fcase(
  level == 1, "primary",
  level == 2, "secondary",
  level == 3, "external"
)]

setnames(curated, "level", "cause_type")

# No longer need ins_index
curated[, ins_index := NULL]

# Load in follow-up information to identify fatal events and remove people with withdrawn
# consent for linkage
followup <- fread("data/curated/ukbiobank/followup/followup.txt")

withdrawn_consent <- followup[lost_to_followup_reason == "Participant has withdrawn consent for future linkage"]
curated <- curated[!withdrawn_consent, on = .(eid)]

# need to withhold events with no date so we can add back after filtering using 
# non-equi join.
#
# As an aside, I'm not sure what to make of these events. They're all (N=21) OPCS-4 
# codes X99.8 which is "No operation performed". Among these, there are N=4 participant
# whose only entry in the hospital records is an operation labelled X99.8.
#
# OPCS-4 code labels are listed in:
# https://biobank.ndph.ox.ac.uk/ukb/coding.cgi?id=240&nl=1
nodate <- curated[is.na(event_date)]

# Need to manually select event date column from the 'curated' table otherwise 
# overwritten with 'earliest_hospital_date' by non-equi join.
valid_follow <- followup[!is.na(earliest_hospital_date), .(eid, earliest_hospital_date, latest_hospital_date)]
valid_follow <- unique(valid_follow)
curated <- curated[valid_follow, 
  on = .(eid, event_date >= earliest_hospital_date, event_date <= latest_hospital_date), nomatch=0,
  .(eid, event_type, event_date=x.event_date, episode_end, hospital_nation, code, code_type, cause_type)] 

# Add back in those without dates
curated <- rbind(curated, nodate)

# Add in hyphen to code type
curated[, code_type := gsub("ICD", "ICD-", code_type)]
curated[, code_type := gsub("OPCS", "OPCS-", code_type)]

# Identify fatal episodes
curated[, fatal_event := FALSE]
curated[followup[(all_cause_mortality)], on = .(eid, event_date = latest_mortality_date), fatal_event := TRUE]

curated[, fatality_in_episode := FALSE]
curated[followup[(all_cause_mortality)], on = .(eid, event_date <= latest_mortality_date, episode_end >= latest_mortality_date), fatality_in_episode := TRUE]

# Impose sensible row order:
curated <- curated[order(code)]
curated <- curated[order(factor(cause_type, levels=c("primary", "secondary", "external")))]
curated <- curated[order(event_type)][order(event_date)][order(eid)]

# Write out 
fwrite(curated, sep="\t", quote=FALSE, file="data/curated/ukbiobank/hospital_records/hospital_records.txt")

# Curate column information
info <- rbind(
  data.table(var="event_type", description="Indicates whether the corresponding entry contains diagnosis codes or operation codes"),
  data.table(var="event_date", description="Date of hospital admission, hospital episode, or operation"),
  data.table(var="episode_end", description="Hospital discharge date or end of hospital episode. NA for operations"),
  data.table(var="hospital_nation", description="Nation of hospital (England, Wales, or Scotland)"),
  data.table(var="code", description="Diagnosis code or operation code, note multiple rows (codes) may be present for each event"),
  data.table(var="code_type", description="Type of code: either ICD-9, ICD-10, OPCS-3, or OPCS-4 indicating source of 'code'"),
  data.table(var="cause_type", description=paste(
    "Whether the diagnosis code is the primary cause, secondary cause, or external cause of hospitalisation",
    "or whether the operation code is the primary or secondary reason for operation"
  )),
  data.table(var="fatal_event", description="TRUE where the person died on the same day as the event/operation"),
  data.table(var="fatality_in_episode", description="TRUE where the person died during the hospital episode (between 'event_date' and 'episode_end')")
)

fwrite(info, sep="\t", quote=FALSE, file="data/curated/ukbiobank/hospital_records/column_info.txt")

