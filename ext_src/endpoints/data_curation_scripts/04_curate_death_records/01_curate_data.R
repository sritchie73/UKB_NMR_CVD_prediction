library(data.table)
library(lubridate)

system("mkdir -p data/curated/ukbiobank/deaths/", wait=TRUE)

# Load death records
deaths <- fread("data/raw/ukbiobank/Death/death_20210727.txt")
deaths[, date_of_death := as.IDate(date_of_death, format="%d/%m/%Y")]

# Load cause of death
death_cause <- fread("data/raw/ukbiobank/Death/death_cause_20210727.txt")

# Get primary cause of death. These are always in arr_index 0. Where a
# participant has entries with ins_index = 1 we take this as the primary
# cause. These entries correspond to a second later issuing of a death
# certificate, e.g. issued after a post-mortem, and thus usually give a 
# more refined cause of death (e.g. instead of R99, unspecified illness, 
# F10.2, alcohol dependence).
# (see https://biobank.ndph.ox.ac.uk/showcase/ukb/docs/DeathLinkage.pdf)
post_mortem <- death_cause[ins_index == 1]
no_post_mortem <- death_cause[!post_mortem, on = .(eid)]
death_cause <- rbind(no_post_mortem, post_mortem)

# organise
death_cause <- death_cause[, .(eid, cause_icd10, cause_number = arr_index + 1, primary_cause = level == 1)]
death_cause <- death_cause[order(cause_number)][order(eid)]

# Add in event date
death_cause[deaths, on = .(eid), date_of_death := i.date_of_death]

# Add in data source
death_cause[deaths, on = .(eid), death_source := fcase(
  i.dsource == "E/W", "NHS Digital records of death within England or Wales",
  i.dsource == "SCOT", "NHS Central Register of deaths within Scotland"
)]

# Filter to events in the curated follow-up table (e.g. removing samples with
# consent withdrawn or lost to follow-up)
followup <- fread("data/curated/ukbiobank/followup/followup.txt")
followup <- followup[(all_cause_mortality), .(eid, date_of_death=latest_mortality_date)]
followup <- unique(followup)

death_cause <- death_cause[followup, on = .(eid, date_of_death), nomatch=0]

# Write out
fwrite(death_cause, sep="\t", quote=FALSE, file="data/curated/ukbiobank/deaths/deaths.txt")

# Curate information
info <- rbind(
  data.table(var="cause_icd10", description="ICD-10 code for cause of death"),
  data.table(var="cause_number", description="Cause number, where 1 = primary cause of death"),
  data.table(var="primary_cause", description="TRUE if the ICD-10 code was recorded as the primary cause of death"),
  data.table(var="date_of_death", description="Date of death"),
  data.table(var="death_source",  description="Either 'NHS Digital records of death within England or Wales' or 'NHS Central Register of deaths within Scotland'")
)

fwrite(info, sep="\t", quote=FALSE, file="data/curated/ukbiobank/deaths/column_info.txt")






