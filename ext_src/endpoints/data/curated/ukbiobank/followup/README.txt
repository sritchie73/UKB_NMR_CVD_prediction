Hospital and death record linkage follow-up
-------------------------------------------
Scott Ritchie, 15 December 2021
sr827@medschl.cam.ac.uk
-------------------------------------------

Data in this folder curates and harmonises the linked electronic health record
data available for hospitalisations and death records for UK Biobank participants.

The main file is 'followup.txt', which contains information for each participant:

  (1) Whether they have had any hospitalisations
  (2) The earliest retrospective date hospital records could possibly be available
  (3) The latest prospective date hospital records could possibly be available
  (4) Whether there is any record that the participant has died
  (5) The latest prospective date death records (or lack thereof) could possibly be available
  (6) Time in years between date of assessment and dates above for survival modelling
  (7) Other basic demographic information, e.g. age and sex.

The data available in 'followup.txt' are summarised in 'column_info.txt'.

These data were harmonized from multiple sources:

  (1) Hospital inpatient records from hospitals in England
  (2) Hospital inpatient records from hospitals in Wales
  (3) Hospital inpatient records from hospitals in Scotland
  (4) Death records from NHS Digital for deaths in England or Wales
  (5) Death records from the National Records of Scotland
  (6) Ongoing UK Biobank records on participants lost to follow-up

Each of these sources of information have different date ranges for which linked follow-up
is available. These are recorded in 'max_followup_by_nation.txt' and change periodically 
as new data is linked and released by the NHS and UK Biobank. As of 15th December 2021, 
these date ranges are:

Hospital records:
     nation min_follow max_follow
1:  England 1993-07-27 2021-03-31
2: Scotland 1980-12-02 2021-05-07
3:    Wales 1991-04-18 2018-03-06

Death records:
           nation max_follow
1:  England/Wales 2021-03-20
2:       Scotland 2021-03-23

Note in particular:

  (1) there's dramatically different retrospective linkage of hospital records depending
      on the nation the hospital is located in.
  (2) prospective follow-up of hospital records is shorter in Wales than England or Scotland

An aside on event dates in the hospital records. These were taken as the start date of the
hospital episode ('epistart') or if missing, whichever was earliest (and non-missing) of 
the hospital admission date ('admidate') or decision to admit to hospital ('elecdate'). The
maximum and minimum follow-up dates were then determined as the earliest and latest event
dates respectively among the hospital records from each nation among all UK Biobank 
participants.

The 'earliest_hospital_date' and 'latest_hospital_date' for each participant represent
hypotheticals, rather than actual hospital events. These are intended for use in downstream
endpoint customisation, i.e. as minimum/maximum follow-up dates for people *without* events.

The 'earliest_hospital_nation' and 'latest_hospital_nation' column indicate nation used to
obtain the hypothetical minimum and maximum follow-up dates, and conseqeuntly also represent
an interpolated (assumed) nation of residence at these cut-off dates. In most (98.5%) of cases
these are identical to the nation of residence at baseline assessment. The remaining cases are
those where there is evidence the participant moved between nations prior to, or after, baseline
assessment, based on hospital inpatient records in different nations. In particular, the 
'earliest_hospital_nation' is the nation the respective hospital is located in for the earliest
hospital record for each participant (prior to baseline), and similarly the 'latest_hospital_nation'
is the nation the respective hospital is located in for the most recent hospital record (after 
baseline) for each participant.

These dates (and nations) were further harmonized against (1) death records, (2) lost to follow-up
information, and (3) information from repeat UK Biobank assessments.

If a death record existed for a participant, then the maximum follow-up date was set to the date
of death. In addition to national death registries, death records could also come from reports
to UK Biobank by relatives in the lost to follow-up reason and date fields #190 and #191. In most
cases, these participants had no death record in the national death registry (e.g. one possibility
is that they may have died abroad), but in cases where a death record existed in national death
registries (N=1, death reported by relative in 2014, death record in national registries in 2021)
we took the earliest date of death as the most conservative approach (e.g. assuming linkage error
and that relative reported death is correct).

If there was a record of the participant being lost to follow-up (e.g. due to leaving the UK) in
fields #190 and #191 then this date was used as the maximum follow-up date in 'latest_hospital_date'.

The exception to this was the lost to follow-up reason "Participant has withdrawn consent for 
future linkage". Although these came with an associated date, we assumed consent was withdrawn at
that date for all linkage, not just for subsequent future linkage, and all fields curated from
linkage were set to have missing data.

The 'latest_hospital_nation' was also harmonized against repeat assessment information. E.g. if a 
person was in England at baseline assessment in 2009, their most recent hospital record was in 
Wales in 2011, then they subsequently attended a repeat assessment in 2014 in England, then their 
nation of residence at maximum follow-up was interpolated to be in England rather than Wales as 
that was their last known location in the data.

The 'lastest_hospital_nation' was also harmonized against the death registry records. E.g. if
a person died in Scotland, then 'latest_hospital_nation' was updated to reflect this (if it was
not already 'Scotland'). In cases where the death record was from England or Wales we set 
'latest_hospital_nation' to "England". There are several reasons for this. First, since the 
maximum follow-up date has been rolled back to date of death so there is no impact on getting
this wrong on 'latest_hospital_date'. Second, all such discrepant events happened well before
maximum follow-up of any nation's hospital records, so getting this wrong for these people should
have no impact on donwstream analysis controlling for differences in follow-up time between nations.
Third, England is the most populous nation, e.g. ~90% of UK Biobank participants are from England, 
so on balance of probabilities these records are far more likely to be from England than Wales.
Fourth, none of the individuals in question were located in Wales at baseline (or repeat) assessment.

The hypothetical maximum follow-up in death records was also similarly curated. 'latest_mortality_date'
gives the date of death, or date of maximum follow-up in the death registry records (for people with
no death records). For people with no death record, we used the 'latest_hospital_nation' to determine
the death registry source to use for the maximum follow-up date. Unlike the hopsital records, this 
currently does not make much difference as there is only a 3 day difference in maximum follow-up date
between sources.

The maximum follow-up date (and death records) were also harmonised to the lost to follow-up information.
Participants with withdrawn consent for linkage were set to have missing data, and the maximum follow-up
date was rolled-back to date lost to follow-up (field #191) if it existed.

We further curated an 'all_cause_mortality' column, indicating for each person whether a death record
exists from any source (and NA where consent for linkage was withdrawn).

Similarly, the 'any_hospitalisation' column indicates whether there are any records of hospitalisations 
(either prospectively or retrospectively) for each participant (or NA where consent for linkage withdrawn).

Follow-up time is computed in years where whole years correspond to calendar years. That is, if the
follow-up time is a whole number, then the corresponding date will always be the same calendar day
in the N years that follow, e.g. 2009-07-27 + 8 years = 2017-07-27. Note this means that rank ordering
based on the follow-up time can in rare cases be incorrect where two people have similar/identical 
follow-up time but there exists a different number of leap-days in their interval dates, e.g. where
the follow-up period for person 1 has 1 leap day, but the follow-up for person 2 has two leap days, they
may have the same follow-up time, but in actuality the number of days of follow-up differs between the
two people.

