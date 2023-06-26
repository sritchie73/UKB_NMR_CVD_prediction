Curate incidence and/or prevalence for a custom endpoint
---------------------------------------------------------
Scott Ritchie, 5th January 2022
sr827@medschl.cam.ac.uk
---------------------------------------------------------

The program 'curate_endpoint.R' curates follow-up time and case status in UK Biobank using a user provided set of ICD 
codes (hospitalisations and deaths), OPCS codes (operations), and UK Biobank self reported medical history field codes 
given in a definition file provided to the --def-file argument.

The program can be run from this directory as follows:

./curate_endpoint.R --def-file <file> --output <directory>

Where <file> is the path to your endpoint definition file, and <directory> is the directory you wish to save the 
curated follow-up time and case status as 'events_and_followup.txt'. The endpoint definition file is also copied to 
this directory.

If --output is not provided, then the 'events_and_followup.txt' file will be saved in the same locations as the
--def-file.

----------------------------------------------------------
Endpoint definition file
----------------------------------------------------------

An endpoint definition file consists of series of pre-defined options and arguments which define the case/control 
status of the endpoint of interest based on ICD, OPSC4 and UK Biobank self-report codes. Example definition files can
be found in sub-directories of the the endpoints/ directory.

For example, a simple endpoint for Coronary Heart Disease would look like:

ICD-10: I20-I25

Here, 'ICD-10' is the pre-defined option, which must be followed by a ':', then the argument is'I20-I25' giving the set 
of code(s) for which to consider an event an endpoint case.

This would identify from the hospital and death records all events with ICD-10 codes within the range of I20 to I25 as 
disease cases for the endpoint.

Codes can be given as ranges, as in I20-I25, single codes, as in I24, or multiple comma-separated codes, e.g. I20, I24,
E40. Disease cases are counted where there is any event with the given codes or sub-codes (e.g. I24 captures events 
with codes I24, I24.0, I24.1, I24.8, or I24.9). Specific sub-codes may also be given, e.g. I20, I24.1, I24.9 (or 
without the "." as I241, I249). 

These codes may be specified separately for incident and prevalent events, e.g. where only one or the other are of 
interest, or where definitions differ for incident and prevalent events:

Incident ICD-10: I20-I25
Prevalent ICD-10: I20, I24

In this case, incident events include hospitalisation and deaths with ICD-10 codes in range of I20 to I25, but for 
prevalent disease cases, only events with ICD-10 codes I20 or I24 (or sub-codes, e.g. I24.8) are counted.

Incident events can be further partitioned into fatal or non-fatal events:

Fatal incident ICD-10: I20-I25
Non-fatal incident ICD-10: I20, I21

Fatal cases are counted using only where the ICD-10 code is listed as a cause of death in the death record data (see 
'data/curated/ukbiobank/deaths/README.txt' for more details), and non-fatal cases are those in the hospital records 
that occur during a non-fatal hospitalisation episode (see 'data/curated/ukbiobank/hospital_records/README.txt' for 
more details)

Finally, the endpoints can be defined using only events where the ICD-10 codes of interest are given as the primary 
cause of hospitalisation (or death):

ICD-10 primary cause only: I21, I22

These can be given in combination, giving the following full set of possible options for defining an endpoint from 
ICD-10 codes:

ICD-10:
Prevalent ICD-10: 
Prevalent ICD-10 primary cause only:
Incident ICD-10:
Incident ICD-10 primary cause only:
Non-fatal incident ICD-10:
Non-fatal incident ICD-10 primary cause only:
Fatal incident ICD-10:
Fatal incident ICD-10 primary cause only:

These options can be given alone, or together, allowing for flexible endpoint definitions. For example, for incident 
coronary heart disease:

Incident ICD-10: I21, I22            # Myocardial infarction
Fatal incident ICD-10: I20-I25       # Fatal CHD

To contribute to a disease case, an event must satisfy at least one of the given conditions. E.g. a fatal event with 
codes I20-I25 or non-fatal event with codes I21 or I22. Note that in this definition, non-fatal events with code I21 or 
I22 contribute to disease cases, even though we specify I20-I25 for fatal events, i.e. the listed options do not 
override each other.

Each of the above options can also be modified with the "Excluding" keyword, enabling one to exclude certain ICD codes,
making it easier to define endpoints which don't cover full ranges:

Incident ICD-10: I60-169                        # Range includes all stroke events plus some extras
Excluding incident ICD-10: 167.1, 167.5, 168.2  # Exclude these event codes (and subcodes therein) from the above range

The excluding keyword can be applied to all options for ICD-10, ICD-9, OPCS-3, and OPCS-4 codes, but not self-report
fields.

The definition file can also include comments starting with `#` which are ignored but can be useful documentation for 
later viewing of definition files.

Options may also be repeated to make the documentation clearer. For example, for incident cardiovascular disease we 
could have the following definition file:

# Incident CVD is essentially incident CHD + Stroke:

# Incident CHD
Incident ICD-10: I21, I22            # Myocardial infarction
Fatal incident ICD-10: I20-I25       # Fatal CHD

# Incident stroke
Incident ICD-10: I63                 # Ischaemic stroke
Incident ICD-10: I61                 # Haemorrhagic stroke
Incident ICD-10: I60                 # Subarachnoid haemorrhage
Incident ICD-10: I64                 # Unclassified stroke
Fatal incident ICD-10: I60-I69, F01  # Fatal cerebrovascular events

Which would be equivalent to a simpler definition file without documentation:

Incident ICD-10: I21, I22, I60, I61, I63, I64
Fatal incident ICD-10: I20-I25, I60-I69, F01

Endpoints may also include hospital operations (e.g. surgeries) which are coded via OPCS-4 codes, with similar 
combinations of options available:

OPCS-4:
Prevalent OPCS-4:
Prevalent OPCS-4 primary cause only:
Incident OPCS-4:
Incident OPCS-4 primary cause only:
Non-fatal incident OPCS-4:
Non-fatal incident OPCS-4 primary cause only:

Note here that there is no option to filter to fatal events: death records include only ICD-10 codes.

When defining prevalent cases, ICD-9 codes and OPCS-3 codes will also need to be given to capture the full range of 
historical hospitalisations and operations. Hospital records up until 1996 are often coded with ICD-9 rather than 
ICD-10 codes. Likewise, hospital operations up until 1988 are often coded with OPCS-3 codes rather than OPCS-4 codes. 
The following options can be used in the definition file to capture these:

Prevalent ICD-9:
Prevalent ICD-9 primary cause only:
Prevalent OPCS-3:
Prevalent OPCS-3 primary cause only:

Note, if you ever need to match to a specific code but not it's subcodes, you can add a "$" to the end of the code. E.g.:

Prevalent ICD-9: 250.0$

Would match codes 250.0, but not 250.00, 250.01 etc.

If you're feeling particularly brave you could substitute in other regular expressions, although please note "." are 
stripped out from the codes, and you will need to double check you're matching the codes you intend.

Prevalent disease events can also be defined using self-reported medical history from the touchscreen survey 
questionairres and verbal interviews with nurse at UK Biobank assessment. These can be specified using one or more of 
the following options:

# Verbal interview with nurse:
Self-report field 20001:  # Cancer-related illnesses
Self-report field 20002:  # Non-cancer related illnesses
Self-report field 20004:  # Medical operations

# Touchscreen survey questions relating to diabetes
Self-report field 2443:   # "Has a doctor ever told you that you have diabetes?"
Self-report field 2986:   # "Did you start insulin within one year of your diagnosis of diabetes?"
Self-report field 4041:   # "Did you only have diabetes during pregnancy?"

# Touchscreen survey questions relating to vascular/heart problems diagnosed by doctor
Self-report field 6150:   # "Has a doctor ever told you that you have had any of the following conditions? 
                          #  (You can select more than one answer)"

# Touchscreen survey questions relating to blood clot, DVT, bronchitis, emphysema, asthma, rhinitis, eczema, allergy 
# diagnosed by doctor
Self-report field 6152:   # "Has a doctor ever told you that you have had any of the following conditions? 
                          #  (You can select more than one answer)"

# Touchscreen survey questions relating to claudication and peripheral artery disease
Self-report field 4728:   # "Do you get a pain in either leg on walking?"
Self-report field 5452:   # "Does this pain ever begin when you are standing still or sitting?"
Self-report field 5463:   # "Do you get this pain in your calf (calves)?"
Self-report field 5474:   # "Do you get pain when you walk uphill or hurry?"
Self-report field 5485:   # "Do you get pain when you walk at an ordinary pace on the level?"
Self-report field 5496:   # "Does the pain you get while walking ever disappear when you continue walking?"
Self-report field 5507:   # "What do you do if you get pain when you are walking?"
Self-report field 5518:   # "What happens to the pain you get while walking if you stand still?"
Self-report field 5529:   # "Have you ever had surgery on the arteries of your legs (other than for varicose veins)?"
Self-report field 5540:   # "Have you ever had surgery to remove any of the following?"

# Touchscreen survey questions relating to bone fractures:
Self-report field 2463:   # "Have you fractured/broken any bones in the last 5 years?" 
Self-report field 3005:   # "Did the fracture result from a simple fall (i.e. from standing height)?"
Self-report field 6151:   # "Which bones did you fracture/break? (You can select more than one answer)"

# Touchscreen survey questions relating to mental health:
Self-report field 20126:  # "Bipolar and major depression status"
Self-report field 20122:  # "Bipolar disorder status"
Self-report field 20127:  # "Neuroticism score"
Self-report field 20124:  # "Probable recurrent major depression (moderate)"
Self-report field 20125:  # "Probable recurrent major depression (severe)"
Self-report field 20123:  # "Single episode of probable major depression"
Self-report field 1920:   # "Mood swings"
Self-report field 1930:   # "Miserableness"
Self-report field 1940:   # "Irritability"
Self-report field 1950:   # "Sensitivity / hurt feelings"
Self-report field 1960:   # "Fed-up feelings"
Self-report field 1970:   # "Nervous feelings"
Self-report field 1980:   # "Worrier / anxious feelings"
Self-report field 1990:   # "Tense / 'highly strung'"
Self-report field 2000:   # "Worry too long after embarrassment"
Self-report field 2010:   # "Suffer from 'nerves'"
Self-report field 2020:   # "Loneliness, isolation"
Self-report field 2030:   # "Guilty feelings"
Self-report field 2040:   # "Risk taking"
Self-report field 4526:   # "Happiness"
Self-report field 4537:   # "Work/job satisfaction"
Self-report field 4548:   # "Health satisfaction"
Self-report field 4559:   # "Family relationship satisfaction"
Self-report field 4570:   # "Friendships satisfaction"
Self-report field 4581:   # "Financial situation satisfaction"
Self-report field 2050:   # "Frequency of depressed mood in last 2 weeks"
Self-report field 2060:   # "Frequency of unenthusiasm / disinterest in last 2 weeks"
Self-report field 2070:   # "Frequency of tenseness / restlessness in last 2 weeks"
Self-report field 2080:   # "Frequency of tiredness / lethargy in last 2 weeks"
Self-report field 2090:   # "Seen doctor (GP) for nerves, anxiety, tension or depression"
Self-report field 2100:   # "Seen a psychiatrist for nerves, anxiety, tension or depression"
Self-report field 4598:   # "Ever depressed for a whole week"
Self-report field 4609:   # "Longest period of depression"
Self-report field 4620:   # "Number of depression episodes"
Self-report field 4631:   # "Ever unenthusiastic/disinterested for a whole week"
Self-report field 5375:   # "Longest period of unenthusiasm / disinterest"
Self-report field 5386:   # "Number of unenthusiastic/disinterested episodes"
Self-report field 4642:   # "Ever manic/hyper for 2 days"
Self-report field 4653:   # "Ever highly irritable/argumentative for 2 days"
Self-report field 6156:   # "Manic/hyper symptoms"
Self-report field 5663:   # "Length of longest manic/irritable episode"
Self-report field 5674:   # "Severity of manic/irritable episodes"
Self-report field 6145:   # "Illness, injury, bereavement, stress in last 2 years"

Codes given must match those used by UK Biobank for the respective field. For example, for history of
other heart disease:

Self-report field 6150: 1 # "Heart attack"
Self-report field 20002: 1066 # heart/cardiac problem
Self-report field 20002: 1075 # heart attack/myocardial infarction

Note that while ICD-10, ICD-9, OPCS-4, and OPCS-3 codes may be given as ranges (e.g. I20-I25) this is not recommended 
for self-report fiels as their codes are not necessarily organised hierarchically like ICD/OPCS. In particular, any
ranges given will be automatically expanded by the program (e.g. 1-10 becomes 1, 2, 3, 4, 5, 6, 7, 8, 9, 10). For self-
report fields codes are matched exactly, rather than treated as the begining of codes (e.g. for self-report data, "90"
would match code "90" but for ICD-9, "90" would match anything in code range 900).

A full list of codes and their labels for each self-report field can be found in 
'data/curated/ukbiobank/self_report/code_labels.txt'

Likewise, a full list of ICD-10, ICD-9, OPCS-4, and OPCS-3 codes and their labels can be found in:

'data/curated/ukbiobank/hospital_records/ICD10_codes.tsv'
'data/curated/ukbiobank/hospital_records/ICD9_codes.tsv'
'data/curated/ukbiobank/hospital_records/OPCS4_codes.tsv'
'data/curated/ukbiobank/hospital_records/OPCS3_codes.tsv'

Although please note in these instances the list of codes and labels gives the full range of possible codes, i.e. there
may be codes that are not actually present in the data (either because they are generally not used or because no such 
events have occurred in UK Biobank participants).

Note some mental health field touchscreen surveys contain a mix of numeric data and codes:

Self-report field 20127:  # "Neuroticism score"
Self-report field 4609:   # "Longest period of depression"
Self-report field 4620:   # "Number of depression episodes"
Self-report field 5375:   # "Longest period of unenthusiasm / disinterest"
Self-report field 5386:   # "Number of unenthusiastic/disinterested episodes"

For these, positive numbers indicate a numeric metric, while negative numbers encode various reasons for missing data
as for other self-report fields (e.g. code -1, "Do not know" and code -3 "Prefer not to answer" are treated as missing
data). For these 5 fields, you may additionally specify thresholds using the ">" or "<" symbols e.g.:

Self-report field 20127: >5 

Will set as a prevalent event anyone with a neuroticism score 6 or higher.

Endpoints may also incorporate UK Biobank's algorithmically defined outcomes https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=42:
This field accepts short variable names as follows:

Algorithmically defined outcome: mi              # UK Biobank's algorithmically defined myocardial infarction endpoint (fields 42000 and 42001)
Algorithmically defined outcome: stemi           # UK Biobank's algorithmically defined STEMI endpoint (fields 42002 and 42003)
Algorithmically defined outcome: nstemi          # UK Biobank's algorithmically defined NSTEMI endpoint (fields 42004 and 42005)
Algorithmically defined outcome: stroke          # UK Biobank's algorithmically defined stroke endpoint (fields 42006 and 42007)
Algorithmically defined outcome: stroke_is       # UK Biobank's algorithmically defined ischaemic stroke endpoint (fields 42008 and 42009)
Algorithmically defined outcome: stroke_ih       # UK Biobank's algorithmically defined intracerebral haemorrhage endpoint (fields 42010 and 42011)
Algorithmically defined outcome: stroke_sh       # UK Biobank's algorithmically defined subarachnoid haemorrhage endpoint (fields 42012 and 42013)
Algorithmically defined outcome: asthma          # UK Biobank's algorithmically defined asthma endpoint (fields 42014 and 42015)
Algorithmically defined outcome: copd            # UK Biobank's algorithmically defined COPD endpoint (fields 42016 and 42017)
Algorithmically defined outcome: dementia        # UK Biobank's algorithmically defined all-cause dementia endpoint (fields 42018 and 42019)
Algorithmically defined outcome: alzheimers      # UK Biobank's algorithmically defined alzheimers endpoint (fields 42020 and 42021)
Algorithmically defined outcome: dementia_vasc   # UK Biobank's algorithmically defined vascular dementia endpoint (fields 42022 and 42023)
Algorithmically defined outcome: dementia_ft     # UK Biobank's algorithmically defined frontotemporal dementia endpoint (fields 42024 and 42025)
Algorithmically defined outcome: esrd            # UK Biobank's algorithmically defined end-stage renal disease endpoint (fields 42026 and 42027)
Algorithmically defined outcome: mnd             # UK Biobank's algorithmically defined motor neurone disease endpoint (fields 42028 and 42029)
Algorithmically defined outcome: parkinsonism    # UK Biobank's algorithmically defined all-cause parkinsonism endpoint (fields 42030 and 42031)
Algorithmically defined outcome: parkinsons      # UK Biobank's algorithmically defined parkinsons endpoint (fields 42032 and 42033)
Algorithmically defined outcome: palsy_ps        # UK Biobank's algorithmically defined progressive supranuclear palsy endpoint (fields 42034 and 42035)
Algorithmically defined outcome: atrophy_ms      # UK Biobank's algorithmically defined multiple system atrophy endpoint (fields 42036 and 42037) 

These can be used in endpoint definition files on their own, i.e. simply curating the follow-up time for the algorithimically
defined outcome (only dates for called events are provided by UK Biobank), or combined or augmented with additional codes from above

It is also possible to restrict the follow-up time when defining incident cases:

max follow years:  # Restricts incident disease follow-up to end within a set time period, e.g. 10 years
max follow date:   # Restricts incident disease follow-up to a given date
max follow age:    # Restricts incident disease follow-up to a maximum age at event for each person

When more than one of these options is given, the most restrictive follow-up time will apply to each person.
For example:

max follow years: 10
max follow date: 2020-03-01
max follow age: 80

Will mean that follow-up is truncated for each person at whichever occurs first: 10 years after UK Biobank assessment, 
the 1st March 2020, or the person turns 80 years of age.

Likewise, retrospective follow-up can restricted with similar options:

min follow years: 
min follow date:
min follow age:

However, please note that (1) hospital records are only available retrospectively until the early 90s (England and 
Wales) or 80s (Scotland) while self-report data goes to birth for each person, and (2) age / date of event is sometimes
not available (could not be inferred by UK Biobank team) for self-report events (see 
'data/curated/ukbiobank/self_report/README.txt' for more details). The output provides separate columns for prevalent 
case/control status with sufficient information for time-from-event ('prevalent_event_with_followup') and prevalent 
case/control status regardless with sufficient information for prevalent case exclusion ('prevalent_event') (see output
section below for more details).

Where the earliest prevalent event date is desired, instead of the most recent, then the following option can be set:

earliest prevalent occurrence: true 

Note in these cases the retrospective follow-up time won't be curated.

Further, one may also set:

prevalent date overrules missing: true 

In which case the most recent/earliest prevalent event with known date will be returned, even when there are missing event
dates in the self-report data requested in the endpoint. 

It is also possible to restrict incident disease events to be considered cases only when there are no prevalent disease
events (with or without follow-up time available) by setting:

Incident only where not prevalent: true

----------------------------------------------
Program output
----------------------------------------------

The output is a tab-separated file named 'events_and_followup.txt' containing the following columns:

Information on UK Biobank participants:

 - 'eid':                            Participant identifier for project 7439
 - 'sex':                            Participant sex (Male or Female); self-reported at baseline assessment.
 - 'visit_index':                    Integer coding the respective assessment visit for the entry (row): 
												              - '0' for baseline assessment (first visit; 2006-2010; N=502,460 participants). 
												              - '1' for first repeat assessment (second visit; 2009-2014; N=20,344 
                                            participants).
												              - '2' for imaging assessment (second/third visit; 2014-2020; N=48,998 
                                            participants, of which 8,284 had measurements at first repeat assessment).
												              - '3' for repeat imaging assessment (third/fourth visit; 2019-; N=4,499 
                                            participants).
 - 'assessment_date':                Date at assessment
 - 'assessment_centre':              Location of assessment centre as coded by UK Biobank Field 54. One of:
                                       "Barts", "Birmingham", "Bristol", "Bristol (imaging)", "Bury", "Cardiff", 
                                       "Cheadle (imaging)", "Cheadle (revisit)", "Croydon", "Edinburgh", "Glasgow", 
                                       "Hounslow", "Leeds", "Liverpool", "Manchester", "Middlesborough", "Newcastle", 
                                       "Newcastle (imaging)", "Nottingham", "Oxford", "Reading", "Reading (imaging)", 
                                       "Sheffield", "Stockport (pilot)", "Stoke", "Swansea", or "Wrexham".
                                     Note some locations are listed separately depending on assessment visit, e.g.
                                     "Bristol" vs. "Bristol (imaging)".
 - 'assessment_nation':              Nation assessment centre is located in: "England", "Wales", or "Scotland".
 - 'age':                            Age in years (decimal) at assessment, computed from birth year and birth month
                                     using the 15th (approximate mid-point of each month) as the birth date (specific
                                     date of birth not provided by UK Biobank for privacy reasons). This provides a 
                                     roughly approximate fractional estimate of age, with maximum error of <0.05 years.
                                     Note age here is rounded at 2 decimal places to reflect this uncertainty.
 - 'ehr_linkage_withdrawn':          TRUE or FALSE: TRUE where the participant has withdrawn consent for electronic
                                     health record linkage (UK Biobank field 190). These individuals will have missing
                                     data in their incident and prevalent disease fields where defined from electronic
                                     health records (i.e. data can only be non missing for self-reported prevalent 
                                     disease events).

Prevalent disease case information:

 - 'earliest_hospital_date':         Earliest date in which it is hypothetically possible to find any hospital records 
                                     for this participant, within any restrictions on minimum follow-up date / time / 
                                     age requested by the user. Availability of retrospective linkage to hospital 
                                     records is dependent on the nation of the reporting hospital. As of the 15th 
                                     December 2020 data release, records go back to 27th July 1993 for hospitals in 
                                     England, 2nd December 1980 for hospitals in Scotland, and 18th April 1991 for 
                                     hospitals in Wales. Here, the date reflects the earliest possible date based on 
                                     the historical most likely nation of residence, stored in 'earliest_hospital_nation'. 
                                     or estimated midpoint in the year of arrival (e.g. midpoint of the year, or midpoint
                                     to earliest hospital record, or from start of hospital record follow-up, if they 
                                     occurred in the same year of arrival). For more details see: 
                                     'data/curated/ukbiobank/followup/README.txt'.
 - 'earliest_hospital_nation':       Historical nation of most likely residence, used to infer the hypothetical date at
                                     which the earliest hospital records could potentially be found for this person 
                                     ('earliest_hospital_date'). This is determined based on the nation of the 
                                     reporting hospital for the earliest available hospital records for each person, or
                                     nation at baseline UK Biobank assessment ('assessment_nation' where 'visit_index' 
                                     == 0) where a person had no hospital records (e.g. for people with no 
                                     hospitalisations). This field contains "Immigrated" where the participant immigrated
                                     to the UK after the start of hospital record follow-up. In 97.1% of cases this field 
                                     is identical to nation at UK Biobank assessment ('assessment_nation'). See also
                                     'data/curated/ukbiobank/followup/README.txt' for additional details.
 - 'prevalent_event':                TRUE, FALSE, or NA. TRUE where the person had one or more event of the prevalent 
                                     ICD-10, ICD-9, OPCS-4, OPCS-3, or self-report codes in the provided endpoint 
                                     definition file (--def-file) prior to the corresponding UK Biobank assessment 
                                     visit (e.g. events occuring between baseline assessment and first repeat 
                                     assessment will be listed as prevalent where 'visit_index' > 0), provided that 
                                     event occurred within any restrictions on minimum follow-up time/date/age (or 
                                     where date/age of event was missing). Missing data ('NA') are given where a 
                                     participant answered 'Do not know' or 'Prefer not to answer' to any of the 
                                     self-report fields listed in the endpoint definition file (--def-file), or where 
                                     'ehr_linkage_withdrawn' is TRUE and the prevalent disease definition includes 
                                     ICD-10, ICD-9, OPCS-4, or OPCS-3 codes and there are no prevalent events in the 
                                     self-report data. FALSE is given only where there were definitively no events for 
                                     any of the requested ICD-10, ICD-9, OPCS-4, OPCS-3, or self-report codes in the 
                                     data.
 - 'prevalent_event_with_followup':  TRUE, FALSE, or NA. TRUE or FALSE where prevalent disease cases (or lack thereof) 
                                     from 'prevalent_event' can be used for time-from-event analysis. NA are given 
                                     where (1) there was a prevalent event in the self-report data but the date of 
                                     event was unknown (could not be inferred by the UK Biobank team; for details see 
                                     'data/curated/ukbiobank/self_report/README.txt'), or (2) the prevalent disease 
                                     definition included hospital records (ICD-10, ICD-9, OPCS-4, or OPCS-3 codes) but
                                     the participant had withdrawn consent for linkage ('ehr_linkage_withdrawn'), e.g.
                                     even where there is a self-reported event there may have been one more recent in
                                     the hospital records we can't know about so it's not possible to do 
                                     time-from-event analysis. Further, if self-report codes are requested in 
                                     conjunction with ICD-10, ICD-9, OPCS-4, or OPCS-3 codes for prevalent case
                                     definition, then the retrospective follow-up time is constrained by the earliest
                                     date in the hospital records ('earliest_hospital_date'), i.e. prevalent events
                                     occuring only in the self-report data prior to these dates do not contribute as 
                                     they occur prior to the minimum follow-up time, and this column will contain FALSE 
                                     (while 'prevalent_event' will be TRUE).
 - 'prevalent_event_followup':       Time in years since the most recent prevalent disease event, or for people without 
                                     prevalent disease cases, the time in years from the minimum possible follow-up
                                     date. The minimum possible follow-up date is the 'earliest_hospital_date' if 
                                     prevalent cases include hospitalisations or operations (ICD-10, ICD-9, OPCS-4, or 
                                     OPCS-3 codes). If prevalent cases include only self-report codes, then the minimum
                                     possible follow-up time is to their (inferred approximate) date of birth (i.e. -1 
                                     times 'age'). Missing data (NA) are given where the date of prevalent disease 
                                     event was not available (could not be inferred by the UK Biobank team; see 
                                     'data/curated/ukbiobank/self_report/README.txt') or where case status could not be 
                                     determined (i.e. where a person answered 'Do not know' or 'Prefer not to answer' 
                                     to a self-report question; NA also in 'prevalent_event'). Note that for 
                                     self-reported events, the date and age at event are approximate (inferred by the
                                     UK Biobank team based on the self-reported age or year of event, i.e. at best
                                     accurate to within half a year; 'data/curated/ukbiobank/self_report/README.txt'),
                                     and have thus been rounded to two decimal places.
 - 'prevalent_event_followup_date':  Date corresponding to 'prevalent_event_followup', i.e. the date of the most recent 
                                     prevalent disease case, or for people without prevalent disease cases the minimum 
                                     possible follow-up date (and NA where 'prevalent_event_followup' is NA). Note that
                                     for self-reported events, the date and age at event are approximate (inferred by
                                     the UK Biobank team based on the self-reported age or year of event, i.e. at best
                                     accurate to within half a year; 'data/curated/ukbiobank/self_report/README.txt').
 - 'prevalent_event_followup_age':   Age of the person at the 'prevalent_event_followup_date', computed as 'age' +
                                     'prevalent_event_followup'. This has been rounded to two decimal places as age 
                                     at assessment is approximate: as noted above in the 'age' field, UK Biobank do not 
                                     make specific date of birth available for privacy reasons, so the fractional age
                                     is inferred based on birth month and accurate up to 0.05 of a year. Further, note 
                                     that for self-reported events the age at event is approximate (inferred by the UK 
                                     Biobank team based on the self-reported age or year of event, i.e. at best 
                                     accurate to within half a year; 'data/curated/ukbiobank/self_report/README.txt').
 - 'prevalent_event_date':           Date of most recent prevalent event case. NA where the person had no prevalent
                                     disease cases or where the date of event could not be inferred. Note this date may
                                     be earlier than 'prevalent_event_followup_date' where the most recent event 
                                     occured prior to the earliest date at which time-from-event analysis can be 
                                     conducted (i.e. where prevalent case definition includes both self-report codes 
                                     and hospitalisation/operation codes and the most recent event was in the 
                                     self-report data prior to the earliest date available in the hospital records).
                                     Note that for self-reported events, the date and age at event are approximate 
                                     (inferred by the UK Biobank team based on the self-reported age or year of event, 
                                     i.e. at best accurate to within half a year; for further details see 
                                     'data/curated/ukbiobank/self_report/README.txt').
 - 'prevalent_event_age':            Age at most recent prevalent event, computed from the 'age' at assessment, 
                                     'assessment_date' and 'prevalent_event_date'. Note that 'prevalent_event_age' may
                                     be earlier than 'prevalent_event_followup_age' where the most recent event occured
                                     prior to the earliest date at which time-from-event analysis can be conducted. As
                                     above, this field has been rounded to 2 decimal places to reflect approximate 
                                     nature of fractional age estimates.
 - 'prevalent_has_missing_dates'     Logical; if the option "prevalent date overrules missing: true" was set in the 
                                     endpoint definition file, then this column contains TRUE where the participant 
                                     also had prevalent events in the requested self-report fields that did not have
                                     an inferred event date, i.e. where its not possible to rule out that those events
                                     occurred more recently than the dated one (or earlier in the case the option 
                                     "earliest prevalent occurrence: true" was set in the endpoint definition file).
 - 'prevalent_event_type':           Type of the most-recent prevalent event case, determined by the data source. One
                                     of "hospitalisation", "operation", or "self-reported". In instances where multiple
                                     event types happened on this same date,the most severe / reliable event is 
                                     reported, e.g. if a hospitalisation occurred this is reported, the operation is
                                     reported if no hospitalisation occurred, and the self-reported event is reported
                                     only if there was no hospitalisation or operation on the date of most recent 
                                     event.
 - 'prevalent_record_source':        Record source for the most-recent prevalent event case. Where 
                                     'prevalent_event_type' is "hospitalisation" or "operation" this will be one of 
                                     "NHS records from hospitals in England", "NHS records from hospitals in Scotland", 
                                     or "NHS records from hospitals in Wales". Where 'prevalent_event_type' is 
                                     "self-reported" this will either be "Verbal interview with nurse" or 
                                     "Touchscreen survey". If the date of most-recent prevalent event included only
                                     self-reported events, but from touchscreen survey and verbal interview, then
                                     the verbal interview with nurse is reported.
 - 'prevalent_cause_type':           One of "primary", "secondary", "external", or "self-reported". For 
                                     hospitalisations, this indicates whether the matched ICD-10 or ICD-9 code was the 
                                     primary cause of hospitalisation, one of the secondary causes of hospitalisation, 
                                     or an external cause (e.g. due to accidental injury). For medical operations, 
                                     this indicates whether the matched OPCS-4 or OPCS-3 code was recorded as the 
                                     primary reason or one of the secondary reasons for medical operations.
 - 'prevalent_code_type':            Gives the type of code for the most-recent prevalent event case. One of "ICD-10", 
                                     "ICD-9", "OPCS-4", "OPCS-3" or "self-report field" followed by the UK Biobank 
                                     field number. Self report fields are ordered as listed above in the endpoint 
                                     definition file documentation, i.e. if the most-recent prevalent event includes
                                     codes matching across multiple self-report fields, the first field with a matching
                                     code is listed (i.e. verbal interview with nurse fields take priority over 
                                     touchscreen survey results, and diseases take priority over other types of medical
                                     history).
 - 'prevalent_code':                 Gives the specific code for the most-recent prevalent event case, e.g. "I20.9". 
                                     Where multiple codes in the endpoint definition file are recorded on the date of 
                                     most-recent prevalent event case, the most relevant code is reported. For example, 
                                     if the code is the primary cause of hospitalisation or this code is listed in 
                                     favour of any other matching codes listed as secondary (or external) cause(s). 
                                     Since codes in the hospital records are ordered as recorded, if the 
                                     'prevalent_cause_type' (above) is "secondary" then the code is the most relevant 
                                     (highest in the ordered list of causes) among the matching codes. In contrast 
                                     self-report fields with multiple codes (e.g. multiple choice touchscreen survey
                                     questions, or the disease lists in the verbal interview fields) are unordered, 
                                     thus the first matching code (by numeric order) is listed. 
 - 'prevalent_code_label':           Gives the label for the respective 'prevalent_code', e.g. "Angina pectoris, 
                                     unspecified" for 'prevalent_code' == "I20.9". 

Incident disease cases and event information:

 - 'lost_to_followup':               TRUE or FALSE. TRUE where the participant has been lost to follow-up at some point
                                     since UK Biobank assessment (see UK Biobank fields 190 and 191).
 - 'lost_to_followup_reason':        Reason lost to follow-up, as reported in UK Biobank field 190. One of "Death 
                                     reported to UK Biobank by a relative", "NHS records indicate they are lost to 
                                     follow-up", "NHS records indicate they have left the UK", "UK Biobank sources 
                                     report they have left the UK", or "Participant has withdrawn consent for future 
                                     linkage" (corresponding to 'ehr_linkage_withdrawn').
 - 'lost_to_followup_date':          Date participant is recorded as being lost to follow-up in UK Biobank field 191,
                                     with the exception of participants with 'ehr_linkage_withdrawn'.
 - 'latest_hospital_date':           Latest date in which it is hypothetically possible to find any hospital records
                                     for this participant, given any restrictions on maximum fullow-up time / date /age
                                     specified by the user. Availability of prospective linkage to hospital records is
                                     dependent on the nation of the reporting hospital. As of the 15th December 2020
                                     data release, records are available until 31st March 2021 for hospitals in 
                                     England, 7th May 2021 for hospitals in Scotland, and 6th March 2018 for hospitals 
                                     in Wales. Here, the latest hospital date reflects the maximum follow-up dates in
                                     hospitals in each nation listed above based on the person's most likely nation
                                     of residence stored in 'latest_hospital_nation', provided they had not died, were 
                                     otherwise lost to follow-up, or other requested restrictions on follow-up time or
                                     date or age occured prior to the maximum follow-up date available in hospitals in
                                     the given nation. See 'data/curated/ukbiobank/followup/README.txt' for further 
                                     details.
 - 'latest_hospital_nation':         Nation of most likely residence at end of available follow-up used to infer the
                                     hypothetical date at which the most recent hospital records might be found for 
                                     this person. This is determined based on (1) the nation of reporting hospital for
                                     the most recent hospital records for that person, (2) the nation from which the 
                                     death record was issued for that person, and (3) nation of assessment centre for
                                     most recent UK Biobank assessment. In 99.0% of cases this field is identical to
                                     nation at UK Biobank assessment ('assessment_nation'). For further information
                                     see: 'data/curated/ukbiobank/followup/README.txt'
 - 'latest_mortality_date':          Latest date in which it is hypothetically possible to find any death records for
                                     this participant, given any restrictions on maximum fullow-up time / date /age
                                     specified by the user. Availability of prospective linkage to death records is
                                     dependent on the nation the death is reported in. As of the 15th December 2020 
                                     data release, records are available until the 20th March 2021 for deaths reported
                                     in England or Wales (the source of death records is the same for both nations) and
                                     23rd March 2021 for deaths reported in Scotland. For further details see:
                                     'data/curated/ukbiobank/followup/README.txt'. Here, the person's latest mortality
                                     date reflects either (1) the date of death, (2) the date lost to follow-up or (3) 
                                     the latest date in the death records we might find the death for this person if it 
                                     has occurred (and occurred in England, Wales, or Scotland) subject to any user 
                                     specified restrictions on maximum follow-up. The 'latest_hospital_nation' is used
                                     to determine the nation in which maximum death record follow-up is applicable for
                                     each person (and has been harmonised as such; for further details see:
                                     'data/curated/ukbiobank/followup/README.txt')
 - 'incident_event':                 TRUE, FALSE, or NA. TRUE where there were any hospitalisations, operations, or
                                     deaths matching the ICD-10 or OPCS-4 codes in the endpoint definition file. 
                                     Missing data (NA) correspond to people with withdrawn consent for record linkage
                                     ('ehr_linkage_withdrawn'), or people with any prevalent events (TRUE in the
                                     'prevalent_event' column) if "Incident only where not prevalent: true" is set in
                                     the endpoind definition file (--def-file).
 - 'incident_event_followup':        Time in years to the first occuring event of interest, or for people without
                                     incident disease cases, the time in years to the maximum possible follow-up
                                     date. The maximum possible follow-up date is determined from the latest date 
                                     records could hypothetically be found for each person. This is typically the 
                                     'latest_hospital_date' (see above for details), unless the endpoint definition 
                                     includes a specific set of criteria for fatal events, in which case the maximum 
                                     follow-up is currently slightly shorter in England (20th March 2021) and Scotland 
                                     (23rd March 2021), or several years longer in Wales (20th March 2021). Data are 
                                     missing (NA) where 'incident_event' is missing. Follow-up is further truncated
                                     where the user requests a maximum follow-up date / time / age.
 - 'incident_event_followup_date':   Date of first occurring event of interest or for people without incident disease,
                                     the date from which the maximum follow up time ('incident_event_followup') was 
                                     derived.
 - 'incident_event_followup_age':    Age of the person at the 'incident_event_followup_date', computed as 'age' +
                                     'incident_event_followup'. This has been rounded to two decimal places as age 
                                     at assessment is approximate: as noted above in the 'age' field, UK Biobank do not 
                                     make specific date of birth available for privacy reasons.
 - 'mortality_at_followup_date':     TRUE, FALSE, or NA. TRUE where the person died from any cause at the maximum 
                                     follow-up date ('incident_event_followup_date'), FALSE if the person was still
                                     alive (or no death record existed in England, Wales, or Scotland and had not 
                                     otherwise been reported to the NHS or UK Biobank) at the maximum follow-up date.
                                     Data are missing (NA) where the person had withdrawn consent for record linkage
                                     ('ehr_linkage_withdrawn). 
 - 'incident_event_type':            Type of the first incident event case, determined by the data source. One
                                     of "death", "hospitalisation", or "operation". In instances where multiple
                                     event types happened on this same date, the most severe  event is reported, e.g.
                                     if a death record exists on the same date as a hospital event, the death is 
                                     reported, and the hospitalisation is reported if there is also matching OPCS-4 
                                     codes for a medical operation on the same date the hospitalisation (or death)
                                     occurred. Note that in some instances the reported event may be a hospitalisation
                                     even where a death occurred on the same date ('mortality_at_followup_date'). This
                                     happens when the codes listed in the endpoint file are not listed as any of the 
                                     causes of death on the person's death record, even though they had a hospital 
                                     event with a matching code occurring on the same date.
 - 'incident_record_source':         Record source for the first incident event case. Where 'incident_event_type' is 
                                     "death" this will be either "NHS Central Register of deaths within Scotland" or
                                     "NHS Digital records of death within England or Wales". For hospitalisations or
                                     operations, this will be one of  "NHS records from hospitals in England", "NHS 
                                     records from hospitals in Scotland", or "NHS records from hospitals in Wales".
 - 'incident_cause_type':            One of "primary", "secondary", or "external". For deaths this indicates that the 
                                     matched ICD-10 code was the primary cause of death, or one of the secondary causes
                                     of death. For hospitalisations, this indicates whether the matched ICD-10 code was 
                                     the primary cause of hospitalisation, one of the secondary causes of 
                                     hospitalisation, or an external cause (e.g. due to accidental injury). For medical 
                                     operations, this indicates whether the matched OPCS-4 code was recorded as the
                                     primary reason or one of the secondary reasons for medical operations.
 - 'incident_code_type':             Gives the type of code for the first incident event case. This may be either 
                                     "ICD-10" for hospitalisations or deaths, or "OPCS-4" for medical operations.
 - 'incident_code':                  Gives the specific code for the first incident event case, e.g. "I20.9". Where 
                                     multiple codes in the endpoint definition file are recorded on the date of event,
                                     the most relevant code is reported. For example, if the code is the primary cause
                                     of death (or hospitalisation etc.) this will be reported instead of any codes 
                                     matching secondary or external causes. Since codes in the hospital and death 
                                     records are ordered, if the 'incident_cause_type' (above) is "secondary" then the
                                     code is the one highest in the list of recorded events (e.g. the "second" cause 
                                     of death/hospitalisation will be reported over the "third" and "forth" and so on
                                     if multiple codes are matched).
 - 'incident_code_label':            Gives the label for the respective 'incident_code', e.g. "Angina pectoris,
                                     unspecified" for 'incident_code' == "I20.9".

--------------------------------
Data limitations
--------------------------------

Note that hospital records pertaining to maternity inpatient episodes and psychiatric inpatient episodes are not 
present in the data for events occurring in Scotland:
https://biobank.ndph.ox.ac.uk/showcase/ukb/docs/HospitalEpisodeStatistics.pdf

Note that the number of causes for any hospital event is limited to different numbers per nation, thus secondary events
may be harded to identify depending on the nation if they tend to be used in conjunction with many others:
https://biobank.ndph.ox.ac.uk/showcase/ukb/docs/HospitalEpisodeStatistics.pdf

Note also that the nation of hospital that an event occurred does not necessarily correspond to the nation of residence 
of a participant at UK Biobank assessment. Participants may have moved before or after assessment, or been admitted in 
a different nation e.g. due to accidents while on holiday within the UK.

Note that follow-up time can dramatically differ depending on nation of assessment, which can impact downstream 
modelling (e.g. people in Wales might appear "healthier" as an artefact of the shorter follow-up time available).

For death records where a second death certificate has been issued (e.g. subsequent to post-mortem, see 
https://biobank.ndph.ox.ac.uk/showcase/ukb/docs/DeathLinkage.pdf for more details) the causes of death from the first 
death certificate have been discarded as the second death certificate typically gives more refined causes of death.

--------------------------------
Other notes
--------------------------------

For further details on input data sources, please see:

data/curated/ukbiobank/followup/README.txt
data/curated/ukbiobank/deaths/README.txt
data/curated/ukbiobank/hospital_records/README.txt
data/curated/ukbiobank/self_report/README.txt

For curated information relating to medication usage, family history of disease, biomarker data, etc. please go up
one directory.

For prevalent and incident diabetes, please go up one directory and see the Eastwood_diabetes/ folder for diabetes
defined using the Eastwood et al. 2016 algorithms.

