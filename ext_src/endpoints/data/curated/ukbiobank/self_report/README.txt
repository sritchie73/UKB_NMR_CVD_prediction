Self-reported medical history
-------------------------------
Scott Ritchie, 14 December 2021
sr827@medschl.cam.ac.uk
-------------------------------

Data in this folder curates the self-reported medical history from touchscreen
and verbal interview questions taken at each UK Biobank assessment into a combined
format that can be used downstream for arbitrary code/field lookup (e.g. where a
user wants to define prevalent disease based on some combination of self-report 
field codes).

Files:
------------

 - README.txt:             You're reading this now!

 - self_report_codes.txt:  Contains long format table of self-reported events,
                           event codes, interpolated age of event, and interpolated
                           date of event (see details below).

 - code_labels.txt:        Contains text labels for each event code in each UK Biobank
                           field curated here (see details below).

 - field_information.txt:  Contains information about the curated fields, the survey
                           question they correspond to, and how age and date of event
                           were interpolated (see details below).

Details:
------------

Data are currently curated from the following UK Biobank fields:

Touchscreen questions relating to medical history, Category 100044
Medical conditions - Health and medical history - Touchscreen - Assessment Centre
https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=100044

 Field ID                                                                    Survey question description
--------- ----------------------------------------------------------------------------------------------
     6150                                                    Vascular/heart problems diagnosed by doctor
     6152  Blood clot, DVT, bronchitis, emphysema, asthma, rhinitis, eczema, allergy diagnosed by doctor
     2443                                                                   Diabetes diagnosed by doctor
     4041                                                                      Gestational diabetes only
     2986                                          Started insulin within one year diagnosis of diabetes
     2463                                                         Fractured/broken bones in last 5 years
     6151                                                                         Fractured bone site(s)
     3005                                                            Fracture resulting from simple fall
--------- ----------------------------------------------------------------------------------------------

Under this category, fields pertaining to the UK Biobank pilot study are omitted 
as they are captured by other fields (i.e. #10844 is captured by #4041 which is
asked of all UK Biobank participants).

Likewise, fields #2453 and #2473 are also omitted as those answering "yes" are 
followed-up in more detail in verbal interview fields #20001 and #20002 curated
below.

Touchscreen questions relating to claudication and peripheral artery disease, Category 100038
Claudication and peripheral artery disease - Health and medical history - Touchscreen - Assessment Centre
https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=100038


 Field ID                                                                    Survey question description
--------- ----------------------------------------------------------------------------------------------
     4728                                                                            Leg pain on walking
     5452                                                        Leg pain when standing still or sitting
     5463                                                                        Leg pain in calf/calves
     5474                                                       Leg pain when walking uphill or hurrying
     5485                                                                 Leg pain when walking normally
     5496                                            Leg pain when walking ever disappears while walking
     5507                                                             Leg pain on walking : action taken
     5518                                                 Leg pain on walking : effect of standing still
     5529                                        Surgery on leg arteries (other than for varicose veins)
     5540                                                               Surgery/amputation of toe or leg
--------- ----------------------------------------------------------------------------------------------

Verbal interview questions relating to cancers, other serious illness, and medical operations:
https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=100074
https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=100076

 Field ID                                                                    Survey question description
--------- ----------------------------------------------------------------------------------------------
    20002                                                                        Non-cancer illness code 
    20001                                                                                    Cancer code 
    20004                                                                                 Operation code 
--------- ----------------------------------------------------------------------------------------------

Touchscreen interview questions relating to mental health:
https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=100060

 Field ID                                                                    Survey question description
--------- ----------------------------------------------------------------------------------------------
    20126                                                            Bipolar and major depression status
    20122                                                                        Bipolar disorder status
    20127                                                                              Neuroticism score
    20124                                                 Probable recurrent major depression (moderate)
    20125                                                   Probable recurrent major depression (severe)
    20123                                                    Single episode of probable major depression
     1920                                                                                    Mood swings
     1930                                                                                  Miserableness
     1940                                                                                   Irritability
     1950                                                                    Sensitivity / hurt feelings
     1960                                                                                Fed-up feelings
     1970                                                                               Nervous feelings
     1980                                                                     Worrier / anxious feelings
     1990                                                                        Tense / 'highly strung'
     2000                                                             Worry too long after embarrassment
     2010                                                                           Suffer from 'nerves'
     2020                                                                          Loneliness, isolation
     2030                                                                                Guilty feelings
     2040                                                                                    Risk taking
     4526                                                                                      Happiness
     4537                                                                          Work/job satisfaction
     4548                                                                            Health satisfaction
     4559                                                               Family relationship satisfaction
     4570                                                                       Friendships satisfaction
     4581                                                               Financial situation satisfaction
     2050                                                    Frequency of depressed mood in last 2 weeks
     2060                                        Frequency of unenthusiasm / disinterest in last 2 weeks
     2070                                          Frequency of tenseness / restlessness in last 2 weeks
     2080                                              Frequency of tiredness / lethargy in last 2 weeks
     2090                                    Seen doctor (GP) for nerves, anxiety, tension or depression
     2100                                 Seen a psychiatrist for nerves, anxiety, tension or depression
     4598                                                                Ever depressed for a whole week
     4609                                                                   Longest period of depression
     4620                                                                  Number of depression episodes
     4631                                             Ever unenthusiastic/disinterested for a whole week
     5375                                                   Longest period of unenthusiasm / disinterest
     5386                                                Number of unenthusiastic/disinterested episodes
     4642                                                                    Ever manic/hyper for 2 days
     4653                                                 Ever highly irritable/argumentative for 2 days
     6156                                                                           Manic/hyper symptoms
     5663                                                      Length of longest manic/irritable episode
     5674                                                           Severity of manic/irritable episodes
     6145                                           Illness, injury, bereavement, stress in last 2 years
--------- ----------------------------------------------------------------------------------------------


In this directory, the file 'self_report_codes.txt' contains a long-format table containing for each
person ('eid') and assessment visit ('visit_index') their answer(s) to these touchscreen and verbal
interview questions in the form of integer codes. The meanings of these codes are detailed in 
'code_labels.txt' (generally 0 = "No", 1 = "Yes", negative numbers = "Missing data", >1 = specific 
illness or operation, and 99999 = "Unclassifiable"). The labels themselves are stored separately to
(1) save storage space, and (2) because the program enabling users to define custom prevalent disease
endpoints will work on these codes (i.e. so we don't have to parse labels or worry about spelling 
mismatches).

For the verbal interview questions (#20002, #20001, #20004) missing data are curated from fields
#2453 (for field #20001), #2473 (for field #20002), and #2415 and #2844 (for field #20004) where 
those fields had answers "Don't know" or "Prefer not to answer" and no codes were recorded in the
respective verbal interview question fields for that person at that assessment visit.

Where the verbal interview questions were not asked (i.e. due to participant reporting no major 
illnesses or operations), the coding "Not applicable" has been added so that downstream code doesn't
interpret these instances as missing data. This has likewise been done for other touchscreen survey
questions which were conditionally asked based on answer(s) to other touchscreen questions.

The file 'self_report_codes.txt' also contains information about the age and date of the event where
available (or can be imputed).

For the verbal interview fields (#20002, #20001, #20004) age and date information was obtained from
their respective fields: #20009, #20007, and #20011 for age of event for fields #20002, #20001, and 
#20004 respectively; and #20008, #20006, and #20010 for date of event for fields #20002, #20001, and
#20004 respectively.

These are interpolated ages and dates, which are not precise. For each event, participants answered
with their age at event or year of event (e.g. both representing a possible one-year timespan). The
following text copied from the UK Biobank showcase describes this in more detail:

 -  If the participant gave a calendar year, then the best-fit time is their age at the mid-point of 
    that year. For example if the year was given as 1970, and the participant was born on 1 April 
    1950, then their age on 1st July 1970 is 20.25 then the value presented is 1970.5.
 -  If the participant gave their age then the value presented is the fractional year corresponding 
    to the mid-point of that age. For example, if the participant said they were 30 years old then 
    the value is 30.5
 -  Interpolated values before the date of birth were truncated forwards to that time.
 -  Interpolated values after the time of data acquisition were truncated back to that time.

We followed this same process to interpolate ages and dates for touchscreen survey questions where
possible. 

For field #6150, "Vascular/heart problems diagnosed by doctor", we obtained age of event information 
from fields #3894 (Heart Attack), #3627 (Angina), #4056 (Stroke), and #2966 (High Blood Pressure). 
For field #6152, "	Blood clot, DVT, bronchitis, emphysema, asthma, rhinitis, eczema, allergy 
diagnosed by doctor", we obtained age of event information from fields #4012 (Deep Vein Thrombosis), 
#4022 (Pulmonary Embolism), #3922 (Emphysema/Chronic Bronchitis), #3786 (Asthma), and #3671. 
For field #2443, "Diabetes diagnosed by doctor", we obtained age of diabetes diagnosis from field #2976.

In each case, these represented an integer age. Following the above, we therefore took the age
of event as the mid-point of that age (e.g. if the record said the participant was 30 years old 
then we set the age of event to 30.5). We then subtracted this from decimal age at assessment 
(see below) to obtain an interpolated date of event. If the interpolated age or date were greater

In the above, decimal age at assessment was computed from the participant's year of birth 
(field #34) and month of birth (field #52). Day of birth is withheld by UK Biobank for 
privacy reasons, so the day of birth in the given month was interpolated as the 15th 
(i.e. the approximate mid-point for each month). This interpolated age of birth was subtracted
from the date of assessment to get an interpolated decimal age at assessment.

For field #2986, "Started insulin within one year diagnosis of diabetes", the age at event 
was interpolated as the mid-point of 1 year after the age of diabetes diagnosis interpolated
from field #2976 (i.e. age recorded in field #2976 + 0.5 years to represent the mid-point of 
age of diabetes diagnosis [the value in the 'age' column for events under field #2443] + 0.5
years to represent the mid-point of the 1 year after interpolated age of diabetes diagnosis).

For field #2463, "Fractured/broken bones in last 5 years", and follow-up questions #6151,
"Fractured bone site(s)", and #3005, ""Fracture resulting from simple fall", the age and date
of event were interpolated as the mid-point of 5 years prior to sample assessment (i.e. age/
date sample assessment minus 2.5 years).

For fields #2050, #2060, #2070, and #2080, each relating to frequency of depressed mood in the
previous two weeks, the age and date were interpolated as the mid-point of 1 week prior to 
sample assessment (i.e. age/date sample assessment minus 7 days).

For field #6145, "Illness, injury, bereavement, stress in last 2 years", the age and date were
interpolated as the mid-point of 1 year prior to sample assessment (i.e. age/date sample 
assessment minus 1 year). 

In all cases, if the interpolated age or date were greater than that at baseline, the age and 
date at event were rolled back to the date of assessment. This included some of the 
interpolated ages and dates recorded for the verbal interview questions which were recorded
as happening the day after assessment. For the touchscreen questions there were no records 
of events within the first year of life, so no ages/dates were rolled forward.

Note for all fields, it was possible for a participant to answer "Don't know" or "Prefer not
to answer" for age or year of event, in which case the corresponding age and date fields 
contain missing values.

Entries are missing entirely where the question was only quantified in a subset of participants,
e.g. only at certain time point. Participants not asked the question are not included in the
data. The exception to this are fields #20124, #20125, and #20123, relating to different ways
of quantifying major depression, in which case the UK Biobank data only contains information for
"yes" answers - here we imputed "no" to be the answer for the rest of the participants at the
time-point (baseline assessment) at which the question was asked.

Ages and dates of events were impossible to interpolate for several touchscreen questions.
These included field #4041, "Gestational diabetes only", which despite being linked to field
#2976, "Age diabetes diagnosed", on the UK Biobank showcase, had no corresponding entries
in field #2976. Touchscreen questions relating to claudication and peripheral artery disease
(fields #4728, #5452, #5463, #5474, #5485, #5496, #5507, #5518, #5529, and #5540) had no 
linked information regarding age of onset of leg pain or age of surgery, and most of the 
mental health questions. 

A summary of all these details are given in the file 'field_information.txt'.

