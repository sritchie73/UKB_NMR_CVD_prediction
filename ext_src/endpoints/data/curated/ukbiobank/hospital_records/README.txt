Curated hospital episode statistics
------------------------------------
Scott Ritchie, 16 December 2021
sr827@medschl.cam.ac.uk
------------------------------------

Data in this folder curates the hospital episode statistics linked to UK Biobank
participants for downstream custom endpoint analysis.

The file 'hospital_records.txt' contains information on hospitalisations and 
hospital operations for each participant. Each entry corresponds to a single
"cause", of which there may be multiple for each event. These are coded using a
cominbation of ICD-9 and ICD-10 codes (hospitalisations) and OPCS-3 and OPCS-4
codes (operations). The primary cause of event (or reason for operation) is 
distinguished from secondary causes/reasons, and external causes of 
hospitalisation. Fatal events, where a death is recorded on the same date as 
the hospitalisation or operation are distinguished, as are hospital events 
occuring within fatal episodes (e.g. where a person is admitted, then dies 
several days later).

Details on the columns available are provided in 'column_info.txt':

                 column                                                        description
1:           event_type       Indicates whether the corresponding entry contains diagnosis 
                                                                  codes or operation codes
2:           event_date         Date of hospital admission, hospital episode, or operation
3:          episode_end         Hospital discharge date or end of hospital episode. NA for 
                                                                                operations
4:      hospital_nation                   Nation of hospital (England, Wales, or Scotland)
5:                 code   Diagnosis code or operation code, note multiple rows (codes) may 
                                                                 be present for each event
6:            code_type       Type of code: either ICD9, ICD10, OPCS3, or OPCS4 indicating 
                                                                          source of 'code'
7:           cause_type  Whether the diagnosis code is the primary cause, secondary cause, 
                             or external cause of hospitalisation or whether the operation 
                                     code is the primary or secondary reason for operation
8:          fatal_event  TRUE where the person died on the same day as the event/operation
9:  fatality_in_episode             TRUE where the person died during the hospital episode 
                                                  (between 'event_date' and 'episode_end')

Data have been filtered to events in the curated follow-up data, notably excluding any
records that exist for people who have withdrawn consent for electronic health record
linkage or otherwise recorded as lost to follow-up by UK Biobank (fields #190 and #191).

Details
--------

Nation of hospital is provided as this impacts the number of causes that can be coded 
for any given event:
https://biobank.ndph.ox.ac.uk/showcase/ukb/docs/HospitalEpisodeStatistics.pdf)

Note there is only one primary cause for each hospitalisation or operation event.

The date of event of hospitalisations is taken as the hospital episode start date
recorded in the UK Biobank HESIN table, or where that information is not provided,
as the hospital admission date, or decision to admit to hospital date (if that 
information is provided and earlier than the actual date of admission).

For operations, the event date is taken as the recorded operation date in the HESIN
operations table, or if not provided the event date for the matching recorded in the 
main HESIN table as described above.

For episode end date is taken as the episode end date recorded in the HESIN table,
or where that information is not provided, the hospital discharge date. If neither
are provided, then the episode end is taken as the event date.

Labels for the ICD-10, ICD-9, OPCS-4, and OPCS-3 codes used by UK Biobank (Data Coding
#19, #87, #240, and #259 respectively) are provided in ICD10_codes.tsv, ICD9_codes.tsv,
OPCS4_codes.tsv, and OPCS4_codes.tsv respectively.

Note there are 21 records for codes without any event date or epsiode end. These are 
all operations with code OPCS-4 X99.8, "No procedure performed", and for N=4 of these
participants this is the only linked hospital record they have. 

Note also that maternity inpatient episodes and psychiatric inpatient episodes are not
available for events occurring in Scotland:
https://biobank.ndph.ox.ac.uk/showcase/ukb/docs/HospitalEpisodeStatistics.pdf

Note also that the nation of hospital that an event occurred does not necessarily
correspond to the nation of residence of a participant at UK Biobank assessment.
Participants may have moved before or after assessment, or been admitted in a different
nation e.g. due to accidents while on holiday within the UK.

