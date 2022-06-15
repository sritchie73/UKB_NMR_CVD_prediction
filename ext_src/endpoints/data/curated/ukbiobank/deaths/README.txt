Curated death records
---------------------------------
Scott Ritchie, 15 December 2021
sr827@medschl.cam.ac.uk
---------------------------------

Data in this folder curates the death registry records linked to UK Biobank participants
for downstream custom endpoint curation and analysis.

The file 'deaths.txt' contains information on causes of death. Each entry corresponds to
a single cause of death, represented by an ICD-10 code, with each death record containing
one or more causes of death. The primary cause of death is distinguished from secondary 
and other causes of death. 

Details on the columns available are provided in 'column_info.txt':

          column                                                          description
1:   cause_icd10                                       ICD-10 code for cause of death
2:  cause_number                       Cause number, where 1 = primary cause of death
3: primary_cause   TRUE if the ICD-10 code was recorded as the primary cause of death
4: date_of_death                                                        Date of death
5:  death_source     Either "NHS Digital records of death within England or Wales" or 
                                     "NHS Central Register of deaths within Scotland"

Data have been filtered to events in the curated follow-up data, notably excluding any
records that exist for people who have withdrawn consent for electronic health record 
linkage or otherwise recorded as lost to follow-up by UK Biobank (fields #190 and #191).

For death records where a second death certificate has been issued (e.g. subsequent to
post mortem, see https://biobank.ndph.ox.ac.uk/showcase/ukb/docs/DeathLinkage.pdf for
more details) the causes of death from the first death certificate have been discarded
as the second death certificate typically gives more refined causes of death.

Labels for the ICD-10 codes used by UK Biobank (Data Coding #19) are provided in
ICD10_codes.tsv
