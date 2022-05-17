QDiabetes risk scores and associated data
---------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
Scott Ritchie, 2nd February 2022
sr827@medschl.cam.ac.uk
------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------

This folder contains information about type 2 diabetes risk in UK Biobank participants computed using the QDiabetes risk
scores (Hippisley-Cox and Coupland. Development and validation of QDiabetes-2018 risk prediction algorithm to estimate
future risk of type 2 diabetes: cohort study. BMJ 359, j5019 (2017)) using the QDiabetes R package, along with
associated risk factors, prevalent and incident diabetes status, as well as potential undiagnosed diabetes and
pre-diabetes.

The key output of interest is stored in 'output/qdiabetes.txt' and a summary of the column headers is given in
'output/column_headers.txt'.

This README file contains extended details on the column contents and their derivation.

It is organised into sections by column type:

 - Row identifiers (participant and visit IDs)
 - QDiabetes risk scores
 - QDiabetes risk factors
 - Prevalent diabetes
 - Possible undiagnosed diabetes
 - Possible prediabetes
 - Incident diabetes
 - Follow-up details

With each section contain details on the following columns:

 - Row identifiers (participant and visit IDs)
   - 'eid'
   - 'visit_index'

 - QDiabetes risk scores
   - 'QDiabetes2018A'
   - 'QDiabetes2018B_fasting'
   - 'QDiabetes2018B_non_fasting'
   - 'QDiabetes2018C'
   - 'QDiabetes2013'

 - QDiabetes risk factors
   - 'ethnicity'
   - 'sex'
   - 'age'
   - 'bmi'
   - 'height'
   - 'weight'
   - 'smoking_status'
   - 'townsend'
   - 'family_history_diabetes'
   - 'history_cvd'
   - 'history_gestational_diabetes'
   - 'history_pcos'
   - 'history_learning_difficulties'
   - 'history_bipolar_schizophrenia'
   - 'hypertension_medication'
   - 'lipid_lowering_medication'
   - 'systematic_corticosteroids'
   - 'atypical_antipsychotics'
   - 'hba1c'
   - 'fasting_glucose'
   - 'fasting_time'
   - 'non_fasting_glucose'

 - Prevalent diabetes
   - 'prevalent_diabetes'
   - 'no_history_any_diabetes'
   - 'type_2_diabetes'
   - 'type_1_diabetes'
   - 'uncertain_diabetes'
   - 'diabetes_history_reason'

 - Possible undiagnosed diabetes
   - 'fasting_undiagnosed_diabetes'
   - 'fasting_undiagnosed_diabetes_reason'
   - 'non_fasting_undiagnosed_diabetes'
   - 'non_fasting_undiagnosed_diabetes_reason'

 - Possible prediabetes
   - 'fasting_prediabetes'
   - 'fasting_prediabetes_reason'
   - 'non_fasting_prediabetes'
   - 'non_fasting_prediabetes_reason'

 - Incident diabetes
   - 'incident_type_2_diabetes'
   - 'incident_type_1_diabetes'
   - 'incident_uncertain_diabetes'
   - 'incident_any_diabetes'
   - 'incident_censor_date'
   - 'incident_censor_years'
   - 'incident_censor_age'
   - 'incident_censor_reason'

 - Follow-up details
   - 'assessment_date'
   - 'assessment_centre'
   - 'assessment_nation'
   - 'death_at_censor_date'
   - 'censor_hospital_nation'
   - 'lost_to_followup'
   - 'lost_to_followup_date'
   - 'lost_to_followup_reason'

------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------
Row identifiers (participant and visit IDs)
---------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
The first set of columns provide sample identifying information:

------
'eid'
------
This column gives the ID for the participant under project 7439.

-------------
'visit_index'
-------------
This column denotes the assessment visit for the corresponding table row for that participant. This index matches the
corresponding integer given when extracting fields from the raw UK Biobank data: '0' for baseline assessment (first
visit; 2006-2010; N=502,460 participants), '1' for first repeat assessment (second visit; 2009-2014; N=20,344
participants), '2' for imaging assessment (second/third visit; 2014-2020; N=48,998 participants, of which 8,284 had
measurements at first repeat assessment), and '3' for repeat imaging assessment (third/fourth visit; 2019-; N=4,499
participants). Currently, the QDiabetes output data is only curated for baseline assessment, so visit_index will always
be 0, but should be used when matching to other curated data under the parent directory.

------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------
QDiabetes risk scores
---------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
The QDiabetes 2018 risk prediction algorithm comprises three models: model A, model B, and model C. Model A predicts 10-
year risk of type 2 diabetes using a combination of lifestyle risk factors, anthropometrics, personal and family medical
history, and medication history. Model B uses these same risk factors, but additionally incorporates information from a
fasting plasma glucose measurement. Model C uses a glycated haemoglobin (HbA1c) measurement instead of fasting glucose
to improve prediction over Model A.

Models A and C are provided as columns in the output file 'QDiabetes2018A' and 'QDiabetes2018C' respectively (see below
for further details).

For Model B, two versions are given: 'QDiabetes2018B_fasting' and 'QDiabetes2018B_non_fasting'. This reflects the
underlying reality of UK Biobank: participants were not fasting when blood samples were taken for biomarker
quantification. The 'QDiabetes2018B_non_fasting' column uses this non-fasting glucose measurement. This will
over-estimate diabetes risk, as glucose measurements will be elevated depending on how recent the participant's last
meal was (see the 'fasting_time' column). The 'QDiabetes2018B_fasting' column uses a glucose measurement which has been
adjusted for fasting time based on the findings of Moebus et al. Impact of time since last caloric intake on blood
glucose levels. Eur. J. Epidemiol. 26, 719–728 (2011) (see below for more details).

The QDiabetes R package can also be used to compute the older QDiabetes risk score from 2013, which uses a subset of the
QDiabetes 2018 model A risk factors. This is provided in the output file as 'QDiabetes2013', with risk computed over a
10-year time horizon to match the time-horizon of the QDiabetes 2018 risk models.

The numeric value in columns 'QDiabetes2018A', 'QDiabetes2018B_fasting', 'QDiabetes2018B_non_fasting', 'QDiabetes2018C',
and 'QDiabetes2013' corresponds to absolute risk (%) of developing (being diagnosed with) type 2 diabetes within 10
years.

Reference:
Hippisley-Cox and Coupland. Development and validation of QDiabetes-2018 risk prediction algorithm to estimate future
risk of type 2 diabetes: cohort study. BMJ 359, j5019 (2017)

----------------
'QDiabetes2018A'
----------------
Absolute risk (%) of developing type 2 diabetes within 10 years predicted by the QDiabetes 2018 model A algorithm.

Predicted absolute risk is missing (N=99,608) where:

  (1) The participant had prevalent diabetes (N=26,445)
      (1) Type 1, Type 2, or Unspecified (but excluding gestational diabetes) see 'prevalent_diabetes' column below
  (2) Prevalent diabetes status was unable to be determined for the participant (N=149)
      (1) No history of any diabetes in self-report touchscreen survey or verbal interview with nurse at UK Biobank
          assessment
      (2) AND, withdrawn consent for linkage to hospital and death records
  (3) The participant was missing data for any of columns used as input to QDiabetes 2018 model A (N=72,836):
      (1)  'ethnicity'                      (N =      0 missing among those without prevalent diabetes)
      (2)  'sex'                            (N =      0 missing among those without prevalent diabetes)
      (3)  'age'                            (N =      0 missing among those without prevalent diabetes)
      (4)  'bmi'                            (N =  2,742 missing among those without prevalent diabetes)
      (5)  'smoking_status'                 (N =  2,638 missing among those without prevalent diabetes)
      (6)  'townsend'                       (N =    585 missing among those without prevalent diabetes)
      (7)  'family_history_diabetes'        (N = 63,206 missing among those without prevalent diabetes)
      (8)  'history_cvd'                    (N =  1,785 missing among those without prevalent diabetes)
      (9)  'history_gestational_diabetes'   (N =      0 missing among those without prevalent diabetes)
      (10) 'history_pcos'                   (N =    571 missing among those without prevalent diabetes)
      (11) 'history_learning_difficulties'  (N =      0 missing among those without prevalent diabetes)
      (12) 'history_bipolar_schizophrenia'  (N =  1,179 missing among those without prevalent diabetes)
      (13) 'hypertension_medication'        (N =  7,716 missing among those without prevalent diabetes)
      (14) 'lipid_lowering_medication'      (N =  7,792 missing among those without prevalent diabetes)
      (15) 'systematic_corticosteroids'     (N =  2,638 missing among those without prevalent diabetes)
      (16) 'atypical_antipsychotics'        (N =  2,638 missing among those without prevalent diabetes)
  (4) The participant's height or weight were outside acceptable ranges for the QDiabetes algorithms (N=178):
      (1) 'weight' <  40 kilograms          (N =     90 among those without missing data or prevalent diabetes)
      (2) 'weight' > 180 kilgorams          (N =     12 among those without missing data or prevalent diabetes)
      (3) 'height' < 1.4 meters             (N =     84 among those without missing data or prevalent diabetes)
      (4) 'height' > 2.1 meters             (N =      0 among those without missing data or prevalent diabetes)

----------------------------
'QDiabetes2018B_non_fasting'
----------------------------
Absolute risk (%) of developing type 2 diabetes within 10 years predicted by the QDiabetes 2018 model B algorithm, here
using non-fasting glucose as proxy for fasting glucose.

Model B additionally incorporates information about (fasting) glucose status in addition to the factors used by model A.
As UK Biobank measurements are non-fasting, we compute two models, one using glucose as-is (non-fasting; the data
captured by this column), and one using glucose adjusted for fasting status/time (see 'QDiabetes2018B_fasting' below).

Note since glucose levels are higher shortly after a meal, the predicted absolute risk in this column will be
overestimated, particularly as 'fasting_time' gets closer to 0.

Additionally, this column contains more missing data due to cut-offs on glucose levels used by the algorithm; risk is
not computed where the input glucose measure > 7 mmol/L as values above this are indicative of type 2 diabetes (when
repeat fasting measurements are taken; see "Possible undiagnosed diabetes" section below).

Here, predicted absolute risk is missing (N=160,116) where:

  (1) The participant had prevalent diabetes (N=26,445)
      (1) Type 1, Type 2, or Unspecified (but excluding gestational diabetes) see 'prevalent_diabetes' column below
  (2) Prevalent diabetes status was unable to be determined for the participant (N=149)
      (1) No history of any diabetes in self-report touchscreen survey or verbal interview with nurse at UK Biobank
          assessment
      (2) AND, withdrawn consent for linkage to hospital and death records
  (3) The participant was missing data for any of columns used as input to QDiabetes 2018 model A (N=127,771):
      (1)  'ethnicity'                        (N =      0 missing among those without prevalent diabetes)
      (2)  'sex'                              (N =      0 missing among those without prevalent diabetes)
      (3)  'age'                              (N =      0 missing among those without prevalent diabetes)
      (4)  'bmi'                              (N =  2,742 missing among those without prevalent diabetes)
      (5)  'smoking_status'                   (N =  2,638 missing among those without prevalent diabetes)
      (6)  'townsend'                         (N =    585 missing among those without prevalent diabetes)
      (7)  'family_history_diabetes'          (N = 63,206 missing among those without prevalent diabetes)
      (8)  'history_cvd'                      (N =  1,785 missing among those without prevalent diabetes)
      (9)  'history_gestational_diabetes'     (N =      0 missing among those without prevalent diabetes)
      (10) 'history_pcos'                     (N =    571 missing among those without prevalent diabetes)
      (11) 'history_learning_difficulties'    (N =      0 missing among those without prevalent diabetes)
      (12) 'history_bipolar_schizophrenia'    (N =  1,179 missing among those without prevalent diabetes)
      (13) 'hypertension_medication'          (N =  7,716 missing among those without prevalent diabetes)
      (14) 'lipid_lowering_medication'        (N =  7,792 missing among those without prevalent diabetes)
      (15) 'systematic_corticosteroids'       (N =  2,638 missing among those without prevalent diabetes)
      (16) 'atypical_antipsychotics'          (N =  2,638 missing among those without prevalent diabetes)
      (17) 'non_fasting_glucose'              (N = 68,876 missing among those without prevalent diabetes)
  (4) The participant's 'non_fasting_glucose' were outside acceptable ranges for fasting glucose for the QDiabetes 2018
      model B algorithm (N = 5,608):
      (1)  'non_fasting_glucose' <  2 mmol/L  (N =      8 among those without missing data or prevalent diabetes)
      (2)  'non_fasting_glucose' >= 7 mmol/L  (N =  5,600 among those without missing data or prevalent diabetes)
  (5) The participant's height or weight were outside acceptable ranges for the QDiabetes algorithms (N=143):
      (1)  'weight' <  40 kilograms           (N =     76 among those not excluded by previous criteria above)
      (2)  'weight' > 180 kilgorams           (N =      5 among those not excluded by previous criteria above)
      (3)  'height' < 1.4 meters              (N =     70 among those not excluded by previous criteria above)
      (4)  'height' > 2.1 meters              (N =      0 among those not excluded by previous criteria above)

------------------------
'QDiabetes2018B_fasting'
------------------------
Absolute risk (%) of developing type 2 diabetes within 10 years predicted by the QDiabetes 2018 model B algorithm, here
using glucose after adjustment for fasting time.

Model B additionally incorporates information about (fasting) glucose status in addition to the factors used by model A.
As UK Biobank measurements are non-fasting, we compute two models, one using glucose as-is (see column
'QDiabetes2018_non_fasting' above), and one using glucose adjusted for fasting status/time (this column).

Here, fasting time was corrected following the findings of Moebus et al. Impact of time since last caloric intake on
blood glucose levels. Eur. J. Epidemiol. 26, 719–728 (2011); that while fasting glucose usually requires a person to
fast for at least 8 hours, glucose levels don't materially change after 3 hours of fasting:

  (1) Median glucose levels were compute for males and females separately at 0 hours fasting, 1 hours fasting, 2 hours
      fasting, and 3 or more hours of fasting.

  (2) Glucose levels were shifted based on the difference in sex-specific medians from those fasting for 3 or more
      hours, i.e. so that the new median glucose would be idetentical across the different fasting time groups within
      each sex (4.903 mmol/L in females, 4.946 mmol/L in males; see below)

This reflects (1) the granularity of the fasting time information, (2) that change in glucose over time is non-linear
within the first 3 hours since last meal, and (3) that the change in glucose levels on average differs between males and
females.

The median glucose and adjustment factors were as follows:

    Fasting time     Sex  Median glucose (mmol/L)  Adjustment factor (mmol/L)
---------------- ------- ------------------------ ---------------------------
 3 hours or more  Female                    4.903                           0
         2 hours  Female                    4.920                      -0.017
         1 hours  Female                    5.111                      -0.208
         0 hours  Female                    5.347                      -0.444
---------------- ------- ------------------------ ---------------------------
 3 hours or more    Male                    4.946                           0
         2 hours    Male                    4.962                      -0.016
         1 hours    Male                    5.088                      -0.142
         0 hours    Male                    5.334                      -0.388
---------------- ------- ------------------------ ---------------------------

Additionally, this column contains more missing data due to cut-offs on glucose levels used by the algorithm; risk is
not computed where the input glucose measure > 7 mmol/L as values above this are indicative of type 2 diabetes (when
repeat fasting measurements are taken; see "Possible undiagnosed diabetes" section below).

Here, predicted absolute risk is missing (N=159,907) where:

  (1) The participant had prevalent diabetes (N=26,445)
      (1) Type 1, Type 2, or Unspecified (but excluding gestational diabetes) see 'prevalent_diabetes' column below
  (2) Prevalent diabetes status was unable to be determined for the participant (N=149)
      (1) No history of any diabetes in self-report touchscreen survey or verbal interview with nurse at UK Biobank
          assessment
      (2) AND, withdrawn consent for linkage to hospital and death records
  (3) The participant was missing data for any of columns used as input to QDiabetes 2018 model A (N=127,771):
      (1)  'ethnicity'                        (N =      0 missing among those without prevalent diabetes)
      (2)  'sex'                              (N =      0 missing among those without prevalent diabetes)
      (3)  'age'                              (N =      0 missing among those without prevalent diabetes)
      (4)  'bmi'                              (N =  2,742 missing among those without prevalent diabetes)
      (5)  'smoking_status'                   (N =  2,638 missing among those without prevalent diabetes)
      (6)  'townsend'                         (N =    585 missing among those without prevalent diabetes)
      (7)  'family_history_diabetes'          (N = 63,206 missing among those without prevalent diabetes)
      (8)  'history_cvd'                      (N =  1,785 missing among those without prevalent diabetes)
      (9)  'history_gestational_diabetes'     (N =      0 missing among those without prevalent diabetes)
      (10) 'history_pcos'                     (N =    571 missing among those without prevalent diabetes)
      (11) 'history_learning_difficulties'    (N =      0 missing among those without prevalent diabetes)
      (12) 'history_bipolar_schizophrenia'    (N =  1,179 missing among those without prevalent diabetes)
      (13) 'hypertension_medication'          (N =  7,716 missing among those without prevalent diabetes)
      (14) 'lipid_lowering_medication'        (N =  7,792 missing among those without prevalent diabetes)
      (15) 'systematic_corticosteroids'       (N =  2,638 missing among those without prevalent diabetes)
      (16) 'atypical_antipsychotics'          (N =  2,638 missing among those without prevalent diabetes)
      (17) 'fasting_glucose'                  (N = 68,887 missing among those without prevalent diabetes)
           (1) Due to missing 'fasting_time'  (N = 11)
  (4) The participant's 'fasting_glucose' were outside acceptable ranges for fasting glucose for the QDiabetes 2018
      model B algorithm (N = 5,608):
      (1)  'fasting_glucose' <  2 mmol/L      (N =      8 among those without missing data or prevalent diabetes)
      (2)  'fasting_glucose' >= 7 mmol/L      (N =  5,380 among those without missing data or prevalent diabetes)
  (5) The participant's height or weight were outside acceptable ranges for the QDiabetes algorithms (N=143):
      (1)  'weight' <  40 kilograms           (N =     76 among those not excluded by previous criteria above)
      (2)  'weight' > 180 kilgorams           (N =      5 among those not excluded by previous criteria above)
      (3)  'height' < 1.4 meters              (N =     70 among those not excluded by previous criteria above)
      (4)  'height' > 2.1 meters              (N =      0 among those not excluded by previous criteria above)

----------------
'QDiabetes2018C'
----------------
Absolute risk (%) of developing type 2 diabetes within 10 years predicted by the QDiabetes 2018 model C algorithm.

Model C additionally incorporates information about glycated haemoglobin (HbA1c) in addition to the factors used by
model A. Unlike glucose, HbA1c is not materially impacted by fasting status.

Additionally, this column contains more missing data due to cut-offs on HbA1c levels used by the algorithm; risk is
not computed where the input HbA1c measure > 48 mmol/mol as values above this are indicative of type 2 diabetes (when
repeat measurements are taken; see "Possible undiagnosed diabetes" section below).

Here, predicted absolute risk is missing (N=127,324) where:

  (1) The participant had prevalent diabetes (N=26,445)
      (1) Type 1, Type 2, or Unspecified (but excluding gestational diabetes) see 'prevalent_diabetes' column below
  (2) Prevalent diabetes status was unable to be determined for the participant (N=149)
      (1) No history of any diabetes in self-report touchscreen survey or verbal interview with nurse at UK Biobank
          assessment
      (2) AND, withdrawn consent for linkage to hospital and death records
  (3) The participant was missing data for any of columns used as input to QDiabetes 2018 model A (N=127,771):
      (1)  'ethnicity'                        (N =      0 missing among those without prevalent diabetes)
      (2)  'sex'                              (N =      0 missing among those without prevalent diabetes)
      (3)  'age'                              (N =      0 missing among those without prevalent diabetes)
      (4)  'bmi'                              (N =  2,742 missing among those without prevalent diabetes)
      (5)  'smoking_status'                   (N =  2,638 missing among those without prevalent diabetes)
      (6)  'townsend'                         (N =    585 missing among those without prevalent diabetes)
      (7)  'family_history_diabetes'          (N = 63,206 missing among those without prevalent diabetes)
      (8)  'history_cvd'                      (N =  1,785 missing among those without prevalent diabetes)
      (9)  'history_gestational_diabetes'     (N =      0 missing among those without prevalent diabetes)
      (10) 'history_pcos'                     (N =    571 missing among those without prevalent diabetes)
      (11) 'history_learning_difficulties'    (N =      0 missing among those without prevalent diabetes)
      (12) 'history_bipolar_schizophrenia'    (N =  1,179 missing among those without prevalent diabetes)
      (13) 'hypertension_medication'          (N =  7,716 missing among those without prevalent diabetes)
      (14) 'lipid_lowering_medication'        (N =  7,792 missing among those without prevalent diabetes)
      (15) 'systematic_corticosteroids'       (N =  2,638 missing among those without prevalent diabetes)
      (16) 'atypical_antipsychotics'          (N =  2,638 missing among those without prevalent diabetes)
      (17) 'hba1c'                            (N = 33,647 missing among those without prevalent diabetes)
  (4) The participant's 'hba1c' were outside acceptable ranges for fasting glucose for the QDiabetes 2018 model C
      algorithm (N = 2,996):
      (1)  'hba1c' <  15 mmol/mol             (N =    167 among those without missing data or prevalent diabetes)
      (2)  'hba1c' >= 48 mmol/mol             (N =  2,829 among those without missing data or prevalent diabetes)
  (5) The participant's height or weight were outside acceptable ranges for the QDiabetes algorithms (N=160):
      (1)  'weight' <  40 kilograms           (N =     81 among those not excluded by previous criteria above)
      (2)  'weight' > 180 kilgorams           (N =      9 among those not excluded by previous criteria above)
      (3)  'height' < 1.4 meters              (N =     77 among those not excluded by previous criteria above)
      (4)  'height' > 2.1 meters              (N =      0 among those not excluded by previous criteria above)

---------------
'QDiabetes2013'
---------------
Absolute risk (%) of developing type 2 diabetes within 10 years predicted by the QDiabetes 2013 algorithm.

The QDiabetes 2013 is an older risk prediction algorithm that uses a subset of the risk factors captured by QDiabetes
2018 model A. While its time-horizon for risk prediction is more flexible (at least in the QDiabetes R package), here
we've computed 10-year risk to match the time-horizon of the QDiabetes 2018 algorithms.

Predicted absolute risk is missing (N=99,316) where:

  (1) The participant had prevalent diabetes (N=26,445)
      (1) Type 1, Type 2, or Unspecified (but excluding gestational diabetes) see 'prevalent_diabetes' column below
  (2) Prevalent diabetes status was unable to be determined for the participant (N=149)
      (1) No history of any diabetes in self-report touchscreen survey or verbal interview with nurse at UK Biobank
          assessment
      (2) AND, withdrawn consent for linkage to hospital and death records
  (3) The participant was missing data for any of columns used as input to QDiabetes 2018 model A (N=72,836):
      (1)  'ethnicity'                      (N =      0 missing among those without prevalent diabetes)
      (2)  'sex'                            (N =      0 missing among those without prevalent diabetes)
      (3)  'age'                            (N =      0 missing among those without prevalent diabetes)
      (4)  'bmi'                            (N =  2,742 missing among those without prevalent diabetes)
      (5)  'smoking_status'                 (N =  2,638 missing among those without prevalent diabetes)
      (6)  'townsend'                       (N =    585 missing among those without prevalent diabetes)
      (7)  'family_history_diabetes'        (N = 63,206 missing among those without prevalent diabetes)
      (8)  'history_cvd'                    (N =  1,785 missing among those without prevalent diabetes)
      (9)  'hypertension_medication'        (N =  7,716 missing among those without prevalent diabetes)
      (10) 'systematic_corticosteroids'     (N =  2,638 missing among those without prevalent diabetes)
  (4) The participant's Townsend deprivation index was outside the acceptable range for the QDiabetes 2013 algorithm:
      (1) 'townsend' > 11                   (N =      1 among those without missing data or prevalent diabetes)
  (4) The participant's height or weight were outside acceptable ranges for the QDiabetes algorithms (N=178):
      (1) 'weight' <  40 kilograms          (N =     90 among those not excluded by previous criteria above)
      (2) 'weight' > 180 kilgorams          (N =     12 among those not excluded by previous criteria above)
      (3) 'height' < 1.4 meters             (N =     84 among those not excluded by previous criteria above)
      (4) 'height' > 2.1 meters             (N =      0 among those not excluded by previous criteria above)

------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------
QDiabetes risk factors
---------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
Columns described in this section are those containing data directly used by the QDiabetes R package and risk models to
compute each participant's 10-year absolute risk of diabetes at the respective UK Biobank assessment visit.

-----------
'ethnicity'
-----------
This column contains the self-reported ethnicity of UK Biobank participants at baseline assessment (UK Biobank field
#21000) mapped to the ethnicity categories used by the QDiabetes risk scores (and those expected by the corresponding R
package).

This mapping is as follows:

    QDiabetes category  R package coding          UK Biobank category  Field #21000 code
---------------------- ----------------- ---------------------------- ------------------
 White or not recorded           WhiteNA                        White                  1
                                                              British               1001
                                                                Irish               1002
                                           Any other white background               1003
                                                          Do not know                 -1
                                                 Prefer not to answer                 -3
                                                                 <NA>               <NA>
                Indian            Indian                       Indian               3001
             Pakistani         Pakistani                    Pakistani               3002
           Bangladeshi       Bangladeshi                  Bangladeshi               3003
           Other Asian        OtherAsian       Asian or Asian British                  3
                                           Any other Asian background               3004
             Caribbean    BlackCaribbean                    Caribbean               4001
         Black African      BlackAfrican                      African               4002
               Chinese           Chinese                      Chinese                  5
                 Other             Other                        Mixed                  2
                                            White and Black Caribbean               2001
                                              White and Black African               2002
                                                      White and Asian               2003
                                           Any other mixed background               2004
                                               Black or Black British                  4
                                           Any other Black background               4003
                                                   Other ethnic group                  6
---------------------- ----------------- ---------------------------- ------------------

-----
'sex'
-----
This column contains self-reported sex at baseline assessment from UK Biobank field #31. Either "Male" or "Female".

-----
'age'
-----
This column contains self-reported age at the assessment visit from UK Biobank field #21003. Note age is in whole
numbers.

-----
'bmi'
-----
This column contains body mass index available in UK Biobank field #21001, derived from standing height measured by UK
Biobank (field #54) and weight recorded by UK Biobank (field #21002). Note data in this column may be missing where
standing height or weight were not recorded (see fields #20047 and #20041 for list of possible reasons, e.g. equipment
failure, participant is an amputee, participant using a wheelchair, participant above max weight, etc.).

Note QDiabetes risk could not be computed for participants with missing BMI (N=3,104).

Note also that the QDiabetes risk algorithms truncate BMI to be between 20 and 40 kg/m^2 for participants within
acceptable weight and height ranges (see column information below). Participants with BMI < 20 kg/m^2 (N=11,738) are set
to have BMI of 20 kg/m^2, and participants with BMI > 40 kg/m^2 (N=9,688) are set to have BMI of 40 kg/m^2.

--------
'height'
--------
This column contains standing height in meters. Derived from UK Biobank field #50, where values were provided in
centimeters.

Note data in this column may be missing where standing height was not recorded (see field #20047 for list of possible
reasons, e.g. equipment failure, participant is an amputee, participant using a wheelchair, etc.) (N=2,538).

Note also that the QDiabetes prediction algorithms are invalid where height < 1.4 meters (N=111) or height > 2.1 meters
(N=0), so people with height outside this range are also missing QDiabetes prediction scores.

--------
'weight'
--------
This column contains weight in kilograms from UK Biobank field #21002.

Note data in this column may be missing where weight was not recorded (see field #20041 for list of possible
reasons, e.g. equipment failure, participant is an amputee, participant using a wheelchair, etc.) (N=2,774).

Note also that the QDiabetes prediction algorithms are invalid where weight < 40 kg (N=119) or weight > 180 kg (N=22),
so people with weight outside this range are also missing QDiabetes prediction scores.

----------------
'smoking_status'
----------------
This column contains details on the participant's smoking status at UK Biobank assessment, mapped to QDiabetes risk
categories.

This mapping is as follows:

   QDiabetes category  R package coding            Current smoking    Cigarettes smoked per day
                                         (UK Biobank field #20116)     (UK Biobank field #3456)
--------------------- ----------------- -------------------------- ----------------------------
           Non-smoker               Non                      Never                         <NA>
        Former smoker                Ex                   Previous                         <NA>
         Light smoker             Light                    Current                          1-9
                                                           Current          Less than one a day
      Moderate smoker          Moderate                    Current                        10-19
                                                           Current                  Do not know
                                                           Current         Prefer not to answer
         Heavy smoker             Heavy                    Current                       20-140
                 <NA>              <NA>       Prefer not to answer                         <NA>
--------------------- ----------------- -------------------------- ----------------------------

Note QDiabetes risk could not be computed for participants with missing smoking status (N=2,948; N=2,057 prefer no to
answer, N=441 with no data available)

----------
'townsend'
----------
This column contains the Townsend Deprivation Index based on the participant's home address at baseline assessment,
obtained from UK Biobank field #189.

Note QDiabetes risk could not be computed for participants with missing Townsend Deprivation Index (N=623), and for the
QDiabetes 2013 risk score, for participants with Townsend Deprivation Index > 11 (N=1).

-------------------------
'family_history_diabetes'
-------------------------
This column curates information about family history of diabetes in first degree relatives (biological father, mother,
or siblings) from UK Biobank fields #20107 (illness of father), #20110 (illness of mother), and #20111 (illness of
siblings).

This column contains 'TRUE' (N=108,580) where any of these three fields included the answer "Diabetes". This column
contains 'FALSE' where all three fields contained "None of the above" as the answer (i.e. definitely no family history
of diabetes).

Data were missing where it was not possible to determine presence OR absence of family history of diabetes (N=67,176).
Data are missing where the three fields contained a combination of "None of the above", "Prefer not to answer", and
"Do not know" answers (N=58,606) or where there were no data on these touchscreen survey questions (N=8,570).

Note QDiabetes risk could not be computed for participants with missing data.

-------------
'history_cvd'
-------------
This column curates information about the participant's history of cardiovascular disease, defined by the QDiabetes
paper as "ischaemic heart disease, stroke, or transient ischaemic attack". Data in this field is curated from a
combination of self-report touchscreen and verbal interview questions as well as retrospective hospital records using
the endpoint definition file 'data/adjudicated_medical_history/CVD/endpoint_definition.txt' in conjuction with the
script 'src/adjudicate_medical_history/curate_endpoint.R' (a symbolic link to the program in '../endpoints/').

Briefly, the column contains 'TRUE' (N=28,966) where:

  (1) The participant is recorded as having a previous hospital episode for reasons relating to:
      (1) Coronary artery disease (ICD-10 codes I21-I24 or ICD-9 codes 410-412)
      (2) OR, chronic ischaemic heart disease (ICD-10 codes I25.1, I25.2, or I25.5-I25.9, or ICD-9 codes 414.0, 414.8,
          or 414.9)
      (3) OR, any type of stroke (ICD-10 codes I60, I61, I63, or I64, or ICD-9 codes 430, 431, 434, or 436)
      (4) OR, coronary surgery (ICD-10 code Z95.1, ICD-9 code V45.81, OPCS-4 codes K40-K46, K49, K50.1, or K75, or
          OPCS-3 codes 309.4 or 884)
  (2) OR, the participant answered with "Heart attack" or "Stroke" on the touchscreen survey question relating to
      "Vascular/heart problems diagnosed by doctor" (UK Biobank field #6150) which asked "Has a doctor ever told you
      that you have had any of the following conditions? (You can select more than one answer)"
  (3) OR, the trained nurse recorded the following after conducting a verbal interview with the participant about
      history of non-cancer related illness (UK Biobank field #20002):
      (1) Heart attack/myocardial infarction (code 1075)
      (2) OR, cardiomyopathy (code 1079)
      (3) OR, stroke (code 1081)
      (4) OR, subarachnoid haemorrhage (code 1086)
      (5) OR, brain haemorrhage (code 1491)
      (6) OR, ischaemic stroke (code 1583)
      (7) OR, transient ischaemic attack (code 1082)
  (4) OR, the trained nurse recorded the following after conducting a verbal interview with the participant about
      history of medical operations (UK Biobank field #20004):
      (1) Coronary angioplasty (ptca) +/- stent (code 1070)
      (2) OR, coronary artery bypass grafts (cabg) (code 1095)
      (3) OR, triple heart bypass (code 1523)

Note that retrospective follow-up in hospital records only goes back so far, depending on the nation of the hospital in
which the hospital episodes occurred. The earliest follow-up dates are currently 1993-07-27 for hospitals in England,
1991-04-18 for hospitals in Wales, and 1980-12-02 for hospitals in Scotland, where baseline UK Biobank assessment
occurred between 2006 and 2010 (and participants were at minimum 37 years old at baseline assessment).

The earliest hypothetical date at which hospital records might be found for each participant, based on inferred nation
of residence at those dates, are coded in the 'earliest_hospital_date' and 'earliest_hospital_nation' columns described
in more detail in the section "Follow-up details" below.

Data were missing where it was not possible to determine presence OR absence of history of CVD (N=1,992), where:

  (1) The person did not satisfy any of the case definitions above,
  (2) AND, they:
      (1) Had withdrawn consent for hospital record linkage (see 'lost_to_follow_up_reason' in section "Follow-up
          details" below)
      (2) OR, had no data in the respective touchscreen survey question (UK Biobank field #6150)
      (3) OR, had no data in the respective nurse interview questions (UK Biobank field #20002 or #20004)
      (4) OR, they:
          (1) Answered "Do not know" or "Prefer not to answer" to the touchscreen survey question (UK Biobank field
              #6150)
          (2) AND, they had no data in the respective follow-up nurse interview question (UK Biobank field #20002), e.g.
              due to a "Prefer not to answer" or "Do no know" with nurse unable to determine any details on medical
              history.

 Note QDiabetes risk could not be computed for participants with missing data.

------------------------------
'history_gestational_diabetes'
------------------------------
This column curates information about the participant's history of gestational diabetes. This data was curated using a
combination of the self-report touchscreen and verbal interview questions as well as retrospective hospital records.

First, gestational diabetes status was curated from self-report touchscreen and verbal interview questions using the
algorithms described in Eastwood et al. Algorithms for the Capture and Adjudication of Prevalent and Incident Diabetes
in UK Biobank. PLoS One 11, e0162388 (2016).

Then, for people classified as "Diabetes Unlikely" from the self-report data, additional gestational diabetes cases were
determined from retrospective hospital records using the endpoint definition file
'data/adjudicated_medical_history/gestational_diabetes/endpoint_definition.txt' in conjuction with the script
'src/adjudicate_medical_history/curate_endpoint.R' (a symbolic link to the program in '../endpoints/').

Briefly, the column contains 'TRUE' (N=844) where:

  (1) They self-reported "Female" sex (see 'sex' field above),
  (2) AND,
      (1) They were classified as a "Possible gestational diabetes" by the Eastwood et al. algorithms:
          (1) They:
              (1) Self-report gestational diabetes in the touchscreen survey (UK Biobank field #4041),
              (2) OR, they self-report gestational diabetes in interview with nurse (UK Biobank field #20002) with age
                  of diagnosis < 50 years age (UK Biobank field #20008)
          (2) AND, they did not self-report concurrent type 1 or type 2 diabetes in interview with nurse (codes 1222 or
              1223 in UK Biobank field #20002)
          (3) AND, they did not self-report any diabetes medication:
              (1) In touchscreen interview questions (UK Biobank fields #6177, #6153, or #2986)
              (2) In interview with nurse (UK Biobank field #20003), see full list of diabetes medications below
      (2) Or, they were classified as "Diabetes Unlikely" or "Uncertain diabetes status" by the Eastwood algorithm
          (1) AND, they had any hospital records prior to UK Biobank assessment for gestational diabetes (ICD-10 code
              O24.4)
          (2) AND, they had no record of other diabetes (type 1, type 2, or uncertain) in the hospital records

In interview with nurse on current medication usage, the following set of medications were considered diabetes
medications for the purposes of the above case classification:

 Medication class  UK Biobank field #20003 code               UK Biobank field #20003 label
----------------- ----------------------------- -------------------------------------------
          Insulin                    1140883066                             Insulin product
        Metformin                    1140884600                                   Metformin
                                     1140874686                     Glucophage 500mg tablet
                                     1141189090  Rosiglitazone 1mg / metformin 500mg tablet
    Sulfonylureas                    1140874718                               Glibenclamide
                                     1140874744                                  Gliclazide
                                     1140874746                       Diamicron 80mg tablet
                                     1141152590                                 Glimepiride
                                     1141156984                           Amaryl 1mg tablet
                                     1140874646                                   Glipizide
                                     1141157284                           Glipizide product
                                     1140874652                       Minodiab 2.5mg tablet
                                     1140874674                                 Tolbutamide
                                     1140874728                       Euglucon 2.5mg tablet
     Meglitinides                    1141173882                                 Nateglinide
                                     1141173786                         Starlix 60mg tablet
                                     1141168660                                 Repaglinide
       Glitazones                    1141171646                                Pioglitazone
                                     1141171652                           Actos 15mg tablet
                                     1141153254                                Troglitazone
                                     1141177600                               Rosiglitazone
                                     1141177606                          Avandia 4mg tablet
            Other                    1140868902                                    Acarbose
                                     1140868908                        Glucobay 50mg tablet
                                     1140857508             Glucotard 5g/sachet mini-tablet
----------------- ----------------------------- -------------------------------------------

Note that retrospective follow-up in hospital records only goes back so far, depending on the nation of the hospital in
which the hospital episodes occurred. The earliest follow-up dates are currently 1993-07-27 for hospitals in England,
1991-04-18 for hospitals in Wales, and 1980-12-02 for hospitals in Scotland, where baseline UK Biobank assessment
occurred between 2006 and 2010 (and participants were at minimum 37 years old at baseline assessment).

The earliest hypothetical date at which hospital records might be found for each participant, based on inferred nation
of residence at those dates, are coded in the 'earliest_hospital_date' and 'earliest_hospital_nation' columns described
in more detail in the section "Follow-up details" below.

Data were missing where it was not possible to determine presence OR absence of history of gestational diabetes
(N=10,114), where:

  (1) They self-reported "Female" sex (see 'sex' field above),
  (2) AND, they were not classified as having possible gestational diabetes above,
  (2) AND, they:
      (1) Were classified as having type 1 diabetes, type 2 diabetes, or uncertain diabetes (either status or type), see
          columns 'type_1_diabetes', 'type_2_diabetes', and 'uncertain_diabetes' in the "Prevalent diabetes" section
          below.
      (2) OR, they had withdrawn consent for hospital record linkage (see 'lost_to_follow_up_reason' in section
          "Follow-up details" below)
      (3) OR, they had any history of unspecified or uncertain diabetes (or complications) during pregnancy in hospital
          records prior to UK Biobank assessment (ICD-10 code O24.3, or ICD-9 codes 648.0 or 648.8).

Note QDiabetes 2018 risk could not be computed for participants with missing data, however, the vast majority of the
missing here data (N=10,017) were females classified as another type of diabetes (either in conjuction with gestational
diabetes, or independently of) which are excluded from the QDiabetes risk score calculations due to presence of
prevalent diabetes.

The 'diabetes_history_reason' contains further details on the adjudication diabetes history for each participant,
including distinguishing between gestational diabetes identified from self-report vs. hospital records. See section
"Prevalent diabetes" below for further details.

For clarification of details on adjudication from self-report data using the Eastwood et al. algorithms, contact Dr.
Sam Lambert (sl925@medschl.cam.ac.uk) who implemented the code in the iPython/Jupyter notebooks in
'data/adjudicated_medical_history/Eastwood_diabetes/' (symbolically links to directory elsewhere on the cluster).

The column always contains 'FALSE' where the participant is self-reported "Male" sex (pregnancy not possible, no
gestational diabetes), even where they were classified as having another type of diabetes.

--------------
'history_pcos'
--------------
This column curates information about the participant's history of polycystic ovary syndrome (PCOS). Data in this field
is curated from a combination of self-report touchscreen and verbal interview questions as well as retrospective
hospital records using the endpoint definition file 'data/adjudicated_medical_history/PCOS/endpoint_definition.txt' in
conjuction with the script 'src/adjudicate_medical_history/curate_endpoint.R' (a symbolic link to the program in
'../endpoints/').

Briefly, the column contains 'TRUE' (N=755) where:

  (1) The participant self-reported "Female" sex (see 'sex' column above)
  (2) AND, they:
      (1) The participant is recorded as having a previous hospital episode for reasons relating to polycystic ovary
          syndrome (ICD-10 code E28.2 or ICD-9 code 256.4)
      (2) OR, The trained nurse recorded "polycystic ovaries/polycystic ovarian syndrome" (code 1350) after conducting a
          verbal interview with the participant about history of non-cancer related illness (UK Biobank field #20002).

Note that retrospective follow-up in hospital records only goes back so far, depending on the nation of the hospital in
which the hospital episodes occurred. The earliest follow-up dates are currently 1993-07-27 for hospitals in England,
1991-04-18 for hospitals in Wales, and 1980-12-02 for hospitals in Scotland, where baseline UK Biobank assessment
occurred between 2006 and 2010 (and participants were at minimum 37 years old at baseline assessment).

The earliest hypothetical date at which hospital records might be found for each participant, based on inferred nation
of residence at those dates, are coded in the 'earliest_hospital_date' and 'earliest_hospital_nation' columns described
in more detail in the section "Follow-up details" below.

Data were missing where it was not possible to determine presence OR absence of history of PCOS (N=621), where:

  (1) The participant self-reported "Female" sex (see 'sex' column above)
  (2) AND they:
      (1) Did not satisfy any of the case definitions above
      (2) AND they:
          (1) Had withdrawn consent for hospital record linkage (see 'lost_to_follow_up_reason' in section "Follow-up
              details" below)
          (2) OR, had no data in the respective nurse interview questions (UK Biobank field #20002), e.g. due to a
              "Prefer not to answer" in the preceeding touchscreen survey question #2473 ("Has a doctor ever told you
              that you have had any other serious medical conditions or disabilities?" or were not interviewed by the
              nurse

Note QDiabetes 2018 risk could not be computed for participants with missing data.

The column always contains 'FALSE' where the participant is self-reported "Male" sex (no ovaries; no PCOS).

-------------------------------
'history_learning_difficulties'
-------------------------------
This column curates any available evidence suggesting a history of learning difficulties/disabilty. Although this is a
risk factor required by the QDiabetes risk score, this data is not something collected by UK Biobank. Here, evidence of
history of learning difficulties/disability is detected via incidental codings in restrospective hospital records using
the endpoint definition file 'data/adjudicated_medical_history/learning_difficulties/endpoint_definition.txt' in
conjuction with the script 'src/adjudicate_medical_history/curate_endpoint.R' (a symbolic link to the program in
'../endpoints/').

Briefly, the column contains 'TRUE' (N=94) where:

  (1) The participant is recorded as having a previous hospital episode which included any coding for:
      (1) "Mental retardation" (ICD-10 codes F70-F79 or ICD-9 codes 317-319)
      (2) OR,  “Specific developmental disorders of scholastic skills” (ICD-10 code F81)
      (3) OR, “Specific delays in development” (ICD-9 code 315)

Note that in addition to hospital records only capturing this data incidentally, the retrospective follow-up in hospital
records only goes back so far, often not covering a participant's childhood, depending on the nation of the hospital in
which the hospital episodes occurred. The earliest follow-up dates are currently 1993-07-27 for hospitals in England,
1991-04-18 for hospitals in Wales, and 1980-12-02 for hospitals in Scotland, where baseline UK Biobank assessment
occurred between 2006 and 2010 (and participants were at minimum 37 years old at baseline assessment).

The earliest hypothetical date at which hospital records might be found for each participant, based on inferred nation
of residence at those dates, are coded in the 'earliest_hospital_date' and 'earliest_hospital_nation' columns described
in more detail in the section "Follow-up details" below.

Data were missing where it was not possible to determine presence OR absence of history of learning difficulties
(N=158), where:

  (1) The participant had withdrawn consent for hospital record linkage (see 'lost_to_follow_up_reason' in section
      "Follow-up details" below)

Note QDiabetes 2018 risk could not be computed for participants with missing data.

-------------------------------
'history_bipolar_schizophrenia'
-------------------------------
This column curates information about the participant's history of bipolar or schizophrenia disorders. Data in this
field is curated from a combination of self-report touchscreen and verbal interview questions as well as retrospective
hospital records using the endpoint definition file
'data/adjudicated_medical_history/bipolar_schizophrenia/endpoint_definition.txt' in conjuction with the script
'src/adjudicate_medical_history/curate_endpoint.R' (a symbolic link to the program in '../endpoints/').

Briefly, the column contains 'TRUE' (N=3,701) where:

  (1) the participant is recorded as having a previous hospital episode for reasons relating to:
      (1) Bipolar disorder (ICD-10 code F31 or ICD-9 codes 296.80 or 296.89)
      (2) OR, Schizophrenia (ICD-10 code F20 or ICD-9 code 295)
      (3) OR, Manic-depressive psychosis (ICD-9 codes 296.0 or 296.4-296.7)
  (2) OR, the trained nurse recorded either of the following after conducting a verbal interview with the participant
      about history of non-cancer related illness (UK Biobank field #20002):
      (1) schizophrenia (code 1289)
      (2) OR, mania/bipolar disorder/manic depression (code 1291)
  (3) OR, the participant was classified as having Bipolar type I or Bipolar type I disorders in UK Biobank fields
      #20122/#20126 (derived from various combinations of mental health and psychosocial factor touchscreen survey
      answers, see https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=158772 for details).

Note that retrospective follow-up in hospital records only goes back so far, depending on the nation of the hospital in
which the hospital episodes occurred. The earliest follow-up dates are currently 1993-07-27 for hospitals in England,
1991-04-18 for hospitals in Wales, and 1980-12-02 for hospitals in Scotland, where baseline UK Biobank assessment
occurred between 2006 and 2010 (and participants were at minimum 37 years old at baseline assessment).

The earliest hypothetical date at which hospital records might be found for each participant, based on inferred nation
of residence at those dates, are coded in the 'earliest_hospital_date' and 'earliest_hospital_nation' columns described
in more detail in the section "Follow-up details" below.

Data were missing where it was not possible to determine presence OR absence of history of bipolar or schizophrenia
disorders (N=1,369), where:

  (1) The person did not satisfy any of the case definitions above,
  (2) AND, they:
      (1) Had withdrawn consent for hospital record linkage (see 'lost_to_follow_up_reason' in section "Follow-up
          details" below)
      (3) OR, had no data in the respective nurse interview question (UK Biobank field #20002), e.g. due to a "Prefer
          not to answer" in the preceeding touchscreen survey question #2473 ("Has a doctor ever told you that you have
          had any other serious medical conditions or disabilities?" or were not interviewed by the nurse
      (4) OR, they:
          (1) Had no data on the mental health touchscreen survey questions (participants attending UK Biobank
              assessment prior to April 2009; https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=100060)
          (2) AND, they had no data in the respective follow-up nurse interview question (UK Biobank field #20002), e.g.
              due to a "Prefer not to answer" or "Do no know" with nurse unable to determine any details on medical
              history.

Note QDiabetes 2018 risk could not be computed for participants with missing data.

-------------------------
'hypertension_medication'
-------------------------
This column curates information about the participant's history of treated hypertension using medication history data
available from touchscreen survey questions (UK Biobank fields #6177 and #6153) and from information curated by the
trained nurse from the participant's current medications (UK Biobank field #20003).

This column contains 'TRUE' (N=109,611) where:

  (1) The participant reported taking medication for blood pressure in the touchscreen survey questions asking "Do you
      regularly take any of the following medications?" (UK Biobank field #6177 for males, field #6153 for females),
  (2) OR, the trained nurse recorded them as taking one of the curated list of drugs for hypertension in UK Biobank
      field #20003.

The list of possible drugs for treating hypertension in UK Biobank field #20003 was curated by me from the British
National Formularly (BNF) chapter 2.5: Hypertension and heart failure (https://openprescribing.net/bnf/0205/). These
were overlaid with the list of medications in field #20003 (https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=4) and
cross-checked against their indications (https://bnf.nice.org.uk/drug/) to ensure that they captured the intended
treatment effect (e.g. here were interested just in drugs used to treat hypertension, not those for heart failure). Note
that not all drugs listed in the BNF are found in UKB, and some need to be looked up by components of their name.

The following were the drugs considered as evidence for current treatment of hypertension:

 UK Biobank Code                                 UK Biobank medication name
---------------- ----------------------------------------------------------
      1141186674                                                   bosentan
      1140888686                                                hydralazine
      1140860532                                                  minoxidil
      1141168936                                                 sildenafil
      1141187810                                                  tadalafil
      1140883468                                                  clonidine
      1140871986                clonidine hydrochloride 25micrograms tablet
      1140860470                                                 methyldopa
      1140860562           methyldopa+hydrochlorothiazide 250mg/15mg tablet
      1140910606                                           alpha methyldopa
      1140928284                                                 moxonidine
      1140888536                                               guanethidine
      1140879778                                                  doxazosin
      1140879782                                                  indoramin
      1141157490                                          indoramin product
      1140879794                                                   prazosin
      1140879798                                                  terazosin
      1141156836                                      candesartan cilexetil
      1140860750                                                  captopril
      1140860764           captopril+hydrochlorothiazide 25mg/12.5mg tablet
      1140860882                                                 cilazapril
      1141181186                             co-zidocapt 25mg/12.5mg tablet
      1140888552                                                  enalapril
      1140860790   enalapril maleate+hydrochlorothiazide 20mg/12.5mg tablet
      1141171336                                                 eprosartan
      1140888556                                                 fosinopril
      1141164148                                    imidapril hydrochloride
      1141152998                                                 irbesartan
      1141172682         irbesartan+hydrochlorothiazide 150mg/12.5mg tablet
      1140860696                                                 lisinopril
      1140864952          lisinopril+hydrochlorothiazide 10mg/12.5mg tablet
      1140916356                                                   losartan
      1141151016  losartan potassium+hydrochlorothiazide 50mg/12.5mg tablet
      1140923712                                                  moexipril
      1141193282                                                 olmesartan
      1140879802                                                 amlodipine
      1140888560                                                perindopril
      1141180592                                     perindopril+indapamide
      1140860728                                                  quinapril
      1140860806                                                   ramipril
      1141165470                                        felodipine+ramipril
      1141166006                                                telmisartan
      1141187788         telmisartan+hydrochlorothiazide 40mg/12.5mg tablet
      1140860904                                               trandolapril
      1141153328                       trandolapril+verapamil hydrochloride
      1140888510                                                  verapamil
      1141145660                                                  valsartan
      1141201038           valsartan+hydrochlorothiazide 80mg/12.5mg tablet
---------------- ----------------------------------------------------------

Data were missing where it was not possible to determine presence OR absence of treatment of hypertension (N=8,080),
where:

  (1) There was no touchscreen survey or verbal interview data for the participant for UK Biobank fields #6177, #6153,
      and #20003
  (2) OR, the participant answered "Do not know" or "Prefer not to answer" to fields #6177 (males) or #6153 (females)
      (1) AND, they had no data in the respective follow-up nurse interview question (UK Biobank field #20003), e.g.
          due to a "Prefer not to answer" or "Do no know" with nurse unable to determine any details on current
          medications,

Note QDiabetes risk could not be computed for participants with missing data.

---------------------------
'lipid_lowering_medication'
---------------------------
This column curates information about the participant's statin usage, using medication history data available from
touchscreen survey questions (UK Biobank fields #6177 and #6153) and from information curated by the trained nurse from
the participant's current medications (UK Biobank field #20003).

This column contains 'TRUE' (N=90,405) where:

  (1) The participant reported taking cholesterol lowering medication in the touchscreen survey questions asking "Do you
      regularly take any of the following medications?" (UK Biobank field #6177 for males, field #6153 for females),
  (2) OR, the trained nurse recorded them as taking one of the curated list of drugs for lipid lowering medication in UK
      Biobank field #20003.

The list of possible drugs for treating hypertension in UK Biobank field #20003 was curated by me from the British
National Formularly (BNF) chapter 2.12: Lipid-regulating drugs (https://openprescribing.net/bnf/0212/). Here, we
considered all drugs that are used to lower lipid levels for CVD prevention, rather than just statins alone. Drugs were
cross-checked against their indications (https://bnf.nice.org.uk/drug/) to ensure they captured the intended treatment
effect (e.g. while ispaghula husk in BNF 2.12, its primary indication is for constipation, with lipid lowering as a
side-effect, so it is omitted), and intersected with the list of medications in UK Biobank field #20003
(https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=4). Note that not all drugs listed in the BNF are found in UKB,
and some need to be looked up by components of their name.

The following were the drugs considered as lipid lowering treatments:

 UK Biobank Code                 UK Biobank medication name
---------------- ------------------------------------------
      1140861892                                   acipimox
      1141146234                               atorvastatin
      1140861924                                bezafibrate
      1141157260                        bezafibrate product                                                                                                                                                                                        1140862026                               ciprofibrate
      1140888590                                 colestipol
      1140909780                              colestyramine
      1141180734                      colestyramine product
      1141180722   colestyramine+aspartame 4g/sachet powder
      1141192736                                  ezetimibe
      1140861954                                fenofibrate
      1140888594                                fluvastatin
      1140861856                                gemfibrozil
      1141157262                        gemfibrozil product
      1140861868                     nicotinic acid product
      1140888648                                pravastatin
      1141192410                               rosuvastatin
      1140861958                                simvastatin
---------------- ------------------------------------------

Data were missing where it was not possible to determine presence OR absence of lipid lowering medication (N=8,085),
where:

  (1) There was no touchscreen survey or verbal interview data for the participant for UK Biobank fields #6177, #6153,
      and #20003
  (2) OR, the participant answered "Do not know" or "Prefer not to answer" to fields #6177 (males) or #6153 (females)
      (1) AND, they had no data in the respective follow-up nurse interview question (UK Biobank field #20003), e.g.
          due to a "Prefer not to answer" or "Do no know" with nurse unable to determine any details on current
          medications,

Note QDiabetes 2018 risk could not be computed for participants with missing data.

----------------------------
'systematic_corticosteroids'
----------------------------
This column curates information about the participant's history of systematic corticosteriod usage (N=5,284 with 'TRUE'
in column), using information curated by the trained nurse from the participant's current medications (UK Biobank field
#20003).

This list is my best-effort to curate the list that most closely matches what is listed in the QDiabetes paper,
Hippisley-Cox and Coupland. Development and validation of QDiabetes-2018 risk prediction algorithm to estimate future
risk of type 2 diabetes: cohort study. BMJ 359, j5019 (2017), using the British National Formulary (BNF) chapter 6.3.2:
Glucocorticoid therapy (https://openprescribing.net/bnf/060302/). Here, the list of drugs is intended to capture
long-term use of (strong) corticosteroids (oral or injected) not transient (e.g. over the counter topical treatments) or
preventive inhalers for asthma.

The following were the drugs considered as systematic corticosteroids:

 UK Biobank Code                 UK Biobank medication name
---------------- ------------------------------------------
      1140874790                              betamethasone
      1141145782                                deflazacort
      1140874816                              dexamethasone
      1140874896                             hydrocortisone
      1140874976                         methylprednisolone
      1140874930                               prednisolone
      1141157402                       prednisolone product
      1140868364                                 prednisone
      1140868426                              triamcinolone
---------------- ------------------------------------------

Data were missing where it was not possible to determine presence OR absence of systematic corticosteroids (N=2,729),
where:

  (1) There was no verbal interview medication data for the participant in field #20003 due to a "Prefer not to answer"
      or "Do no know" with nurse unable to determine any details on current medications (e.g. distinct from people with
      no data in field #20003 due to no current medications, which are set to 'FALSE')

Note QDiabetes risk could not be computed for participants with missing data.

-------------------------
'atypical_antipsychotics'
-------------------------
This column curates information about the participant's history of second generation atypical antipsychotics usage
(N=1,346 with 'TRUE' in column) using information curated by the trained nurse from the participant's current
medications (UK Biobank field #20003).

This list of drugs was directly obtained from the QDiabetes paper, Hippisley-Cox and Coupland. Development and
validation of QDiabetes-2018 risk prediction algorithm to estimate future risk of type 2 diabetes: cohort study. BMJ
359, j5019 (2017):

 UK Biobank Code                 UK Biobank medication name
---------------- ------------------------------------------
      1141153490                                amisulpride
      1141195974                               aripiprazole
      1140867420                                  clozapine
      1140928916                                 olanzapine
      1141152848                                 quetiapine
      1140867444                                risperidone
      1140927956                                 sertindole
      1141169714                                   zotepine
---------------- ------------------------------------------

Data were missing where it was not possible to determine presence OR absence of atypical antipsychotics (N=2,729),
where:

  (1) There was no verbal interview medication data for the participant in field #20003 due to a "Prefer not to answer"
      or "Do no know" with nurse unable to determine any details on current medications (e.g. distinct from people with
      no data in field #20003 due to no current medications, which are set to 'FALSE')

Note QDiabetes 2018 risk could not be computed for participants with missing data.

-------
'hba1c'
-------
This column contains glycated haemoglobin (HbA1c) measurements (UK Biobank field #30750). Values reported as missing due
to being above or below detection limits (see missingness reason from UK Biobank field #30756) were set to the maximum
and minimum values present in the data (with a +/- 0.0001 offset to distinguish from values within detection limits).

Values below detection limit were those with "Reportable at assay but not reportable after any corrections (too low)" or
"Not reportable at assay (too low)" in UK Biobank field #30576, and conversely, values above detection limit were those
with "Reportable at assay but not reportable after any corrections (too high)" or "Not reportable at assay (too high)"
in field #30576.

Note that QDiabetes 2018 model C risk could not be computed for participants missing HbA1c measurement (N=35,802) or
where HbA1c >= 48 mmol/mol (N=17,608).

-----------------
'fasting_glucose'
-----------------
This column contains glucose measurements (UK Biobank field #30740) corrected for fasting time (UK Biobank field #74) as
described above. Values reported as missing due to being above or below detection limits (see missingness reason from UK
Biobank field #30746) were set to the maximum and minimum values present in the data prior to correction for fasting
time (with a +/- 0.0001 offset to distinguish from values within detection limits).

Values below detection limit were those with "Reportable at assay but not reportable after any corrections (too low)" or
"Not reportable at assay (too low)" in UK Biobank field #30746, and conversely, values above detection limit were those
with "Reportable at assay but not reportable after any corrections (too high)" or "Not reportable at assay (too high)"
in field #30746.

Here, fasting time was corrected following the findings of Moebus et al. Impact of time since last caloric intake on
blood glucose levels. Eur. J. Epidemiol. 26, 719–728 (2011); that while fasting glucose usually requires a person to
fast for at least 8 hours, glucose levels don't materially change after 3 hours of fasting:

  (1) Median glucose levels were compute for males and females separately at 0 hours fasting, 1 hours fasting, 2 hours
      fasting, and 3 or more hours of fasting.

  (2) Glucose levels were shifted based on the difference in sex-specific medians from those fasting for 3 or more
      hours, i.e. so that the new median glucose would be idetentical across the different fasting time groups within
      each sex (4.903 mmol/L in females, 4.946 mmol/L in males; see below)

This reflects (1) the granularity of the fasting time information, (2) that change in glucose over time is non-linear
within the first 3 hours since last meal, and (3) that the change in glucose levels on average differs between males and
females.

The median glucose and adjustment factors were as follows:

    Fasting time     Sex  Median glucose (mmol/L)  Adjustment factor (mmol/L)
---------------- ------- ------------------------ ---------------------------
 3 hours or more  Female                    4.903                           0
         2 hours  Female                    4.920                      -0.017
         1 hours  Female                    5.111                      -0.208
         0 hours  Female                    5.347                      -0.444
---------------- ------- ------------------------ ---------------------------
 3 hours or more    Male                    4.946                           0
         2 hours    Male                    4.962                      -0.016
         1 hours    Male                    5.088                      -0.142
         0 hours    Male                    5.334                      -0.388
---------------- ------- ------------------------ ---------------------------

Note that QDiabetes 2018 model B risk could not be computed using the 'fasting_glucose' column
('QDiabetes2018B_fasting_glucose' in this instance) where participants missing glucose measurement (N=72,922), fasting
time ('fasting_glucose' could not be determined; N=11 participants with non-missing 'non_fasting_glucose'), or where
'fasting_glucose' >= 7 mmol/L (N=15,863).

--------------
'fasting_time'
--------------
This column contains the number of hours since consumption of last food or drink at time of sample collection, from
UK Biobank field #74.

Data are missing for N=4,100 participants, but most of these (N=4,089) are also missing glucose measurements (for
unrelated reasons).

---------------------
'non_fasting_glucose'
---------------------
This column contains glucose measurements (UK Biobank field #30740). Values reported as missing due to being above or
below detection limits (see missingness reason from UK Biobank field #30746) were set to the maximum and minimum values
present in the data (with a +/- 0.0001 offset to distinguish from values within detection limits).

Values below detection limit were those with "Reportable at assay but not reportable after any corrections (too low)" or
"Not reportable at assay (too low)" in UK Biobank field #30746, and conversely, values above detection limit were those
with "Reportable at assay but not reportable after any corrections (too high)" or "Not reportable at assay (too high)"
in field #30746.

Here, no correction for fasting time has taken place.

Note that QDiabetes 2018 model B risk could not be computed using the 'non_fasting_glucose' column
('QDiabetes2018B_non_fasting_glucose' in this instance) where participants missing glucose measurement (N=72,922), or
where 'non_fasting_glucose' >= 7 mmol/L (N=16,136).

------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------
Prevalent diabetes
---------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
Columns described in this section contain details on the adjudication of diabetes history for each participant.

--------------------
'prevalent_diabetes'
--------------------
This column curates information about whether the participant has prevalent diabetes (N=26,445 with 'TRUE' in column),
for whom QDiabetes risk prediction and incident disease modelling and incident diabetes analysis cannot be performed.
This includes anyone adjudicated to have prevalent type 2 diabetes (see 'type_2_diabetes' below), type 1 diabetes (see
'type_1_diabetes' below), or uncertain diabetes status (or type) (see 'uncertain_diabetes' below).

Participants with history of gestational diabetes (see 'history_gestational_diabetes' above) but no other history of
diabetes are not counted as prevalent cases (i.e. are 'FALSE' in this column).

Data is missing for N=149 participants where:

  (1) The participants had withdrawn consent for hospital record linakge,
  (2) AND, they were not adjudicated to have type 1 diabetes, type 2 diabetes, or uncertain diabetes from the
      self-report touchscreen survey or verbal interview data.

Note that QDiabetes risk predictions and incident disease data are only available where 'prevalent_diabetes' is FALSE.

-------------------------
'no_history_any_diabetes'
-------------------------
This column curates information about whether the participant has history of any diabetes, including gestational
diabetes.

I.e., this column contains 'TRUE' (N=475,021) where:

  (1) 'prevalent_diabetes' is 'FALSE'
  (2) AND, 'history_gestational_diabetes' is 'FALSE'

Data is missing for N=149 participants where:

  (1) The participants had withdrawn consent for hospital record linakge,
  (2) AND, they were adjudicated as "Diabetes unlikely" in the self-report data.

-----------------
'type_2_diabetes'
-----------------
This column curates information about the participant's type 2 diabetes status at time of assessment. This data was
curated using a combination of the self-report touchscreen and verbal interview questions as well as retrospective
hospital records.

First, type 2 diabetes status was curated from self-report touchscreen and verbal interview questions using the
algorithms described in Eastwood et al. Algorithms for the Capture and Adjudication of Prevalent and Incident Diabetes
in UK Biobank. PLoS One 11, e0162388 (2016).

Then, additional type 2 diabetes cases were determined from retrospective hospital records using the endpoint definition
file 'data/adjudicated_medical_history/type_2_diabetes/endpoint_definition.txt' in conjuction with the script
'src/adjudicate_medical_history/curate_endpoint.R' (a symbolic link to the program in '../endpoints/').

Briefly, this column contains 'TRUE' (N=24,381) where:

  (1) They reported previous type 2 diabetes diagnosis in interview with trained nurse (UK Biobank field #20002, code
      1223) and age diagnosis > 30 years of age for participants of self-reported South Asian or Africa/Caribbean
      ethnicity, or age diagnosis > 36 years of age for participants of self-reported European or Mixed/Other ethnicity
      (see table below for ethnicity mappings used here).
  (2) OR, they reported prescription of metformin medication in interview with trained nurse (UK Biobank field #20003),
      (1) AND, they did not report prescription of insulin medication
      (2) AND, they reported any diabetes in either the touchscreen survey (UK Biobank field #2443) or in interview with
          trained nurse (UK Biobank field #20002, codes 1220 [diabetes] or 1223 [type 2 diabetes]).
  (3) OR, they reported taking any current oral diabetes medication other than metformin (or insulin) in interview with
      trained nurse (UK Biobank field #20003) (see table below for list of medications)
  (4) OR, they had any hospital records prior to UK Biobank assessment for type 2 diabetes (ICD-10 codes E11 or O24.1,
      or ICD-9 codes 250.00 or 250.10)
      (1) AND, they did not have any hospital records indicative of uncertain diabetes type or status (ICD-10 codes E13
          or E14, or ICD-9 codes 250.0, 250.09, 250.19, 250.29, or 250.99)
      (2) AND, they did not have any hospital records indicative of type 1 diabetes (ICD-10 codes E10 or O24.0, or ICD-9
          codes 250.01 or 250.11)

Ethnicity mapping used for determining age-at-diagnosis cut-offs for type 2 diabetes adjudication by the Eastwood et al.
algorithm as described above:

 Ethnicity group for Eastwood algorithms          UK Biobank category  Field #21000 code
---------------------------------------- ---------------------------- ------------------
                                European                        White                  1
                                                              British               1001
                                                                Irish               1002
                                           Any other white background               1003
                             South Asian       Asian or Asian British                  3
                                                               Indian               3001
                                                            Pakistani               3002
                                                          Bangladeshi               3003
                       African/Caribbean       Black or Black British                  4
                                                            Caribbean               4001
                                                              African               4002
                                           Any other Black background               4003
                             Mixed/Other                        Mixed                  2
                                            White and Black Caribbean               2001
                                              White and Black African               2002
                                                      White and Asian               2003
                                           Any other mixed background               2004
                                           Any other Asian background               3004
                                                              Chinese                  5
                                                   Other ethnic group                  6
                                    <NA>                  Do not know                 -1
                                                 Prefer not to answer                 -3
                                                                 <NA>               <NA>
---------------------- ----------------- ---------------------------- ------------------

Note, rule 1 was not applied (returns 'FALSE') for people in the missing data <NA> ethnicity group.

The full list of diabetes medication codes (insulin/metformin/non-metformin) is given in the table below:

 Medication class  UK Biobank field #20003 code               UK Biobank field #20003 label
----------------- ----------------------------- -------------------------------------------
          Insulin                    1140883066                             Insulin product
        Metformin                    1140884600                                   Metformin
                                     1140874686                     Glucophage 500mg tablet
                                     1141189090  Rosiglitazone 1mg / metformin 500mg tablet
    Sulfonylureas                    1140874718                               Glibenclamide
                                     1140874744                                  Gliclazide
                                     1140874746                       Diamicron 80mg tablet
                                     1141152590                                 Glimepiride
                                     1141156984                           Amaryl 1mg tablet
                                     1140874646                                   Glipizide
                                     1141157284                           Glipizide product
                                     1140874652                       Minodiab 2.5mg tablet
                                     1140874674                                 Tolbutamide
                                     1140874728                       Euglucon 2.5mg tablet
     Meglitinides                    1141173882                                 Nateglinide
                                     1141173786                         Starlix 60mg tablet
                                     1141168660                                 Repaglinide
       Glitazones                    1141171646                                Pioglitazone
                                     1141171652                           Actos 15mg tablet
                                     1141153254                                Troglitazone
                                     1141177600                               Rosiglitazone
                                     1141177606                          Avandia 4mg tablet
            Other                    1140868902                                    Acarbose
                                     1140868908                        Glucobay 50mg tablet
                                     1140857508             Glucotard 5g/sachet mini-tablet
----------------- ----------------------------- -------------------------------------------

Note that retrospective follow-up in hospital records only goes back so far, depending on the nation of the hospital in
which the hospital episodes occurred. The earliest follow-up dates are currently 1993-07-27 for hospitals in England,
1991-04-18 for hospitals in Wales, and 1980-12-02 for hospitals in Scotland, where baseline UK Biobank assessment
occurred between 2006 and 2010 (and participants were at minimum 37 years old at baseline assessment).

The earliest hypothetical date at which hospital records might be found for each participant, based on inferred nation
of residence at those dates, are coded in the 'earliest_hospital_date' and 'earliest_hospital_nation' columns described
in more detail in the section "Follow-up details" below.

Data were missing where it was not possible to determine case OR control status for type 2 diabetes at UK Biobank
assessment (N=344), where:

  (1) They were classified as having uncertain diabetes status or type (see 'uncertain_diabetes' below)
  (2) OR, they had withdrawn consent for hospital record linkage (see 'lost_to_follow_up_reason' in section "Follow-up
      details" below)

Note that QDiabetes risk predictions and incident disease data are only available where 'type_2_diabetes' is FALSE.

The 'diabetes_history_reason' below contains further details on the adjudication diabetes history for each participant,
including distinguishing between type 2 diabetes identified from self-report vs. hospital records, and classification of
type 2 diabetes self-reprot adjudication into to "Possible" and "Probable" type 2 diabetes.

For clarification of details on adjudication from self-report data using the Eastwood et al. algorithms, contact Dr.
Sam Lambert (sl925@medschl.cam.ac.uk) who implemented the code in the iPython/Jupyter notebooks in
'data/adjudicated_medical_history/Eastwood_diabetes/' (symbolically links to directory elsewhere on the cluster).

-----------------
'type_1_diabetes'
-----------------
This column curates information about the participant's type 1 diabetes status at time of assessment. This data was
curated using a combination of the self-report touchscreen and verbal interview questions as well as retrospective
hospital records.

First, type 1 diabetes status was curated from self-report touchscreen and verbal interview questions using the
algorithms described in Eastwood et al. Algorithms for the Capture and Adjudication of Prevalent and Incident Diabetes
in UK Biobank. PLoS One 11, e0162388 (2016).

Then, additional type 1 diabetes cases were determined from retrospective hospital records using the endpoint definition
file 'data/adjudicated_medical_history/type_1_diabetes/endpoint_definition.txt' in conjuction with the script
'src/adjudicate_medical_history/curate_endpoint.R' (a symbolic link to the program in '../endpoints/').

Briefly, this column contains 'TRUE' (N=1,869) where:

 (1) They reported previous type 1 diabetes diagnosis in interview with trained nurse (UK Biobank field #20002, code
     1222),
 (2) OR, they reported prescription of current insulin in interview with trained nurse (UK Biobank field #20003, code
     1140883066)
 (3) OR, they reported starting insulin treatment within 12 months of diabetes diagnosis in the touchscreen survey
     question (UK Biobank field #2986) asked to those reporting previous diabetes diagnosis by doctor (UK Biobank field
     #2443)
  (4) OR, they had any hospital records prior to UK Biobank assessment for type 1 diabetes (ICD-10 codes E10 or O24.0,
      or ICD-9 codes 250.01 or 250.11)
      (1) AND, they did not have any hospital records indicative of uncertain diabetes type or status (ICD-10 codes E13
          or E14, or ICD-9 codes 250.0, 250.09, 250.19, 250.29, or 250.99)
      (2) AND, they did not have any hospital records indicative of type 2 diabetes (ICD-10 codes E11 or O24.1, or ICD-9
          codes 250.00 or 250.10)


Note that retrospective follow-up in hospital records only goes back so far, depending on the nation of the hospital in
which the hospital episodes occurred. The earliest follow-up dates are currently 1993-07-27 for hospitals in England,
1991-04-18 for hospitals in Wales, and 1980-12-02 for hospitals in Scotland, where baseline UK Biobank assessment
occurred between 2006 and 2010 (and participants were at minimum 37 years old at baseline assessment).

The earliest hypothetical date at which hospital records might be found for each participant, based on inferred nation
of residence at those dates, are coded in the 'earliest_hospital_date' and 'earliest_hospital_nation' columns described
in more detail in the section "Follow-up details" below.

Data were missing where it was not possible to determine case OR control status for type 1 diabetes at UK Biobank
assessment (N=344), where:

  (1) They were classified as having uncertain diabetes status or type (see 'uncertain_diabetes' below)
  (2) OR, they had withdrawn consent for hospital record linkage (see 'lost_to_follow_up_reason' in section "Follow-up
      details" below)

Note that QDiabetes risk predictions and incident disease data are only available where 'type_1_diabetes' is FALSE.

The 'diabetes_history_reason' below contains further details on the adjudication diabetes history for each participant,
including distinguishing between type 1 diabetes identified from self-report vs. hospital records, and classification of
type 1 diabetes self-reprot adjudication into to "Possible" and "Probable" type 1 diabetes.

For clarification of details on adjudication from self-report data using the Eastwood et al. algorithms, contact Dr.
Sam Lambert (sl925@medschl.cam.ac.uk) who implemented the code in the iPython/Jupyter notebooks in
'data/adjudicated_medical_history/Eastwood_diabetes/' (symbolically links to directory elsewhere on the cluster).

---------------------------------
'uncertain_diabetes'
---------------------------------
This column curates information about cases where it was not possible to ascertain a participant's diabetes status,
or type of diabetes from either the self-report touchscreen questions, self-report verbal interview questions, or
retrospective hospital records.

Briefly, this column contains 'TRUE' (N=195) where:

  (1) The trained nurse recorded "gestational diabetes" in the non-cancer related illness field (UK Biobank field
      #20002) for participants self-reporting as male (these most likely represent coding errors by the nurse, as the
      genetic sex determined from genotype data is consistent with self-reported sex in these instances).
      (1) AND, the participant was not adjudicated as having type 1 diabetes or type 2 diabetes from their hospital
          records prior to UK Biobank assessment
  (2) OR, the participant had a history of uncertain diabetes or diabetes complications in their hospital records prior
      to UK Biobank assessment (ICD-10 codes E13, E14, G59.0, G63.2, H28.0, H36.0, M14.2, N08.3, or ICD-9 codes 250.0,
      or 250.1-250.90)
      (1) AND, the participant was not adjudicated as having type 1 diabetes or type 2 diabetes from either their
          self-report data or hospital records above.
  (3) OR, the participant had a history of unspecified diabetes/complications in pregnancy (ICD-10 code O24.3 or ICD-9
      codes 648.0 or 648.8)
      (1) AND, the participant was not adjudicated as having gestational diabetes from either their self-report data or
          hospital records (see 'history_gestational_diabetes' above).

Note that retrospective follow-up in hospital records only goes back so far, depending on the nation of the hospital in
which the hospital episodes occurred. The earliest follow-up dates are currently 1993-07-27 for hospitals in England,
1991-04-18 for hospitals in Wales, and 1980-12-02 for hospitals in Scotland, where baseline UK Biobank assessment
occurred between 2006 and 2010 (and participants were at minimum 37 years old at baseline assessment).

The earliest hypothetical date at which hospital records might be found for each participant, based on inferred nation
of residence at those dates, are coded in the 'earliest_hospital_date' and 'earliest_hospital_nation' columns described
in more detail in the section "Follow-up details" below.

Data were set to missing where (N=149):

  (1) The participant was not classified as having type 1 diabetes or type 2 diabetes from self-report data as above,
  (2) AND, they had withdrawn consent for hospital record linkage (see 'lost_to_follow_up_reason' in section "Follow-up
      details" below)

The 'diabetes_history_reason' below contains further details on the adjudication diabetes history for each participant.

-------------------------
'diabetes_history_reason'
-------------------------
This column curates information about the adjudication and stratification of participants by diabetes case/control
status and type.

In general, these rules can be summarised as:

  (1) Adjudication in self-report takes precendence over hospital records
  (2) But, adjudication of Type 1 or Type 2 diabetes takes precedence over gestational diabetes or uncertain diabetes
      (1) e.g. if self-report adjudication results in "Uncertain diabetes" or "Possible gestational diabetes" but
          hospital records support type 2 diabetes, then the adjudication status becomes "Possible Type 2 diabetes
          (hospital records)
  (3) Gestational diabetes is only adjudicated when it can specifically be determined in the absence of any other
      diabetes (including subsequent diagnosis of type 2 diabetes years after preganancy)
  (4) Uncertain diabetes is only adjudicated where self-report or hospital record data indicates presence of diabetes,
      but no call can be made as to type 1, or type 2 (or gestational diabetes specifically for unspecified diabetes
      in pregnancy).
  (5) Withdrawn consent for hospital linkage occurs only where self-report adjudicates "Diabetes unlikely" but hospital
      record data is absent due to withdrawn consent for linkage.

When adjudicating type 2 or type 1 diabetes from self-report data, the Eastwood et al. algorithms further stratify
classification into "Probable" or "Possible" type 2 or type 1 diabetes, depending on the confidence of information.

For type 2 diabetes, participants are downgraded from "Probable Type 2 diabetes" to "Possible Type 2 diabetes" where:

  (1) The trained nurse notes a self-report of either any diabetes or type 2 diabetes, but not type 1 diabetes, in
      verbal interview (UK Biobank field #20002, code 1220 or 1223, but not 1222)
  (2) AND, they self-report taking insulin medication in trained interview with nurse (UK Biobank field #20003, code
      1140883066)
  (3) AND, they report age diagnosis > 30 years of age for participants of self-reported South Asian or Africa/Caribbean
      ethnicity, or age diagnosis > 36 years of age for participants of self-reported European or Mixed/Other ethnicity
      (see table above in 'type_2_diabetes' for ethnicity mappings used here).

For type 1 diabetes, participants are downgraded from "Probable Type 1 diabetes" to "Possible Type 1 diabetes" where:

  (1) They self-report any diabetes in verbal interview with nurse (UK Biobank field #20002, code 1220)
  (2) AND, they *do not* self-report type 1 or type 2 diabetes in interview with nurse (UK Biobank field #20002, codes
      1222 or 1223)
  (3) BUT, they *do* report taking current insulin medication and/or starting insulin medication within 1 year of
      diagnosis (UK Biobank field #20003 code 1140883066 or "Yes" in UK Biobank field #2986)
  (3) AND, they *do not* self-report any current non-metformin diabetes medication in interview with nurse (UK Biobank
      field #20003, see table above in 'type_2_diabetes' for full list of medication codes)
  (4) AND, they self-report age diagnosis < 30 years of age for participants of self-reported South Asian or
      Africa/Caribbean ethnicity, or age diagnosis < 36 years of age for participants of self-reported European or
      Mixed/Other ethnicity (see table above in 'type_2_diabetes' for ethnicity mappings used here).

Overall, this column contains the following possible values:

                     Diabetes adjudication status   Number of participants
------------------------------------------------- ------------------------
           Probable Type 2 diabetes (self-report)                   20,564
           Possible Type 2 diabetes (self-report)                    3,277
      Possible Type 2 diabetes (hospital records)                      540
           Probable Type 1 diabetes (self-report)                    1,398
           Possible Type 1 diabetes (self-report)                      433
      Possible Type 1 diabetes (hospital records)                       38
          Uncertain diabetes status (self-report)                       43
     Uncertain diabetes status (hospital records)                      152
      Possible gestational diabetes (self-report)                      779
 Possible gestational diabetes (hospital records)                       65
                                Diabetes unlikely                  475,021
    Withdrawn consent for hospital record linkage                      149
------------------------------------------------- ------------------------

This adjudication status determines the values of the other columns above as follows:

                     Diabetes adjudication status  prevalent_diabetes  no_history_any_diabetes
------------------------------------------------- ------------------- ------------------------
           Probable Type 2 diabetes (self-report)                TRUE                    FALSE
           Possible Type 2 diabetes (self-report)                TRUE                    FALSE
      Possible Type 2 diabetes (hospital records)                TRUE                    FALSE
           Probable Type 1 diabetes (self-report)                TRUE                    FALSE
           Possible Type 1 diabetes (self-report)                TRUE                    FALSE
      Possible Type 1 diabetes (hospital records)                TRUE                    FALSE
          Uncertain diabetes status (self-report)                TRUE                    FALSE
     Uncertain diabetes status (hospital records)                TRUE                    FALSE
      Possible gestational diabetes (self-report)               FALSE                    FALSE
 Possible gestational diabetes (hospital records)               FALSE                    FALSE
                                Diabetes unlikely               FALSE                     TRUE
    Withdrawn consent for hospital record linkage                <NA>                     <NA>
------------------------------------------------- ------------------- ------------------------

                     Diabetes adjudication status  type_2_diabetes  type_1_diabetes  diabetes_uncertain
------------------------------------------------- ---------------- ---------------- -------------------
           Probable Type 2 diabetes (self-report)             TRUE            FALSE               FALSE
           Possible Type 2 diabetes (self-report)             TRUE            FALSE               FALSE
      Possible Type 2 diabetes (hospital records)             TRUE            FALSE               FALSE
           Probable Type 1 diabetes (self-report)            FALSE             TRUE               FALSE
           Possible Type 1 diabetes (self-report)            FALSE             TRUE               FALSE
      Possible Type 1 diabetes (hospital records)            FALSE             TRUE               FALSE
          Uncertain diabetes status (self-report)             <NA>             <NA>                TRUE
     Uncertain diabetes status (hospital records)             <NA>             <NA>                TRUE
      Possible gestational diabetes (self-report)            FALSE            FALSE               FALSE
 Possible gestational diabetes (hospital records)            FALSE            FALSE               FALSE
                                Diabetes unlikely            FALSE            FALSE               FALSE
    Withdrawn consent for hospital record linkage             <NA>             <NA>                <NA>
------------------------------------------------- ---------------- ---------------- -------------------

                     Diabetes adjudication status     sex  history_gestational_diabetes
------------------------------------------------- ------- -----------------------------
           Probable Type 2 diabetes (self-report)  Female                          <NA>
           Possible Type 2 diabetes (self-report)  Female                          <NA>
      Possible Type 2 diabetes (hospital records)  Female                          <NA>
           Probable Type 1 diabetes (self-report)  Female                          <NA>
           Possible Type 1 diabetes (self-report)  Female                          <NA>
      Possible Type 1 diabetes (hospital records)  Female                          <NA>
          Uncertain diabetes status (self-report)  Female                          <NA>
     Uncertain diabetes status (hospital records)  Female                          <NA>
      Possible gestational diabetes (self-report)  Female                          TRUE
 Possible gestational diabetes (hospital records)  Female                          TRUE
                                Diabetes unlikely  Female                         FALSE
    Withdrawn consent for hospital record linkage  Female                          <NA>
                                                -    Male                         FALSE
------------------------------------------------- ------- -----------------------------

------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------
Possible undiagnosed diabetes
---------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------

The following columns flag possible cases of undiagnosed diabetes based on thresholds of HbA1c, fasting glucose, or
random (non-fasting) glucose outlined in the UK NICE guidelines for type 2 diabetes diagnosis:
(https://cks.nice.org.uk/topics/diabetes-type-2/diagnosis/diagnosis-in-adults/)

 - HbA1c >= 48 mmol/mol
 - Fasting glucose >= 7 mmol/L
 - Random (non-fasting) glucose >= 11.1 mmol/L

Two sets of columns are provided: one based on the non-fasting glucose measurement taken by UK Biobank, and one based on
the fasting glucose measurement we obtain by adjusting for fasting time as described above (see 'fasting_glucose'
column).

Columns based on HbA1c and fasting glucose measurements:

 - 'fasting_undiagnosed_diabetes'
 - 'fasting_undiagnosed_diabetes_reason'

Columns based on HbA1c and random (non-fasting) glucose measurements:

 - 'non_fasting_undiagnosed_diabetes'
 - 'non_fasting_undiagnosed_diabetes_reason'

In both cases, the first column  ('fasting_undiagnosed_diabetes' or 'non_fasting_undiagnosed_diabetes') give 'TRUE' or
'FALSE' values indicating whether the HbA1c or glucose levels suggest the possibility the person has undiagnosed type 2
diabetes at UK Biobank assessment, and the second column ('fasting_undiagnosed_diabetes_reason' or
'non_fasting_undiagnosed_diabetes_reason') gives the reason for each TRUE or FALSE value respectively, organised by
tiers of evidence based on the concordance/discordance of HbA1c and glucose levels and/or presence of missing data.

Importantly, note that the NICE guidelines state that a single measurement (in time) exceeding these given thresholds is
not sufficient for type 2 diabetes diagnosis, especially in the absence of symptoms, so people with 'TRUE' in these
columns should not be treated as type 2 diabetes cases, but could be excluded from controls.

------------------------------
'fasting_undiagnosed_diabetes'
------------------------------
This column contains 'TRUE' or 'FALSE' values indicating whether HbA1c ('hba1c') or fasting glucose ('fasting_glucose')
levels suggest the possibility the person has undiagnosed type 2 diabetes at UK Biobank assessment. The reason for each
'TRUE' or 'FALSE' value are given in the 'fasting_undiagnosed_diabetes_reason' column (see below for more details).

Importantly, note that the NICE guidelines state that a single measurement (in time) exceeding these given thresholds is
not sufficient for type 2 diabetes diagnosis, especially in the absence of symptoms, so people with 'TRUE' in the
'fasting_undiagnosed_diabetes' should not be treated as type 2 diabetes cases, but could be excluded from controls.

-------------------------------------
'fasting_undiagnosed_diabetes_reason'
-------------------------------------
Reason for each 'TRUE' or 'FALSE' value in the 'fasting_undiagnosed_diabetes' column, organised into tiers of evidence
based on agreement (or disagreement) of HbA1c and fasting glucose levels and/or presence of missing data in either.

Where 'fasting_undiagnosed_diabetes' is 'TRUE' (N=10,057), values in this column include:

 - "Possible undiagnosed diabetes (tier 1) | HbA1c >= 48 mmol/mol and fasting glucose >= 7 mmol/L"
 - "Possible undiagnosed diabetes (tier 2) | HbA1c >= 48 mmol/mol, missing fasting glucose measurement"
 - "Possible undiagnosed diabetes (tier 2) | Fasting glucose >= 7 mmol/L, missing HbA1c measurement"
 - "Possible undiagnosed diabetes (tier 3) | HbA1c >= 48 mmol/mol but fasting glucose < 7 mmol/L"
 - "Possible undiagnosed diabetes (tier 3) | Fasting glucose >= 7 mmol/L but HbA1c < 48 mmol/mol"

Where 'fasting_undiagnosed_diabetes' is 'FALSE', values in this column include:

 - "Adjudicated T2D case | (Eastwood et al. algorithm)"
 - "Adjudicated T2D case | (hospital records)"
 - "No evidence of diabetes (tier 1) | HbA1c < 48 mmol/mol and fasting glucose < 7 mmol/L"
 - "No evidence of diabetes (tier 2) | Fasting glucose < 7 mmol/L but missing HbA1c measurement"
 - "No evidence of diabetes (tier 2) | HbA1c < 48 mmol/mol but missing fasting glucose measurement"
 - "No evidence of diabetes (tier 3) | Missing HbA1c and fasting glucose measurements"

The thresholds for HbA1c and fasting glucose used here are based on the UK NICE guidelines for type 2 diabetes diagnosis
(https://cks.nice.org.uk/topics/diabetes-type-2/diagnosis/diagnosis-in-adults/).

Importantly, note that the NICE guidelines state that a single measurement (in time) exceeding these given thresholds is
not sufficient for type 2 diabetes diagnosis, especially in the absence of symptoms, so people with 'TRUE' in the
'fasting_undiagnosed_diabetes' should not be treated as type 2 diabetes cases, but could be excluded from controls.

-----------------------------------
'non_fasting_undiagnosed_diabetes'
-----------------------------------
This column contains 'TRUE' or 'FALSE' values indicating whether HbA1c ('hba1c') or non-fasting glucose
('non_fasting_glucose') levels suggest the possibility the person has undiagnosed type 2 diabetes at UK Biobank
assessment. The reason for each 'TRUE' or 'FALSE' value are given in the 'non_fasting_undiagnosed_diabetes_reason'
column (see below for more details).

Importantly, note that the NICE guidelines state that a single measurement (in time) exceeding these given thresholds is
not sufficient for type 2 diabetes diagnosis, especially in the absence of symptoms, so people with 'TRUE' in this
column should not be treated as type 2 diabetes cases, but could be excluded from controls.

------------------------------------------
'non_fasting_undiagnosed_diabetes_reason'
-----------------------------------------
Reason for each 'TRUE' or 'FALSE' value in the 'non_fasting_undiagnosed_diabetes' column, organised into tiers of
evidence based on agreement (or disagreement) of HbA1c and fasting glucose levels and/or presence of missing data in
either.

Where 'non_fasting_undiagnosed_diabetes' is 'TRUE' (N=5,157), values in this column include:

 - "Possible undiagnosed diabetes (tier 1) | HbA1c >= 48 mmol/mol and non-fasting glucose >= 11.1 mmol/L"
 - "Possible undiagnosed diabetes (tier 2) | HbA1c >= 48 mmol/mol, missing non-fasting glucose measurement"
 - "Possible undiagnosed diabetes (tier 2) | Non-fasting glucose >= 11.1 mmol/L, missing HbA1c measurement"
 - "Possible undiagnosed diabetes (tier 3) | HbA1c >= 48 mmol/mol but non-fasting glucose < 11.1 mmol/L"
 - "Possible undiagnosed diabetes (tier 3) | Non-fasting glucose >= 11.1 mmol/L but HbA1c < 48 mmol/mol"

Where 'non_fasting_undiagnosed_diabetes' is 'FALSE', values in this column include:

 - "Adjudicated T2D case | (Eastwood et al. algorithm)"
 - "Adjudicated T2D case | (hospital records)"
 - "No evidence of diabetes (tier 1) | HbA1c < 48 mmol/mol and non-fasting glucose < 11.1 mmol/L"
 - "No evidence of diabetes (tier 2) | Non-fasting glucose < 11.1 mmol/L but missing HbA1c measurement"
 - "No evidence of diabetes (tier 2) | HbA1c < 48 mmol/mol but missing non-fasting glucose measurement"
 - "No evidence of diabetes (tier 3) | Missing HbA1c and non-fasting glucose measurements"

The thresholds for HbA1c and non-fasting glucose used here are based on the UK NICE guidelines for type 2 diabetes
diagnosis (https://cks.nice.org.uk/topics/diabetes-type-2/diagnosis/diagnosis-in-adults/).

Importantly, note that the NICE guidelines state that a single measurement (in time) exceeding these given thresholds is
not sufficient for type 2 diabetes diagnosis, especially in the absence of symptoms, so people with 'TRUE' in the
'non_fasting_undiagnosed_diabetes' should not be treated as type 2 diabetes cases, but could be excluded from controls.

------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------
Possible prediabetes
---------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------

Similarly to the possilbe undiagnosed diabetes columns above, the following set of columns flag possible cases of
pre-diabetes based on HbA1c and fasting glucose levels:

 - HbA1c >= 42 mmol/mol
 - Fasting glucose >= 5.6 mmol/L

Two sets of columns are provided: one based on HbA1c and the fasting glucose measurement we obtain by adjusting for
fasting time as described above (see 'fasting_glucose' column), and one based on HbA1c alone (since I was unable to find
any suggested thresholds for random/non-fasting glucose to flag pre-diabetes/impaired glucose control).

Columns based on HbA1c and fasting glucose measurements:

 - 'fasting_prediabetes'
 - 'fasting_prediabetes_reason'

Columns based on HbA1c measurement alone:

 - 'non_fasting_prediabetes'
 - 'non_fasting_prediabetes_reason'

In both cases, the first column  ('fasting_prediabetes' or 'non_fasting_prediabetes') give 'TRUE' or 'FALSE' values
indicating whether the HbA1c and/or fasting glucose levels suggest the possibility the person has impaired glucose
control / pre-diabetes, and the second column ('fasting_undiagnosed_diabetes_reason' or
'non_fasting_prediabetes_reason') gives the reason for each TRUE or FALSE value respectively, organised by
tiers of evidence based on the concordance/discordance of HbA1c and/or fasting glucose levels and/or presence of missing
data.

Following the caveats of flagging possible undiagnosed diabetics above, I would also suggest that a single measurement
in time of HbA1c/fasting-glucose exceeding these given thresholds should not be sufficient for pre-diabetes diagnosis,
but could be used as additional criteria for excluding a participant from non-diabetic controls.

---------------------
'fasting_prediabetes'
---------------------
This column contains 'TRUE' or 'FALSE' values indicating whether HbA1c ('hba1c') or fasting glucose ('fasting_glucose')
levels suggest the possibility the person has impaired glucose control/pre-diabetes at UK Biobank assessment. The reason
for each 'TRUE' or 'FALSE' value are given in the 'fasting_prediabetes_reason' column (see below for more details).

----------------------------
'fasting_prediabetes_reason'
----------------------------
Reason for each 'TRUE' or 'FALSE' value in the 'fasting_prediabetes' column, organised into tiers of evidence based on
agreement (or disagreement) of HbA1c and fasting glucose levels and/or presence of missing data in either.

Where 'fasting_prediabetes' is 'TRUE' (N=52,842), values in this column include:

 - "Possible pre-diabetes (tier 1) | HbA1c >= 42 mmol/mol and fasting glucose >= 5.6 mmol/L"
 - "Possible pre-diabetes (tier 2) | HbA1c >= 42 mmol/mol, missing fasting glucose measurement"
 - "Possible pre-diabetes (tier 2) | Fasting glucose >= 5.6 mmol/L, missing HbA1c measurement"
 - "Possible pre-diabetes (tier 3) | HbA1c >= 42 mmol/mol but fasting glucose < 5.6 mmol/L"
 - "Possible pre-diabetes (tier 3) | Fasting glucose >= 5.6 mmol/L but HbA1c < 42 mmol/mol"

Where 'fasting_prediabetes' is 'FALSE', values in this column include:

 - "Adjudicated T2D case | (Eastwood et al. algorithm)"
 - "Adjudicated T2D case | (hospital records)"
 - "Possible undiagnosed diabetes (tier 1) | HbA1c >= 48 mmol/mol and fasting glucose >= 7 mmol/L"
 - "Possible undiagnosed diabetes (tier 2) | HbA1c >= 48 mmol/mol, missing fasting glucose measurement"
 - "Possible undiagnosed diabetes (tier 2) | Fasting glucose >= 7 mmol/L, missing HbA1c measurement"
 - "Possible undiagnosed diabetes (tier 3) | HbA1c >= 48 mmol/mol but fasting glucose < 7 mmol/L"
 - "Possible undiagnosed diabetes (tier 3) | Fasting glucose >= 7 mmol/L but HbA1c < 48 mmol/mol"
 - "No evidence of pre-diabetes (tier 1) | HbA1c < 42 mmol/mol and fasting glucose < 5.6 mmol/L"
 - "No evidence of pre-diabetes (tier 2) | Fasting glucose < 5.6 mmol/L but missing HbA1c measurement"
 - "No evidence of pre-diabetes (tier 2) | HbA1c < 42 mmol/mol but missing fasting glucose measurement"
 - "No evidence of pre-diabetes (tier 3) | Missing HbA1c and fasting glucose measurements"

-------------------------
'non_fasting_prediabetes'
-------------------------
This column contains 'TRUE' or 'FALSE' values indicating whether HbA1c ('hba1c') levels alone suggest the possibility
the person has impaired glucose control/pre-diabetes at UK Biobank assessment. The reason for each 'TRUE' or 'FALSE'
value are given in the 'fasting_prediabetes_reason' column (see below for more details).

--------------------------------
'non_fasting_prediabetes_reason'
--------------------------------
Reason for each 'TRUE' or 'FALSE' value in the 'non_fasting_prediabetes' column, organised into tiers of evidence.

Where 'fasting_prediabetes' is 'TRUE' (N=16,220), values in this column include:

 - "Possible pre-diabetes (tier 1) | HbA1c >= 42 mmol/mol"

Where 'fasting_prediabetes' is 'FALSE', values in this column include:

 - "Adjudicated T2D case | (Eastwood et al. algorithm)"
 - "Adjudicated T2D case | (hospital records)"
 - "Possible undiagnosed diabetes (tier 1) | HbA1c >= 48 mmol/mol and non-fasting glucose >= 11.1 mmol/L"
 - "Possible undiagnosed diabetes (tier 2) | HbA1c >= 48 mmol/mol, missing non-fasting glucose measurement"
 - "Possible undiagnosed diabetes (tier 2) | Non-fasting glucose >= 11.1 mmol/L, missing HbA1c measurement"
 - "Possible undiagnosed diabetes (tier 3) | HbA1c >= 48 mmol/mol but non-fasting glucose < 11.1 mmol/L"
 - "Possible undiagnosed diabetes (tier 3) | Non-fasting glucose >= 11.1 mmol/L but HbA1c < 48 mmol/mol"
 - "No evidence of pre-diabetes (tier 1) | HbA1c < 42 mmol/mol"
 - "No evidence of pre-diabetes (tier 3) | Missing HbA1c measurement"

------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------
Incident diabetes
---------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
Columns described in this section contain details on incident diabetes - whether, when, and which type of diabetes each
participant over 10 years of follow-up.

Incident diabetes status and estimated date of onset were adjudicated from hospital records and death records using the
algorithms described in Eastwood et al. Algorithms for the Capture and Adjudication of Prevalent and Incident Diabetes
in UK Biobank. PLoS One 11, e0162388 (2016).

Briefly;

  (1) A participant was considered to develop diabetes if, at any point in the follow-up data available in UK Biobank
      (up until 1st February 2020, to prevent confounding from the COVID19 pandemic, e.g. due to changes in behaviour
      or environment due to lockdowns, symptomatic COVID-19, or SARS-CoV2 exposure), they had a hospital record or death
      record with ICD-10 codes E10, E11, E13, or E14 (or any subcodes) as primary or secondary outcomes.
  (2) The date of onset was estimated as the mid-point between that record and the last diabetes free record (or UK
      Biobank assessment if the first occuring incident record had diabetes ICD-10 codes).
  (3) Type of incident diabetes was partitioned based on the first occuring ICD-10 code:
      (1) Type 1 diabetes for E10 (or subcodes)
      (2) Type 2 diabetes for E11 (or subcodes)
      (3) Or uncertain diabetes for E13 or E14 (or subcodes)

Data are missing (N=26,594) for incident diabetes and follow-up where:

  (1) The participant has withdrawn consent for hospital record linkage,
  (2) OR, the participant is adjudicated as having prevalent diabetes (excluding gestational) at UK Biobank assessment
      ('prevalent_diabetes') i.e., type 2 diabetes, type 1 diabetes, or uncertain diabetes.

Note that the incident diabetes definition here only captures diabetes incidentally picked up in hospital visits,
whereas diagnosis typically occurs in a primary care setting (e.g. for type 2 diabetes), i.e. the cases being identified
here are more likely to be those that are severe (or unmanaged) leading to downstream complications or being picked up
later in life as frailty increases (along with risk of hospitalisation for age-related diseases).

For clarification of details on adjudication of incident diabetes status and inference of date of onset using the
Eastwood et al. algorithms, contact Dr. Sam Lambert (sl925@medschl.cam.ac.uk) who implemented the code in the
iPython/Jupyter notebooks in 'data/adjudicated_medical_history/Eastwood_diabetes/' (symbolically links to directory
elsewhere on the cluster).

Here, I have further restricted the follow-up time to 10 years to match the time-horizon of the QDiabetes 2018 risk
scores, i.e. we only want to considered incident events where the estimated date of onset was within 10 calendar years
of UK Biobank assessment.

Here, I also preserve the human notion of calendar years when restricting the time horizon and computing follow-up time.
I.e. <N> years after a specific date always gives the same calendar day and month in the <N> following years regardless
of the difference in number of days in each year due to leap days. E.g., 10 years after 2009-05-27 will be 2019-05-27.
When UK Biobank assessment occurred on a leap day, e.g. 2008-02-29, the <N> years after is set as the 28th February if
the corresponding year is not a leap year.

Follow-up has also been updated where the participant is known by UK Biobank to be lost to follow-up, e.g. due to
participant leaving the UK or death reported by relative (see 'lost_to_followup_reason' below).

Note also that maximum follow-up time available is dependent on the nation of the reporting hospital (or death record):

Death records:

   Nation                         Record source  Maximum follow-up available
--------- ------------------------------------- ----------------------------
  England  Death registry for England and Wales                   2021-03-20
    Wales  Death registry for England and Wales                   2021-03-20
 Scotland           Death registry for Scotland                   2021-03-23
--------- ------------------------------------- ----------------------------

Hospital records:

   Nation  Maximum follow-up available
--------- ----------------------------
  England                   2021-03-31
 Scotland                   2021-05-07
    Wales                   2018-03-06
--------- ----------------------------

Note in particular that participants with hospital events in Wales have less follow-up than those in other nations,
and that this follow-up is often earlier than the 10-year time horizon (baseline assessment occurred between 2006 and
2010), so you will need to adjust any incident disease analyses for yes/no nation of hospital at date of censor in Wales
(see 'censor_hospital_nation' below), which if not accounted for will appear as a protective effect for 10-year diabetes
risk. Note that participants residing in Wales at maximum follow-up make up a small fraction of participants (N=19,728;
4.3%).

Here, I have updated the maximum follow-up to be consistent with the inferred nation of residence at maximum follow-up
(as in ../endpoints/), which has been harmonised based on: (1) the nation of reporting hospital for the most recent
hospital records for that participant, (2) the nation from which the death record was issued for that participant (if
there is one), and (3) nation of assessment centre for most recent UK Biobank assessment visit for that participant.

This inferred nation of residence at maximum follow-up can be found in 'censor_hospital_nation', and does not always
correspond to the nation at baseline assessment:

 Nation at baseline assessment  Inferred nation of residence at follow-up  Samples
         ('assessment_nation')                 ('censor_hospital_nation')
------------------------------ ------------------------------------------ --------
                      Scotland                                    England    2,558
                       England                                      Wales      964
                         Wales                                    England      859
                       England                                   Scotland      436
                         Wales                                   Scotland       17
------------------------------ ------------------------------------------ --------

Note also that baseline assessment tended to occur earlier for participants in Scotland at baseline assessment, which
means they are more likely to have more follow-up (i.e. participants are less likely than other nations to have attended
baseline assessment after 1st February 2010, so are less likely to have their follow-up truncated by the maximum
follow-up date of 1st February 2020 chosen to avoid confounding from the COVID19 pandemic), so you will need to adjust
any incident disease analyses for yes/no nation at baseline assessment in Scotland (see 'assessment_nation' below),
which if not accounted for will appear as a risk factor for 10-year diabetes risk. Note that participants residing in
Scotland at maximum follow-up make up a small fraction of participants (N=34,159; 7.7%).

--------------------------
'incident_type_2_diabetes'
--------------------------
This column curates information about whether each participant developed type 2 diabetes over the 10 years of follow-up.

This column contains 'TRUE' (N=13,555) where:

  (1) The participant had any hospital record or death record with ICD-10 code E11 as primary or secondary outcome in
      any of the available follow-up data.
  (2) AND, that record occurred before any other record coding for type 1 diabetes (ICD-10 code E10) or uncertain
      diabetes (ICD-10 codes E13 or E14)
  (3) The date of onset, estimated as the mid-point between that record and the last diabetes free record (or UK Biobank
      assessment if the first occuring incident record had type 2 diabetes code E11), was within 10 years of follow-up,
      and occurred prior to the 1st Februrary 2020.

This column contains missing data (N=27,773) where:

  (1) The participant has withdrawn consent for hospital record linkage (N=158)
  (2) OR, the participant is adjudicated as having prevalent diabetes (excluding gestational) at UK Biobank assessment
      ('prevalent_diabetes') i.e., type 2 diabetes, type 1 diabetes, or unspecified diabetes (N=26,445)
  (3) OR, the first occurring incident diabetes in hospital (or death) records was diabetes of uncertain type (ICD-10
      codes E13 or E14) (N=1,179)

Notably, this column contains FALSE where the participant's first occuring incident diabetes record was for type 1
diabetes (N=221), in addition to those that were diabetes free within the available follow-up (N=460,910).

--------------------------
'incident_type_1_diabetes'
--------------------------
This column curates information about whether each participant developed type 1 diabetes over the 10 years of follow-up.

This column contains 'TRUE' (N=221) where:

  (1) The participant had any hospital record or death record with ICD-10 code E10 as primary or secondary outcome in
      any of the available follow-up data.
  (2) AND, that record occurred before any other record coding for type 2 diabetes (ICD-10 code E11) or uncertain
      diabetes (ICD-10 codes E13 or E14)
  (3) The date of onset, estimated as the mid-point between that record and the last diabetes free record (or UK Biobank
      assessment if the first occuring incident record had type 1 diabetes code E10), was within 10 years of follow-up,
      and occurred prior to the 1st Februrary 2020.

This column contains missing data (N=27,773) where:

  (1) The participant has withdrawn consent for hospital record linkage (N=158)
  (2) OR, the participant is adjudicated as having prevalent diabetes (excluding gestational) at UK Biobank assessment
      ('prevalent_diabetes') i.e., type 2 diabetes, type 1 diabetes, or unspecified diabetes (N=26,445)
  (3) OR, the first occurring incident diabetes in hospital (or death) records was diabetes of uncertain type (ICD-10
      codes E13 or E14) (N=1,179)

Notably, this column contains FALSE where the participant's first occuring incident diabetes record was for type 2
diabetes (N=13,555), in addition to those that were diabetes free within the available follow-up (N=460,910).

-----------------------------
'incident_uncertain_diabetes'
-----------------------------
This column curates information about whether each participant developed diabetes of uncertain type over the 10 years of
follow-up.

This column contains 'TRUE' (N=1,179) where:

  (1) The participant had any hospital record or death record with ICD-10 codes E13 or E14 as primary or secondary
      outcome in any of the available follow-up data.
  (2) AND, that record occurred before any other record coding for type 1 diabetes (ICD-10 code E10) or type 2 diabetes
      (ICD-10 code E11)
  (3) The date of onset, estimated as the mid-point between that record and the last diabetes free record (or UK Biobank
      assessment if the first occuring incident record had uncertain diabetes codes E13 or E14), was within 10 years of
      follow-up and occurred prior to the 1st Februrary 2020.

This column contains missing data (N=26,594) where:

  (1) The participant has withdrawn consent for hospital record linkage (N=158)
  (2) OR, the participant is adjudicated as having prevalent diabetes (excluding gestational) at UK Biobank assessment
      ('prevalent_diabetes') i.e., type 2 diabetes, type 1 diabetes, or unspecified diabetes (N=26,445)

Notably, this column contains FALSE where the participant's first occuring incident diabetes record was for type 1
diabetes (N=221), or type 2 diabetes (N=13,555), in addition to those that were diabetes free within the available
follow-up (N=460,910).

-----------------------
'incident_any_diabetes'
-----------------------
This column curates information about whether each participant developed any form of diabetes over the 10 years of
follow-up.

This column contains 'TRUE' (N=14,955) where:

  (1) The participant had any hospital record or death record with ICD-10 codes E10, E11, E13, or E14 as primary or
      secondary outcome in any of the available follow-up data.
  (2) The date of onset, estimated as the mid-point between that record and the last diabetes free record (or UK Biobank
      assessment if the first occuring incident record had diabetes codes E10, E11, E13, or E14), was within 10 years of
      follow-up and occurred prior to the 1st Februrary 2020.

This column contains missing data (N=26,594) where:

  (1) The participant has withdrawn consent for hospital record linkage (N=158)
  (2) OR, the participant is adjudicated as having prevalent diabetes (excluding gestational) at UK Biobank assessment
      ('prevalent_diabetes') i.e., type 2 diabetes, type 1 diabetes, or unspecified diabetes (N=26,445)

----------------------
'incident_censor_date'
----------------------
This column contains the date at maximum follow-up for each participant. This is either:

  (1) the date of diabetes onset
  (2) the date the participant is lost to follow-up:
      (1) Due to death unrelated to diabetes, see 'death_at_censor_date' below
      (2) OR, the participant is lost to follow-up, e.g. due to leaving the UK, see 'lost_to_followup' and
          'lost_to_followup_reason' below.
  (3) the 1st February 2020
  (4) the date at 10 years after UK Biobank assessment, at which the participant remains diabetes free
  (5) the date at maximum follow-up available (e.g. < 10 years for participants in Wales; see 'censor_hospital_nation'),
      at which the participant remains diabetes free

This column contains missing data (N=26,594) where:

  (1) The participant has withdrawn consent for hospital record linkage (N=158)
  (2) OR, the participant is adjudicated as having prevalent diabetes (excluding gestational) at UK Biobank assessment
      ('prevalent_diabetes') i.e., type 2 diabetes, type 1 diabetes, or unspecified diabetes (N=26,445)

-----------------------
'incident_censor_years'
-----------------------
This column contains the follow-up time in years between the 'assessment_date' and 'incident_censor_date'. This is the
time in calendar years between UK Biobank assessment and either:

  (1) the date of type 2 diabetes onset
  (2) the date the participant is lost to follow-up:
      (1) Due to death unrelated to diabetes, see 'death_at_censor_date' below
      (2) OR, the participant is lost to follow-up, e.g. due to leaving the UK, see 'lost_to_followup' and
          'lost_to_followup_reason' below.
  (3) the 1st February 2020
  (4) the date at 10 years after UK Biobank assessment, at which the participant remains diabetes free
  (5) the date at maximum follow-up available (e.g. < 10 years for participants in Wales; see 'censor_hospital_nation'),
      at which the participant remains diabetes free

The numeric value corresponds to fractional calendar years. I.e., a whole number always indicates the same calendar
day + month in the corresponding year. E.g., 3 years after 2009-05-27 will be 2011-05-27. When UK Biobank assessment
occurred on a leap day, e.g. 2008-02-29, the <N> years after is set as the 28th February if the corresponding year is
not a leap year.

This column contains missing data (N=26,594) where:

  (1) The participant has withdrawn consent for hospital record linkage (N=158)
  (2) OR, the participant is adjudicated as having prevalent diabetes (excluding gestational) at UK Biobank assessment
      ('prevalent_diabetes') i.e., type 2 diabetes, type 1 diabetes, or unspecified diabetes (N=26,445)

-------------------------
'incident_t2d_censor_age'
-------------------------
This column contains the age at incident disease censoring, as a whole number.

Here, the age calculation is based on the follow up time ('incident_censor_years') added to the *fractional* age of the
participant at UK Biobank assessment, estimated from year of birth (UK Biobank field #34) and birth month (UK Biobank
field #52) with day of birth set at the midpoint (15th) of the month (specific date of birth not made available by UK
Biobank for privacy reasons; fractional age therefore accurate to with ~15 days, or 0.05 of a year).

----------------------------
'incident_censor_reason'
----------------------------
This column contains the reason for censoring each participant at the given 'incident_censor_date': whether they reached
the maximum follow-up date (1st February 2020, or 6th March 2018 for hospitals in Wales) diabetes free, whether they
reached the maximum follow-up time (10 years) diabetes free, whether they were lost to follow-up, or whether they
developed diabetes (and which type) within the follow-up period. Where the participant developed diabetes, the reason
distinguished between cases identified from hospital records (HES) or death records. Cases where there were no hospital
records for the participant (e.g. due to no hospital admissions) are also flagged.

The possible values in this column, organised by type, are:

                                                            incident_censor_reason  Samples
---------------------------------------------------------------------------------- --------
No evidence of diabetes within available follow-up:
---------------------------------------------------------------------------------- --------
                        No evidence of diabetes (HES) | Max follow time (10 years)  316,911
                             No in-patient data (HES) | Max follow time (10 years)   49,368
    No evidence of diabetes (HES) | Max follow-date (2020-02-01; pandemic cut-off)   56,221
         No in-patient data (HES) | Max follow-date (2020-02-01; pandemic cut-off)    9,956
      No evidence of diabetes (HES) | Max follow-up available (Hospitals in Wales)    6,595
           No in-patient data (HES) | Max follow-up available (Hospitals in Wales)    1,315
 No evidence of diabetes (HES) | Lost to follow-up (see 'lost_to_followup_reason')      716
      No in-patient data (HES) | Lost to follow-up (see 'lost_to_followup_reason')      370
---------------------------------------------------------------------------------- --------
No evidence of diabetes, death within follow-up due to other causes:
---------------------------------------------------------------------------------- --------
                        No evidence of diabetes (HES) | Death Record (no diabetes)   19,009
                             No in-patient data (HES) | Death Record (no diabetes)      449
---------------------------------------------------------------------------------- --------
Incident type 2 diabetes
---------------------------------------------------------------------------------- --------
                                                                Incident T2D (HES)   13,534
                       No evidence of diabetes (HES) | Incident T2D (Death Record)       20
                            No in-patient data (HES) | Incident T2D (Death Record)        1
---------------------------------------------------------------------------------- --------
Incident type 1 diabetes
---------------------------------------------------------------------------------- --------
                                                 Incident Uncertain Diabetes (HES)   13,534
        No evidence of diabetes (HES) | Incident Uncertain Diabetes (Death Record)        0
             No in-patient data (HES) | Incident Uncertain Diabetes (Death Record)        0
---------------------------------------------------------------------------------- --------
Incident diabetes of uncertain type
---------------------------------------------------------------------------------- --------
                                                                Incident T1D (HES)    1,148
                       No evidence of diabetes (HES) | Incident T1D (Death Record)       29
                            No in-patient data (HES) | Incident T1D (Death Record)        2
---------------------------------------------------------------------------------- --------

This column contains missing data (N=26,594) where:

  (1) The participant has withdrawn consent for hospital record linkage (N=158)
  (2) OR, the participant is adjudicated as having prevalent diabetes (excluding gestational) at UK Biobank assessment
      ('prevalent_diabetes') i.e., type 2 diabetes, type 1 diabetes, or unspecified diabetes (N=26,445)

------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------
Follow-up details
---------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
The final set of columns provide additional details on the prospective and retrospective follow-up available in hospital
and death records.

-----------------
'assessment_date'
-----------------
This column gives the date the participant attended the UK Biobank assessment centre for the corresponding
'visit_index', extracted from UK Biobank field #53.

-------------------
'assessment_centre'
-------------------
This column gives the location of the assessment centre attended, extracted from UK Biobank field #54. One of: "Barts",
"Birmingham", "Bristol", "Bristol (imaging)", "Bury", "Cardiff", "Cheadle (imaging)", "Cheadle (revisit)", "Croydon",
"Edinburgh", "Glasgow", "Hounslow", "Leeds", "Liverpool", "Manchester", "Middlesborough", "Newcastle",
"Newcastle (imaging)", "Nottingham", "Oxford", "Reading", "Reading (imaging)", "Sheffield", "Stockport (pilot)",
"Stoke", "Swansea", or "Wrexham". Note some locations are listed separately depending on assessment visit, e.g.
"Bristol" vs. "Bristol (imaging)".

-------------------
'assessment_nation'
-------------------
This column gives the nation within the UK the assessment centre is located, curated by Google Maps lookup of the
corresponding 'assessment_centre'. One of "England", "Wales", or "Scotland". Note the vast majority of participants
are located in England at baseline assessment (N=445,814; 88.8%).

----------------------
'death_at_censor_date'
----------------------
This column contains 'TRUE' where the maximum date of follow-up 'incident_censor_date' is due to participant death
for any cause (i.e. all-cause mortality).

Note that this column does not capture deaths due to diabetes (either as primary or secondary cause), as their date of
onset in 'incident_censor_date' will be earlier than the date of death due to the between records mid-point estimation
of date of onset (see above).

This column contains missing data (N=26,594) where:

  (1) The participant has withdrawn consent for hospital record linkage (N=158)
  (2) OR, the participant is adjudicated as having prevalent diabetes (excluding gestational) at UK Biobank assessment
      ('prevalent_diabetes') i.e., type 2 diabetes, type 1 diabetes, or unspecified diabetes (N=26,445)

------------------------
'censor_hospital_nation'
------------------------
Nation of most likely residence at end of available follow-up, based on (1) the nation of reporting hospital for the
most recent hospital records for that participant, (2) the nation from which the death record was issued for that
participant (if there is one), and (3) nation of assessment centre for most recent UK Biobank assessment visit for that
participant.

This information is important, as the follow-up time available for hospital records differs dramaticaly depending on the
nation of the reporting hospital (and there are currently minor differences when it comes to death records too):

Hospital records:

   Nation  Maximum follow-up available
--------- ----------------------------
  England                   2021-03-31
 Scotland                   2021-05-07
    Wales                   2018-03-06
--------- ----------------------------

Death records:

   Nation                         Record source  Maximum follow-up available
--------- ------------------------------------- ----------------------------
  England  Death registry for England and Wales                   2021-03-20
    Wales  Death registry for England and Wales                   2021-03-20
 Scotland           Death registry for Scotland                   2021-03-23
--------- ------------------------------------- ----------------------------

I.e. this means that for incident disease analysis you will need to adjust for 'censor_hospital_nation == "Wales"' to
account for the shorter follow-up time available (which if not accounted for appears as a protective effect on type 2
diabetes risk).

Note in 99.0% of cases this field is identical to nation at UK Biobank assessment ('assessment_nation'):

 Nation at baseline assessment  Inferred nation of residence at follow-up  Samples
         ('assessment_nation')                 ('censor_hospital_nation')
------------------------------ ------------------------------------------ --------
                      Scotland                                    England    2,558
                       England                                      Wales      964
                         Wales                                    England      859
                       England                                   Scotland      436
                         Wales                                   Scotland       17
------------------------------ ------------------------------------------ --------

This column contains missing data (N=26,594) where:

  (1) The participant has withdrawn consent for hospital record linkage (N=158)
  (2) OR, the participant is adjudicated as having prevalent diabetes (excluding gestational) at UK Biobank assessment
      ('prevalent_diabetes') i.e., type 2 diabetes, type 1 diabetes, or unspecified diabetes (N=26,445)

------------------
'lost_to_followup'
------------------
This column indicates ('TRUE') whether the participant is known to be lost to follow-up by UK Biobank; as recorded in
UK Biobank fields #190 and #191.

Note this column contains 'TRUE' where participants have withdrawn consent for linkage to hospital records (N=158) (see
'lost_to_followup_reason' below)

This column contains missing data (N=26,436) where the participant is adjudicated as having prevalent diabetes
(excluding gestational) at UK Biobank assessment ('prevalent_diabetes') i.e., type 2 diabetes, type 1 diabetes, or
unspecified diabetes (N=26,445), but has not withdrawn consent for linkage to hospital records.

-----------------------
'lost_to_followup_date'
-----------------------
Date participant is recorded as being lost to follow-up by UK Biobank, obtained from UK Biobank field #191 (with the
exception of participants with withdrawn consent for hospital and death record linkage, which we assume applies to all
such data, not just those occurring after the recorded date in field #191).

Data is missing (N=501,210) where:

 (1) the participant is not recorded as lost to follow-up ('FALSE' in 'lost_to_followup'),
 (2) has withdrawn consent for hospital and death record linkage (see 'lost_to_followup_reason' below),
 (3) or where there was no incident follow-up due to prevalent diabetes at UK Biobank assessment (see
     'prevalent_diabetes' above)

-------------------------
'lost_to_followup_reason'
-------------------------
Reason recorded by UK Biobank that the participant is lost to follow-up. One of:

                            'lost_to_followup_reason'  Samples
----------------------------------------------------- --------
 Participant has withdrawn consent for future linkage      158
           NHS records indicate they have left the UK      565
      UK Biobank sources report they have left the UK      491
           Death reported to UK Biobank by a relative       32
      NHS records indicate they are lost to follow-up        3
----------------------------------------------------- --------

Data is missing (N=501,210) where:

 (1) the participant is not recorded as lost to follow-up ('FALSE' in 'lost_to_followup'),
 (2) has withdrawn consent for hospital and death record linkage (see 'lost_to_followup_reason' below),
 (3) or where there was no incident follow-up due to prevalent diabetes at UK Biobank assessment (see
     'prevalent_diabetes' above)

------------------------
'earliest_hospital_date'
------------------------
As with incident disease follow-up in hospital records, retrospective follow-up time available also differs depending on
the nation of the reporting hospital:

   Nation  Minumum follow-up available
--------- ----------------------------
  England                   1993-07-27
 Scotland                   1980-12-02
    Wales                   1991-04-18
--------- ----------------------------

I.e. this determines how far back hospital records are available when defining prevalent diabetes and disease history
(columns 'history_cvd', 'history_gestational_diabetes', 'history_pcos', 'history_learning_difficulties',
'history_bipolar_schizophrenia', 'prevalent_diabetes', 'no_history_any_diabetes', 'type_2_diabetes', 'type_1_diabetes',
'uncertain_diabetes', and 'diabetes_history_reason').

--------------------------
'earliest_hospital_nation'
--------------------------
Nation of most likely residence at start available follow-up in hospital records, based on (1) the nation of reporting
hospital for the earliest hospital records for that participant, and (2) the nation of assessment centre for at baseline
assessment for that participant.

Note in 99.5% of cases this field is identical to nation at UK Biobank assessment ('assessment_nation'):

 Nation at baseline assessment  Inferred nation of residence at earliest hospital records  Samples
         ('assessment_nation')                               ('earliest_hospital_nation')
------------------------------ ---------------------------------------------------------- --------
                       England                                                      Wales    1,118
                       England                                                   Scotland      838
                         Wales                                                    England      549
                      Scotland                                                    England       46
                         Wales                                                   Scotland       23
------------------------------ ---------------------------------------------------------- --------
