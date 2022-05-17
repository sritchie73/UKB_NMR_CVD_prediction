Medication data from touchscreen survey and verbal interview
-------------------------------------------------------------

Data here are curated from the following UK Biobank fields:

    Field                                                                      Description
--------- --------------------------------------------------------------------------------
     6177                           Medication for cholesterol, blood pressure or diabetes
     6153 Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones
     2492                                            Taking other prescription medications
      137                                           Number of treatments/medications taken
    20003                                                        Treatment/medication code
     6671                                     Number of antibiotics taken in last 3 months
    20199                                               Antibiotic codes for last 3 months
--------- --------------------------------------------------------------------------------

These are split into several files:

- output/detailed_medications_field_20003.txt
- output/detailed_antibiotics_field_20199.txt
- output/medications_simple.txt
- output/medications_simple_info.txt
- output/detailed_medications_summariesed.txt

The file 'output/detailed_medications_field_20003.txt' contains a list of medications either 
currently being taken, or prescribed to, each participant at each assessment visit. This data were
obtained and curated by a trained nurse in verbal interview with the participant, and made available
in field 20003. The data curated here are in long format, with each row corresponding to a single 
medication for a single participant at a single assessment visit. A person may be taking multiple
medications, in which case there are multiple rows. The 'medication_code' column contains the 
medication code given in field 20003, while the 'medication_name' column contains the corresponding
label (https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=4). Rows with 'medication_code == -2',
and 'medication_name == "Not applicable"' indicate cases where the person was not taking any 
medications at that assessment visit ('0' in field #137, Number of treatments/medications taken).

The file 'output/detailed_antibiotics_field_20199.txt' is similarly formated, containing a list of 
antibiotics prescribed/taken in the 3 months prior to assessment. This data were obtained and 
curated by a trained nurse in verbal interview with the participant, and made available in field 
20199. However, unlike general medications in field 20003, this question was only asked at the 
imaging assessments ('visit_index' '2' or '3'). Rows with 'antibiotic_code == -2', and 
'antibiotic_name == "Not applicable"' likewise indicate cases where the person was not taking any
antiobitics in 3 months prior to that assessment visit ('0' in field #6671, Number of antibiotics 
taken in last 3 months)

The file 'output/medications_simple.txt' contains data curated from the other fields, and unlike the
previous two files there is only one row per participant and assessment visit. Here, the 
'cholesterol_medication', 'blood_pressure_medication', 'insulin_medication', 
'hormone_replacement_therapy', and 'oral_contraceptive' correspond to possible answers to the 
multiple-choice, multiple-answer, touchscreen survey questions 6153 and 6177. These columns contain
TRUE/FALSE values indicating whether the respective option was selected, or missing where the 
participant answered "Do not know" or "Prefer not to answer" to fields 6153 or 6177. The 
'hormone_replacement_therapy' and 'oral_contraceptive' are always FALSE for men, as these options
were not listed in the general medication touchscreen question presented to self-identifying males
(field 6177; field 6153 was the survey question presented to self-identifying females). The 
'other_prescription_meds' is TRUE where the participant answered "Yes" to the "Taking other 
prescription medications" question (field 2492) (in which case these may be found in 
'output/detailed_medications_field_20003.txt'). The 'num_current_meds' column gives the number of 
medications being currently taken (field 137); note this may be > 0 where 'other_prescription_meds'
is FALSE, e.g. where the participant is taking non-prescription/over-the-counter medication(s) (e.g.
painkillers) which are also listed in 'output/detailed_medications_field_20003.txt'. Finally, the 
column 'num_antibiotics_last_3mo' gives the number of antibiotics taken in the 3 months prior to 
assessment - although note this data is not available at baseline assessment (visit_index == 0) or
first repeat assessment (visit_index == 1) as the corresponding verbal interview questions were not
asked (see above). Data in these columns are also missing (NA) where a participant answered 
"Do not know" or "Prefer not to answer".

Descriptions of the columns in 'output/medications_simple.txt' can also be viewed in 
'output/medications_simple_info.txt'

In all datasets, here and in other folders, the 'eid' column contains the participant identifier in
project 7439, and the 'visit_index' column indicates the assessment visit: '0' for baseline 
assessment (first visit; 2006-2010; N=502,460 participants), '1' for first repeat assessment (second 
visit; 2009-2014; N=20,344 participants), '2' for imaging assessment (second/third visit; 2014-2020; 
N=48,998 participants, of which 8,284 had measurements at first repeat assessment), and '3' for 
for repeat imaging assessment (third/fourth visit; 2019-; N=4,499 participants).

The file 'output/detailed_medications_summarised.txt' contains columns summarising presence/absence
of various lists of medications found in 'output/detailed_medications_field_20003.txt'. 

Currently, these are designed to match the columns in 'output/medications_simple.txt':

 - 'lipid_lowering_medication' gives TRUE or FALSE where the participant is taking any type of lipid 
   lowering medication, analogous to the 'cholesterol_medication' column in 'output/medications_simple.txt'

 - 'hypertension_medication' gives TRUE or FALSE where the participant is taking any medications to 
   treat hypertension, analogous to the 'blood_pressure_medication' column in 'output/medications_simple.txt'

 - 'insulin_medication' gives TRUE or FALSE where the participant is taking any insulin products, analogous
   to the 'insulin_medication' column in 'output/medications_simple.txt'

The lists of medications contributing to each column have been curated based on the respective pages
in the British National Formulary (BNF) at https://openprescribing.net/bnf/. 

---------------------------
Lipid Lowering Medications
---------------------------

List of drugs from BNF chapter 2.12: Lipid-regulating drugs https://openprescribing.net/bnf/0212/
that overlap with those in the UK Biobank medication code list https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=4

Note that not all drugs listed in the BNF are found in UKB, and some need to be looked up by components of their
name. All drugs have been cross-checked against their indications using https://bnf.nice.org.uk/drug/ to ensure 
that they capture the intended treatment effect (in this instance, lipid lowering to treat hypercholesterolemia to 
prevent CVD)

For example, while Ispaghula husk is listed, its primary indication is for constipation, not lipid lowering, so its 
not included. https://bnf.nice.org.uk/drug/ispaghula-husk.html

Conversely, colestyramine is listed, even though it has multiple indications, because one of its main indications is 
for primary prevention of coronary heart disease. https://bnf.nice.org.uk/drug/colestyramine.html

 UK Biobank Code                 UK Biobank medication name 
---------------- ------------------------------------------
      1140861892                                   acipimox
      1141146234                               atorvastatin
      1140861924                                bezafibrate
      1141157260                        bezafibrate product
      1140862026                               ciprofibrate
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

---------------------------
Hypertension Medications
---------------------------

List of drugs from BNF chapter 2.5: Hypertension and heart failure https://openprescribing.net/bnf/0205/
that overlap with those in the UK Biobank medication code list https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=4

Note that not all drugs listed in the BNF are found in UKB, and some need to be looked up by components of their
name. All drugs have been cross-checked against their indications using https://bnf.nice.org.uk/drug/ to ensure 
that they capture the intended treatment effect (in this instance, to treat hypertension as opposed to heart failure) 


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

-------------------
Insulin Medication
-------------------

List of drugs from BNF chapter 6.1.1: Insulin: https://openprescribing.net/bnf/060101/
that overlap with those in the UK Biobank medication code list https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=4

Note that not all drugs listed in the BNF are found in UKB, and some need to be looked up by components of their
name. All drugs have been cross-checked against their indications using https://bnf.nice.org.uk/drug/ to ensure 
that they capture the intended treatment effect (in this instance, insulin to treat type 1 diabetes)

 UK Biobank Code                 UK Biobank medication name 
---------------- ------------------------------------------
      1140883066                            insulin product
---------------- ------------------------------------------
