History of family illness
--------------------------

Data curated from fields #20107 (illness of father), #20110 (illness of mother), 
and #20111 (illness of siblings).

Individual diseases have been split out into separate columns based on available 
answers: (see https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=1010)

prostate_cancer: "Prostate cancer"
depression: "Severe depression"
parkinsons: "Parkinson's disease"
alzheimers: "Alzheimer's disease/dementia"
diabetes: "Diabetes"
hypertension: "High blood pressure"
copd: "Chronic bronchitis/emphysema"
breast_cancer: "Breast cancer"
bowel_cancer: "Bowel cancer"
lung_cancer: "Lung cancer"
stroke: "Stroke"
heart_disease: "Heart disease"

Note "breast_cancer" is always FALSE for illness of father and "prostate_cancer" is
always false for illness of mother as these options were not presented on the 
respective touchscreen questions.

Additional combined diseases have also been defined:

cardiovascular_disease: answered yes to either Heart disease or Stroke
any_cancer: answered yes to any of the cancers 

Data are missing where the participant answered either "Do not know" or 
"Prefer not to answer". These answers were available by two disease 
groups (see Notes tab for fields above):

Group 1 : Heart disease, Stroke, High blood pressure, Chronic bronchitis/emphysema, 
          Alzheimer's disease/dementia, Diabetes.

Group 2 : Parkinson's disease, Severe Depression, Lung cancer, Bowel cancer, 
          Prostate cancer.

E.g. a "prefer not to answer (group 2)" only resulted in missing data for diseases
listed in group 2 above, not those in group 1.

Another combined field, 'no_disease' is also provided, containing TRUE where a 
participant answered "None of the above" to both group 1 and group 2 diseases.

I've also curated files for "illness of parents" and "illness of first degree relatives",
which aggregate yes/no/NA answers across father+mother, and father+mother+siblings fields
respectively.


