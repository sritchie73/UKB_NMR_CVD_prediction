This folder houses copies of scripts used in symlinked folders under data/, e.g. those used to extract UK Biobank data 
variables, curate them, and derive new variables (e.g. QDiabetes risk scores). 

The reason these are separate is that each of the symlinked folders is its own self-contained "project" - i.e. with its
own set of (git-tracked) source files, data directories (again with symbolic links to raw UK biobank data, or other 
curated UK Biobank datasets). This is so they can be re-used across multiple projects, e.g. symlinked under a data 

Any calls to 'ukbconv' or other similar tools stored in a src/ukbtools/ folder (or folder
folder as I have done here.

Folders tracked here:
----------------------


 - anthropometrics/: scripts and information used to extract and curate basic anthropometrics data for UK Biobank 
   participants, e.g. age, sex, bmi, and other related information such as assessment visit data.

 - blood_pressure/: scripts and information used to extract blood pressure measurements and harmonize to a single
   measure per UK Biobank participant.

 - medication/: scripts and information used to extract and curate medication history from UK Biobank touchscreen survey
   and verbal interview questions.

 - smoking/: scripts and information used to extract and curate smoking status and related variables from UK Biobank
   touchscreen survey data.

 - family_history/: scripts and information used to extract and curated family history of disease from UK Biobank survey
   questions. 

 - QDiabetes/: scripts and information used to extract and curate data in UK Biobank relevant to QDiabetes risk scores,
   including prevalent diabetes status. Note that the data the QDiabetes folder imports includes several other extracted 
   + curated UK Biobank datasets, which have been copied into separate folders here (described below), including: 
   'anthropometrics/', 'biomarkers/', 'Eastwood_Diabetes/', 'endpoints/', 'family_history/', 'medication/', and 
   'smoking/'.

 - Eastwood_Diabetes/: IPython notebooks created by Dr. Sam Lambert (sl925@medschl.cam.ac.uk) to curate diabetes status
   from self-report touchscreen survey and verbal interview data and incident diabetes from hospital records following
   the algorithms described in Eastwood et al. Algorithms for the Capture and Adjudication of Prevalent and Incident 
   Diabetes in UK Biobank. PLoS One 11, e0162388 (2016). Used by the 'QDiabetes/' folder above (with additional 
   modifications, e.g. adjudication of prevalent diabetes additionally including hospital records, and incident diabetes
   restricted to 10 years of follow-up to match QDiabetes risk horizon).

 - endpoints/: program used to define custom endpoints from ICD-10 codes, ICD-9 codes, OPCS-4 codes, OPCS-3 codes, 
   touchscreen survey medical and mental health history, and verbal interview questions on medical history, along with
   scripts used to extract and curate said underlying data. Includes definition of incident CVD and prevalent vascular
   disease. Note that there are also scripts under this folder in 'data_curation_scripts/' for curating self-report
   data, hospital episode statistics, death records, and per-participant follow-up information. 

 - NMR_metabolomics/: scripts and information used to extract the NMR metabolomics biomarkers and adjust them for 
   technical variation using the ukbnmr package.

 - biomarkers/: scripts and information used to extract and curate clinical biochemistry biomarker data for UK Biobank
   participant.

 - PGS_resources/: tool for downloading PGS from the PGS catalog and general purpose pipeline I've developed for 
   computing polygenic score levels from genotype data and one or more files of varaint weights. The scripts have 
   several hard coded aspects (to make them runnable by others in the location on the cluster they live in).

 - ukb_geno_conversion/: scripts used to convert UKB imputed genotype data (BGEN format) to plink pgen format, with 
   basic QC filtering of duplicate variants (by pos and alleles), and linkage to specific project sample identifiers.

