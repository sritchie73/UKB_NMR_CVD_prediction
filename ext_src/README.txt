This folder houses copies of scripts used in symlinked folders under data/, e.g. those used to extract UK Biobank data 
variables, curate them, and derive new variables (e.g. QDiabetes risk scores). 

The reason these are separate is that each of the symlinked folders is its own self-contained "project" - i.e. with its
own set of (git-tracked) source files, data directories (again with symbolic links to raw UK biobank data, or other 
curated UK Biobank datasets). This is so they can be re-used across multiple projects, e.g. symlinked under a data 

Any calls to 'ukbconv' or other similar tools stored in a src/ukbtools/ folder (or folder
folder as I have done here.

Folders tracked here:
----------------------

 - algorithmically_defined_outcomes/: scripts and information used to extract the UK Biobank algorithmically defined
   outcomes (see respective README.md)

 - anthropometrics/: scripts and information used to extract and curate basic anthropometrics data for UK Biobank 
   participants, e.g. age, sex, bmi, and other related information such as assessment visit data.

 - biomarkers/: scripts and information used to extract and curate clinical biochemistry biomarker data for UK Biobank
   participant.

 - blood_pressure/: scripts and information used to extract blood pressure measurements and harmonize to a single
   measure per UK Biobank participant.

 - Eastwood_Diabetes/: scripts used to curate diabetes status from self-report touchscreen survey and verbal interview 
   data and incident diabetes from hospital records following the algorithms described in Eastwood et al. Algorithms 
   for the Capture and Adjudication of Prevalent and Incident Diabetes in UK Biobank. PLoS One 11, e0162388 (2016).

 - endpoints/: program used to define custom endpoints from ICD-10 codes, ICD-9 codes, OPCS-4 codes, OPCS-3 codes, 
   touchscreen survey medical and mental health history, and verbal interview questions on medical history, along with
   scripts used to extract and curate said underlying data. Includes definition of incident CVD and prevalent vascular
   disease. Note that there are also scripts under this folder in 'data_curation_scripts/' for curating self-report
   data, hospital episode statistics, death records, and per-participant follow-up information. 

 - medication/: scripts and information used to extract and curate medication history from UK Biobank touchscreen survey
   and verbal interview questions.

 - NMR_metabolomics/: scripts and information used to extract the NMR metabolomics biomarkers and adjust them for 
   technical variation using the ukbnmr package.

 - PGS_resources/: tool for downloading PGS from the PGS catalog and general purpose pipeline I've developed for 
   computing polygenic score levels from genotype data and one or more files of varaint weights. The scripts have 
   several hard coded aspects (to make them runnable by others in the location on the cluster they live in).

 - smoking/: scripts and information used to extract and curate smoking status and related variables from UK Biobank
   touchscreen survey data.

 - ukb_geno_conversion/: scripts used to convert UKB imputed genotype data (BGEN format) to plink pgen format, with 
   basic QC filtering of duplicate variants (by pos and alleles), and linkage to specific project sample identifiers.

