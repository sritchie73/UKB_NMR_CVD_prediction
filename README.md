# Cardiovascular risk prediction using metabolomic biomarkers and polygenic risk scores: A cohort study and modelling analyses

This repository houses and documents the code used to generate the results in the study Ritchie SC et al. Cardiovascular risk prediction using metabolomic biomarkers and polygenic risk scores: A cohort study and modelling analyses. medRxiv, doi: XXXXXX.

## Disclaimer

This code has not been designed to regenerate the results as-is for third-parties. It contains hard-coded filepaths pointing towards UK Biobank data from a specific UK Biobank project number that has been extracted and processed in a particular way (more details below). To some extent, it also relies on the specific compute configuration of the HPC cluster at the University of Cambridge. 

## Project organisation

The source code written for this project is contained in the `src/` folder in this repository. The scripts contained in this folder contain numerous hard-coded filepaths to the `data/` folder, which is expected to contain the UK Biobank used for this project. The `data/` directory is not provided in this repository as we cannot make this publicly available. All data described are available through UK Biobank subject to approval from the UK Biobank access committee. See https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access for further details. 

## Source code organisation

Source code is provided in two folders: `src/` and `ext_src/`. The `src/` folder contains scripts specifically written for this project. The `ext_src/` folder contains copies of scripts used to 
