# Cardiovascular risk prediction using metabolomic biomarkers and polygenic risk scores: A cohort study and modelling analyses

This repository houses and documents the code used to generate the results in the study Ritchie SC et al. Cardiovascular risk prediction using metabolomic biomarkers and polygenic risk scores: A cohort study and modelling analyses. medRxiv, doi: 10.1101/2023.10.31.23297859.

## Repository information

The purpose of this repository is to provide a public record of the source code (i.e. methods) used for the entitled manuscript. The source code provided in this repository has not been designed to regenerate the results as-is for third-parties. The scripts herein contain numerous hard-coded filepaths and rely on data that cannot be made publicly available through this repository. 

The analysis scripts used to generate the results in our paper are contained in the `src/` folder. The `ext_src/` folder contains additional scripts primarily used for extraction and curation of UK Biobank and other data used for this project. This delineation has been made for pragmatic purposes: the scripts stored in `ext_src/` are copies of scripts located elsewhere on our HPC cluster as they have been written for cross-project purposes (see `ext_src/README.txt` for details).

The data underlying this project are primarily from UK Biobank, obtained as part of UK Biobank project 30418. Access to UK Biobank data is subject to approval from the UK Biobank access committee. See https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access for further details.

## Software and versions used

The following software and versions were used to run these scripts:

- CentOS Linux release 7.9.2009 (Core) (HPC operating system)
- Slurm version 22.05.10 (HPC queue manager and job submission system)
- GNU bash version 4.2.46(2) (shell environment)
- docopts v0.6.3-rc1 commit 13f0bbc (part of some of the commandline tools in `ext_src/`) 
- PLINK v2.00a3LM AVX2 Intel (2 Mar 2021) (www.cog-genomics.org/plink/2.0/), aliased as plink2 (used to curate UKB genotype data and calculate PRS)
- UKBiobank ukbconv_lx (c) CTSU. Compiled Mar 14 2018 (Used for converting UK Biobank data to csv)
- R version 4.2.0 (2022-04-22), along with the R packages:
  - Data wrangling:
    - data.table version 1.14.8
    - bit64 version 4.0.5
    - forcats version 1.0.0
    - lubridate version 1.9.2
    - openxlsx version 4.2.5
  - Data quality control:
    - MASS version 7.3-60
    - ukbnmr version 2.2.1
  - Programming:
    - docopt version 0.7.1
    - foreach version 1.5.2
    - doMC version 1.3.8
  - Statistics:
    - survival version 3.3-1
    - boot version 1.3-28
    - impute version 1.70.0
    - glmnet version 4.1-6
    - caret version 6.0-92
    - nricens version 1.6
  - Visualisation:
    - ggplot2 version 3.4.2
    - ggh4x version 0.2.4
    - ggpp version 0.5.2
    - ggrastr version 1.0.1
    - ggstance version 0.3.6
    - ggthemes version 4.2.4
    - scales version 1.2.0
    - RColorBrewer version 1.1-3
    - cowplot version 1.1.1

Inkscape version 1.2 was used to layout and annotate figures from the figure components generated within the R scripts. Microsoft Office for Mac (2019 edition) was used to draft the manuscript (Microsoft Word) and curate supplemental tables (Microsoft Excel) on MacOS Ventura 13.6
