library(data.table)
library(impute)

system("mkdir -p data/imputed/")

dat <- fread("data/cleaned/analysis_cohort.txt")

# Set missing smoking status (prefer not to answer, or otherwise NA) to FALSE.
# This is consistent with SCORE2, which categorizes smoking status into "current" or "other"
dat[is.na(smoking), smoking := FALSE]

# extract the non-derived nmr biomarkers for which we'll impute missing data
nmr_info <- fread("data/ukb/NMR_metabolomics/biomarker_information.txt")
nmr_dat <- dat[,.SD,.SDcols=c("eid", nmr_info[Type == "Non-derived", Biomarker])]

# Drop all nmr biomarkers (ratios included) from the main dataset
dat[, c(nmr_info[,Biomarker]) := NULL]

# Impute the missing data using K-nearest neighbours
nmr_dat <- as.matrix(nmr_dat, rownames="eid")
nmr_dat <- impute::impute.knn(t(nmr_dat)) # Expects samples as columns
nmr_dat <- t(nmr_dat$data)
nmr_dat <- as.data.table(nmr_dat, keep.rownames="eid")
nmr_dat <- ukbnmr::recompute_derived_biomarkers(nmr_dat)
nmr_dat[, eid := as.integer(eid)]
dat <- dat[nmr_dat, on = .(eid)]

# Write out
fwrite(dat, sep="\t", quote=FALSE, file="data/imputed/analysis_cohort.txt")

