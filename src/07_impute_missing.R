library(data.table)
library(impute)

system("mkdir -p data/imputed/")

dat <- fread("data/cleaned/analysis_cohort.txt")

# Set missing smoking status (prefer not to answer, or otherwise NA) to FALSE.
# This is consistent with SCORE2, which categorizes smoking status into "current" or "other"
dat[is.na(smoking), smoking := FALSE]

# extract the non-derived nmr biomarkers for which we'll impute missing data
nmr_info <- fread("data/ukb/NMR_metabolomics/biomarker_information.txt")
non_derived <- nmr_info[Type == "Non-derived", Biomarker]
non_derived <- setdiff(non_derived, "Clinical_LDL_C")
nmr_dat <- dat[,.SD,.SDcols=c("eid", non_derived)]

# Drop all nmr biomarkers (ratios included) from the main dataset
dat[, c(nmr_info[,Biomarker]) := NULL]

# Convert to matrix
nmr_dat <- as.matrix(nmr_dat, rownames="eid")

# Do first pass imputation with default KNN so we can run PCA to determine better number for K
pca_dat <- impute::impute.knn(t(nmr_dat)) # Expects samples as columns
pcs <- prcomp(t(pca_dat$data), scale=TRUE)
pc_info <- as.data.table(t(summary(pcs)$importance), keep.rownames="PC")
setnames(pc_info, c("PC", "SD", "Prop.Var", "Cumul.Prop.Var"))
fwrite(pc_info, sep="\t", quote=FALSE, file="data/imputed/pca.txt")

# Set K to number of PCs explaining >95% of the variation (K=20)
K <- pc_info[,which(Cumul.Prop.Var > 0.95)[1]]
nmr_dat <- impute::impute.knn(t(nmr_dat), k=20)
nmr_dat <- t(nmr_dat$data)
nmr_dat <- as.data.table(nmr_dat, keep.rownames="eid")
nmr_dat <- ukbnmr::recompute_derived_biomarkers(nmr_dat)
nmr_dat[, eid := as.integer(eid)]
dat <- dat[nmr_dat, on = .(eid)]

# Write out
fwrite(dat, sep="\t", quote=FALSE, file="data/imputed/analysis_cohort.txt")

