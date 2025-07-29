library(data.table)
library(impute)

system("mkdir -p data/imputed/")

dat <- fread("data/cleaned/analysis_cohort.txt")

# Set missing smoking status (prefer not to answer, or otherwise NA) to FALSE.
# This is consistent with SCORE2, which categorizes smoking status into "current" or "other"
dat[is.na(smoking), smoking := FALSE]

# extract the assays for which we'll impute missing data
assay_info <- fread('data/ukb/biomarkers/output/biomarker_info.txt')
assays <- assay_info[sample_type != "Urine" & !is.na(UKB.Field.ID) & var != "tchol" & var != "hdl", var]
assay_dat <- dat[,.SD,.SDcols=c("eid", assays)]

# Convert to matrix
assay_dat <- as.matrix(assay_dat, rownames="eid")

# Do first pass imputation with default KNN so we can run PCA to determine better number for K
pca_dat <- impute::impute.knn(t(assay_dat)) # Expects samples as columns
pcs <- prcomp(t(pca_dat$data), scale=TRUE)
pc_info <- as.data.table(t(summary(pcs)$importance), keep.rownames="PC")
setnames(pc_info, c("PC", "SD", "Prop.Var", "Cumul.Prop.Var"))
fwrite(pc_info, sep="\t", quote=FALSE, file="data/imputed/pca_assays.txt")

# Set K to number of PCs explaining >95% of the variation (K=23)
K <- pc_info[,which(Cumul.Prop.Var > 0.95)[1]]
assay_dat <- impute::impute.knn(t(assay_dat), k=23)
assay_dat <- t(assay_dat$data)
assay_dat <- as.data.table(assay_dat, keep.rownames="eid")
assay_dat <- ukbassay::recompute_derived_biomarkers(assay_dat)
assay_dat[, eid := as.integer(eid)]

# Write out
fwrite(assay_dat, sep="\t", quote=FALSE, file="data/imputed/assays.txt")

