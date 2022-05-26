library(data.table)
library(NetRep)
source("src/utils/factor_by_size.R")

# Make output directory
system("mkdir -p analyses/test/feature_correlations/")

# Load test data
test <- fread("data/processed/test/processed_test_data.txt")

# Code factors, using non-risk/lower-risk group as reference
test[, sex := factor(sex, levels=c("Female", "Male"))]
test[, diabetes := factor(diabetes, levels=c("FALSE", "TRUE"))]
test[, smoking := factor(smoking, levels=c("FALSE", "TRUE"))]
test[, family_history_cvd := factor(family_history_cvd, levels=c("FALSE", "TRUE"))]
test[, assessment_centre := factor_by_size(assessment_centre)]
test[, earliest_hospital_nation := factor_by_size(earliest_hospital_nation)]
test[, latest_hospital_nation := factor_by_size(latest_hospital_nation)]

# Load hazard ratios so we know what features have been selected
hrs <- fread("analyses/test/hazard_ratios.txt")

# Load model info
model_info <- fread("analyses/test/model_fit_information.txt")

# Plot correlations between model coefficients
for (midx in model_info[,.I]) {
  this_model <- model_info[midx]

  this_hrs <- hrs[long_name == this_model$long_name]

	# Split out multi-level factors to multiple columns, then recode as integer
	d1 <- model.matrix(~ 1 + ., test[,.SD,.SDcols=unique(this_hrs$var)])
	d1 <- d1[,-1] # drop intercept - we just used it because we don't want one-hot coding here.
  c1 <- cor(d1, method="spearman", use="pairwise.complete.obs")
  n1 <- abs(c1)^6
  m1 <- rep(1, ncol(c1))
  names(m1) <- colnames(c1)

  pdf(width=14, height=6, file=sprintf("analyses/test/feature_correlations/%s%s%s.pdf",
    tolower(gsub(" \\+? ?", "_", this_model[, name])),
    ifelse(this_model[,PGS], "_PGS", ""),
    ifelse(this_model[,lambda] == "", "", paste0("_", this_model[,lambda]))
  ))
  par(mar=c(15, 7, 2, 2))
  plotCorrelation(data=NULL, correlation=c1, network=n1, moduleAssignments=m1)
  dev.off()
}
