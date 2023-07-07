library(data.table)
library(ukbnmr)
library(foreach)
library(survival)
source("src/utils/cv_coxph.R")
source("src/utils/score_cindex.R")
source("src/utils/cox_test.R")
registerDoMC(10)

# Create output directory
system("mkdir -p analyses/univariate")

# Load biomarker information sheets
bio_info <- fread("data/ukb/biomarkers/output/biomarker_info.txt")
nmr_info <- ukbnmr::nmr_info

# Load analysis cohort
dat <- fread("data/cleaned/analysis_cohort.txt")

# Code factors, using non-risk/lower-risk group as reference
dat[, sex := factor(sex, levels=c("Female", "Male"))]
dat[, smoking := factor(smoking, levels=c("FALSE", "TRUE"))]

# Build set of models to test
test_nmr <- nmr_info[Type == "Non-derived", Biomarker]
test_assay <- bio_info[!is.na(UKB.Field.ID) & sample_type != "Urine" & var != "tchol" & var != "hdl", var]
models <- rbind(
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ scale(age)", type="demographics", name="age"),
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ scale(age) + smoking", type="risk_factors", name="smoking"),
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ scale(age) + scale(sbp)", type="risk_factors", name="sbp"),
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ scale(age) + scale(tchol)", type="risk_factors", name="tchol"),
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ scale(age) + scale(hdl)", type="risk_factors", name="hdl"),
  data.table(formula=sprintf("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + scale(%s)", test_nmr), type="NMR", name=test_nmr),
  data.table(formula=sprintf("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + scale(%s)", test_assay), type="assays", name=test_assay)
)

# Compute C-index in 10-fold cross-validation
foldid <- fread("analyses/CVD_score_weighting/cross_validation_fold_allocation.txt")
dat <- dat[foldid, on = .(eid)]

# Takes about an hour to run on 10 cores
## cinds <- foreach(this_sex=c("Male", "Female"), .combine=rbind) %:% 
##   foreach(mIdx = models[,.I], .combine=rbind) %do% {
##     this_model <- models[mIdx]
##     cbind(this_model, sex=paste0(this_sex, "s"), cv.cindex(this_model$formula, dat[sex == this_sex], "coxph_cv_foldid"))
## }
## fwrite(cinds, sep="\t", quote=FALSE, file="analyses/univariate/cindices.txt")

# Estimate hazards ratios
cox_list <- foreach(this_sex = c("Male", "Female")) %:% 
  foreach(mIdx = models[,.I]) %dopar% {
    this_model <- models[mIdx]
    cox.test(this_model$formula, "incident_cvd", dat[sex == this_sex])
}
saveRDS(cox_list, file="analyses/test/univariate_cox_list.rds")





