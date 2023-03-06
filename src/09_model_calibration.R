library(data.table)
library(ukbnmr)
library(foreach)
library(survival)
library(ggplot2)
library(ggthemes)
source("src/utils/calibration.R")
source("src/utils/factor_by_size.R")

# Create output directory
system("mkdir -p analyses/test/calibration")

# Load biomarker information sheets
bio_info <- fread("data/ukb/biomarkers/output/biomarker_info.txt")
nmr_info <- ukbnmr::nmr_info

# Load test dataset
test <- fread("data/processed/test/processed_test_data.txt")

# Code factors, using non-risk/lower-risk group as reference
test[, sex := factor(sex, levels=c("Female", "Male"))]
test[, diabetes := factor(diabetes, levels=c("FALSE", "TRUE"))]
test[, smoking := factor(smoking, levels=c("FALSE", "TRUE"))]
test[, family_history_cvd := factor(family_history_cvd, levels=c("FALSE", "TRUE"))]
test[, assessment_centre := factor_by_size(assessment_centre)]
test[, earliest_hospital_nation := factor_by_size(earliest_hospital_nation)]
test[, latest_hospital_nation := factor_by_size(latest_hospital_nation)]

# Load information about models
model_info <- fread("analyses/test/model_fit_information.txt")

# Check model calibration
for (midx in model_info[,.I]) {
  this_model <- model_info[midx]
  ggdt <- calibration.fit(this_model$formula, test, 10)

  g <- ggplot(ggdt) + 
    aes(x=obs.risk, xmin=obs.L95, xmax=obs.U95,
        y=pred.risk, ymin=pred.L95, ymax=pred.U95,
        color=factor(risk_decile)) +
    geom_abline(linetype=2, slope=1, intercept=0) +
    geom_errorbarh(height=0) + geom_errorbar(width=0) +
    geom_point(shape=19) + 
    scale_color_tableau(name="Risk decile") + 
    facet_wrap(~ sex, scales="free") +
    xlab("Observed absolute risk") +
    ylab("Predicted absolute risk") +
    theme_bw()
  ggsave(g, width=6, height=3, file=sprintf("analyses/test/calibration/%s%s%s.pdf",
    gsub(" \\+? ?", "_", this_model$name), ifelse(this_model$PGS, "_with_PGS", ""),
    ifelse(this_model$lambda == "", "", paste0("_", this_model$lambda))))
   
}

