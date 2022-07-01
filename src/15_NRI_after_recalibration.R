library(data.table)
library(nricens)
source("src/utils/aki_absrisk.R")
source("src/utils/factor_by_size.R")
source("src/utils/risk_recalibration.R")

# Make output directory
system("mkdir -p analyses/public_health_modelling/net_reclassification")

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

# Reapply filtering used for risk recalibration - we need to redo this
# for pairwise comparisons after fitting Cox models in the subset of
# samples with non-missing data for each pair of models being compared
# in NRI analysis
test[, age_group := age %/% 5 * 5]
test <- test[age >= 40 & age < 70]
test <- test[(incident_cvd) | incident_followup == 10]

# Load model information
model_info <- fread("analyses/test/model_fit_information.txt")

# Wrapper function to compute continuous/categorical NRI 
nri.test <- function(base_model, new_model) {
  if (base_model == new_model) {
    # Nothing to do here, both models are identical
    return(NULL)
  }

  base_model <- as.formula(base_model)
  new_model <- as.formula(new_model)

  # Fit Cox proportional hazards models - in the first instance this is
  # just to find the set of samples that have non-missing data for both
  # models.
  base_cph <- coxph(base_model, data=test, x=TRUE)
  new_cph <- coxph(new_model, data=test, x=TRUE)

  # Now we re-fit the cox models in the samples with non-missing data for
  # both models so that the results are directly comparable for NRI analysis
  shared_rows <- as.integer(intersect(rownames(base_cph$x), rownames(new_cph$x)))
  base_cph <- coxph(base_model, data=test[shared_rows], x=TRUE)
  new_cph <- coxph(new_model, data=test[shared_rows], x=TRUE)

  # Extract predicted 10-year risk
  pred_risk <- test[shared_rows, .(eid, age, sex, incident_cvd, incident_followup)]
  pred_risk[as.integer(rownames(base_cph$x)), base_cph_risk := Coxar(base_cph, 10)]
  pred_risk[as.integer(rownames(new_cph$x)), new_cph_risk := Coxar(new_cph, 10)]

  # Recalibrate risk to sex-specific 5-year age group incidence rates from CPRD
  pred_risk[, base_cph_recalibrated := recalibrate_risk(base_cph_risk, eid, age, sex, male="Male")$recalibrated_risk]
  pred_risk[, new_cph_recalibrated := recalibrate_risk(new_cph_risk, eid, age, sex, male="Male")$recalibrated_risk]

  # Calculate continuous NRI analysis
  contNRI <- nricens(event = pred_risk$incident_cvd, time = pred_risk$incident_followup,
                     p.std = pred_risk$base_cph_recalibrated, p.new = pred_risk$new_cph_recalibrated,
                     updown = "diff", cut = 0, t0 = 10, niter = 1000)

  # Calculate categorical NRI using NICE 2014 cutoffs
  niceNRI <- nricens(event = pred_risk$incident_cvd, time = pred_risk$incident_followup,
                     p.std = pred_risk$base_cph_recalibrated, p.new = pred_risk$new_cph_recalibrated,
                     updown = "category", cut = c(0.05, 0.1), t0 = 10, niter = 1000)

  # Calculate cotegorical NRI using AHA/ACC 2019 cutoffs
  ahaaccNRI <- nricens(event = pred_risk$incident_cvd, time = pred_risk$incident_followup,
                       p.std = pred_risk$base_cph_recalibrated, p.new = pred_risk$new_cph_recalibrated,
                       updown = "category", cut = c(0.05, 0.075), t0 = 10, niter = 1000)

  # Add in sample size and case numbers
  contNRI$n <- niceNRI$n <- ahaaccNRI$n <- pred_risk[,.N]
  contNRI$nevent <- niceNRI$nevent <- ahaaccNRI$nevent <- pred_risk[, sum(incident_cvd)]

  list("Continuous NRI"=contNRI, "NICE 2014 Categorical NRI"=niceNRI, "ACC/AHA 2019 Categorical NRI"=ahaaccNRI)
}

# Run NRI analysis 
nri_list <- list(
  "Base to PGS" = nri.test(model_info[!(PGS) & name == "Conventional RF", formula], model_info[(PGS) & name == "Conventional RF", formula]),
  "Base to CRP" = nri.test(model_info[!(PGS) & name == "Conventional RF", formula], model_info[!(PGS) & name == "CRP", formula]),
  "Base to PGS + CRP" = nri.test(model_info[!(PGS) & name == "Conventional RF", formula], model_info[(PGS) & name == "CRP", formula]),
  "Base to Assays (1se)" = nri.test(model_info[!(PGS) & name == "Conventional RF", formula], model_info[!(PGS) & name == "Assays" & lambda == "lambda.1se", formula]),
  "Base to NMR (1se)" = nri.test(model_info[!(PGS) & name == "Conventional RF", formula], model_info[!(PGS) & name == "NMR" & lambda == "lambda.1se", formula]),
  "Base to NMR + Assays (1se)" = nri.test(model_info[!(PGS) & name == "Conventional RF", formula], model_info[!(PGS) & name == "NMR + Assays" & lambda == "lambda.1se", formula]),
  "Base to PGS + Assays (1se)" = nri.test(model_info[!(PGS) & name == "Conventional RF", formula], model_info[(PGS) & name == "Assays" & lambda == "lambda.1se", formula]),
  "Base to PGS + NMR (1se)" = nri.test(model_info[!(PGS) & name == "Conventional RF", formula], model_info[(PGS) & name == "NMR" & lambda == "lambda.1se", formula]),
  "Base to PGS + NMR + Assays (1se)" = nri.test(model_info[!(PGS) & name == "Conventional RF", formula], model_info[(PGS) & name == "NMR + Assays" & lambda == "lambda.1se", formula]),
  "Base to Assays (min)" = nri.test(model_info[!(PGS) & name == "Conventional RF", formula], model_info[!(PGS) & name == "Assays" & lambda == "lambda.min", formula]),
  "Base to NMR (min)" = nri.test(model_info[!(PGS) & name == "Conventional RF", formula], model_info[!(PGS) & name == "NMR" & lambda == "lambda.min", formula]),
  "Base to NMR + Assays (min)" = nri.test(model_info[!(PGS) & name == "Conventional RF", formula], model_info[!(PGS) & name == "NMR + Assays" & lambda == "lambda.min", formula]),
  "Base to PGS + Assays (min)" = nri.test(model_info[!(PGS) & name == "Conventional RF", formula], model_info[(PGS) & name == "Assays" & lambda == "lambda.min", formula]),
  "Base to PGS + NMR (min)" = nri.test(model_info[!(PGS) & name == "Conventional RF", formula], model_info[(PGS) & name == "NMR" & lambda == "lambda.min", formula]),
  "Base to PGS + NMR + Assays (min)" = nri.test(model_info[!(PGS) & name == "Conventional RF", formula], model_info[(PGS) & name == "NMR + Assays" & lambda == "lambda.min", formula])
)
saveRDS(nri_list, file="analyses/public_health_modelling/net_reclassification/nri.rds")

# Extract tables of estimates
nri_estimates <- rbindlist(idcol="model", fill=TRUE, lapply(nri_list, function(l1) {
  rbindlist(idcol="nri_type", fill=TRUE, lapply(l1, function(l2) {
    cbind(samples=l2$n, cases=l2$nevent, as.data.table(l2$nri, keep.rownames="metric"))
  }))
}))

# Extract tables of reclassifications for categorical nris
reclassified_cases <- rbindlist(idcol="model", fill=TRUE, lapply(nri_list, function(l1) {
  rbindlist(idcol="nri_type", fill=TRUE, lapply(l1, function(l2) {
    as.data.table(l2[["rtab.case"]])
  }))
}))

reclassified <- rbindlist(idcol="model", fill=TRUE, lapply(nri_list, function(l1) {
  rbindlist(idcol="nri_type", fill=TRUE, lapply(l1, function(l2) {
    as.data.table(l2[["rtab"]])
  }))
}))

# Merge
setnames(reclassified, "N", "All")
reclassified[reclassified_cases, on = .(model, nri_type, Standard, New), Cases := N]
setnames(reclassified, "Standard", "Old")

# Add in feature selection lambda information
nri_estimates[, lambda := fcase(
  model %like% "(min)", "Best model",
  model %like% "(1se)", "Best model with fewest features",
  default = "No feature selection")]

reclassified[, lambda := fcase(
  model %like% "(min)", "Best model",
  model %like% "(1se)", "Best model with fewest features",
  default = "No feature selection")]

# Add in total sample size and total cases to reclassified table
reclassified[nri_estimates, on = .(model), Total_Samples := i.samples]
reclassified[nri_estimates, on = .(model), Total_Cases := i.cases]

fwrite(nri_estimates, sep="\t", quote=FALSE, file="analyses/public_health_modelling/net_reclassification/nri_estimates.txt")
fwrite(reclassified, sep="\t", quote=FALSE, file="analyses/public_health_modelling/net_reclassification/nri_reclassified.txt")
