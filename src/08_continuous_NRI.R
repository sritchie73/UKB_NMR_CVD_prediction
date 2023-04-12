library(data.table)
library(nricens)
source("src/utils/aki_absrisk.R")
source("src/utils/factor_by_size.R")

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

# Load model information
model_info <- fread("analyses/test/model_fit_information.txt")

# Filter to complete data so all model comparisons use the same samples
test <- test[(complete_data)]

# Wrapper function for continuous NRI test
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
  pred_risk <- test[shared_rows, .(eid, incident_cvd, incident_followup)]
  pred_risk[as.integer(rownames(base_cph$x)), base_cph_risk := Coxar(base_cph, 10)]
  pred_risk[as.integer(rownames(new_cph$x)), new_cph_risk := Coxar(new_cph, 10)]

  # Run NRI analysis
  # Originally I tried passing the fitted Cox models directly here, but a tonne
  # of warnings were generated and the bootstrap confidence intervals were very
  # weird (e.g. looking like:  -o------ or ------o instead of ---o---). 
  contNRI <- nricens(event = pred_risk$incident_cvd, time = pred_risk$incident_followup,
                     p.std = pred_risk$base_cph_risk, p.new = pred_risk$new_cph_risk,
                     updown = "diff", cut = 0, t0 = 10, niter = 1000) 

  # Add in sample size and case numbers
  contNRI$n <- pred_risk[,.N]
  contNRI$nevent <- pred_risk[, sum(incident_cvd)]
  
  return(contNRI)
}

# Run NRI analysis with 1000 bootstraps:
nri_list <- list(
  "Base to NMR" = nri.test(model_info[!(PGS) & name == "Conventional RF", formula], model_info[!(PGS) & name == "NMR" & lambda == "lambda.min", formula]),
  "Base to PGS" = nri.test(model_info[!(PGS) & name == "Conventional RF", formula], model_info[(PGS) & name == "Conventional RF", formula]),
  "Base to CRP" = nri.test(model_info[!(PGS) & name == "Conventional RF", formula], model_info[!(PGS) & name == "CRP", formula]),
  "Base to PGS + NMR" = nri.test(model_info[!(PGS) & name == "Conventional RF", formula], model_info[(PGS) & name == "NMR" & lambda == "lambda.min", formula]),
  "Base to PGS + CRP" = nri.test(model_info[!(PGS) & name == "Conventional RF", formula], model_info[(PGS) & name == "CRP", formula])
)

saveRDS(nri_list, file="analyses/test/nri.rds")

# Extract tables of estimates
nri_estimates <- rbindlist(idcol="model", fill=TRUE, lapply(nri_list, function(l1) {
  cbind(samples=l1$n, cases=l1$nevent, as.data.table(l1$nri, keep.rownames="metric"))
}))

# Add in feature selection lambda information
nri_estimates[, lambda := fcase(
  model %like% "(min)", "Best model",
  default = "No feature selection")]

fwrite(nri_estimates, sep="\t", quote=FALSE, file="analyses/test/nri_estimates.txt")
