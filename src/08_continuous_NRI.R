library(data.table)
library(nricens)
source("src/utils/aki_absrisk.R")
source("src/utils/factor_by_size.R")

# Make output directory
system("mkdir -p analyses/test/net_reclassification")

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

# Wrapper function for continuous NRI test
nri.test <- function(base_model, new_model) {
  base_model <- as.formula(base_model)
  new_model <- as.formula(new_model)

  # Fit Cox proportional hazards models
  base_cph <- coxph(base_model, data=test, x=TRUE)
  new_cph <- coxph(new_model, data=test, x=TRUE)

  # Extract predicted 10-year risk
  pred_risk <- test[, .(eid, incident_cvd, incident_followup)]
  pred_risk[as.integer(rownames(base_cph$x)), base_cph_risk := Coxar(base_cph, 10)]
  pred_risk[as.integer(rownames(new_cph$x)), new_cph_risk := Coxar(new_cph, 10)]
  
  # Filter to samples that both models could be fit for
  pred_risk <- pred_risk[!is.na(base_cph_risk) & !is.na(new_cph_risk)]

  # Run NRI analysis
  # Originally I tried passing the fitted Cox models directly here, but a tonne
  # of warnings were generated and the bootstrap confidence intervals were very
  # weird (e.g. looking like:  -o------ or ------o instead of ---o---). 
  nricens(event = pred_risk$incident_cvd, time = pred_risk$incident_followup,
          p.std = pred_risk$base_cph_risk, p.new = pred_risk$new_cph_risk,
          updown = "diff", cut = 0, t0 = 10, niter = 1000) 
}

# Run NRI analysis with 1000 bootstraps:
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

saveRDS(nri_list, file="analyses/test/net_reclassification/nri.rds")

# Extract tables of estimates
nri_estimates <- rbindlist(idcol="model", fill=TRUE, lapply(nri_list, function(l1) {
  cbind(samples=l1$mdl.std$n, cases=l1$mdl.std$nevent, as.data.table(l1$nri, keep.rownames="metric"))
}))

# Add in feature selection lambda information
nri_estimates[, lambda := fcase(
  model %like% "(min)", "Best model",
  model %like% "(1se)", "Best model with fewest features",
  default = "No feature selection")]

fwrite(nri_estimates, sep="\t", quote=FALSE, file="analyses/test/net_reclassification/nri_estimates.txt")
