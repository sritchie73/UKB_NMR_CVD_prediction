library(data.table)
source("src/utils/nri.test.R")

# Make output directory
system("mkdir -p analyses/public_health_modelling/net_reclassification")

# Load model information
model_info <- fread("analyses/test/model_fit_information.txt")

# Load recalibrated predicted CVD risk
risk <- fread("analyses/public_health_modelling/risk_recalibration/absolute_risks.txt")

# Helper function to aid NRI tests:
cvd.nri.test <- function(m1_risk, m2_risk) {
  shared_eid <- data.table(eid=intersect(m1_risk$eid, m2_risk$eid))

  m1_risk <- m1_risk[shared_eid, on = .(eid)]
  m2_risk <- m2_risk[shared_eid, on = .(eid)]

  nri.test(m1_risk$recalibrated_risk, m2_risk$recalibrated_risk, 10, NULL, m1_risk$incident_cvd, m1_risk$incident_followup, 1000)
}

# Run NRI analysis with 1000 bootstraps:
nri_list <- list(
	"Base to PGS" = cvd.nri.test(risk[!(PGS) & name == "Conventional RF"], risk[(PGS) & name == "Conventional RF"]),
	"Base to Blood (1se)" = cvd.nri.test(risk[!(PGS) & name == "Conventional RF"], risk[!(PGS) & name == "Blood" & lambda == "lambda.1se"]),
	"Base to Nightingale (1se)" = cvd.nri.test(risk[!(PGS) & name == "Conventional RF"], risk[!(PGS) & name == "Nightingale" & lambda == "lambda.1se"]),
	"Base to Blood + Nightingale (1se)" = cvd.nri.test(risk[!(PGS) & name == "Conventional RF"], risk[!(PGS) & name == "Blood + Nightingale" & lambda == "lambda.1se"]),
	"Base + PGS to PGS + Blood (1se)" = cvd.nri.test(risk[(PGS) & name == "Conventional RF"], risk[(PGS) & name == "Blood" & lambda == "lambda.1se"]),
	"Base + PGS to PGS + Nightingale (1se)" = cvd.nri.test(risk[(PGS) & name == "Conventional RF"], risk[(PGS) & name == "Nightingale" & lambda == "lambda.1se"]),
	"Base + PGS to PGS + Blood + Nightingale (1se)" = cvd.nri.test(risk[(PGS) & name == "Conventional RF"], risk[(PGS) & name == "Blood + Nightingale" & lambda == "lambda.1se"]),
	"Base to PGS + Blood (1se)" = cvd.nri.test(risk[!(PGS) & name == "Conventional RF"], risk[(PGS) & name == "Blood" & lambda == "lambda.1se"]),
	"Base to PGS + Nightingale (1se)" = cvd.nri.test(risk[!(PGS) & name == "Conventional RF"], risk[(PGS) & name == "Nightingale" & lambda == "lambda.1se"]),
	"Base to PGS + Blood + Nightingale (1se)" = cvd.nri.test(risk[!(PGS) & name == "Conventional RF"], risk[(PGS) & name == "Blood + Nightingale" & lambda == "lambda.1se"]),
	"Base to Blood (min)" = cvd.nri.test(risk[!(PGS) & name == "Conventional RF"], risk[!(PGS) & name == "Blood" & lambda == "lambda.min"]),
	"Base to Nightingale (min)" = cvd.nri.test(risk[!(PGS) & name == "Conventional RF"], risk[!(PGS) & name == "Nightingale" & lambda == "lambda.min"]),
	"Base to Blood + Nightingale (min)" = cvd.nri.test(risk[!(PGS) & name == "Conventional RF"], risk[!(PGS) & name == "Blood + Nightingale" & lambda == "lambda.min"]),
	"Base + PGS to PGS + Blood (min)" = cvd.nri.test(risk[(PGS) & name == "Conventional RF"], risk[(PGS) & name == "Blood" & lambda == "lambda.min"]),
	"Base + PGS to PGS + Nightingale (min)" = cvd.nri.test(risk[(PGS) & name == "Conventional RF"], risk[(PGS) & name == "Nightingale" & lambda == "lambda.min"]),
	"Base + PGS to PGS + Blood + Nightingale (min)" = cvd.nri.test(risk[(PGS) & name == "Conventional RF"], risk[(PGS) & name == "Blood + Nightingale" & lambda == "lambda.min"]),
	"Base to PGS + Blood (min)" = cvd.nri.test(risk[!(PGS) & name == "Conventional RF"], risk[(PGS) & name == "Blood" & lambda == "lambda.min"]),
	"Base to PGS + Nightingale (min)" = cvd.nri.test(risk[!(PGS) & name == "Conventional RF"], risk[(PGS) & name == "Nightingale" & lambda == "lambda.min"]),
	"Base to PGS + Blood + Nightingale (min)" = cvd.nri.test(risk[!(PGS) & name == "Conventional RF"], risk[(PGS) & name == "Blood + Nightingale" & lambda == "lambda.min"]),
  "Base to CRP" = cvd.nri.test(risk[!(PGS) & name == "Conventional RF"], risk[!(PGS) & name == "CRP"]),
  "Base to GlycA" = cvd.nri.test(risk[!(PGS) & name == "Conventional RF"], risk[!(PGS) & name == "GlycA"]),
  "Base to PGS + CRP" = cvd.nri.test(risk[!(PGS) & name == "Conventional RF"], risk[(PGS) & name == "CRP"]),
  "Base to PGS + GlycA" = cvd.nri.test(risk[!(PGS) & name == "Conventional RF"], risk[(PGS) & name == "GlycA"]),
  "Base + PGS to PGS + CRP" = cvd.nri.test(risk[(PGS) & name == "Conventional RF"], risk[(PGS) & name == "CRP"]),
  "Base + PGS to PGS + GlycA" = cvd.nri.test(risk[(PGS) & name == "Conventional RF"], risk[(PGS) & name == "GlycA"])
)
saveRDS(nri_list, file="analyses/public_health_modelling/net_reclassification/nri.rds")

# Extract tables of estimates
nri_estimates <- rbindlist(idcol="model", fill=TRUE, lapply(nri_list, function(l1) {
  rbindlist(idcol="nri_type", fill=TRUE, lapply(l1, function(l2) {
    as.data.table(l2[["nri"]], keep.rownames="metric")
  }))
}))

# Extract tables of reclassifications for categorical nris
reclassified_cases <- rbindlist(idcol="model", fill=TRUE, lapply(nri_list, function(l1) {
  rbindlist(idcol="nri_type", fill=TRUE, lapply(l1, function(l2) {
    as.data.table(l2[["rtab.case"]])
  }))
}))

reclassified_controls <- rbindlist(idcol="model", fill=TRUE, lapply(nri_list, function(l1) {
  rbindlist(idcol="nri_type", fill=TRUE, lapply(l1, function(l2) {
    as.data.table(l2[["rtab.ctrl"]])
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
reclassified[reclassified_controls, on = .(model, nri_type, Standard, New), Controls := N]
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

fwrite(nri_estimates, sep="\t", quote=FALSE, file="analyses/public_health_modelling/net_reclassification/nri_estimates.txt")
fwrite(reclassified, sep="\t", quote=FALSE, file="analyses/public_health_modelling/net_reclassification/nri_reclassified.txt")
