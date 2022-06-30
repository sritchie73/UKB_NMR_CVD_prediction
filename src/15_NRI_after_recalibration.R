library(data.table)
library(nricens)

# Make output directory
system("mkdir -p analyses/public_health_modelling/net_reclassification")

# Load model information
model_info <- fread("analyses/test/model_fit_information.txt")

# Load recalibrated predicted CVD risk
risk <- fread("analyses/public_health_modelling/risk_recalibration/absolute_risks.txt")

# Wrapper function to compute continuous/categorical NRI from recalibrated risks
nri.test <- function(base_risk, new_risk) {
  shared_eid <- data.table(eid=intersect(base_risk$eid, new_risk$eid))

  base_risk <- base_risk[shared_eid, on = .(eid)]
  new_risk <- new_risk[shared_eid, on = .(eid)]

  # Calculate continuous NRI
  contNRI <- nricens(event = base_risk$incident_cvd, time = base_risk$incident_followup,
                     p.std = base_risk$recalibrated_risk, p.new = new_risk$recalibrated_risk,
                     updown = "diff", cut = 0, t0 = 10, niter = 1000)

  # Calculate NRI for NICE 2014
  niceNRI <- nricens(event = base_risk$incident_cvd, time = base_risk$incident_followup,
                     p.std = base_risk$recalibrated_risk, p.new = new_risk$recalibrated_risk,
                     updown = "category", cut = c(0.05, 0.1), t0 = 10, niter = 1000)

  # Calculate NRI for AHA/ACC 2019
  ahaaccNRI <- nricens(event = base_risk$incident_cvd, time = base_risk$incident_followup,
                       p.std = base_risk$recalibrated_risk, p.new = new_risk$recalibrated_risk,
                       updown = "category", cut = c(0.05, 0.075), t0 = 10, niter = 1000)

  # Add in sample size and case numbers
  contNRI$n <- niceNRI$n <- ahaaccNRI$n <- base_risk[,.N]
  contNRI$nevent <- niceNRI$nevent <- ahaaccNRI$nevent <- base_risk[, sum(incident_cvd)]

  list("Continuous NRI"=contNRI, "NICE 2014 Categorical NRI"=niceNRI, "ACC/AHA 2019 Categorical NRI"=ahaaccNRI)
}

# Run NRI analysis 
nri_list <- list(
	"Base to PGS" = nri.test(risk[!(PGS) & name == "Conventional RF"], risk[(PGS) & name == "Conventional RF"]),
  "Base to CRP" = nri.test(risk[!(PGS) & name == "Conventional RF"], risk[!(PGS) & name == "CRP"]),
  "Base to PGS + CRP" = nri.test(risk[!(PGS) & name == "Conventional RF"], risk[(PGS) & name == "CRP"]),
	"Base to Assays (1se)" = nri.test(risk[!(PGS) & name == "Conventional RF"], risk[!(PGS) & name == "Assays" & lambda == "lambda.1se"]),
	"Base to NMR (1se)" = nri.test(risk[!(PGS) & name == "Conventional RF"], risk[!(PGS) & name == "NMR" & lambda == "lambda.1se"]),
	"Base to NMR + Assays (1se)" = nri.test(risk[!(PGS) & name == "Conventional RF"], risk[!(PGS) & name == "NMR + Assays" & lambda == "lambda.1se"]),
	"Base to PGS + Assays (1se)" = nri.test(risk[!(PGS) & name == "Conventional RF"], risk[(PGS) & name == "Assays" & lambda == "lambda.1se"]),
	"Base to PGS + NMR (1se)" = nri.test(risk[!(PGS) & name == "Conventional RF"], risk[(PGS) & name == "NMR" & lambda == "lambda.1se"]),
	"Base to PGS + NMR + Assays (1se)" = nri.test(risk[!(PGS) & name == "Conventional RF"], risk[(PGS) & name == "NMR + Assays" & lambda == "lambda.1se"]),
	"Base to Assays (min)" = nri.test(risk[!(PGS) & name == "Conventional RF"], risk[!(PGS) & name == "Assays" & lambda == "lambda.min"]),
	"Base to NMR (min)" = nri.test(risk[!(PGS) & name == "Conventional RF"], risk[!(PGS) & name == "NMR" & lambda == "lambda.min"]),
	"Base to NMR + Assays (min)" = nri.test(risk[!(PGS) & name == "Conventional RF"], risk[!(PGS) & name == "NMR + Assays" & lambda == "lambda.min"]),
	"Base to PGS + Assays (min)" = nri.test(risk[!(PGS) & name == "Conventional RF"], risk[(PGS) & name == "Assays" & lambda == "lambda.min"]),
	"Base to PGS + NMR (min)" = nri.test(risk[!(PGS) & name == "Conventional RF"], risk[(PGS) & name == "NMR" & lambda == "lambda.min"]),
	"Base to PGS + NMR + Assays (min)" = nri.test(risk[!(PGS) & name == "Conventional RF"], risk[(PGS) & name == "NMR + Assays" & lambda == "lambda.min"])
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
