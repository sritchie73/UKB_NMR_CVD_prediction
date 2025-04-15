library(data.table)
library(foreach)
library(doMC)
library(nricens)

# request > 10 cores when running, otherwise not enough memory per core
# Ran successfully in ~55 minutes using 15 cores on cclake partition
registerDoMC(10) 
setDTthreads(10)

dat <- rbind(
  fread('analyses/CVD_weight_training/CVD_linear_predictors_and_risk.txt'),
  fread('analyses/CVD_weight_training/phase3_CVD_linear_predictors_and_risk.txt')
)

# Get unique models to fit
model_info <- unique(dat[model_type != "", .(endpoint, score, model_type, model)])
model_info <- rbind(idcol="model_sex", "Sex-stratified"=model_info, "Males"=model_info, "Females"=model_info)

# Run NRI analysis
nri_lists <- foreach(modelIdx = model_info[,.I]) %dopar% {
  # Extract pair of models to test
  this_model <- model_info[modelIdx]
  comp_risk <- dat[endpoint == this_model$endpoint & score == this_model$score]
  comp_risk <- comp_risk[model_type %in% c("", this_model$model_type)]
  comp_risk[, comp_model := ifelse(model_type == "", "base", "new")]
  if (this_model$model_sex == "Males") {
    comp_risk <- comp_risk[sex == "Male"]
  } else if (this_model$model_sex == "Females") {
    comp_risk <- comp_risk[sex == "Female"]
  }
  comp_risk <- dcast(comp_risk, eid + age + incident_cvd_followup + incident_cvd ~ comp_model, value.var="uk_calibrated_risk")
  comp_risk <- comp_risk[!is.na(base) & !is.na(new)]
  
  # Set risk thresholds depending on the score
  if (this_model$score == "QRISK3") {
    # NICE guidelines simply use a 10% risk threshold for QRISK3
    risk_thresholds = 0.1
  } else {
    # The ESC 2021 guidelines use different risk thresolds for different age groups:
    #   <2.5%, 2.5-7.5%, and >7.5% in people under 50
    #   <5%, 5-10%, and >10% in people over 50
    # To compute a single NRI across all ages, we can therefor simply add 2.5% risk to everyone under 50 so
    # we can then use the same thresholds in everyone
    comp_risk[age < 50, base := base + 0.025]
    comp_risk[age < 50, new := new + 0.025]
    risk_thresholds = c(0.05, 0.1)
  }

	# Run NRI analysis
	NRI <- nricens(event = comp_risk$incident_cvd, time = comp_risk$incident_cvd_followup,
								 p.std = comp_risk$base, p.new = comp_risk$new,
								 updown = "category", cut = risk_thresholds, 
                 t0 = 10, niter = 1000)

  # Add in sample size, case numbers, and model_info
  NRI$n <- comp_risk[,.N]
  NRI$nevent <- comp_risk[, sum(incident_cvd)]
  NRI$model_info <- this_model

  return(NRI)
}
saveRDS(nri_lists, file="analyses/test/pooled_categorical_nri.rds")

# Extract tables of estimates
nri_estimates <- foreach(modelIdx = model_info[,.I], .combine=rbind) %do% {
  cbind(
    nri_lists[[modelIdx]][["model_info"]],
    shared.samples=nri_lists[[modelIdx]][["n"]],
    shared.cases=nri_lists[[modelIdx]][["nevent"]],
    as.data.table(nri_lists[[modelIdx]][["nri"]], keep.rownames="metric")
  )
}

# Compute bootstrap standard errors
nri_bsse <- foreach(modelIdx = model_info[,.I], .combine=rbind) %do% {
  se_vec <- apply(nri_lists[[modelIdx]][["bootstrapsample"]], 2, sd)
  cbind(
    nri_lists[[modelIdx]][["model_info"]],
    data.table(metric=names(se_vec), SE=se_vec)
  )
}
nri_estimates[nri_bsse, on = .(model_sex, endpoint, score, model_type, model, metric), SE := i.SE]

# Compute 95% CI and P-value from BSSE 
# (Reassuringly, 95% CI are very very similar to percentile method)
nri_estimates[, L95 := Estimate - qnorm(1-(0.05/2))*SE]
nri_estimates[, U95 := Estimate + qnorm(1-(0.05/2))*SE]
nri_estimates[, Pval := pmin(1, pnorm(abs(Estimate/SE), lower.tail=FALSE)*2)]

# Extract tables of reclassifications for categorical nris
reclassified_cases <- foreach(modelIdx = model_info[,.I], .combine=rbind) %do% {
  cbind(
    nri_lists[[modelIdx]][["model_info"]],
    as.data.table(nri_lists[[modelIdx]][["rtab.case"]])
  )
}

reclassified <- foreach(modelIdx = model_info[,.I], .combine=rbind) %do% {
  cbind(
    nri_lists[[modelIdx]][["model_info"]],
    as.data.table(nri_lists[[modelIdx]][["rtab"]])
  )
}

# Merge
setnames(reclassified, "N", "All")
reclassified[reclassified_cases, on = .(model_sex, endpoint, score, model_type, model, Standard, New), Cases := N]
setnames(reclassified, "Standard", "Old")

# Add in total sample size and total cases to reclassified table
reclassified[nri_estimates, on = .(model_sex, endpoint, score, model_type), Total_Samples := i.shared.samples]
reclassified[nri_estimates, on = .(model_sex, endpoint, score, model_type), Total_Cases := i.shared.cases]

fwrite(nri_estimates, sep="\t", quote=FALSE, file="analyses/test/pooled_categorical_nri_estimates.txt")
fwrite(reclassified, sep="\t", quote=FALSE, file="analyses/test/pooled_categorical_nri_reclassified.txt")

