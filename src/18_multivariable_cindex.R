library(data.table)
library(foreach)
source("src/utils/Cindex.R")

# Create output directory
system("mkdir -p analyses/rest")

# Load in pre-computed linear predictors
pred_scores <- fread("analyses/CVD_weight_training/CVD_linear_predictors_and_risk.txt")

# Extract set of models to test
model_info <- unique(pred_scores[, .(endpoint, score, model_type, model)])
model_info <- model_info[model_type != ""]
model_info <- rbind(idcol="model_sex", "Males"=model_info, "Females"=model_info, "Sex-stratified"=model_info)

# Fit models, extract C-indices, and compute delta C-indices
fits <- foreach(modelIdx = model_info[,.I], .combine=rbind) %do% {
  this_model <- model_info[modelIdx]

  # Extract relevant linear predictors
  new_lp <- pred_scores[this_model, on = .(endpoint, score, model_type, model)]
  old_lp <- pred_scores[model_type == ""][this_model, on = .(endpoint, score)]
  this_lp <- merge(
    old_lp[,.(eid, sex, incident_cvd, incident_cvd_followup, old=linear_predictor)],
    new_lp[,.(eid, sex, new=linear_predictor)],
    by=c("eid", "sex"), all=FALSE
  )

  # Filter to relevant sex
  if (this_model$model_sex == "Males") {
    this_lp <- this_lp[sex == "Male"]
  } else if (this_model$model_sex == "Females") {
    this_lp <- this_lp[sex == "Female"]
  }

	# Compute C-index for both scores and and delta-C index between old and new score
	if (this_model$model_sex == "Sex-stratified") {
		dc_res <- delta_cindex(this_lp$old, this_lp$new, this_lp[,Surv(incident_cvd_followup, incident_cvd)], this_lp$sex)
	} else {
		dc_res <- delta_cindex(this_lp$old, this_lp$new, this_lp[,Surv(incident_cvd_followup, incident_cvd)])
	}

	# Extract results
	res <- data.table(
		shared.samples=dc_res$shared.n, shared.cases=dc_res$shared.nevent, 
		old.C.index=dc_res$m1.C.index, old.C.SE=dc_res$m1.SE, old.C.L95=dc_res$m1.L95, old.C.U95=dc_res$m1.U95,
		new.C.index=dc_res$m2.C.index, new.C.SE=dc_res$m2.SE, new.C.L95=dc_res$m2.L95, new.C.U95=dc_res$m2.U95,
		deltaC=dc_res$deltaC, deltaC.SE=dc_res$deltaC.SE, deltaC.L95=dc_res$deltaC.L95, deltaC.U95=dc_res$deltaC.U95, deltaC.pval=dc_res$deltaC.pval
	)

	res <- cbind(this_model, res)
	return(res)
}

# Write out
fwrite(fits, sep="\t", quote=FALSE, file="analyses/test/discovery_cindices.txt")

