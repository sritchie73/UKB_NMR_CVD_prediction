library(data.table)
library(foreach)
library(doMC)
library(survival)
source("src/utils/Cindex.R")

registerDoMC(10) # request > 10 cores when running, otherwise not enough memory per core
setDTthreads(10)

# Create output directory
system("mkdir -p analyses/univariate")

# Load analysis cohort
dat <- rbind(fill=TRUE,
  fread("data/cleaned/analysis_cohort.txt"),
  fread("data/cleaned/phase3_analysis_cohort.txt")
)

# Load in results from discovery cohort
disc <- fread("analyses/univariate/discovery_analysis.txt")
disc <- disc[,.SD[1], by=.(endpoint, sex, score, model_type)]

# extract unique model information
model_info <- disc[,.(endpoint, sex, score, model, biomarker, model_type)]

# Fit models, extract C-indices, and compute delta C-indices
fits <- foreach(modelIdx = model_info[,.I], .combine=rbind) %dopar% {
  this_model <- model_info[modelIdx]

  # Filter to relevant sex
  if (this_model$sex == "Males") {
    this_dat <- dat[sex == "Male"]
  } else if (this_model$sex == "Females") {
    this_dat <- dat[sex == "Female"]
  } else {
    this_dat <- copy(dat)
  }

  # Extract relevant endpoint
  if (this_model$endpoint == "cvd") {
    this_y <- Surv(this_dat[["incident_cvd_followup"]], this_dat[["incident_cvd"]])
  } else if (this_model$endpoint == "cvd_narrow") {
    this_y <- Surv(this_dat[["incident_cvd2_followup"]], this_dat[["incident_cvd2"]])
  }

  # Extract relevant risk score
  this_score <- this_dat[[this_model$score]]

  # Extract strata grouping
  this_strata <- this_dat[["sex"]]

  # Compute C-index for the risk score on its own
	if (this_model$sex == "Sex-stratified") {
		rs_res <- cindex(this_score, this_y, this_strata)
	} else {
		rs_res <- cindex(this_score, this_y)
	}

  # If the models are only the risk scores, we can just return the relevant info now
  if (is.na(this_model$biomarker) || this_model$biomarker == "") {
    # Extract results
    res <- data.table(
      old.samples=rs_res$n, old.cases=rs_res$nevent, old.C.index=rs_res$C.index, old.C.SE=rs_res$SE, old.C.L95=rs_res$L95, old.C.U95=rs_res$U95,
      new.samples=NA_real_, new.cases=NA_real_, new.C.index=NA_real_, new.C.SE=NA_real_, new.C.L95=NA_real_, new.C.U95=NA_real_,
      shared.samples=NA_real_, shared.cases=NA_real_, deltaC=NA_real_, deltaC.SE=NA_real_, deltaC.L95=NA_real_, deltaC.U95=NA_real_, deltaC.pval=NA_real_, deltaC.fdr=NA_real_,
      biomarker.mean=NA_real_, biomarker.sd=NA_real_, biomarker.HR=NA_real_, biomarker.HR.SE=NA_real_, biomarker.HR.L95=NA_real_, biomarker.HR.U95=NA_real_, biomarker.HR.pval=NA_real_, biomarker.HR.fdr=NA_real_
    )
    res <- cbind(this_model, res)
    return(res)
  } else { 
    # For models including a biomarker, use the hazard ratio from the discovery cohort to build a new linear predictor
    bio_vec <- this_dat[[this_model$biomarker]]
    bio_vec <- scale(bio_vec) # key assumption is that mean and sd in discovery cohort were representative and generalizable when applied to other cohorts
    bio_score <- this_score + disc[modelIdx, log(biomarker.HR)] * as.vector(bio_vec)

    # Compute mean and sd of biomarker for comparison to discovery cohort if needed
    bio_vec <- this_dat[[this_model$biomarker]]
    bio_mean <- mean(bio_vec, na.rm=TRUE)
    bio_sd <- sd(bio_vec, na.rm=TRUE)

    # Compute C-index directly from new risk score
    if (this_model$sex == "Sex-stratified") {
      cind <- cindex(bio_score, this_y, this_strata)
    } else {
      cind <- cindex(bio_score, this_y)
    }
  
    # Compute delta-C index from target score
    if (this_model$sex == "Sex-stratified") {
      dc_res <- delta_cindex(this_score, bio_score, this_y, this_strata)
    } else {
      dc_res <- delta_cindex(this_score, bio_score, this_y)
    }

    # Extract results
    res <- data.table(
      old.samples=rs_res$n, old.cases=rs_res$nevent, old.C.index=rs_res$C.index, old.C.SE=rs_res$SE, old.C.L95=rs_res$L95, old.C.U95=rs_res$U95,
      new.samples=cind$n, new.cases=cind$nevent, new.C.index=cind$C.index, new.C.SE=cind$SE, new.C.L95=cind$L95, new.C.U95=cind$U95,
      shared.samples=dc_res$shared.n, shared.cases=dc_res$shared.nevent, deltaC=dc_res$deltaC, deltaC.SE=dc_res$deltaC.SE, deltaC.L95=dc_res$deltaC.L95, deltaC.U95=dc_res$deltaC.U95, deltaC.pval=dc_res$deltaC.pval, deltaC.fdr=NA_real_,
      biomarker.mean=bio_mean, biomarker.sd=bio_sd, biomarker.HR=NA_real_, biomarker.HR.SE=NA_real_, biomarker.HR.L95=NA_real_, biomarker.HR.U95=NA_real_, biomarker.HR.pval=NA_real_, biomarker.HR.fdr=NA_real_
    )
    res <- cbind(this_model, res)
    return(res)
  }
}

# add FDR correction
fits[model_type != "risk score", deltaC.fdr := p.adjust(deltaC.pval, method="fdr"), by=.(endpoint, sex, score, model_type)]

# Write out
fwrite(fits, sep="\t", quote=FALSE, file="analyses/univariate/pooled_analysis.txt")
