library(data.table)
library(foreach)
library(doMC)
library(survival)
source("src/utils/Cindex.R")

registerDoMC(10)
setDTthreads(10)

# Create output directory
system("mkdir -p analyses/univariate")

# Load biomarker information sheet
nmr_info <- fread("data/ukb/NMR_metabolomics/biomarker_information.txt")
assay_info <- fread("data/ukb/biomarkers/output/biomarker_info.txt")

# Load analysis cohort
dat <- fread("data/cleaned/analysis_cohort.txt")

# Get list of biomarkers to test
test_nmr <- nmr_info[(Nightingale), Biomarker]
test_assay <- assay_info[sample_type != "Urine" & !is.na(UKB.Field.ID) & var != "tchol" & var != "hdl", var]

# Build set of models to test
model_info <- foreach(this_endpoint=c("cvd", "cvd_narrow"), .combine=rbind) %:%
  foreach(this_score=c("SCORE2", "SCORE2_excl_UKB", "QRISK3"), .combine=rbind) %:%
    foreach(this_sex=c("Males", "Females", "Sex-stratified"), .combine=rbind) %:%
      foreach(this_model=c(this_score, paste(this_score, "+", test_nmr), paste(this_score, "+", test_assay)), .combine=rbind) %do% { 
        data.table(endpoint=this_endpoint, sex=this_sex, score=this_score, model=this_model)
}
model_info[score == "SCORE2" & model != "SCORE2", biomarker := gsub("SCORE2 \\+ ", "", model)]
model_info[score == "SCORE2_excl_UKB" & model != "SCORE2_excl_UKB", biomarker := gsub("SCORE2_excl_UKB \\+ ", "", model)]
model_info[score == "QRISK3" & model != "QRISK3", biomarker := gsub("QRISK3 \\+ ", "", model)]
model_info[score == "SCORE2_excl_UKB", model := gsub("SCORE2_excl_UKB", "SCORE2", model)]
model_info[, model_type := fcase(
  model %in% c("SCORE2", "SCORE2_excl_UKB", "QRISK3"), "risk score",
  biomarker %in% test_assay, "assays",
  biomarker %in% test_nmr, "NMR"
)]

# Convert sex to integer for strata
strata_num <- dat[,.GRP, by=sex]
dat[strata_num, on = .(sex), sex_int := i.GRP]

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
  if (is.na(this_model$biomarker)) {
    # Extract results
    res <- data.table(
      old.samples=rs_res$n, old.cases=rs_res$nevent, old.C.index=rs_res$C.index, old.C.SE=rs_res$SE, old.C.L95=rs_res$L95, old.C.U95=rs_res$U95,
      new.samples=NA_real_, new.cases=NA_real_, new.C.index=NA_real_, new.C.SE=NA_real_, new.C.L95=NA_real_, new.C.U95=NA_real_,
      shared.samples=NA_real_, shared.cases=NA_real_, deltaC=NA_real_, deltaC.SE=NA_real_, deltaC.L95=NA_real_, deltaC.U95=NA_real_, deltaC.pval=NA_real_, deltaC.fdr=NA_real_,
      biomarker.HR=NA_real_, biomarker.HR.SE=NA_real_, biomarker.HR.L95=NA_real_, biomarker.HR.U95=NA_real_, biomarker.HR.pval=NA_real_, biomarker.HR.fdr=NA_real_
    )
    res <- cbind(this_model, res)
    return(res)
  } else { 
    # Now handle case where we have to fit a Cox prortional hazards model

    # Build the model formula based on model options
    if (this_model$endpoint == "cvd") {
      mf <- "Surv(incident_cvd_followup, incident_cvd) ~"
    } else if (this_endpoint == "cvd_narrow") {
      mf <- "Surv(incident_cvd2_followup, incident_cvd2) ~"
    }
    if (this_model$sex == "Sex-stratified") {
      mf <- paste(mf, "strata(sex) +")
    }
    mf <- paste(mf, sprintf("offset(%s) +", this_model$score))
    mf <- paste(mf, sprintf("scale(%s)", this_model$biomarker)) 

    # Fit Cox proportional hazards model
		cx <- coxph(as.formula(mf), data=this_dat, x=TRUE)

    # Extra hazard ratios
		cf <- coef(summary(cx))
		ci <- confint(cx)

    HR <- cf[,2]
    HR.SE <- cf[,3]
    HR.L95 <- exp(ci[,1])
    HR.U95 <- exp(ci[,2])
    HR.pval <- cf[,5]

    # Extract C-index information
    cind <- cindex(cx)
  
    # Compute delta-C index from target score
    if (this_model$sex == "Sex-stratified") {
      dc_res <- delta_cindex(this_score, cx, this_y, this_strata)
    } else {
      dc_res <- delta_cindex(this_score, cx, this_y)
    }

    # Extract results
    res <- data.table(
      old.samples=rs_res$n, old.cases=rs_res$nevent, old.C.index=rs_res$C.index, old.C.SE=rs_res$SE, old.C.L95=rs_res$L95, old.C.U95=rs_res$U95,
      new.samples=cind$n, new.cases=cind$nevent, new.C.index=cind$C.index, new.C.SE=cind$SE, new.C.L95=cind$L95, new.C.U95=cind$U95,
      shared.samples=dc_res$shared.n, shared.cases=dc_res$shared.nevent, deltaC=dc_res$deltaC, deltaC.SE=dc_res$deltaC.SE, deltaC.L95=dc_res$deltaC.L95, deltaC.U95=dc_res$deltaC.U95, deltaC.pval=dc_res$deltaC.pval, deltaC.fdr=NA_real_,
      biomarker.HR=HR, biomarker.HR.SE=HR.SE, biomarker.HR.L95=HR.L95, biomarker.HR.U95=HR.U95, biomarker.HR.pval=HR.pval, biomarker.HR.fdr=NA_real_
    )
    res <- cbind(this_model, res)
    return(res)
  }
}

# add FDR correction
fits[model_type != "risk score", deltaC.fdr := p.adjust(deltaC.pval, method="fdr"), by=.(endpoint, sex, score, model_type)]
fits[model_type != "risk score", biomarker.HR.fdr := p.adjust(biomarker.HR.pval, method="fdr"), by=.(endpoint, sex, score, model_type)]

# Add human friendly display name
fits[model_type == "NMR", model := gsub("_pct", " %", model)]
fits[model_type == "NMR", model := gsub("_by_", " / ", model)]
fits[model_type == "NMR", model := gsub("_", "-", model)]
fits[model_type == "NMR", model := gsub("Clinical-", "Clinical ", model)]
fits[model_type == "NMR", model := gsub("Total-", "Total ", model)]
fits[model_type == "NMR", model := gsub("Remnant-", "Remnant ", model)]
fits[model_type == "NMR", model := gsub("-size", " size", model)]
fits[model_type == "assays", model := paste(score, "+", fcase(
  biomarker == "alb", "Albumin",
  biomarker == "alt", "ALT",
  biomarker == "alp", "ALP",
  biomarker == "apoa1", "ApoA1",
  biomarker == "apob", "ApoB",
  biomarker == "asp", "AST",
  biomarker == "dbili", "Bilirubin (direct)",
  biomarker == "tbili", "Bilirubin (total)",
  biomarker == "calcium", "Calcium",
  biomarker == "creat", "Creatinine",
  biomarker == "crp", "CRP",
  biomarker == "cyst", "Cystatin-C",
  biomarker == "ggt", "GGT",
  biomarker == "glucose", "Glucose",
  biomarker == "hba1c", "HbA1c",
  biomarker == "igf1", "IGF-1",
  biomarker == "lpa", "Lp(a)",
  biomarker == "ldl", "LDL cholesterol",
  biomarker == "oest", "Oestradiol",
  biomarker == "phos", "Phosphate",
  biomarker == "rheuf", "RF",
  biomarker == "shbg", "SHBG",
  biomarker == "testos", "Testosterone",
  biomarker == "protein", "Total protein",
  biomarker == "trig", "Triglycerides",
  biomarker == "uric", "Urate",
  biomarker == "urea", "Urea",
  biomarker == "vitd25", "Vitamin D"
))]

# Write out
fwrite(fits, sep="\t", quote=FALSE, file="analyses/univariate/discovery_analysis.txt")
