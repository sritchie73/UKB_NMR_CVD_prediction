# Requires at least 1:10:0 to run on compute nodes, recommend running with 20 cores
library(data.table)
library(foreach)
library(survival)
source("src/utils/score_cindex.R")
registerDoMC(10)

# Create output directory
system("mkdir -p analyses/univariate")

# Load biomarker information sheets
bio_info <- fread("data/ukb/biomarkers/output/biomarker_info.txt")
nmr_info <- fread("data/ukb/NMR_metabolomics/biomarker_information.txt")

# Load analysis cohort
dat <- fread("data/cleaned/analysis_cohort.txt")

# Code risk factors as per SCORE2
dat[, age := (age - 60)/5]
dat[is.na(smoking), smoking := FALSE]
dat[, smoking := factor(smoking, levels=c(FALSE, TRUE))]
dat[, sbp := (sbp - 120)/20]
dat[, tchol := (tchol - 6)/1]
dat[, hdl := (hdl - 1.3)/0.5]

# Code factors, using non-risk/lower-risk group as reference
dat[, sex := factor(sex, levels=c("Female", "Male"))]

# Get list of biomarkers to test
test_nmr <- nmr_info[(Nightingale), Biomarker]
test_assay <- bio_info[!is.na(UKB.Field.ID) & sample_type != "Urine" & var != "tchol" & var != "hdl", var]

# Build set of models to test
models <- rbind(
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ sex", type="demographics", name="sex"),
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ age", type="demographics", name="age"),
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ age + smoking", type="risk_factors", name="smoking"),
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ age + sbp", type="risk_factors", name="sbp"),
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ age + tchol", type="risk_factors", name="tchol"),
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ age + hdl", type="risk_factors", name="hdl"),
  data.table(formula=sprintf("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2_excl_UKB) + scale(%s)", test_nmr), type="NMR", name=test_nmr),
  data.table(formula=sprintf("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2_excl_UKB) + scale(%s)", test_assay), type="assays", name=test_assay)
)

# Compute C-indices
cinds <- foreach(this_sex=c("Males", "Females", "Sex-stratified"), .combine=rbind) %do% { 

  # Extract subset of data needed
  if (this_sex == "Males") {
    this_dat <- dat[sex == "Male"]
  } else if (this_sex == "Females") {
    this_dat <- dat[sex == "Female"]
  } else {
    this_dat <- dat
  }

  # Fit individual models
  foreach(midx = models[,.I], .combine=rbind) %dopar% {
    # Extract model information
    this_model <- models[midx]

    # Set up formula
    mf <- this_model$formula
    if (this_model$name == "sex") {
      if (this_sex != "Sex-stratified") {
        return(NULL) 
      }
    } else if (this_sex == "Sex-stratified") {
      mf <- paste(mf, "+ strata(sex)")
    }

    # Fit Cox proportional hazards model
    cph <- coxph(as.formula(mf), data=this_dat)

    # Get C-index and its 95% CI 
    cindex <- summary(cph)$concordance[1]
    cindex.se <- summary(cph)$concordance[2]
    cindex.l95 <- cindex - qnorm(0.975)*cindex.se
    cindex.u95 <- cindex + qnorm(0.975)*cindex.se
     
    # Return results
    cbind(this_model, sex=this_sex, samples=cph$n, cases=cph$nevent, 
      C.index=cindex, L95=cindex.l95, U95=cindex.u95)
  }
}

# Load in and add SCORE2 C-index
score2_cindex <- fread("analyses/test/cindex_by_SCORE2_method.txt")
score2_cindex <- score2_cindex[SCORE2_method == "Weights derived excluding UK Biobank"]
score2_cindex[, formula := "Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2_excl_UKB)"]
score2_cindex[, type := "SCORE2"]
score2_cindex[, name := "SCORE2"]
score2_cindex[sex == "Males", c("samples", "cases") := dat[sex == "Male", .(.N, sum(incident_cvd))]]
score2_cindex[sex == "Females", c("samples", "cases") := dat[sex == "Female", .(.N, sum(incident_cvd))]]
score2_cindex[sex == "Females", c("samples", "cases") := dat[sex == "Female", .(.N, sum(incident_cvd))]]
score2_cindex[sex == "Sex-stratified", c("samples", "cases") := dat[, .(.N, sum(incident_cvd))]]
score2_cindex[, c("SCORE2_method", "SE") := NULL]
cinds <- rbind(cinds, score2_cindex)

# Add in missing strata term to documented formula
cinds[sex == "Sex-stratified" & name != "sex", formula := paste(formula, "+ strata(sex)")]

# Add in delta change compared to reference
sex_ref <- cinds[name == "sex"]
sex_ref[, name := "age"]
cinds[sex_ref, on = .(name, sex), c("deltaC", "deltaC.L95", "deltaC.U95") := .(C.index - i.C.index, L95 - i.C.index, U95 - i.C.index)]

age_ref <- cinds[name == "age"]
age_ref[, type := "risk_factors"]
cinds[age_ref, on = .(type, sex), c("deltaC", "deltaC.L95", "deltaC.U95") := .(C.index - i.C.index, L95 - i.C.index, U95 - i.C.index)]

score2_ref <- cinds[name == "SCORE2"]
score2_ref[, type := NULL]
score2_ref <- rbind(idcol="type", "NMR"=score2_ref, "assays"=score2_ref)
cinds[score2_ref, on = .(type, sex), c("deltaC", "deltaC.L95", "deltaC.U95") := .(C.index - i.C.index, L95 - i.C.index, U95 - i.C.index)]

# Add human friendly display name
cinds[, display_name := name]
cinds[name == "age", display_name := "Age"]
cinds[name == "sex", display_name := "Sex"]
cinds[name == "smoking", display_name := "Smoker"]
cinds[name == "sbp", display_name := "SBP"]
cinds[name == "hdl", display_name := "HDL-C"]
cinds[name == "tchol", display_name := "Total-C"]
cinds[nmr_info[Group %like% "Amino" & Type == "Non-derived"], on = .(name=Biomarker), display_name := Description]
cinds[, display_name := gsub("_", "-", display_name)]
cinds[, display_name := gsub("-pct", " %", display_name)]
cinds[, display_name := gsub("-by-", " / ", display_name)]
cinds[name == "Clinical_LDL_C", display_name := "Clinical LDL-C"]
cinds[name %in% test_assay, display_name := fcase(
  name == "alt", "ALT",
  name == "alb", "Albumin",
  name == "alp", "ALP",
  name == "apoa1", "ApoA1",
  name == "apob", "ApoB",
  name == "asp", "ASP",
  name == "calcium", "Calcium",
  name == "creat", "Creatinine",
  name == "crp", "CRP",
  name == "cyst", "CST3",
  name == "dbili", "Direct bili",
  name == "ggt", "GGT",
  name == "glucose", "Glucose",
  name == "hba1c", "HbA1c",
  name == "igf1", "IGF-1",
  name == "lpa", "Lp(a)",
  name == "ldl", "LDL-C",
  name == "hdl", "HDL-C",
  name == "oest", "Oestradiol",
  name == "phos", "Phosphate",
  name == "rheuf", "RF",
  name == "shbg", "SHBG",
  name == "testos", "Testosterone",
  name == "tbili", "Total bili",
  name == "protein", "Total prot",
  name == "tchol", "Total-C",
  name == "trig", "Triglyceride",
  name == "uric", "Urate",
  name == "urea", "Urea",
  name == "vitd25", "Vitamin D"
)]

# Write out
fwrite(cinds, sep="\t", quote=FALSE, file="analyses/univariate/cindices.txt")

