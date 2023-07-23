# Requires at least 1:10:0 to run on compute nodes, recommend running with 20 cores
library(data.table)
library(foreach)
library(survival)
source("src/utils/cv_coxph.R")
source("src/utils/score_cindex.R")
source("src/utils/cox_test.R")
registerDoMC(10)

# Create output directory
system("mkdir -p analyses/univariate")

# Load biomarker information sheets
bio_info <- fread("data/ukb/biomarkers/output/biomarker_info.txt")
nmr_info <- fread("data/ukb/NMR_metabolomics/biomarker_information.txt")

# Load analysis cohort
dat <- fread("data/cleaned/analysis_cohort.txt")

# Code factors, using non-risk/lower-risk group as reference
dat[, sex := factor(sex, levels=c("Female", "Male"))]
dat[, smoking := factor(smoking, levels=c("FALSE", "TRUE"))]

# Build set of models to test
test_nmr <- nmr_info[, Biomarker]
test_assay <- bio_info[!is.na(UKB.Field.ID) & sample_type != "Urine" & var != "tchol" & var != "hdl", var]
models <- rbind(
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ scale(age)", type="demographics", name="age"),
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ scale(age) + smoking", type="risk_factors", name="smoking"),
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ scale(age) + scale(sbp)", type="risk_factors", name="sbp"),
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ scale(age) + scale(tchol)", type="risk_factors", name="tchol"),
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ scale(age) + scale(hdl)", type="risk_factors", name="hdl"),
  data.table(formula=sprintf("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + scale(%s)", test_nmr), type="NMR", name=test_nmr),
  data.table(formula=sprintf("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + scale(%s)", test_assay), type="assays", name=test_assay)
)

# Takes about an hour to run on 10 cores
cinds <- foreach(this_sex=c("Male", "Female"), .combine=rbind) %:% 
  foreach(mIdx = models[,.I], .combine=rbind) %do% {
    this_model <- models[mIdx]
    cbind(this_model, sex=paste0(this_sex, "s"), cv.cindex(this_model$formula, dat[sex == this_sex], "cvd_prediction_foldid"))
}

# Load in and add SCORE2 C-index
score2_cindex <- fread("analyses/test/cindices.txt")
cinds <- rbind(cinds, 
  score2_cindex[model == "SCORE2" & sex == "Males", 
   .(formula="Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2)", type="SCORE2", name="SCORE2",
     sex, C.index, SE, L95, U95, Samples=dat[sex == "Male", .N], Cases=dat[sex == "Male", sum(incident_cvd)], Missing=0, Missing.Cases=0)],
  score2_cindex[model == "SCORE2" & sex == "Females", 
   .(formula="Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2)", type="SCORE2", name="SCORE2",
     sex, C.index, SE, L95, U95, Samples=dat[sex == "Female", .N], Cases=dat[sex == "Male", sum(incident_cvd)], Missing=0, Missing.Cases=0)])

# Add in delta change compared to reference
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

# Order columns
cinds <- cinds[, .(type, name, display_name, sex, Samples, Cases, Missing, Missing.Cases, C.index, SE, L95, U95, deltaC, deltaC.L95, deltaC.U95)]

# Write out
fwrite(cinds, sep="\t", quote=FALSE, file="analyses/univariate/cindices.txt")

# Estimate hazards ratios (takes about a minute to run on 10 cores)
cox_list <- foreach(this_sex = c("Male", "Female")) %:% 
  foreach(mIdx = models[,.I]) %dopar% {
    this_model <- models[mIdx]
    cox.test(this_model$formula, "incident_cvd", dat[sex == this_sex])
}
names(cox_list) <- c("Male", "Female")
saveRDS(cox_list, file="analyses/test/univariate_cox_list.rds")

# Extract hazard ratios
hrs <- foreach(this_sex=c("Male", "Female"), .combine=rbind) %:%
  foreach(mIdx = models[,.I], .combine=rbind) %do% {
    this_hrs <- cox_list[[this_sex]][[mIdx]][["coefficients"]]
    cbind(models[mIdx], sex=paste0(this_sex, "s"), this_hrs)
}

# Extract variables corresponding to coefficients
hrs[, var := gsub("scale\\(", "", coefficient)]
hrs[, var := gsub("\\)", "", var)]
hrs[coefficient == "smokingTRUE", var := "smoking"]

# Add in human friendly coefficient names
hrs[coefficient == "smokingTRUE", coef_name := "Smoker"]
hrs[var == "age", coef_name := "Age"]
hrs[var == "sbp", coef_name := "SBP"]
hrs[var == "tchol", coef_name := "Total Cholesterol"]
hrs[var == "hdl", coef_name := "HDL Cholesterol"]
hrs[nmr_info[Group %like% "Amino" & Type == "Non-derived"], on = .(var=Biomarker), coef_name := Description]
hrs[is.na(coef_name), coef_name := var]
hrs[, coef_name := gsub("_", "-", coef_name)]
hrs[var == "Clinical_LDL_C", coef_name := "Clinical LDL-C"]
hrs[coef_name %in% test_assay, coef_name := fcase(
  var == "alt", "ALT",
  var == "alb", "Albumin",
  var == "alp", "ALP",
  var == "apoa1", "ApoA1",
  var == "apob", "ApoB",
  var == "asp", "ASP",
  var == "calcium", "Calcium",
  var == "creat", "Creatinine",
  var == "crp", "CRP",
  var == "cyst", "CST3",
  var == "dbili", "Direct bili",
  var == "ggt", "GGT",
  var == "glucose", "Glucose",
  var == "hba1c", "HbA1c",
  var == "igf1", "IGF-1",
  var == "lpa", "Lp(a)",
  var == "ldl", "LDL-C",
  var == "hdl", "HDL-C",
  var == "oest", "Oestradiol",
  var == "phos", "Phosphate",
  var == "rheuf", "RF",
  var == "shbg", "SHBG",
  var == "testos", "Testosterone",
  var == "tbili", "Total bili",
  var == "protein", "Total prot",
  var == "tchol", "Total-C",
  var == "trig", "Triglyceride",
  var == "uric", "Urate",
  var == "urea", "Urea",
  var == "vitd25", "Vitamin D"
)]

# Reorder columns
hrs <- hrs[,.(type, model=name, sex, coefficient, coef_name, var,
              logHR, SE, HR, L95, U95, Pvalue, Proportionality.chisq, Proportionality.df, Proportionality.Pvalue)]

# Write out
fwrite(hrs, sep="\t", quote=FALSE, file="analyses/univariate/hazard_ratios.txt")

