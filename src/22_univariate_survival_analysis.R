library(data.table)
library(ukbnmr)
library(foreach)
library(survival)
source("src/utils/cox_test.R")
source("src/utils/factor_by_size.R")

# Create output directory
system("mkdir -p analyses/test/univariate")

# Load biomarker information sheets
bio_info <- fread("data/ukb/biomarkers/output/biomarker_info.txt")
nmr_info <- ukbnmr::nmr_info

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

# Build set of models
models <- rbind(
  data.table(formula="Surv(incident_followup, incident_cvd) ~ sex", type="demographics", name="sex"),
  data.table(formula="Surv(incident_followup, incident_cvd) ~ scale(age)", type="demographics", name="age"),
  data.table(formula="Surv(incident_followup, incident_cvd) ~ strata(sex) + scale(age)", type="demographics", name="age+sex"),
  data.table(formula="Surv(incident_followup, incident_cvd) ~ strata(sex) + scale(age) + diabetes", type="risk_factors", name="diabetes"),
  data.table(formula="Surv(incident_followup, incident_cvd) ~ strata(sex) + scale(age) + smoking", type="risk_factors", name="smoking"),
  data.table(formula="Surv(incident_followup, incident_cvd) ~ strata(sex) + scale(age) + family_history_cvd", type="risk_factors", name="family_history_cvd"),
  data.table(formula="Surv(incident_followup, incident_cvd) ~ strata(sex) + scale(age) + scale(sbp)", type="risk_factors", name="sbp"),
  data.table(formula="Surv(incident_followup, incident_cvd) ~ strata(sex) + scale(age) + scale(tchol)", type="risk_factors", name="tchol"),
  data.table(formula="Surv(incident_followup, incident_cvd) ~ strata(sex) + scale(age) + scale(hdl)", type="risk_factors", name="hdl")
)
conv_rf <- "Surv(incident_followup, incident_cvd) ~ strata(sex) + scale(age) + scale(sbp) + diabetes + smoking + family_history_cvd + scale(tchol) + scale(hdl)"
clin_chem <- bio_info[!is.na(UKB.Field.ID) & sample_type != "Urine" & var != "tchol" & var != "hdl", var]
models <- rbind(models,
  data.table(formula=sprintf("%s + scale(%s)", conv_rf, c("CAD_metaGRS", "Stroke_metaGRS")), type="PRS", name=c("CAD_metaGRS", "Stroke_metaGRS")),
  data.table(formula=sprintf("%s + scale(CAD_metaGRS) + scale(Stroke_metaGRS)", conv_rf), type="PRS", name="PRS"),
  data.table(formula=sprintf("%s + scale(%s)", conv_rf, clin_chem), type="assays", name=clin_chem),
  data.table(formula=sprintf("%s + scale(%s)", conv_rf, nmr_info$Biomarker), type="NMR", name=nmr_info$Biomarker)
)

# Fit cox proportional hazards models
cox_list <- foreach(mf = models$formula) %dopar% {
  cox.test(mf, "incident_cvd", test)
}
saveRDS(cox_list, file="analyses/test/univariate_cox_list.rds")
stop()

# Extract hazard ratios
hrs <- rbindlist(lapply(cox_list, `[[`, 1), idcol="model_idx")

# Add in key model information
models[, model_idx := .I]
hrs <- models[hrs, on = .(model_idx)]
hrs[, c("model_idx", "formula") := NULL]

# Add variable names corresponding to each coefficient
hrs[, var := gsub("scale\\(", "", coefficient)]
hrs[, var := gsub("\\)", "", var)]
hrs[coefficient == "diabetesTRUE", var := "diabetes"]
hrs[coefficient == "smokingTRUE", var := "smoking"]
hrs[coefficient %like% "sex", var := "sex"]
hrs[coefficient == "family_history_cvdTRUE", var := "family_history_cvd"]

# Add in more human friendly names
hrs[coefficient == "diabetesTRUE", coef_name := "Diabetic"]
hrs[coefficient == "smokingTRUE", coef_name := "Smoker"]
hrs[coefficient == "family_history_cvdTRUE", coef_name := "Family history of CVD in first-degree relatives"]
hrs[var %like% "metaGRS", coef_name := gsub("_", " ", var)]
hrs[var == "age", coef_name := "Age"]
hrs[var == "sbp", coef_name := "SBP"]
hrs[bio_info, on = .(var), coef_name := i.biomarker]
hrs[var == "tchol", coef_name := "Total Cholesterol"]
hrs[var == "hdl", coef_name := "HDL Cholesterol"]
hrs[nmr_info[Group %like% "Amino" & Type == "Non-derived"], on = .(var=Biomarker), coef_name := Description]
hrs[is.na(coef_name), coef_name := var]
hrs[, coef_name := gsub("_by_", " / ", coef_name)]
hrs[, coef_name := gsub("_pct_", " % of ", coef_name)]
hrs[, coef_name := gsub("_pct$", " %", coef_name)]
hrs[, coef_name := gsub("Total_", "Total ", coef_name)]
hrs[, coef_name := gsub("_", "-", coef_name)]
hrs[coef_name %in% clin_chem, coef_name := fcase(
  coef_name == "alt", "ALT",
  coef_name == "alb", "Albumin",
  coef_name == "alp", "ALP",
  coef_name == "apoa1", "ApoA1",
  coef_name == "apob", "ApoB",
  coef_name == "asp", "ASP",
  coef_name == "calcium", "Calcium", 
  coef_name == "creat", "Creatinine",
  coef_name == "crp", "CRP",
  coef_name == "cyst", "CST3",
  coef_name == "dbili", "Direct bili",
  coef_name == "ggt", "GGT", 
  coef_name == "glucose", "Glucose",
  coef_name == "hba1c", "HbA1c",
  coef_name == "igf1", "IGF-1", 
  coef_name == "lpa", "Lp(a)",
  coef_name == "ldl", "LDL-C",
  coef_name == "hdl", "HDL-C",
  coef_name == "oest", "Oestradiol",
  coef_name == "phos", "Phosphate",
  coef_name == "rheuf", "RF",
  coef_name == "shbg", "SHBG",
  coef_name == "testos", "Testosterone",
  coef_name == "tbili", "Total bili",
  coef_name == "protein", "Total prot",
  coef_name == "tchol", "Total-C",
  coef_name == "trig", "Triglyceride",
  coef_name == "uric", "Urate",
  coef_name == "urea", "Urea",
  coef_name == "vitd25", "Vitamin D"
)]

# Reorder columns
hrs <- hrs[,.(type, var, name, coefficient, coef_name, 
              logHR, SE, HR, L95, U95, Pvalue, Proportionality.chisq, Proportionality.df, Proportionality.Pvalue)]

# Write out
fwrite(hrs, sep="\t", quote=FALSE, file="analyses/test/univariate/hazard_ratios.txt")

# Extract C-indices
cinds <- rbindlist(lapply(cox_list, `[[`, 2), idcol="model_idx")

# Add in key model information
cinds <- models[cinds, on = .(model_idx)]
cinds[, c("model_idx", "formula") := NULL]

# Add in more human friendly names
cinds[, display_name := name]
cinds[name == "age", display_name := "Age"]
cinds[name == "sex", display_name := "Sex"]
cinds[name == "age+sex", display_name := "Age + Sex"]
cinds[name == "diabetes", display_name := "Diabetic"]
cinds[name == "smoking", display_name := "Smoker"]
cinds[name == "family_history_cvd", display_name := "Family history"]
cinds[name == "sbp", display_name := "SBP"]
cinds[name == "hdl", display_name := "HDL-C"]
cinds[name == "tchol", display_name := "Total-C"]
cinds[display_name %like% "metaGRS", display_name := gsub("_", " ", display_name)]
cinds[nmr_info[Group %like% "Amino" & Type == "Non-derived"], on = .(name=Biomarker), display_name := Description]
cinds[, display_name := gsub("_by_", " / ", display_name)]
cinds[, display_name := gsub("_pct_", " % of ", display_name)]
cinds[, display_name := gsub("_pct$", " %", display_name)]
cinds[, display_name := gsub("Total_", "Total ", display_name)]
cinds[, display_name := gsub("_", "-", display_name)]
cinds[name %in% clin_chem, display_name := fcase(
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
  name == "oest", "Oestradiol",
  name == "phos", "Phosphate",
  name == "rheuf", "RF",
  name == "shbg", "SHBG",
  name == "testos", "Testosterone",
  name == "tbili", "Total bili",
  name == "protein", "Total prot",
  name == "trig", "Triglyceride",
  name == "uric", "Urate",
  name == "urea", "Urea",
  name == "vitd25", "Vitamin D"
)]

# Get delta C-index
cinds[type == "risk_factors", deltaC := C.index - cinds[name == "age+sex", C.index]]
cinds[type == "risk_factors", deltaC.L95 := L95 - cinds[name == "age+sex", C.index]]
cinds[type == "risk_factors", deltaC.U95 := U95 - cinds[name == "age+sex", C.index]]

conv_rf_c <- fread("analyses/test/C_indices.txt")

cinds[type != "risk_factors" & type != "demographics", deltaC := C.index - conv_rf_c[!(PGS) & name == "Conventional RF", C.index]]
cinds[type != "risk_factors" & type != "demographics", deltaC.L95 := L95 - conv_rf_c[!(PGS) & name == "Conventional RF", C.index]] 
cinds[type != "risk_factors" & type != "demographics", deltaC.U95 := U95 - conv_rf_c[!(PGS) & name == "Conventional RF", C.index]]

# Reorder columns
cinds <- cinds[,.(type, name, display_name, Samples, Cases, Missing, C.index, SE, L95, U95, deltaC, deltaC.L95, deltaC.U95)]

# Write out
fwrite(cinds, sep="\t", quote=FALSE, file="analyses/test/univariate/C_indices.txt")

# Extract other model fit information
model_info <- rbindlist(lapply(cox_list, `[[`, 2), idcol="model_idx")

# Add in key model information
model_info <- models[model_info, on = .(model_idx)]
model_info[, model_idx := NULL]

# Write out 
fwrite(model_info, sep="\t", quote=FALSE, file="analyses/test/univariate/model_fit_information.txt")

