library(data.table)
library(foreach)
library(survival)
source("src/utils/cox_test.R")
source("src/utils/factor_by_size.R")

# Create output directory
system("mkdir -p analyses/test")

# Load biomarker information sheets
bio_info <- fread("data/ukb/biomarkers/output/biomarker_info.txt")
nmr_info <- fread("data/ukb/NMR_metabolomics/output/biomarker_information.txt")

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

# Load information about models from feature selection:
models <- fread("analyses/train/cox_lasso_models.txt")

# Fit cox proportional hazards models using selected features
cox_list <- lapply(models$formula, cox.test, "incident_cvd", test)

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
hrs[coefficient %like% "assessment_centre", var := "assessment_centre"]
hrs[coefficient %like% "earliest_hospital_nation", var := "earliest_hospital_nation"]
hrs[coefficient %like% "latest_hospital_nation", var := "latest_hospital_nation"]

# Add in more human friendly names
hrs[coefficient == "diabetesTRUE", coef_name := "Diabetic"]
hrs[coefficient == "smokingTRUE", coef_name := "Smoker"]
hrs[coefficient == "family_history_cvdTRUE", coef_name := "Family history of CVD in first-degree relatives"]
hrs[coefficient %like% "assessment_centre", coef_name := sprintf("Assessment centre: %s vs. Leeds", gsub("factor_by_size\\(assessment_centre\\)", "", coefficient))]
hrs[coefficient %like% "earliest_hospital_nation", coef_name := sprintf("Retrospective hospital nation: %s vs. England", gsub("factor_by_size\\(earliest_hospital_nation\\)", "", coefficient))]
hrs[coefficient %like% "latest_hospital_nation", coef_name := sprintf("Follow-up hospital nation: %s vs. England", gsub("factor_by_size\\(latest_hospital_nation\\)", "", coefficient))]
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

# Add in information on biomarker type
hrs[, coef_type := "Conventional RF"]
hrs[var %in% nmr_info$Biomarker, coef_type := "NMR Metabolomics"]
hrs[var %in% bio_info$var, coef_type := "Clinical Biochemistry"]
hrs[var %in% c("tchol", "hdl"), coef_type := "Conventional RF"]
hrs[var %in% c("CAD_metaGRS", "Stroke_metaGRS"), coef_type := "Polygenic Risk Score"]
hrs[var %in% c("assessment_centre", "earliest_hospital_nation", "latest_hospital_nation"), coef_type := "Dataset-specific covariate"]

# Add in information on measurement platform
hrs[var %in% c("CAD_metaGRS", "Stroke_metaGRS"), platform := "Genetics"]
hrs[var %in% c("age", "diabetes", "smoking", "sbp", "sex", "family_history_cvd"), platform := "Clinician"]
hrs[var %in% c("assessment_centre", "earliest_hospital_nation", "latest_hospital_nation"), platform := "Dataset-specific covariate"]
hrs[var %in% nmr_info$Biomarker, platform := "NMR spectroscopy"]
hrs[bio_info, on = .(var), platform := sprintf("%s method on %s instrument from %s", analysis_method, instrumentation, supplier)]
hrs[var == "nonhdl", platform := "CHO-POD and Enzyme immunoinhibition methods on AU5800 instrument from Beckman Coulter"]
hrs[var == "apobapoa1", platform := "Immunoturbidimetric method on AU5800 instrument from Beckman Coulter"]

# Reorder columns
hrs <- hrs[,.(name, lambda, PGS, long_name, coefficient, var, coef_name, coef_type, platform, 
              logHR, SE, HR, L95, U95, Pvalue, Proportionality.chisq, Proportionality.df, Proportionality.Pvalue)]

# Write out
fwrite(hrs, sep="\t", quote=FALSE, file="analyses/test/hazard_ratios.txt")

# Extract C-indices
cinds <- rbindlist(lapply(cox_list, `[[`, 2), idcol="model_idx")

# Add in key model information
cinds <- models[cinds, on = .(model_idx)]
cinds[, c("model_idx", "formula") := NULL]

# Get delta C-index relative to conventional risk factors alone
base_cind <- cinds[!(PGS) & name == "Conventional RF", C.index]
cinds[, deltaC := C.index - base_cind]
cinds[, deltaC.L95 := L95 - base_cind]
cinds[, deltaC.U95 := U95 - base_cind]

# Get detla C-index relative to conventional risk factors + PGS
base_cind <- cinds[(PGS) & name == "Conventional RF", C.index]
cinds[, deltaC.PGS := C.index - base_cind]
cinds[, deltaC.PGS.L95 := L95 - base_cind]
cinds[, deltaC.PGS.U95 := U95 - base_cind]

# Reorder columns
cinds <- cinds[,.(Samples, Cases, Missing, name, lambda, PGS, long_name, C.index, SE, L95, U95, deltaC, deltaC.L95, 
                  deltaC.U95, deltaC.PGS, deltaC.PGS.L95, deltaC.PGS.U95)]

# Write out
fwrite(cinds, sep="\t", quote=FALSE, file="analyses/test/C_indices.txt")


# Extract other model fit information
model_info <- rbindlist(lapply(cox_list, `[[`, 2), idcol="model_idx")

# Add in key model information
model_info <- models[model_info, on = .(model_idx)]
model_info[, model_idx := NULL]

# Reorder columns
model_info <- model_info[, .(formula, name, lambda, PGS, long_name, Samples, Cases, Missing, Missing.Cases,
                             Proportionality.chisq, Proportionality.df, Proportionality.Pvalue)]

# Write out 
fwrite(model_info, sep="\t", quote=FALSE, file="analyses/test/model_fit_information.txt")

