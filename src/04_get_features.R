library(data.table)
library(foreach)
library(survival)
library(caret)
library(glmnet)
source("src/utils/factor_by_size.R")

# Load preprocessed training dataset
train <- fread("data/processed/training/processed_training_data.txt")

# Load biomarker information sheets
bio_info <- fread("data/ukb/biomarkers/output/biomarker_info.txt")
nmr_info <- fread("data/ukb/NMR_metabolomics/output/biomarker_information.txt")

# Load selected features in lasso
active <- fread("analyses/train/cox_lasso_features.txt")

# Collate information about each model
models <- foreach(with_pgs = c(TRUE, FALSE), .combine=rbind) %:%
  foreach(model_name = c("Conventional RF", "CRP", "GlycA", "Blood", "Nightingale", "Blood + Nightingale"), .combine=rbind) %:%
      foreach(this_lambda = c("lambda.min", "lambda.1se"), .combine=rbind) %do% {

        # Extract model coefficients
        mcf <- active[coef_type %in% c("blood", "nightingale") & model == tolower(gsub(" ", "", model_name)) & PGS == with_pgs & lambda.fit == this_lambda]

        # Build model formula
        mf <- "Surv(incident_followup, incident_cvd) ~ strata(sex) + age + tchol + hdl + sbp + diabetes + smoking + family_history_cvd"
        # mf <- paste(mf, "+ factor_by_size(assessment_centre) + factor_by_size(earliest_hospital_nation) + factor_by_size(latest_hospital_nation)")
        if (with_pgs) mf <- sprintf("%s + CAD_metaGRS + Stroke_metaGRS", mf)
        if (model_name == "CRP") mf <- mf <- sprintf("%s + crp", mf)
        if (model_name == "GlycA") mf <- mf <- sprintf("%s + GlycA", mf)
        if (nrow(mcf) > 0) mf <- sprintf("%s + %s", mf, paste(mcf[, coef], collapse=" + "))

        # Build long form model name
        n_bio <- mcf[,.N,by=coef_type]
        n_bio[, coef_type := paste0(toupper(substr(coef_type, 1, 1)), substr(coef_type, 2, nchar(coef_type)))]
        if (model_name %in% c("Blood", "Blood + Nightingale") & n_bio[coef_type == "Blood", .N == 0]) {
          n_bio <- rbind(n_bio, data.table(coef_type = "Blood", N = 0))
        } else if (model_name %in% c("Nightingale", "Blood + Nightingale") & n_bio[coef_type == "Nightingale", .N == 0]) {
          n_bio <- rbind(n_bio, data.table(coef_type = "Nightingale", N = 0))
        }
        n_bio[, coef_type := factor(coef_type, levels=c("Blood", "Nightingale"))]
        n_bio <- n_bio[order(n_bio)]
        bio_text <- n_bio[, paste(sprintf("%s %s biomarkers", N, coef_type), collapse=" + ")]
        if (model_name == "CRP") bio_text <- "CRP"
        if (model_name == "GlycA") bio_text <- "GlycA"

        if (with_pgs) {
           long_name <- sprintf("Conventional RF + PGS + %s", bio_text)
        } else {
           long_name <- sprintf("Conventional RF + %s", bio_text)
        }
        long_name <- gsub(" \\+ $", "", long_name) # if just conventional risk factors

				# Scale variables so HRs are per SD
				mf <- strsplit(mf, split=" \\+ ")[[1]]
				needs_scale <- !grepl("^Surv", mf) & !grepl("^strata", mf) & !grepl("^scale", mf) & !grepl("^factor", mf) & !(mf %in% c("diabetes", "smoking", "family_history_cvd"))
				mf[needs_scale] <- paste0("scale(", mf[needs_scale], ")")
				mf <- paste(mf, collapse = " + ")

        # Return model information and formula
        data.table(endpoint = "CVD", PGS = with_pgs, lambda = this_lambda, name = model_name, long_name = long_name, formula = mf)
}
models[name %in% c("Conventional RF", "CRP", "GlycA"), lambda := NA]
models <- unique(models)
fwrite(models, sep="\t", quote=FALSE, file="analyses/train/cox_lasso_models.txt")

# Extract and label model coefficients
setnames(active, "coef", "coefficient")

# Add variable names corresponding to each coefficient
active[, var := gsub("scale\\(", "", coefficient)]
active[, var := gsub("\\)", "", var)]
active[coefficient == "diabetesTRUE", var := "diabetes"]
active[coefficient == "smokingTRUE", var := "smoking"]
active[coefficient %like% "sex", var := "sex"] 
active[coefficient == "family_history_cvdTRUE", var := "family_history_cvd"]
active[coefficient %like% "assessment_centre", var := "assessment_centre"]
active[coefficient %like% "earliest_hospital_nation", var := "earliest_hospital_nation"]
active[coefficient %like% "latest_hospital_nation", var := "latest_hospital_nation"]

# Add in more human friendly names
active[coefficient == "diabetesTRUE", coef_name := "Diabetic"]
active[coefficient == "smokingTRUE", coef_name := "Smoker"]
active[coefficient == "sexMale", coef_name := "Sex (Males vs. Females)"]
active[coefficient == "family_history_cvdTRUE", coef_name := "Family history of CVD in first-degree relatives"]
active[coefficient %like% "assessment_centre", coef_name := sprintf("Assessment centre: %s (one-hot coding)", gsub("assessment_centre", "", coefficient))]
active[coefficient %like% "earliest_hospital_nation", coef_name := sprintf("Retrospective hospital nation: %s vs. England", gsub("earliest_hospital_nation", "", coefficient))]
active[coefficient %like% "latest_hospital_nation", coef_name := sprintf("Follow-up hospital nation: %s vs. England", gsub("latest_hospital_nation", "", coefficient))]
active[var %like% "metaGRS", coef_name := gsub("_", " ", var)]
active[var == "age", coef_name := "Age"]
active[var == "sbp", coef_name := "SBP"]
active[bio_info, on = .(var), coef_name := i.biomarker]
active[var == "tchol", coef_name := "Total Cholesterol"]
active[var == "hdl", coef_name := "HDL Cholesterol"]
active[nmr_info[Group %like% "Amino" & Type == "Non-derived"], on = .(var=Biomarker), coef_name := Description]
active[is.na(coef_name), coef_name := var]
active[, coef_name := gsub("_by_", " / ", coef_name)]
active[, coef_name := gsub("_pct_", " % of ", coef_name)]
active[, coef_name := gsub("_pct$", " %", coef_name)]
active[, coef_name := gsub("Total_", "Total ", coef_name)]
active[, coef_name := gsub("_", "-", coef_name)]

# Add in information on biomarker type
active[, coef_type := "Conventional RF"]
active[var %in% nmr_info$Biomarker, coef_type := "NMR Metabolomics"]
active[var %in% bio_info$var, coef_type := "Clinical Biochemistry"]
active[var %in% c("tchol", "hdl"), coef_type := "Conventional RF"]
active[var %in% c("CAD_metaGRS", "Stroke_metaGRS"), coef_type := "Polygenic Risk Score"]
active[var %in% c("assessment_centre", "earliest_hospital_nation", "latest_hospital_nation"), coef_type := "Dataset-specific covariate"]

# Add in information on measurement platform
active[var %in% c("CAD_metaGRS", "Stroke_metaGRS"), platform := "Genetics"]
active[var %in% c("age", "diabetes", "smoking", "sbp", "sex", "family_history_cvd"), platform := "Clinician"]
active[var %in% c("assessment_centre", "earliest_hospital_nation", "latest_hospital_nation"), platform := "Dataset-specific covariate"]
active[var %in% nmr_info$Biomarker, platform := "NMR spectroscopy"]
active[bio_info, on = .(var), platform := sprintf("%s method on %s instrument from %s", analysis_method, instrumentation, supplier)]
active[var == "nonhdl", platform := "CHO-POD and Enzyme immunoinhibition methods on AU5800 instrument from Beckman Coulter"]
active[var == "apobapoa1", platform := "Immunoturbidimetric method on AU5800 instrument from Beckman Coulter"]

# Add in model information
active[, model := fcase(
  model == "blood", "Blood",
  model == "nightingale", "Nightingale",
  model == "blood+nightingale", "Blood + Nightingale"
)]
active[models, on = .(lambda.fit=lambda, PGS, model=name), long_name := i.long_name]

# Reorder columns
active <- active[,.(name=model, lambda=lambda.fit, PGS, coefficient, var, coef_name, coef_type, platform, beta)]

# Write out
fwrite(active, sep="\t", quote=FALSE, file="analyses/train/lasso_coefficients.txt")



