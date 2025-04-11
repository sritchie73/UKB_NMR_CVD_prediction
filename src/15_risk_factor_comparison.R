library(data.table)
library(foreach)
library(survival)
library(ggplot2)
library(patchwork)

# Make output directory
system("mkdir -p analyses/test", wait=TRUE)

# Load discovery data
dat <- fread("data/cleaned/analysis_cohort.txt")
dat2 <- fread("data/cleaned/sensitivity_analysis_extended_data.txt")

# Load and add predicted NMR scores in the discovery data
nmr_scores <- fread("analyses/nmr_score_training/aggregate_test_non_derived_NMR_scores.txt")
dat <- dat[nmr_scores, on = .(eid)]

# Load FDR-significant clinical chemistry biomarkers
bio_res <- fread("analyses/univariate/cindices_sensitivity_analysis.txt")
bio_res <- bio_res[cohort == "pooled" & sex == "Sex-stratified" & endpoint == "ASCVD" & score == "SCORE2" & model_type == "Clinical biochemistry assay" & deltaC.fdr < 0.05]

# Fit joint models with risk factors
rf_hrs <- foreach(this_sex = c("Males", "Females", "Sex-stratified"), .combine=rbind) %:% 
  foreach(this_score = c("SCORE2", "QRISK3"), .combine=rbind) %:% 
    foreach(this_model = c("biochemistry", "NMR"), .combine=rbind) %do% {
      # Extract relevant model data
      if (this_model == "NMR") {
         this_dat <- dat[,.(eid, CAD_metaGRS, Stroke_metaGRS, CAD_NMR_score, Stroke_NMR_score)]
      } else {
         this_dat <- dat[,.SD,.SDcols=c("eid", "CAD_metaGRS", "Stroke_metaGRS", bio_res$biomarker)]
      }

      # Add common model data
      this_dat <- this_dat[dat[,.(eid, sex, incident_cvd, incident_cvd_followup)], on = .(eid)]

      # Add in score-specific risk factors
      if (this_score == "SCORE2") {
         this_dat <- this_dat[dat[,.(eid, age, sbp, tchol, hdl, smoking)], on = .(eid)]
      } else {
         this_dat <- this_dat[dat2[,.(eid, age, sbp, sd_sbp, tchol_hdl_ratio, smoke, townsend, QRisk_ethnicity,
           weight, height, atrial_fibrillation, chronic_kidney_disease, severe_mental_illness, migraine,
           systemic_lupus_erythematosis, erectile_dysfunction, rheumatoid_arthritis, prevalent_t1d, prevalent_t2d,
           atypical_antipsychotics, systematic_corticosteroids, blood_pressure_treatment)], on = .(eid)]
      }

      # Filter to relevant sex
      if (this_sex == "Males") {
        this_dat <- this_dat[sex == "Male"]
      } else if (this_sex == "Females") {
        this_dat <- this_dat[sex == "Female"]
      }

      # Apply standardisations etc
      this_dat[, CAD_metaGRS := scale(CAD_metaGRS)]
      this_dat[, Stroke_metaGRS := scale(Stroke_metaGRS)]
      this_dat[, sbp := scale(sbp)]
     
      if (this_model == "NMR") {
        this_dat[, CAD_NMR_score := scale(CAD_NMR_score)]
        this_dat[, Stroke_NMR_score := scale(Stroke_NMR_score)]
      } else {
        for (bio_var in bio_res$biomarker) {
          this_dat[, c(bio_var) := as.vector(scale(this_dat[[bio_var]]))]
        }
      }

      if (this_score == "SCORE2") {
        this_dat[is.na(smoking), smoking := FALSE]
        this_dat[, age := scale(age), by=sex]
        this_dat[, hdl := scale(hdl), by=sex]
        this_dat[, tchol := scale(tchol), by=sex]
      } else {
        this_dat[, sd_sbp := scale(sd_sbp), by=sex]
        this_dat[, tchol_hdl_ratio := scale(tchol_hdl_ratio), by=sex]
        this_dat[, townsend := scale(townsend), by=sex]

        this_dat[, smoke := factor(smoke, levels=c("non-smoker", "ex-smoker", "light smoker", "moderate smoker", "heavy smoker"))]
        this_dat[, QRisk_ethnicity := factor(QRisk_ethnicity, levels=c(
          "White or not stated", "Indian", "Pakistani", "Bangladeshi", "Other Asian",
          "Black Caribbean", "Black African", "Chinese", "Other ethnic group"
        ))]

        this_dat[, age_1 := ifelse(sex == "Female", age^-2, age^-1)]
        this_dat[, age_2 := ifelse(sex == "Female", age, age^3)]
        this_dat[, age_1 := scale(age_1), by=sex]
        this_dat[, age_2 := scale(age_2), by=sex]

        this_dat[, height := height/100]
        this_dat[, bmi := weight / height^2]
			  this_dat[, bmi_1 := bmi^-2]
			  this_dat[, bmi_2 := bmi^-2 * log(bmi)]
        this_dat[, bmi_1 := scale(bmi_1), by=sex]
        this_dat[, bmi_2 := scale(bmi_2), by=sex]
      }

      # Build model formula
      mf <- "Surv(incident_cvd_followup, incident_cvd) ~"
 
      if (this_sex == "Sex-stratified") { 
        mf <- paste(mf, "strata(sex) +")
      }
    
      mf <- paste(mf, "CAD_metaGRS + Stroke_metaGRS +")

      if (this_model == "NMR") {
        mf <- paste(mf, "CAD_NMR_score + Stroke_NMR_score +")
      } else {
        mf <- paste(mf, paste(bio_res$biomarker, collapse=" + "), "+")
      }

      if (this_score == "SCORE2") {
        mf <- paste(mf, "age*smoking + age*sbp + age*tchol + age*hdl")
      } else {
        if (this_sex != "Females") {
          mf <- paste(mf, "age_1*erectile_dysfunction + age_2*erectile_dysfunction +")
        }
        mf <- paste(mf, 
          "age_1*smoke + age_1*atrial_fibrillation + age_1*systematic_corticosteroids + age_1*migraine +",
          "age_1*chronic_kidney_disease + age_1*systemic_lupus_erythematosis + age_1*blood_pressure_treatment +",
          "age_1*bmi_1 + age_1*bmi_2 + age_1*sbp + age_1*townsend +",
          "age_2*smoke + age_2*atrial_fibrillation + age_2*systematic_corticosteroids + age_2*migraine +",
          "age_2*chronic_kidney_disease + age_2*systemic_lupus_erythematosis + age_2*blood_pressure_treatment +",
          "age_2*bmi_1 + age_2*bmi_2 + age_2*sbp + age_2*townsend + QRisk_ethnicity +",
          "tchol_hdl_ratio + sd_sbp + atypical_antipsychotics + rheumatoid_arthritis + severe_mental_illness")
      }

      # Fit cox proportional hazards model
      cx <- coxph(as.formula(mf), data=this_dat)

		  # Extra hazard ratios
	    cf <- coef(summary(cx))
		  ci <- confint(cx)

		  HR <- cf[,2]
		  HR.SE <- cf[,3]
		  HR.L95 <- exp(ci[,1])
		  HR.U95 <- exp(ci[,2])
		  HR.pval <- cf[,5]

			# Extract relevant info and return
			data.table(
        score = this_score, model_sex = this_sex, model_type = this_model,
				samples=cx$n, events=cx$nevent, coefficient = names(HR), HR, HR.SE, HR.L95, HR.U95, HR.pval
			)
}

# Give readable names to coefficients
rf_hrs[, coefficient := fcase(
  coefficient == "CAD_metaGRS", "CHD PRS (PGS000018)",
  coefficient == "Stroke_metaGRS", "Stroke PRS (PGS000039)",
  coefficient == "CAD_NMR_score", "CHD NMR score",
  coefficient == "Stroke_NMR_score", "IS NMR score",
  coefficient == "cyst", "Cystatin-C",
  coefficient == "crp", "C-reactive protein",
  coefficient == "alp", "Alkaline phosphatase",
  coefficient == "alb", "Albumin",
  coefficient == "ggt", "Gamma glutamyltransferase",
  coefficient == "lpa", "Lipoprotein(a)",
  coefficient == "asp", "Aspartate aminotransferase",
  coefficient == "hba1c", "Glycated haemoglobin (HbA1c)",
  coefficient == "uric", "Urate",
  coefficient == "vitd25", "Vitamin D",
  coefficient == "apoa1", "Apolipoprotein A1",
  coefficient == "age", "Age",
  coefficient == "smokingTRUE", "Smoker",
  coefficient == "sbp", "SBP",
  coefficient == "tchol", "Total cholesterol",
  coefficient == "hdl", "HDL cholesterol",
  coefficient == "age:smokingTRUE", "Age x smoker",
  coefficient == "age:sbp", "Age x SBP",
  coefficient == "age:tchol", "Age x Total cholesterol",
  coefficient == "age:hdl", "Age x HDL cholesterol",
  coefficient == "age_1", "Age P1",
  coefficient == "age_2", "Age P2", 
  coefficient == "erectile_dysfunctionTRUE", "Erectile dysfunction",
  coefficient == "smokeex-smoker", "Ex-smoker",
  coefficient == "smokelight smoker", "Light smoker",
  coefficient == "smokemoderate smoker", "Moderate smoker",
  coefficient == "smokeheavy smoker", "Heavy smoker",
  coefficient == "atrial_fibrillationTRUE", "Atrial fibrillation",
  coefficient == "systematic_corticosteroidsTRUE", "Corticosteroid usage",
  coefficient == "migraineTRUE", "Migraines",
  coefficient == "chronic_kidney_diseaseTRUE", "Chronic kidney disease",
  coefficient == "systemic_lupus_erythematosisTRUE", "Systemic lupus erythematosis",
  coefficient == "blood_pressure_treatmentTRUE", "Blood pressure medication",
  coefficient == "bmi_1", "BMI P1",
  coefficient == "bmi_2", "BMI P2",
  coefficient == "townsend", "Townsend deprivation index",
  coefficient == "tchol_hdl_ratio", "Ratio of total to HDL cholesterol",
  coefficient == "sd_sbp", "Standard deviation of SBP",
  coefficient == "atypical_antipsychoticsTRUE", "Atypical antipsychotic use",
  coefficient == "rheumatoid_arthritisTRUE", "Rheumatoid arthritis",
  coefficient == "severe_mental_illnessTRUE", "Severe mental illness",
  coefficient == "age_1:erectile_dysfunctionTRUE", "Age P1 x Erectile dysfunction",
  coefficient == "erectile_dysfunctionTRUE:age_2", "Age P2 x Erectile dysfunction",
  coefficient == "age_1:smokeex-smoker", "Age P1 x Ex-smoker",
  coefficient == "age_1:smokelight smoker", "Age P1 x Light smoker",
  coefficient == "age_1:smokemoderate smoker", "Age P1 x Moderate smoker",
  coefficient == "age_1:smokeheavy smoker", "Age P1 x Heavy smoker",
  coefficient == "age_1:atrial_fibrillationTRUE", "Age P1 x Atrial fibrillation",
  coefficient == "age_1:systematic_corticosteroidsTRUE", "Age P1 x Corticosteroid usage",
  coefficient == "age_1:migraineTRUE", "Age P1 x Migraine",
  coefficient == "age_1:chronic_kidney_diseaseTRUE", "Age P1 x Chronic kideny disease",
  coefficient == "age_1:systemic_lupus_erythematosisTRUE", "Age P1 x Systemic lupus erythematosis",
  coefficient == "age_1:blood_pressure_treatmentTRUE", "Age P1 x Blood pressure medication",
  coefficient == "age_1:bmi_1", "Age P1 x BMI P1",
  coefficient == "age_1:bmi_2", "Age P1 x BMI P2",
  coefficient == "age_1:sbp", "Age P1 x SBP",
  coefficient == "age_1:townsend", "Age P1 x Townsend deprivation index",
  coefficient == "age_2:smokeex-smoker", "Age P2 x Ex-smoker",
  coefficient == "age_2:smokelight smoker", "Age P2 x Light smoker",
  coefficient == "age_2:smokemoderate smoker", "Age P2 x Moderate smoker",
  coefficient == "age_2:smokeheavy smoker", "Age P2 x Heavy smoker",
  coefficient == "age_2:atrial_fibrillationTRUE", "Age P2 x Atrial fibrillation",
  coefficient == "age_2:systematic_corticosteroidsTRUE", "Age P2 x Corticosteroid usage",
  coefficient == "age_2:migraineTRUE", "Age P2 x Migraine",
  coefficient == "age_2:chronic_kidney_diseaseTRUE", "Age P2 x Chronic kidney disease",
  coefficient == "age_2:systemic_lupus_erythematosisTRUE", "Age P2 x Systemic lupus erythematosis",
  coefficient == "age_2:blood_pressure_treatmentTRUE", "Age P2 x Blood pressure medication",
  coefficient == "age_2:bmi_1", "Age P2 x BMI P1",
  coefficient == "age_2:bmi_2", "Age P2 x BMI P2",
  coefficient == "age_2:sbp", "Age P2 x SBP",
  coefficient == "age_2:townsend", "Age P2 x Townsend deprivation index",
  coefficient == "smokeex-smoker:age_2", "Age P2 x Ex-smoker",
  coefficient == "smokelight smoker:age_2", "Age P2 x Light smoker",
  coefficient == "smokemoderate smoker:age_2", "Age P2 x Moderate smoker",
  coefficient == "smokeheavy smoker:age_2", "Age P2 x Heavy smoker",
  coefficient == "atrial_fibrillationTRUE:age_2", "Age P2 x Atrial fibrillation",
  coefficient == "systematic_corticosteroidsTRUE:age_2", "Age P2 x Corticosteroid usage",
  coefficient == "migraineTRUE:age_2", "Age P2 x Migraine",
  coefficient == "chronic_kidney_diseaseTRUE:age_2", "Age P2 x Chronic kidney disease",
  coefficient == "systemic_lupus_erythematosisTRUE:age_2", "Age P2 x Systemic lupus erythematosis",
  coefficient == "blood_pressure_treatmentTRUE:age_2", "Age P2 x Blood pressure medication",
  coefficient == "bmi_1:age_2", "Age P2 x BMI P1",
  coefficient == "bmi_2:age_2", "Age P2 x BMI P2",
  coefficient == "sbp:age_2", "Age P2 x SBP",
  coefficient == "townsend:age_2", "Age P2 x Townsend",
  coefficient == "QRisk_ethnicityIndian", "Indian ethnicity",
  coefficient == "QRisk_ethnicityPakistani", "Pakistani ethnicity",
  coefficient == "QRisk_ethnicityBangladeshi", "Bangladeshi ethnicity",
  coefficient == "QRisk_ethnicityOther Asian", "Other Asian ethnicity",
  coefficient == "QRisk_ethnicityBlack Caribbean", "Black Caribbean ethnicity",
  coefficient == "QRisk_ethnicityBlack African", "Black African ethnicity",
  coefficient == "QRisk_ethnicityChinese", "Chinese ethnicity",
  coefficient == "QRisk_ethnicityOther ethnic group", "Other ethnicity"
)]

# write out supp table
rf_hrs[, model_sex := factor(model_sex, levels=c("Sex-stratified", "Males", "Females"))]
rf_hrs[, score := factor(score, levels=c("SCORE2", "QRISK3"))]
rf_hrs[, model_type := fcase(
  model_type == "NMR", "Risk score + NMR scores + PRSs",
  model_type == "biochemistry", "Risk score + clinical biomarkers + PRSs"
)]
rf_hrs[, model_type := factor(model_type, levels=c("Risk score + NMR scores + PRSs", "Risk score + clinical biomarkers + PRSs"))]
rf_hrs <- rf_hrs[order(model_sex)]
rfo <- rf_hrs[model_sex == "Sex-stratified"][order(HR.pval)][order(model_type)][order(score)][,.(score, model_type, coefficient)]
rf_hrs <- rf_hrs[rfo, on = .(score, model_type, coefficient)]
fwrite(rf_hrs, sep="\t", quote=FALSE, file="analyses/test/multivariable_risk_factor_associations.txt")

# Create sex-stratified plot
hr_plot <- function(ggdt) {
  ggdt[, xorder := factor(1:.N), by=model_sex]

	ggplot(ggdt) +
		aes(x=xorder, y=HR, ymin=HR.L95, ymax=HR.U95, color=model_sex) +
		geom_hline(yintercept=1, linetype=2) +
		geom_errorbar(width=0, position=position_dodge(width=0.6)) +
		geom_point(shape=23, fill="white", position=position_dodge(width=0.6)) +
		scale_color_manual("Sex", values=c("Males"="#e41a1c", "Females"="#377eb8", "Sex-stratified"="#006d2c")) +
    scale_x_discrete(labels=structure(unique(ggdt$coefficient), names=unique(as.character(ggdt$xorder)))) +
		ylab("HR per SD (95% CI)") +
		theme_bw() +
		theme(
			axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5), axis.text.y=element_text(size=6),
			axis.title.x=element_blank(), axis.title.y=element_text(size=8),
			panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
			legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=8)
		)
}

g1 <- hr_plot(rf_hrs[score == "SCORE2" & model_type == "Risk score + NMR scores + PRSs"])
g2 <- hr_plot(rf_hrs[score == "SCORE2" & model_type == "Risk score + clinical biomarkers + PRSs"])
g <- g1 / g2
ggsave(g, width=7.2, height=8, file="analyses/test/multivariable_risk_factor_associations.pdf")

