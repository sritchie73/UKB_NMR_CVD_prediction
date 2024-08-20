library(data.table)
library(foreach)

# Make output directory
system("mkdir -p analyses/replication/")

# Load analysis cohort
pheno <- fread("data/cleaned/analysis_cohort.txt")

# Get list of biomarkers to compute densities for
nmr_weights <- fread("analyses/nmr_score_training/avg_best_coefs.txt")
biomarkers <- sort(unique(nmr_weights$coef))

# Compute densities and save as nested list
pheno_long <- melt(pheno, id.vars=c("eid", "sex"), measure.vars=biomarkers, variable.name="biomarker", value.name="concentration")
densities <- list(
  "Male"=foreach(this_biomarker=biomarkers) %do% {
    pheno_long[sex == "Male" & biomarker == this_biomarker, density(concentration, na.rm=TRUE)]
  },
  "Female"=foreach(this_biomarker=biomarkers) %do% {
    pheno_long[sex == "Female" & biomarker == this_biomarker, density(concentration, na.rm=TRUE)]
  }
)
names(densities[[1]]) <- biomarkers
names(densities[[2]]) <- biomarkers

saveRDS(densities, file="analyses/replication/UKB_NMR_biomarker_densities.rds")

# Compute pairwise correlations between biomarkers
male_cc <- cor(as.matrix(pheno[sex == "Male", .SD, .SDcols=biomarkers]), use="pairwise.complete.obs")
female_cc <- cor(as.matrix(pheno[sex == "Female", .SD, .SDcols=biomarkers]), use="pairwise.complete.obs")
cc <- male_cc
cc[lower.tri(cc)] <- female_cc[lower.tri(female_cc)]
saveRDS(cc, file="analyses/replication/UKB_NMR_biomarker_correlations.rds")

# Compute densities for clinical chemistry biomarkers as well
biomarkers <- fread("data/ukb/biomarkers/output/biomarker_info.txt")
biomarkers <- intersect(biomarkers$var, names(pheno))

# Compute densities and save as nested list
pheno_long <- melt(pheno, id.vars=c("eid", "sex"), measure.vars=biomarkers, variable.name="biomarker", value.name="concentration")
densities <- list(
  "Male"=foreach(this_biomarker=biomarkers) %do% {
    pheno_long[sex == "Male" & biomarker == this_biomarker, density(concentration, na.rm=TRUE)]
  },
  "Female"=foreach(this_biomarker=biomarkers) %do% {
    pheno_long[sex == "Female" & biomarker == this_biomarker, density(concentration, na.rm=TRUE)]
  }
)
names(densities[[1]]) <- biomarkers
names(densities[[2]]) <- biomarkers

saveRDS(densities, file="analyses/replication/UKB_assay_biomarker_densities.rds")

# Compute pairwise correlations between biomarkers
male_cc <- cor(as.matrix(pheno[sex == "Male", .SD, .SDcols=biomarkers]), use="pairwise.complete.obs")
female_cc <- cor(as.matrix(pheno[sex == "Female", .SD, .SDcols=biomarkers]), use="pairwise.complete.obs")
cc <- male_cc
cc[lower.tri(cc)] <- female_cc[lower.tri(female_cc)]
saveRDS(cc, file="analyses/replication/UKB_assay_biomarker_correlations.rds")

# Also other covariates, e.g. age, bmi, SBP
covars <- c("age", "bmi", "sbp")
pheno_long <- melt(pheno, id.vars=c("eid", "sex"), measure.vars=covars, variable.name="covar", value.name="concentration")
densities <- list(
  "Male"=foreach(this_covar=covars) %do% {
    pheno_long[sex == "Male" & covar == this_covar, density(concentration, na.rm=TRUE)]
  },
  "Female"=foreach(this_covar=covars) %do% {
    pheno_long[sex == "Female" & covar == this_covar, density(concentration, na.rm=TRUE)]
  }
)
names(densities[[1]]) <- covars
names(densities[[2]]) <- covars

saveRDS(densities, file="analyses/replication/UKB_covariate_densities.rds")

# Also compute densities for SCORE2, NMR scores and PRSs for checking
nmr_scores <- fread("analyses/replication/aggregate_test_non_derived_NMR_scores.txt")
pheno <- pheno[nmr_scores, on = .(eid)]

densities <- list(
  "SCORE2"=list(
     "Male"=pheno[sex == "Male", density(SCORE2, na.rm=TRUE)],
     "Female"=pheno[sex == "Female", density(SCORE2, na.rm=TRUE)]
  ),
  "SCORE2_excl_UKB"=list(   
     "Male"=pheno[sex == "Male", density(SCORE2_excl_UKB, na.rm=TRUE)], 
     "Female"=pheno[sex == "Female", density(SCORE2_excl_UKB, na.rm=TRUE)]
  ),
  "CAD_metaGRS"=list(
    "Male"=pheno[sex == "Male", density(CAD_metaGRS, na.rm=TRUE)],
    "Female"=pheno[sex == "Female", density(CAD_metaGRS, na.rm=TRUE)]
  ),
  "Stroke_metaGRS"=list(
    "Male"=pheno[sex == "Male", density(Stroke_metaGRS, na.rm=TRUE)],
    "Female"=pheno[sex == "Female", density(Stroke_metaGRS, na.rm=TRUE)]
  ),
  "CAD_NMR_score"=list(
    "Male"=pheno[sex == "Male", density(CAD_NMR_score, na.rm=TRUE)],
    "Female"=pheno[sex == "Female", density(CAD_NMR_score, na.rm=TRUE)]
  ),
  "Stroke_NMR_score"=list(
    "Male"=pheno[sex == "Male", density(Stroke_NMR_score, na.rm=TRUE)],
    "Female"=pheno[sex == "Female", density(Stroke_NMR_score, na.rm=TRUE)]
  )
)
saveRDS(densities, file="analyses/replication/UKB_NMR_score_densities.rds")

# Get densities of overall combined scores
comb_scores <- fread("analyses/CVD_weight_training/CVD_linear_predictors_and_risk.txt")
comb_scores <- comb_scores[score_type == "non-derived" & model != "SCORE2"]
densities <- foreach(this_model = unique(comb_scores$model)) %do% {
  foreach(this_sex = c("Male", "Female")) %do% {
    comb_scores[sex == this_sex & model == this_model, density(linear_predictor)]
  }
}
names(densities) <- unique(comb_scores$model)
names(densities[[1]]) <- c("Male", "Female")
names(densities[[2]]) <- c("Male", "Female")
names(densities[[3]]) <- c("Male", "Female")
saveRDS(densities, file="analyses/replication/UKB_combined_score_densities.rds")

####################################################################################
# As a pre-check for the predictor patch approach, generate weightings to adjust 
# each biomarker for SCORE2 risk factors and recalculate densities and correlation
####################################################################################
biomarkers <- sort(unique(nmr_weights$coef))
pheno_long <- melt(pheno, id.vars=c("eid", "sex", "age", "smoking", "sbp", "tchol", "hdl"), 
  measure.vars=biomarkers, variable.name="biomarker", value.name="concentration")

# Apply risk factor standardisations from SCORE2
pheno_long[, age := (age - 60)/5]
pheno_long[, smoking := as.integer(smoking)]
pheno_long[is.na(smoking), smoking := 0]
pheno_long[, sbp := (sbp - 120)/20]
pheno_long[, tchol := (tchol - 6)/1]
pheno_long[, hdl := (hdl - 1.3)/0.5]

# Adjust each biomarker for risk factors, and get table of coefficients from regression models
rf_adj_coef <- foreach(this_sex = c("Male", "Female"), .combine=rbind) %do% {
  foreach(this_biomarker = biomarkers, .combine=rbind) %do% {
    l1 <- lm(concentration ~ age + factor(smoking) + sbp + tchol + hdl, data=pheno_long[sex == this_sex & biomarker == this_biomarker])
    pheno_long[sex == this_sex & biomarker == this_biomarker & !is.na(concentration), rf_adjusted := l1$residuals]
    cbind(sex=this_sex, biomarker=this_biomarker, as.data.table(as.list(l1$coef)))
  }
}
setnames(rf_adj_coef, "(Intercept)", "intercept")
setnames(rf_adj_coef, "factor(smoking)1", "smoking")
fwrite(rf_adj_coef, sep="\t", quote=FALSE, file="analyses/replication/UKB_NMR_biomarker_rf_adjustment_coefficients.txt")

# Compute densities and save as nested list
densities <- list(
  "Male"=foreach(this_biomarker=biomarkers) %do% {
    pheno_long[sex == "Male" & biomarker == this_biomarker, density(rf_adjusted, na.rm=TRUE)]
  },
  "Female"=foreach(this_biomarker=biomarkers) %do% {
    pheno_long[sex == "Female" & biomarker == this_biomarker, density(rf_adjusted, na.rm=TRUE)]
  }
)
names(densities[[1]]) <- biomarkers
names(densities[[2]]) <- biomarkers

saveRDS(densities, file="analyses/replication/UKB_NMR_biomarker_rf_adjusted_densities.rds")

# Compute pairwise correlations between biomarkers
adj_wide <- dcast(pheno_long, eid + sex ~ biomarker, value.var="rf_adjusted")
male_cc <- cor(as.matrix(adj_wide[sex == "Male", .SD, .SDcols=biomarkers]), use="pairwise.complete.obs")
female_cc <- cor(as.matrix(adj_wide[sex == "Female", .SD, .SDcols=biomarkers]), use="pairwise.complete.obs")
cc <- male_cc
cc[lower.tri(cc)] <- female_cc[lower.tri(female_cc)]
saveRDS(cc, file="analyses/replication/UKB_NMR_biomarker_rf_adjusted_correlations.rds")

