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




