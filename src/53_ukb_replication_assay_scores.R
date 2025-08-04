library(data.table)

# Load analaysis cohort
pheno <- fread("data/cleaned/phase3_analysis_cohort.txt")

# Load sex-specific biomarker weightings used to calculate the assay scores
assay_weights <- fread("analyses/assay_score_training/supp_table_coefficients.txt")
assay_weights <- rbind(idcol="sex",
  "Male"=assay_weights[,.(coef, scaling_mean=scaling_mean_Male, scaling_sd=scaling_sd_Male, coef_CAD=CAD_Male, coef_Stroke=Stroke_Male)],
  "Female"=assay_weights[,.(coef, scaling_mean=scaling_mean_Female, scaling_sd=scaling_sd_Female, coef_CAD=CAD_Female, coef_Stroke=Stroke_Female)]
)

# Extract assay data, filter to people with complete data, and transform to long format
assay <- pheno[,.SD,.SDcols=c("eid", "sex", unique(assay_weights$coef))]
assay <- melt(assay, id.vars=c("eid", "sex"), variable.name="biomarker", value.name="concentration")

# Standardise concentrations in males and females separately
assay[, scaled := scale(concentration), by=.(biomarker, sex)]

# Calculate assay biomarker scores 
assay[assay_weights, on = .(biomarker=coef, sex), CAD_weighted_scaled := scaled * i.coef_CAD]
assay[assay_weights, on = .(biomarker=coef, sex), Stroke_weighted_scaled := scaled * i.coef_Stroke]

assay_scores <- assay[, .(
  CAD_assay_score=sum(CAD_weighted_scaled, na.rm=TRUE),
  Stroke_assay_score=sum(Stroke_weighted_scaled, na.rm=TRUE)
), by=.(eid, sex)]

# Write out
fwrite(assay_scores, sep="\t", quote=FALSE, file="analyses/assay_score_training/phase3_assay_scores.txt")

