library(data.table)

# Load analaysis cohort
pheno <- fread("data/cleaned/phase3_analysis_cohort.txt")

# Load sex-specific biomarker weightings used to calculate the NMR scores
nmr_weights <- fread("analyses/nmr_score_training/supp_table_coefficients.txt")
nmr_weights <- rbind(idcol="sex",
  "Male"=nmr_weights[,.(coef, scaling_mean=scaling_mean_Male, scaling_sd=scaling_sd_Male, coef_CAD=CAD_Male, coef_Stroke=Stroke_Male)],
  "Female"=nmr_weights[,.(coef, scaling_mean=scaling_mean_Female, scaling_sd=scaling_sd_Female, coef_CAD=CAD_Female, coef_Stroke=Stroke_Female)]
)

# Extract nmr data, filter to people with complete data, and transform to long format
nmr <- pheno[,.SD,.SDcols=c("eid", "sex", unique(nmr_weights$coef))]
nmr <- melt(nmr, id.vars=c("eid", "sex"), variable.name="biomarker", value.name="concentration")

# Standardise concentrations in males and females separately
nmr[, scaled := scale(concentration), by=.(biomarker, sex)]

# Calculate NMR biomarker scores 
nmr[nmr_weights, on = .(biomarker=coef, sex), CAD_weighted_scaled := scaled * i.coef_CAD]
nmr[nmr_weights, on = .(biomarker=coef, sex), Stroke_weighted_scaled := scaled * i.coef_Stroke]

nmr_scores <- nmr[, .(
  CAD_NMR_score=sum(CAD_weighted_scaled, na.rm=TRUE),
  Stroke_NMR_score=sum(Stroke_weighted_scaled, na.rm=TRUE)
), by=.(eid, sex)]

# Write out
fwrite(nmr_scores, sep="\t", quote=FALSE, file="analyses/nmr_score_training/phase3_nmr_scores.txt")

