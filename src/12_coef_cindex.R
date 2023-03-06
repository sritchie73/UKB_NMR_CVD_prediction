library(data.table)
library(ukbnmr)
library(foreach)
library(survival)
source("src/utils/cox_test.R")
source("src/utils/factor_by_size.R")

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

# Load information about model coefficients
coef <- fread("analyses/train/lasso_coefficients.txt")
coef <- coef[, .(coefficient, var, coef_name, coef_type)]
coef <- unique(coef)
coef <- coef[coef_type != "Dataset-specific covariate"] # not of interest here
coef <- coef[!(var %in% c("age", "sex"))] # baseline model to compare to

# Add in other risk factors not in models by default to compare to
coef <- rbind(coef, use.names=FALSE,
  data.table("CAD_metaGRS", "CAD_metaGRS", "CAD metaGRS", "PRS"),
  data.table("Stroke_metaGRS", "Stroke_metaGRS", "Stroke metaGRS", "PRS")
)

# Build models
mf1 <- "Surv(incident_followup, incident_cvd) ~ strata(sex) + age"
mf2 <- "Surv(incident_followup, incident_cvd) ~ strata(sex) + age + %s"

# Get C-indices
cind1 <- cbind(model="reference", cox.test(mf1, "incident_cvd", test)$model_fit)
cinds <- foreach(this_var = coef[,unique(var)], .combine=rbind, .init=cind1) %do% {
  # Scale variable if not a factor
  varterm <- if (coef[var == this_var, var != coefficient]) { this_var } else { sprintf("scale(%s)", this_var) }

  # Fit cox model
  cph <- cox.test(sprintf(mf2, varterm), "incident_cvd", test) 

  # Extract cindex
  cbind(model=this_var, cph$model_fit)
}

# add in extra details
cinds[coef, on = .(model=var), model_name := coef_name]
cinds[coef, on = .(model=var), model_type := coef_type]
cinds[model == "reference", model_name := "Reference"]
cinds[model == "reference", model_type := "Reference"]

# Apply ordering
cinds[, model_type := factor(model_type, levels=c("Reference", "Conventional RF", "PRS", "NMR Metabolomics", "Clinical Biochemistry"))]
cinds <- cinds[order(C.index)][order(model_type)]
cinds[, model_name := factor(model_name, levels=model_name)]

# Write out
fwrite(cinds, sep="\t", quote=FALSE, file="analyses/test/age_sex_feature_cindices.txt")

# Generate plots of C-indices
g <- ggplot(cinds) +
	aes(x=model_name, y = C.index, ymin = L95, ymax = U95) +
	geom_errorbar(width=0) +
	geom_point(size=2, shape=18) +
	geom_hline(yintercept = cinds[model == "reference", C.index], linetype = 2) +
	facet_grid(. ~ model_type, scales="free_x", space="free_x") +
	ylab("C-index (95% CI)") +
	xlab("") +
	theme_bw() +
	theme(legend.position="bottom", legend.box="vertical", axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

ggsave(g, width=7.2, height=5, units="in", file="analyses/test/age_sex_feature_cindices.pdf")


