library(data.table)
library(foreach)
library(survival)
source("src/utils/cox_test.R")
source("src/utils/factor_by_size.R")

# Create output directory
system("mkdir -p analyses/test/feature_cindices/")

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

# Load information about model coefficients
coef <- fread("analyses/train/lasso_coefficients.txt")
coef <- coef[coef_type != "Dataset-specific covariate"] # not of interest here
coef <- coef[!(var %in% c("age", "sex"))] # baseline model to compare to

# Build models
mf1 <- "Surv(incident_followup, incident_cvd) ~ strata(sex) + age"
mf2 <- "Surv(incident_followup, incident_cvd) ~ strata(sex) + age + %s"

# Get C-indices
cind1 <- cbind(model="reference", cox.test(mf1, "incident_cvd", test)$model_fit)
cinds <- foreach(this_var = coef[,unique(var)], .combine=rbind, .init=cind1) %do% {
  # Scale variable if not a factor
  varterm <- if (coef[var == this_var, var != coefficient][1]) { this_var } else { sprintf("scale(%s)", this_var) }

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
cinds[, model_type := factor(model_type, levels=c("Reference", "Conventional RF", "Polygenic Risk Score", "NMR Metabolomics", "Clinical Biochemistry"))]
cinds <- cinds[order(C.index)][order(model_type)]
cinds[, model_name := factor(model_name, levels=model_name)]

# Write out
fwrite(cinds, sep="\t", quote=FALSE, file="analyses/test/feature_cindices/cindex.txt")

# Generate plots of C-indices per model
models <- unique(coef[,.(name, lambda, PGS)])
foreach(mIdx = models[,.I]) %do% {
  this_model <- models[mIdx]
  this_coef <- coef[this_model, on = .(name, lambda, PGS)]
  this_cind <- cinds[model == "reference" | model %in% this_coef$var]

  g <- ggplot(this_cind) +
    aes(x=model_name, y = C.index, ymin = L95, ymax = U95) +
    geom_errorbar(width=0) +
    geom_point(size=2, shape=18) +
    geom_hline(yintercept = this_cind[model == "reference", C.index], linetype = 2) +
    facet_grid(. ~ model_type, scales="free_x", space="free_x") +
    ylab("C-index (95% CI)") +
    xlab("") +
    theme_bw() +
    theme(legend.position="bottom", legend.box="vertical", axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

  ggsave(g, width=13, height=5, units="in", file=sprintf("analyses/test/feature_cindices/%s%s%s.pdf",
    tolower(gsub(" \\+? ?", "_", this_model[, name])),
    ifelse(this_model[,PGS], "_PGS", ""),
    ifelse(this_model[,lambda] == "", "", paste0("_", this_model[,lambda]))
  ))  
}
 
