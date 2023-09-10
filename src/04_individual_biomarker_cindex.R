# Requires at least 1:10:0 to run on compute nodes, recommend running with 20 cores
library(data.table)
library(foreach)
library(survival)
library(ggplot2)
library(ggh4x)
library(forcats)
source("src/utils/score_cindex.R")
registerDoMC(10)

# Create output directory
system("mkdir -p analyses/univariate")

# Load biomarker information sheets
bio_info <- fread("data/ukb/biomarkers/output/biomarker_info.txt")
nmr_info <- fread("data/ukb/NMR_metabolomics/biomarker_information.txt")

# Load analysis cohort
dat <- fread("data/cleaned/analysis_cohort.txt")

# Code risk factors as per SCORE2
dat[, age := (age - 60)/5]
dat[is.na(smoking), smoking := FALSE]
dat[, smoking := factor(smoking, levels=c(FALSE, TRUE))]
dat[, sbp := (sbp - 120)/20]
dat[, tchol := (tchol - 6)/1]
dat[, hdl := (hdl - 1.3)/0.5]

# Code factors, using non-risk/lower-risk group as reference
dat[, sex := factor(sex, levels=c("Female", "Male"))]

# Get list of biomarkers to test
test_nmr <- nmr_info[(Nightingale), Biomarker]
test_assay <- bio_info[!is.na(UKB.Field.ID) & sample_type != "Urine" & var != "tchol" & var != "hdl", var]

# Build set of models to test
models <- rbind(
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ age", type="risk_factors", name="age"),
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ smoking", type="risk_factors", name="smoking"),
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ sbp", type="risk_factors", name="sbp"),
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ tchol", type="risk_factors", name="tchol"),
  data.table(formula="Surv(incident_cvd_followup, incident_cvd) ~ hdl", type="risk_factors", name="hdl"),
  data.table(formula=sprintf("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2_excl_UKB) + scale(%s)", test_nmr), type="NMR", name=test_nmr),
  data.table(formula=sprintf("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2_excl_UKB) + scale(%s)", test_assay), type="assays", name=test_assay)
)

# Compute C-indices
cinds <- foreach(this_sex=c("Males", "Females", "Sex-stratified"), .combine=rbind) %do% { 

  # Extract subset of data needed
  if (this_sex == "Males") {
    this_dat <- dat[sex == "Male"]
  } else if (this_sex == "Females") {
    this_dat <- dat[sex == "Female"]
  } else {
    this_dat <- dat
  }

  # Fit individual models
  foreach(midx = models[,.I], .combine=rbind) %dopar% {
    # Extract model information
    this_model <- models[midx]

    # Set up formula
    mf <- this_model$formula
    if (this_sex == "Sex-stratified") {
      mf <- paste(mf, "+ strata(sex)")
    }

    # Fit Cox proportional hazards model
    cph <- coxph(as.formula(mf), data=this_dat)

    # Get C-index and its 95% CI 
    cindex <- summary(cph)$concordance[1]
    cindex.se <- summary(cph)$concordance[2]
    cindex.l95 <- cindex - qnorm(0.975)*cindex.se
    cindex.u95 <- cindex + qnorm(0.975)*cindex.se
     
    # Return results
    cbind(this_model, sex=this_sex, samples=cph$n, cases=cph$nevent, 
      C.index=cindex, L95=cindex.l95, U95=cindex.u95)
  }
}

# Load in and add SCORE2 C-index
score2_cindex <- fread("analyses/test/cindex_by_SCORE2_method.txt")
score2_cindex <- score2_cindex[SCORE2_method == "Weights derived excluding UK Biobank"]
score2_cindex[, formula := "Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2_excl_UKB)"]
score2_cindex[, type := "SCORE2"]
score2_cindex[, name := "SCORE2"]
score2_cindex[sex == "Males", c("samples", "cases") := dat[sex == "Male", .(.N, sum(incident_cvd))]]
score2_cindex[sex == "Females", c("samples", "cases") := dat[sex == "Female", .(.N, sum(incident_cvd))]]
score2_cindex[sex == "Females", c("samples", "cases") := dat[sex == "Female", .(.N, sum(incident_cvd))]]
score2_cindex[sex == "Sex-stratified", c("samples", "cases") := dat[, .(.N, sum(incident_cvd))]]
score2_cindex[, c("SCORE2_method", "SE") := NULL]
cinds <- rbind(cinds, score2_cindex)

# Add in missing strata term to documented formula
cinds[sex == "Sex-stratified" & name != "sex", formula := paste(formula, "+ strata(sex)")]

# Add in delta change compared to SCORE2
score2_ref <- cinds[name == "SCORE2"]
score2_ref[, type := NULL]
score2_ref <- rbind(idcol="type", "NMR"=score2_ref, "assays"=score2_ref)
cinds[score2_ref, on = .(type, sex), c("deltaC", "deltaC.L95", "deltaC.U95") := .(C.index - i.C.index, L95 - i.C.index, U95 - i.C.index)]

# Add human friendly display name
cinds[, display_name := name]
cinds[name == "age", display_name := "Age"]
cinds[name == "sex", display_name := "Sex"]
cinds[name == "smoking", display_name := "Smoker"]
cinds[name == "sbp", display_name := "SBP"]
cinds[nmr_info[Group %like% "Amino" & Type == "Non-derived"], on = .(name=Biomarker), display_name := Description]
cinds[, display_name := gsub("_", "-", display_name)]
cinds[, display_name := gsub("-pct", " %", display_name)]
cinds[, display_name := gsub("-by-", " / ", display_name)]
cinds[name == "Clinical_LDL_C", display_name := "Clinical LDL-C"]
cinds[name == "hdl", display_name := "HDL cholesterol"]
cinds[name == "tchol", display_name := "Total cholesterol"]
cinds[name %in% test_assay, display_name := fcase(
  name == "alb", "Albumin",
  name == "alt", "ALT",
  name == "alp", "ALP",
  name == "apoa1", "ApoA1",
  name == "apob", "ApoB",
  name == "asp", "AST",
  name == "dbili", "Bilirubin (direct)",
  name == "tbili", "Bilirubin (total)",
  name == "calcium", "Calcium",
  name == "creat", "Creatinine",
  name == "crp", "CRP",
  name == "cyst", "Cystatin-C",
  name == "ggt", "GGT",
  name == "glucose", "Glucose",
  name == "hba1c", "HbA1c",
  name == "igf1", "IGF-1",
  name == "lpa", "Lp(a)",
  name == "ldl", "LDL cholesterol",
  name == "oest", "Oestradiol",
  name == "phos", "Phosphate",
  name == "rheuf", "RF",
  name == "shbg", "SHBG",
  name == "testos", "Testosterone",
  name == "protein", "Total protein",
  name == "trig", "Triglycerides",
  name == "uric", "Urate",
  name == "urea", "Urea",
  name == "vitd25", "Vitamin D"
)]

# Write out
fwrite(cinds, sep="\t", quote=FALSE, file="analyses/univariate/cindices.txt")

# Plot top NMR biomarkers + SCORE2
ggdt <- rbind(
  cinds[type == "SCORE2"],
  cinds[type == "NMR"][order(-C.index)][,.SD[1:10], by=sex]
)
ggdt[type != "SCORE2", display_name := paste("SCORE2 +", display_name)]
ggdt[, sex := factor(sex, c("Males", "Females", "Sex-stratified"))]
ggdt[, type := factor(type, levels=c("SCORE2", "NMR"))]
ggdt <- ggdt[order(sex)]
ggdt[, rank := factor(.I)]

g <- ggplot(ggdt) + 
  aes(y=fct_rev(rank), x=C.index, xmin=L95, xmax=U95, color=type) +
  facet_grid2(. ~ sex, scales="free", independent="all") +
  geom_errorbarh(height=0) +
  geom_point(shape=23, fill="white", size=1.2) +
  scale_color_manual(values=c("SCORE2"="#2166ac", "NMR"="black")) +
  geom_vline(data=ggdt[type == "SCORE2"], aes(xintercept=C.index), linetype=2, color="#4393c3") +
  scale_y_discrete(labels=structure(ggdt$display_name, names=as.character(ggdt$rank))) +
  xlab("C-index (95% CI)") +
  theme_bw() +
  theme(
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    axis.text.y=element_text(size=6, color="black"), axis.title.y=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    strip.background=element_blank(), strip.text=element_text(size=8, face="bold"),
    legend.position="none"
  )
ggsave(g, width=7.2, height=2, file="analyses/univariate/top10_NMR.pdf")

# Format table for output
dt <- cinds[type %in% c("SCORE2", "NMR")]
dt <- dcast(dt, display_name ~ sex, value.var=c("samples", "cases", "C.index", "L95", "U95", "deltaC", "deltaC.L95", "deltaC.U95"))
dt <- dt[,.SD, .SDcols=c("display_name",
  sapply(c("_Males", "_Females", "_Sex-stratified"), function(f) {
    paste0(c("samples", "cases", "C.index", "L95", "U95", "deltaC", "deltaC.L95", "deltaC.U95"), f)
  })
)]

dt <- dt[order(-C.index_Males)]
dt <- rbind(dt[display_name == "SCORE2"], dt[display_name != "SCORE2"])
dt[display_name != "SCORE2", display_name := paste("SCORE2 +", display_name)]
fwrite(dt, sep="\t", quote=FALSE, file="analyses/univariate/NMR_cindex_supp.txt")

