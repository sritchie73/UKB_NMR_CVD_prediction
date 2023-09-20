library(data.table)
library(foreach)
library(doMC)
library(boot)
library(survival)
library(ggplot2)
library(ggh4x)
library(forcats)
registerDoMC(30) # Takes ~2 hours?

# Create output directory
system("mkdir -p analyses/univariate")

# Load biomarker information sheets
bio_info <- fread("data/ukb/biomarkers/output/biomarker_info.txt")

# Load dataset
dat <- fread("data/cleaned/full_UKB_analysis_cohort.txt")

# Get list of biomarkers to test
test_assay <- bio_info[!is.na(UKB.Field.ID) & sample_type != "Urine" & var != "tchol" & var != "hdl", var]

# Build set of models to test
model_info <- foreach(this_sex=c("Males", "Females", "Sex-stratified"), .combine=rbind) %:%
  foreach(this_model = c("SCORE2", "SCORE2_excl_UKB", paste("SCORE2 +", test_assay)), .combine=rbind) %do% {
    data.table(sex=this_sex, model=this_model)
}
model_info[model != "SCORE2" & model != "SCORE2_excl_UKB", biomarker := gsub("SCORE2 \\+ ", "", model)]
model_info[, model_var := paste0("model_", 1:.N), by=sex]

# Convert sex to integer for strata
strata_num <- dat[,.GRP, by=sex]
dat[strata_num, on = .(sex), sex_int := i.GRP]

# Bootstrap package rejecs data with missing values, so we have to do some wrangling here.
# The function we provide to the bootstrap procedure automatically handles missing data, so
# here we set these to an indicator value to get past the censboot checks, before filtering 
# in the function passed to the bootstrap procedure
dat <- melt(dat, id.vars=c("eid", "sex", "sex_int", "incident_cvd_followup", "incident_cvd", "SCORE2", "SCORE2_excl_UKB"), measure.vars=test_assay)
dat[is.na(value), value := -1]
dat <- dcast(dat, eid + sex + sex_int + incident_cvd_followup + incident_cvd + SCORE2 + SCORE2_excl_UKB ~ variable, value.var="value")

# Create bootstrap function that: 
# (1) calculates the absolute C-index for each model
# (2) calculates the C-index standard error for each model
# (3) calculates the difference in C-index from SCORE2 for each model
# (4) calculates the % change in C-index from SCORE2 for each model
# For each entry in the model_info table
boot.fun <- function(dt) {
  model_cinds <- foreach(midx = model_info[,.I], .combine=rbind) %dopar% {
    this_minfo <- model_info[midx]
    if (this_minfo$model != "SCORE2" & this_minfo$model != "SCORE2_excl_UKB") {
      # Drop missing data for this biomarker (masked with -1 to pass censboot checks)
      this_dt <- dt[dt[[this_minfo$biomarker]] >= 0]

      # Get C-index and SE from cox model
      if (this_minfo$sex == "Sex-stratified") {
        mf <- sprintf("Surv(incident_cvd_followup, incident_cvd) ~ strata(sex) + offset(SCORE2_excl_UKB) + %s", this_minfo$biomarker)
      } else {
        this_dt <- dt[sex == gsub("s$", "", this_minfo$sex)]
        mf <- sprintf("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2_excl_UKB) + %s", this_minfo$biomarker)
      }
      cx <- coxph(as.formula(mf), data=this_dt)
      C.index <- summary(cx)$concordance[1]
      SE <- summary(cx)$concordance[2]

      # Compute C-index for SCORE2 in samples with non-missing data directly from linear predictor
      this_y <- Surv(this_dt[["incident_cvd_followup"]], this_dt[["incident_cvd"]])
      this_x <- this_dt[["SCORE2_excl_UKB"]]
      if (this_minfo$sex == "Sex-stratified") {
        cf <- survival::concordancefit(y = this_y, x = this_x, strata = this_dt[["sex_int"]], reverse = TRUE, timefix = TRUE)
      } else {
        cf <- survival::concordancefit(y = this_y, x = this_x, reverse = TRUE, timefix = TRUE)
      }   
      delta.C <- C.index - cf$concordance
      pct_change <- delta.C/(C.index - 0.5)*100

      return(cbind(this_minfo, C.index, SE, delta.C, pct_change))
    } else {
      # Compute C-index from linear predictor directly for SCORE2
      if (this_minfo$sex == "Sex-stratified") {
        this_y <- Surv(dt[["incident_cvd_followup"]], dt[["incident_cvd"]])
        this_x <- dt[[this_minfo$model]]
        cf <- survival::concordancefit(y = this_y, x = this_x, strata = dt[["sex_int"]], reverse = TRUE, timefix = TRUE)
      } else {
        this_dt <- dt[sex == gsub("s$", "", this_minfo$sex)]
        this_y <- Surv(this_dt[["incident_cvd_followup"]], this_dt[["incident_cvd"]])
        this_x <- this_dt[[this_minfo$model]]
        cf <- survival::concordancefit(y = this_y, x = this_x, reverse = TRUE, timefix = TRUE)
      }    
      return(cbind(this_minfo, C.index=cf$concordance, SE=sqrt(cf$var), delta.C=0, pct_change=0)) # jackknife SE computed same way as in coxph()
    } 
  } 
  # Need to return as flat vector
  res <- melt(model_cinds, measure.vars=c("C.index", "SE", "delta.C", "pct_change"))
  res$value
}

# Run bootstrap analysis
surv_cols_idx <- match(c("incident_cvd_followup", "incident_cvd"), names(dat))
boot_res <- censboot(dat, boot.fun, 1000, index=surv_cols_idx)
saveRDS(boot_res, file="analyses/univariate/full_UKB_cindices_bootstraps.rds")

# Now we need to extract and collate the results
cinds <- foreach(this_metric = c("C.index", "SE", "delta.C", "pct_change"), .combine=rbind) %do% {
  cbind(model_info, metric = this_metric)
}
cinds[, estimate := boot_res$t0]
cinds[, L95 := apply(boot_res$t, 2, function(v) {  sort(v)[25] })]
cinds[, U95 := apply(boot_res$t, 2, function(v) {  sort(v)[975] })]
cinds[, pval := apply(boot_res$t, 2, function(v) {
  nulls <- pmin(sum(v <= 0), sum(v >= 0))
  n <- nrow(boot_res$t)
  pmin(2*(nulls + 1)/(n + 1), 1)
})]

# Cast to wide
cinds <- dcast(cinds, sex + model + biomarker ~ metric, value.var=c("estimate", "L95", "U95", "pval"))
cinds <- cinds[,.(sex, model, biomarker,
  C.index = estimate_C.index, C.L95 = L95_C.index, C.U95 = U95_C.index,
  SE = estimate_SE, SE.L95 = L95_SE, SE.U95 = U95_SE,
  deltaC = estimate_delta.C, deltaC.L95 = L95_delta.C, deltaC.U95 = U95_delta.C, deltaC.pval = pval_delta.C,
  pct_change = estimate_pct_change, pct.L95 = L95_pct_change, pct.U95 = U95_pct_change, pct.pval = pval_pct_change
)]

# Extract comparison of SCORE2 with and without UKB weights in full dataset
score2_cinds <- cinds[model %in% c("SCORE2", "SCORE2_excl_UKB")]
score2_cinds <- score2_cinds[, .(sex, SCORE2_method = ifelse(model == "SCORE2", "Weights derived from all datasets", "Weights derived excluding UK Biobank"),
  C.index, C.L95, C.U95, SE, SE.L95, SE.U95)]
fwrite(score2_cinds, sep="\t", quote=FALSE, file="analyses/test/full_UKB_cindex_by_SCORE2_method.txt")

score2_cinds[, sex := factor(sex, levels=c("Sex-stratified", "Males", "Females"))]
score2_cinds[, SCORE2_method := factor(SCORE2_method, levels=c("Weights derived excluding UK Biobank", "Weights derived from all datasets"))]
g <- ggplot(score2_cinds) +
  aes(x=C.index, xmin=C.L95, xmax=C.U95, y=SCORE2_method, color=SCORE2_method) +
  facet_wrap(~ sex, scales="free_x") +
  geom_errorbarh(height=0) +
  geom_point(shape=18) +
  scale_color_manual(values=c("Weights derived excluding UK Biobank"="#2166ac", "Weights derived from all datasets"="#b2182b")) +
  xlab("C-index (95% CI)") +
  guides(color=guide_legend(title="SCORE2 model", reverse=TRUE)) +
  theme_bw() +
  theme(
    axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    strip.background=element_blank(), strip.text=element_text(size=8, face="bold"),
    legend.text=element_text(size=6), legend.title=element_text(size=8)
  )
ggsave(g, width=7.2, height=1.5, file="analyses/test/full_UKB_cindex_by_SCORE2_method.pdf")

# Now extract the per-biomarker C-indices
cinds <- cinds[model != "SCORE2"]
cinds[model == "SCORE2_excl_UKB", model := "SCORE2"]

# Add human friendly display name
cinds[model != "SCORE2", model := paste("SCORE2 +", fcase(
  biomarker == "alb", "Albumin",
  biomarker == "alt", "ALT",
  biomarker == "alp", "ALP",
  biomarker == "apoa1", "ApoA1",
  biomarker == "apob", "ApoB",
  biomarker == "asp", "AST",
  biomarker == "dbili", "Bilirubin (direct)",
  biomarker == "tbili", "Bilirubin (total)",
  biomarker == "calcium", "Calcium",
  biomarker == "creat", "Creatinine",
  biomarker == "crp", "CRP",
  biomarker == "cyst", "Cystatin-C",
  biomarker == "ggt", "GGT",
  biomarker == "glucose", "Glucose",
  biomarker == "hba1c", "HbA1c",
  biomarker == "igf1", "IGF-1",
  biomarker == "lpa", "Lp(a)",
  biomarker == "ldl", "LDL cholesterol",
  biomarker == "oest", "Oestradiol",
  biomarker == "phos", "Phosphate",
  biomarker == "rheuf", "RF",
  biomarker == "shbg", "SHBG",
  biomarker == "testos", "Testosterone",
  biomarker == "protein", "Total protein",
  biomarker == "trig", "Triglycerides",
  biomarker == "uric", "Urate",
  biomarker == "urea", "Urea",
  biomarker == "vitd25", "Vitamin D"
))]

# Compute FDR- adjusted P-value
cinds[model != "SCORE2", deltaC.fdr := p.adjust(deltaC.pval, method="fdr"), by=sex]

# Write out
fwrite(cinds, sep="\t", quote=FALSE, file="analyses/univariate/full_UKB_cindex.txt")

# Plot delta C-index for top-10 biomarkers in Sex-stratified results
ggdt <- rbind(
  cinds[model != "SCORE2" & sex == "Sex-stratified"][order(-deltaC)][1:10]
)
ggdt[, model := factor(model, levels=rev(unique(model)))]

g <- ggplot(ggdt) +
  aes(x=deltaC, xmin=deltaC.L95, xmax=deltaC.U95, y=model) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, size=2, fill="white") +
  scale_x_continuous("delta C-index (95% CI)", limits=c(0, 0.01)) +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )
ggsave(g, width=3.5, height=2, file="analyses/univariate/full_UKB_cindex_sex_stratified.pdf")

# Create formatted table for manuscript
dt <- cinds[sex == "Sex-stratified"]
dt[model == "SCORE2", c("deltaC", "deltaC.L95", "deltaC.U95", "deltaC.pval", "pct_change", "pct.L95", "pct.U95") := NA]
dt[, pct.pval := NULL]
dt[, pct_change := pct_change / 100]
dt[, pct.L95 := pct.L95 / 100]
dt[, pct.U95 := pct.U95 / 100]
dt <- dt[order(-deltaC, na.last=FALSE)][order(deltaC.pval, na.last=FALSE)]

# Need to also tabulate complete data
dat <- fread("data/cleaned/full_UKB_analysis_cohort.txt")
dat <- melt(dat, id.vars=c("eid", "sex", "incident_cvd"), measure.vars=c("SCORE2", test_assay), variable.name="biomarker")
dat <- dat[,.(samples=sum(!is.na(value)), cases=sum(!is.na(value) & incident_cvd)), by=biomarker]
dat <- dat[biomarker == "SCORE2", biomarker := NA]
dt <- dat[dt, on = .(biomarker)]

# Reorganize columns and write out
dt <- dt[,.(model, samples, cases, C.index, C.L95, C.U95, SE, SE.L95, SE.U95, deltaC, deltaC.L95, deltaC.U95,
  deltaC.pval, deltaC.fdr, pct_change, pct.L95, pct.U95)]
fwrite(dt, sep="\t", quote=FALSE, file="analyses/univariate/full_UKB_cindex_for_supp.txt")




