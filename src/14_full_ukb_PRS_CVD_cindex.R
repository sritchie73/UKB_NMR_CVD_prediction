library(data.table)
library(caret)
library(survival)
library(foreach)
source('src/utils/mean_pvalue.R')
source('src/utils/format_pval.R')
source('src/utils/score_cindex.R')

# First load in and subset the full UKB cohort to those where PRS can be modelled
dat <- fread('data/cleaned/full_UKB_analysis_cohort.txt')
dat <- dat[!(no_genetics) & !(prs_training_samples)]

# Partition into 5 folds for estimating C-index
dat[, cvd_prediction_foldid := createFolds(paste(incident_cvd, sex), k=5, list=FALSE)]

# Write out test folds
fwrite(dat[,.(eid, cvd_prediction_foldid)], sep="\t", quote=FALSE, file="analyses/CVD_weight_training/full_UKB_PRS_cross_validation_fold_allocation.txt")

# Standardise PRSs
dat[, CAD_metaGRS := scale(CAD_metaGRS)]
dat[, Stroke_metaGRS := scale(Stroke_metaGRS)]

# Estimate relative contribution of PRSs weights using Cox proportional hazards models
cvd_weights <- foreach(this_sex = c("Male", "Female"), .combine=rbind) %:%
  foreach(this_test_fold = 1:5, .combine=rbind) %do% {
    # extract training data
    this_dat <- dat[cvd_prediction_foldid != this_test_fold & sex == this_sex]

    # Build model formula
    mf <- "Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2_excl_UKB) + CAD_metaGRS + Stroke_metaGRS"

    # Fit cox model
    cx <- coxph(as.formula(mf), data=this_dat)

    # Extract relevant coefficients
    ci <- confint(cx)
    cf <- as.data.table(coef(summary(cx)), keep.rownames="score")
    cf[, score := gsub("_", " ", score)]
    cf <- cf[, .(score, weight=coef, L95=ci[,1], U95=ci[,2], pval=`Pr(>|z|)`)]

    # add information
    cbind(test_fold = this_test_fold, sex = this_sex, model =  "SCORE2 + PRSs", cf)
}

# Compute average weights
avg_cvd_weights <- cvd_weights[, .(
  weight=round(mean(weight), digits=3),
  L95=round(mean(L95), digits=3),
  U95=round(mean(U95), digits=3),
  pval=format_pval(mean_pvalue(pval, weight))
), by=.(sex, model, score)]

# Write out
fwrite(cvd_weights, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/full_UKB_PRS_cvd_score_weights_in_5fold_cross_validation.txt")
fwrite(avg_cvd_weights, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/full_UKB_PRS_average_cvd_score_weights.txt")

# Extract linear predictors in each test fold for downstream analyses
pred_scores <- foreach(this_model = c("SCORE2", "SCORE2 + PRSs"), .combine=rbind) %:%
  foreach(this_sex = c("Male", "Female"), .combine=rbind) %:%
    foreach(this_test_fold = 1:5, .combine=rbind) %do% {
      # extract test data
      this_dat <- dat[cvd_prediction_foldid == this_test_fold & sex == this_sex]

      # Melt columns of interest to long format for summing
      lp_cols <- "SCORE2_excl_UKB"
      if (this_model == "SCORE2 + PRSs") {
        lp_cols <- c(lp_cols, "CAD_metaGRS", "Stroke_metaGRS")
      }
      this_dat <- melt(this_dat, id.vars=c("eid", "sex", "age", "incident_cvd", "incident_cvd_followup", "cvd_prediction_foldid"), measure.vars=lp_cols, variable.name="score")

      # Multiply scores by weights
      this_weights <- cvd_weights[test_fold == this_test_fold & sex == this_sex & model == this_model]
      this_weights[, score := gsub(" ", "_", score)]
      this_dat[this_weights, on = .(score), value := value * i.weight]

      # Sum to obtain linear predictor
      this_dat <- this_dat[, .(linear_predictor=sum(value)), by=.(eid, sex, age, incident_cvd, incident_cvd_followup, cvd_prediction_foldid)]

      # add model info
      cbind(model = this_model, this_dat)
}

# Write out
fwrite(pred_scores, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/full_UKB_PRS_CVD_linear_predictors_and_risk.txt")

# Compute C-indices
cinds <- foreach(this_sex=c("Males", "Females", "Sex-stratified"), .combine=rbind) %:%
  foreach(this_model=c("SCORE2", "SCORE2 + PRSs"), .combine=rbind) %do% {
    this_dat <- pred_scores[model == this_model]
    this_info <- data.table(sex=this_sex, model=this_model)
    if (this_sex == "Sex-stratified") {
      this_cind <- score_cindex("Surv(incident_cvd_followup, incident_cvd) ~ strata(sex) + linear_predictor", this_dat)
    } else {
      this_cind <- score_cindex("Surv(incident_cvd_followup, incident_cvd) ~ linear_predictor", this_dat[sex == gsub("s$", "", this_sex)])
    }
    cbind(this_info, this_cind)
}

# Compute delta C-index from SCORE2
ref <- cinds[model == "SCORE2"]
cinds[ref, on = .(sex), c("deltaC", "deltaC.L95", "deltaC.U95") := .(C.index - i.C.index, L95 - i.C.index, U95 - i.C.index)]
cinds[model == "SCORE2", c("deltaC", "deltaC.L95", "deltaC.U95") := NA]

# Compute % improvement over SCORE2
cinds[ref, on = .(sex), c("pct_change", "pct.L95", "pct.U95") := .(deltaC/(i.C.index - 0.5)*100, deltaC.L95/(i.C.index - 0.5)*100, deltaC.U95/(i.C.index - 0.5)*100)]
cinds[model == "SCORE2", c("pct_change", "pct.L95", "pct.U95") := NA]

# Write out
fwrite(cinds, sep="\t", quote=FALSE, file="analyses/test/full_UKB_PRS_cindices.txt")

# Encode factors for plotting
cinds[, model := factor(model, levels=rev(c("SCORE2", "SCORE2 + PRSs")))]
cinds[, sex := factor(sex, levels=c("Males", "Females", "Sex-stratified"))]
ref[, model := NULL]
ref[, sex := factor(sex, levels=c("Males", "Females", "Sex-stratified"))]

# Plot full biomarker scores
g <- ggplot(cinds) +
  aes(x=C.index, xmin=L95, xmax=U95, y=model, color=sex) +
  facet_wrap(~ sex, nrow=1, scales="free_x") +
  geom_vline(data=ref, aes(xintercept=C.index), linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, size=2, fill="white") +
  scale_color_manual("Sex", values=c("Males"="#e41a1c", "Females"="#377eb8", "Sex-stratified"="#006d2c")) +
  ylab("") + xlab("C-index (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )
ggsave(g, width=7.2, height=1.5, file="analyses/test/full_UKB_cindices_non_derived.pdf")

# Plot change in C-index
g <- ggplot(cinds[model != "SCORE2"]) +
  aes(x=deltaC, xmin=deltaC.L95, xmax=deltaC.U95, y=model, color=sex) +
  facet_wrap(~ sex, nrow=1) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, size=2, fill="white") +
  scale_color_manual("Sex", values=c("Males"="#e41a1c", "Females"="#377eb8", "Sex-stratified"="#006d2c")) +
  ylab("") + xlab("ΔC-index over SCORE2 (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )
ggsave(g, width=7.2, height=1.3, file="analyses/test/full_UKB_delta_cindices_non_derived.pdf", device=cairo_pdf)

# Plot % improvement
g <- ggplot(cinds[model != "SCORE2"]) +
  aes(x=pct_change, xmin=pct.L95, xmax=pct.U95, y=model, color=sex) +
  facet_wrap(~ sex, nrow=1) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, size=2, fill="white") +
  scale_color_manual("Sex", values=c("Males"="#e41a1c", "Females"="#377eb8", "Sex-stratified"="#006d2c")) +
  ylab("") + xlab("% improvement in C-index over SCORE2 (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )
ggsave(g, width=7.2, height=1.3, file="analyses/test/full_UKB_pct_improvement_cindices_non_derived.pdf")

# Create formatted table for manuscript
dt <- copy(cinds)
dt[, SE := NULL]
dt[, pct_change := pct_change / 100]
dt[, pct.L95 := pct.L95 / 100]
dt[, pct.U95 := pct.U95 / 100]
dt[sex == "Males", c("samples", "cases") := .(dat[sex == "Male", .N], dat[sex == "Male", sum(incident_cvd)])]
dt[sex == "Females", c("samples", "cases") := .(dat[sex == "Female", .N], dat[sex == "Female", sum(incident_cvd)])]
dt[sex == "Sex-stratified", c("samples", "cases") := .(dat[,.N], dat[, sum(incident_cvd)])]
fwrite(dt, sep="\t", quote=FALSE, file="analyses/test/full_UKB_cindices_for_supp.txt")









