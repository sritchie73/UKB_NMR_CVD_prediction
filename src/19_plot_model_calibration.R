library(data.table)
library(survival)
library(foreach)
library(ggplot2)
library(ggh4x)
library(scales)
library(ggthemes)

# Make output directory
system("mkdir -p analyses/test/")

###
# Calibration by risk decile
###

# Stratify by risk decile
pred_risk <- fread("analyses/CVD_weight_training/CVD_linear_predictors_and_risk.txt")
pred_risk <- pred_risk[score_type == "non-derived"]
pred_risk <- pred_risk[order(uk_calibrated_risk)][order(model)][order(sex)]
pred_risk[, risk_decile := floor((seq_len(.N)-1)/(.N/10))+1, by=.(model, sex)]

# Compute mean risk and 95% CI in each decile
risk_comp <- pred_risk[,.(mean_risk=mean(uk_calibrated_risk), L95=quantile(uk_calibrated_risk, 0.025), U95=quantile(uk_calibrated_risk, 0.975)), by=.(sex, model, risk_decile)]

# Get Kaplan Meier estimates of observed risk in UK Biobank
dat <- fread("data/cleaned/analysis_cohort.txt", select=c("eid", "sex", "age", "incident_cvd_followup", "incident_cvd"))
dat <- foreach(this_model = unique(pred_risk$model), .combine=rbind) %do% {
  cbind(model=this_model, dat)
}
dat[pred_risk, on = .(eid, sex, model), risk_decile := i.risk_decile]

obs_risk <- foreach(this_sex = c("Male", "Female"), .combine=rbind) %:% 
  foreach(this_model = unique(pred_risk$model), .combine=rbind) %:% 
    foreach(this_risk_decile = 1:10, .combine=rbind) %do% {
      this_dat <- dat[sex == this_sex & model == this_model & risk_decile == this_risk_decile]
      km <- survfit(Surv(incident_cvd_followup, incident_cvd) ~ 1, data=this_dat)
      kms <- summary(km, time=10)
      data.table(sex=this_sex, model=this_model, risk_decile=this_risk_decile, mean_risk=1-kms$surv, L95=1-kms$upper, U95=1-kms$lower)
}

risk_comp <- rbind(idcol="risk_model", observed=obs_risk, predicted=risk_comp)

# Write out risk table
fwrite(risk_comp, sep="\t", quote=FALSE, file="analyses/test/predicted_vs_observed_by_decile.txt")

# Plot model calibration by risk decile
risk_comp <- dcast(risk_comp, model + sex + risk_decile ~ risk_model, value.var=c("mean_risk", "L95", "U95"))
risk_comp[, risk_decile := sprintf("%s-%s%%", (risk_decile - 1)*10, risk_decile*10)]
risk_comp[, sex := factor(paste0(sex, "s"), levels=c("Males", "Females"))]
risk_comp[, model := factor(model, levels=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]

narrow_ci_y <- risk_comp[
  risk_decile %in% sprintf("%s-%s%%", (1:7 - 1)*10, (1:7)*10)
]

narrow_ci_x <- risk_comp[
  risk_decile %in% sprintf("%s-%s%%", (1:4 - 1)*10, (1:4)*10)
]

g <- ggplot(risk_comp) +
  aes(x=mean_risk_observed, xmin=L95_observed, xmax=U95_observed, y=mean_risk_predicted, ymin=L95_predicted, ymax=U95_predicted, color=risk_decile) +
  facet_grid2(sex ~ model, scales="free", independent="all") +
  geom_abline(intercept=0, slope=1, linetype=2) +
  geom_errorbarh(height=0) +
  geom_errorbar(width=0) +
  geom_point(shape=23, fill="white", size=1.2) +
  geom_errorbarh(data=narrow_ci_x, height=0, show.legend=FALSE) +
  geom_errorbar(data=narrow_ci_y, width=0, show.legend=FALSE) +
  scale_x_continuous("Observed risk (95% CI) in UK Biobank", labels=percent, limits=c(0, NA)) +
  scale_y_continuous("Predicted risk (95% CI)", labels=percent, limits=c(0, NA)) +
  scale_color_tableau(name="Decile of predicted risk") +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8),
    strip.background=element_blank(), strip.text.x=element_text(size=6, face="bold"), strip.text.y=element_text(size=8, face="bold"),
    legend.title=element_text(size=8), legend.text=element_text(size=6), legend.box.margin=margin(-0.5, 0, 0, -1, unit="cm"),
    legend.position="bottom", legend.justification="left", legend.box.background=element_blank(), legend.background=element_blank()
  ) +
  guides(color=guide_legend(title.position="top", nrow=1))
ggsave(g, width=7.2, height=4, file="analyses/test/model_calibration.pdf")

###
# Calibration by age-group
###

# Compute mean risk and 95% CI in each age-group
pred_risk <- fread("analyses/CVD_weight_training/CVD_linear_predictors_and_risk.txt")
pred_risk <- pred_risk[score_type == "non-derived"]
risk_comp <- pred_risk[,.(mean_risk=mean(uk_calibrated_risk), L95=quantile(uk_calibrated_risk, 0.025), U95=quantile(uk_calibrated_risk, 0.975)), by=.(sex, age_group, model)]

# Get Kaplan Meier estimates of observed risk in UK Biobank
dat <- fread("data/cleaned/analysis_cohort.txt", select=c("eid", "sex", "age", "incident_cvd_followup", "incident_cvd"))
dat <- foreach(this_model = unique(pred_risk$model), .combine=rbind) %do% {
  cbind(model=this_model, dat)
}
dat[pred_risk, on = .(eid, sex, model), age_group := i.age_group]

obs_risk <- foreach(this_sex = c("Male", "Female"), .combine=rbind) %:% 
  foreach(this_age_group = sort(unique(pred_risk$age_group)), .combine=rbind) %:% 
    foreach(this_model = unique(pred_risk$model), .combine=rbind) %do% {
      this_dat <- dat[sex == this_sex & model == this_model & age_group == this_age_group]
      km <- survfit(Surv(incident_cvd_followup, incident_cvd) ~ 1, data=this_dat)
      kms <- summary(km, time=10)
      data.table(sex=this_sex, age_group=this_age_group, model=this_model, mean_risk=1-kms$surv, L95=1-kms$upper, U95=1-kms$lower)
}

risk_comp <- rbind(idcol="risk_model", observed=obs_risk, predicted=risk_comp)

# Write out risk table
fwrite(risk_comp, sep="\t", quote=FALSE, file="analyses/test/predicted_vs_observed_by_age_group.txt")

# Plot model calibration by risk decile
risk_comp <- dcast(risk_comp, model + sex + age_group ~ risk_model, value.var=c("mean_risk", "L95", "U95"))
risk_comp[, sex := factor(paste0(sex, "s"), levels=c("Males", "Females"))]
risk_comp[, model := factor(model, levels=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]

g <- ggplot(risk_comp) +
  aes(x=mean_risk_observed, xmin=L95_observed, xmax=U95_observed, y=mean_risk_predicted, ymin=L95_predicted, ymax=U95_predicted, color=age_group) +
  facet_grid2(sex ~ model, scales="free", independent="all") +
  geom_abline(intercept=0, slope=1, linetype=2) +
  geom_errorbarh(height=0) +
  geom_errorbar(width=0) +
  geom_point(shape=23, fill="white", size=1.2) +
  scale_x_continuous("Observed risk (95% CI) in UK Biobank", labels=percent, limits=c(0, NA)) +
  scale_y_continuous("Predicted risk (95% CI)", labels=percent, limits=c(0, NA)) +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8),
    strip.background=element_blank(), strip.text.x=element_text(size=6, face="bold"), strip.text.y=element_text(size=8, face="bold"),
    legend.title=element_text(size=8), legend.text=element_text(size=6), legend.box.margin=margin(-0.5, 0, 0, -1, unit="cm"),
    legend.position="bottom", legend.justification="left", legend.box.background=element_blank(), legend.background=element_blank()
  ) +
  guides(color=guide_legend(title="Age (years)", title.position="top", nrow=1))
ggsave(g, width=7.2, height=4, file="analyses/test/model_calibration_by_age.pdf")


