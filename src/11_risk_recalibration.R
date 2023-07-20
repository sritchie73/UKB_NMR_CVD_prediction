library(data.table)
library(foreach)
library(ggplot2)
library(ggh4x)
library(scales)
library(ggthemes)

# Load predicted risk
pred_risk <- fread("analyses/CVD_score_weighting/CVD_linear_predictors_and_risk.txt")

# Add in 5-year age group
pred_risk[, age_group := sprintf("%s-%s", age %/% 5 * 5, age %/% 5 * 5 + 4)]

# Compute risk quantiles in each age group for each model and risk calculation
risk_comp <- melt(pred_risk, measure.vars=c("ukb_absrisk", "score2_absrisk", "uk_risk"), variable.name="risk_model", value.name="risk")
risk_comp <- risk_comp[,.(mean_risk=mean(risk), L95=quantile(risk, 0.025), U95=quantile(risk, 0.975)),
  by=.(sex, age_group, model, risk_model)]
risk_comp[, risk_model := fcase(
  risk_model == "ukb_absrisk", "Risk predicted from Cox models fit in UK Biobank",
  risk_model == "score2_absrisk", "Risk predicted using baseline hazards estimated across 44 cohorts (SCORE2)",
  risk_model == "uk_risk", "Risk recalibrated to European CVD low-risk region (SCORE2)"
)]

# Add in CPRD-dervied estimates of CVD incidence published in Sun et al. 2021
CPRD = data.table(
  age_at_risk_start = seq(40, 75, by=5),
  age_at_risk_end = seq(44, 79, by=5),
  male_rate_per_1000py = c(1.578, 2.964, 5.141, 7.561, 10.848, 13.731, 19.271, 25.667),
  female_rate_per_1000py = c(0.726, 1.309, 2.165, 2.969, 4.560, 7.007, 10.840, 17.007)
)

# Calculate annual CVD incidence in each age-at-risk group from the
# 1000 person-year CVD rates:
CPRD[, male_annual_incidence := male_rate_per_1000py / 1000]
CPRD[, female_annual_incidence := female_rate_per_1000py / 1000]

# The expected 10-year CVD risk is calculated for each age group
# based on the mid-point of the next interval ahead, e.g. for the
# 40-44 year age-group the expected 10-year risk is calculated based
# on the annual incidence rates in the 45-49 year olds in CPRD:
CPRD[, age_group_start := age_at_risk_start - 5]
CPRD[, age_group_end := age_at_risk_end - 5]

# The expected risk is calculated assuming exponential survival (i.e.
# constant hazard) in the 10-years ahead.
CPRD[, male_expected_risk := 1 - exp(-male_annual_incidence * 10)]
CPRD[, female_expected_risk := 1 - exp(-female_annual_incidence * 10)]

# Add to risk_comp table
CPRD[, age_group := sprintf("%s-%s", age_group_start, age_group_end)]
CPRD <- CPRD[age_group_start >= 40 & age_group_end < 70]
setnames(CPRD, c("male_expected_risk", "female_expected_risk"), c("Male", "Female"))
CPRD <- melt(CPRD, id.vars="age_group", measure.vars=c("Male", "Female"), variable.name="sex", value.name="mean_risk")
CPRD <- foreach(this_model = c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %do% {
  CPRD[, .(sex, age_group, model=this_model, risk_model="CVD incidence rates in CPRD (Sun et al. 2021)", mean_risk)]
}
risk_comp <- rbind(risk_comp, CPRD, fill=TRUE)

# Get Kaplan Meier estimates of observed risk in UK Biobank
dat <- fread("data/cleaned/analysis_cohort.txt", select=c("eid", "sex", "age", "incident_cvd_followup", "incident_cvd"))
dat[, age_group := sprintf("%s-%s", age %/% 5 * 5, age %/% 5 * 5 + 4)]
obs_risk <- foreach(this_sex = c("Male", "Female"), .combine=rbind) %:% 
  foreach(this_age_group = c("40-44", "45-49", "50-54", "55-59", "60-64", "65-69"), .combine=rbind) %do% {
    this_dat <- dat[sex == this_sex & age_group == this_age_group]
    km <- survfit(Surv(incident_cvd_followup, incident_cvd) ~ 1, data=this_dat)
    kms <- summary(km, time=10)
    data.table(sex=this_sex, age_group=this_age_group, risk_model="Kaplan-Meier estimate of observed risk in UK Biobank", mean_risk=1-kms$surv, L95=1-kms$upper, U95=1-kms$lower) 
}
obs_risk <- foreach(this_model = c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %do% {
  cbind(model=this_model, obs_risk)
}
risk_comp <- rbind(risk_comp, obs_risk)

# Write out risk tables
fwrite(risk_comp, sep="\t", quote=FALSE, file="analyses/test/risk_table_by_sex_age_model.txt")

# Make plot sanity checking absolute risk calculations
risk_comp[, sex := factor(paste0(sex, "s"), levels=c("Males", "Females"))]
risk_comp[, age_group := factor(age_group, levels=c("40-44", "45-49", "50-54", "55-59", "60-64", "65-69"))]
risk_comp[, model := factor(model, levels=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
risk_comp[, risk_model := factor(risk_model, levels=c(
  "Kaplan-Meier estimate of observed risk in UK Biobank",
  "Risk predicted from Cox models fit in UK Biobank",
  "Risk predicted using baseline hazards estimated across 44 cohorts (SCORE2)",
  "Risk recalibrated to European CVD low-risk region (SCORE2)",
  "CVD incidence rates in CPRD (Sun et al. 2021)"
))]

g <- ggplot(risk_comp) +
  aes(x=age_group, y=mean_risk, ymin=L95, ymax=U95, color=risk_model) +
  facet_grid(sex ~ model) +
  geom_errorbar(width=0, position=position_dodge(width=0.75)) +
  geom_point(shape=23, fill="white", size=1.2, position=position_dodge(width=0.75)) +
  geom_errorbar(data=risk_comp[risk_model == "Kaplan-Meier estimate of observed risk in UK Biobank"], width=0, position=position_nudge(x=-0.3), show.legend=FALSE) +
  scale_color_manual(values=c(
    "Kaplan-Meier estimate of observed risk in UK Biobank"="#41ab5d",
    "Risk predicted from Cox models fit in UK Biobank"="#006837",
    "Risk predicted using baseline hazards estimated across 44 cohorts (SCORE2)"="#1d91c0",
    "Risk recalibrated to European CVD low-risk region (SCORE2)"="#253494",
    "CVD incidence rates in CPRD (Sun et al. 2021)"="#b30000"
  )) +
  scale_y_continuous(name="10-year CVD risk (95% CI)", labels=percent) +
  xlab("Age-group") + 
  theme_bw() +
  theme(
    axis.text.y=element_text(size=6), axis.title.y=element_text(size=8),
    axis.text.x=element_text(size=6, angle=45, hjust=1), axis.title.x=element_text(size=8),
    strip.text.x=element_text(size=7, face="bold"), strip.text.y=element_text(size=8, face="bold"),
    strip.background=element_blank(),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
    legend.position="bottom", legend.justification="left", legend.box.margin=margin(-0.3, 0, 0, 0, unit="cm"),
    legend.title=element_blank(), legend.text=element_text(size=6),
    legend.spacing.y=unit(-0.1, 'cm')
  ) +
  guides(color=guide_legend(order=1, nrow=3, byrow=TRUE))
ggsave(g, width=7.2, height=5, file="analyses/test/risk_computation_check.pdf")

# Make a plot without the old absrisk and incidence
risk_comp <- risk_comp[risk_model != "Risk predicted from Cox models fit in UK Biobank"]
risk_comp <- risk_comp[risk_model != "CVD incidence rates in CPRD (Sun et al. 2021)"]
g <- ggplot(risk_comp) +
  aes(x=age_group, y=mean_risk, ymin=L95, ymax=U95, color=risk_model) +
  facet_grid(sex ~ model) +
  geom_errorbar(width=0, position=position_dodge(width=0.75)) +
  geom_point(shape=23, fill="white", size=1.2, position=position_dodge(width=0.75)) +
  geom_errorbar(data=risk_comp[risk_model == "Kaplan-Meier estimate of observed risk in UK Biobank"], width=0, position=position_nudge(x=-0.25), show.legend=FALSE) +
  scale_color_manual(values=c(
    "Kaplan-Meier estimate of observed risk in UK Biobank"="#225ea8",
    "Risk predicted using baseline hazards estimated across 44 cohorts (SCORE2)"="#238443",
    "Risk recalibrated to European CVD low-risk region (SCORE2)"="#dd3497"
  )) +
  scale_y_continuous(name="10-year CVD risk (95% CI)", labels=percent) +
  xlab("Age-group") + 
  theme_bw() +
  theme(
    axis.text.y=element_text(size=6), axis.title.y=element_text(size=8),
    axis.text.x=element_text(size=6, angle=45, hjust=1), axis.title.x=element_text(size=8),
    strip.text.x=element_text(size=7, face="bold"), strip.text.y=element_text(size=8, face="bold"),
    strip.background=element_blank(),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
    legend.position="bottom", legend.justification="left", legend.box.margin=margin(-0.3, 0, 0, 0, unit="cm"),
    legend.title=element_blank(), legend.text=element_text(size=6),
    legend.spacing.y=unit(-0.1, 'cm')
  ) +
  guides(color=guide_legend(order=1, nrow=3, byrow=TRUE))
ggsave(g, width=7.2, height=5, file="analyses/test/risk_recalibration.pdf")

# Show model calibration
risk_comp <- risk_comp[risk_model != "Risk recalibrated to European CVD low-risk region (SCORE2)"]
risk_comp[, risk_model := fcase(
  risk_model == "Kaplan-Meier estimate of observed risk in UK Biobank", "observed",
  risk_model == "Risk predicted using baseline hazards estimated across 44 cohorts (SCORE2)", "predicted"
)]
risk_comp <- dcast(risk_comp, sex + age_group + model ~ risk_model, value.var=c("mean_risk", "L95", "U95"))

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

# Stratify by risk decile
pred_risk <- fread("analyses/CVD_score_weighting/CVD_linear_predictors_and_risk.txt")
pred_risk <- pred_risk[order(score2_absrisk)][order(model)][order(sex)]
pred_risk[, risk_decile := floor((seq_len(.N)-1)/(.N/10))+1, by=.(model, sex)]

# Compute mean risk and 95% CI in each decile
risk_comp <- pred_risk[,.(mean_risk=mean(score2_absrisk), L95=quantile(score2_absrisk, 0.025), U95=quantile(score2_absrisk, 0.975)), by=.(sex, model, risk_decile)]

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


  
  





