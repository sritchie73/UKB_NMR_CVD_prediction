library(data.table)
library(survival)
library(foreach)
library(ggplot2)
library(forcats)
library(ggh4x)
library(scales)
library(ggthemes)

# Make output directory
system("mkdir -p analyses/risk_recalibration/")

# We want to compute the scaling factors needed to recalibrate risk so that the group as a whole
# has risk distributions similar to the general UK population.

# First, here are the per-1000 year CVD incidence rates from CPRD published by Sun et al. 2021:
CPRD = data.table(
  age_at_risk_start = seq(40, 75, by=5),
  age_at_risk_end = seq(44, 79, by=5),
  male_rate_per_1000py = c(1.578, 2.964, 5.141, 7.561, 10.848, 13.731, 19.271, 25.667),
  female_rate_per_1000py = c(0.726, 1.309, 2.165, 2.969, 4.560, 7.007, 10.840, 17.007)
)

# Melt table to better format
CPRD <- melt(CPRD, id.vars=c("age_at_risk_start", "age_at_risk_end"), variable.name="sex", value.name="rate_per_1000py")
CPRD[, sex := ifelse(sex == "male_rate_per_1000py", "Male", "Female")]

# The expected 10-year CVD risk is calculated for each age group
# based on the mid-point of the next interval ahead, e.g. for the 
# 40-44 year age-group the expected 10-year risk is calculated based
# on the annual incidence rates in the 45-49 year olds in CPRD:
CPRD[, age_group := sprintf("%s-%s", age_at_risk_start - 5, age_at_risk_end - 5)]
CPRD[, age_at_risk := sprintf("%s-%s", age_at_risk_start, age_at_risk_end)]

# Filter to eligible age groups for SCORE2 screening
CPRD <- CPRD[age_at_risk_start >= 45 & age_at_risk_end < 75]

# Reorganize
CPRD <- CPRD[, .(sex, age_group, age_at_risk, rate_per_1000py)][order(age_group)][order(-sex)]

# Calculate annual CVD incidence in each age-at-risk group from the 1000 person-year CVD rates:
CPRD[, annual_incidence := rate_per_1000py / 1000]

# The expected risk is calculated assuming exponential survival (i.e.
# constant hazard) in the 10-years ahead.
CPRD[, expected_risk := 1 - exp(-annual_incidence * 10)]

# Write out
fwrite(CPRD, sep="\t", quote=FALSE, file="analyses/risk_recalibration/CPRD_incidence_and_risk.txt")

# Load predicted risk and add 5-year age groups
pred_risk <- fread("analyses/CVD_score_weighting/CVD_linear_predictors_and_risk.txt")
pred_risk[, age_group := sprintf("%s-%s", age %/% 5 * 5, age %/% 5 * 5 + 4)]

# Build comparison table of expected 10-year risk so we can fit recalibration models
comp_risk <- pred_risk[, .(predicted=mean(uk_risk)), by=.(model, sex, age_group)]
comp_risk[CPRD, on = .(sex, age_group), CPRD := i.expected_risk]

# Transform with link function
comp_risk[, CPRD := log(-log(1 - CPRD))]
comp_risk[, predicted := log(-log(1 - predicted))]

# Fit linear models to get scaling factors
recalibration_fit <- comp_risk[,as.list(coef(lm(CPRD ~ predicted))), by=.(model, sex)]
setnames(recalibration_fit, c("model", "sex", "scale1", "scale2"))
fwrite(recalibration_fit, sep="\t", quote=FALSE, file="analyses/risk_recalibration/recalibration_scaling_factors.txt")

# Compute recalibrated risk - we do this from the UK-calibrated absolute risk: i.e. we're
# transforming the "true" risk of UK Biobank participants to the risk distribution that we would expect
# in age-groups of those numbers if we had instead sampled from the general UK population
pred_risk[recalibration_fit, on = .(model, sex), 
  recalibrated_risk := 1 - exp(-exp(scale1 + scale2 * log(-log(1 - uk_risk)))),
]

# Write out 
fwrite(pred_risk, sep="\t", quote=FALSE, file="analyses/risk_recalibration/CVD_linear_predictors_and_risk.txt")

# Load incident CVD data and add 5-year age groups
dat <- fread("data/cleaned/analysis_cohort.txt", select=c("eid", "sex", "age", "incident_cvd_followup", "incident_cvd", "mortality_at_cvd_followup"))
dat[, age_group := sprintf("%s-%s", age %/% 5 * 5, age %/% 5 * 5 + 4)]

# Now compute the observed risk in UK Biobank using Kaplan Meier estimates
ukb_inci <- foreach(this_sex = c("Male", "Female"), .combine=rbind) %:%
  foreach(this_age_group = sort(unique(dat$age_group)), .combine=rbind) %do% {
    km <- survfit(Surv(incident_cvd_followup, incident_cvd) ~ 1, data=dat[sex == this_sex & age_group == this_age_group])
    kms <- summary(km, time=10)
    
    m1 <- data.table(estimate=kms$surv, L95=kms$lower, U95=kms$upper)                      # Kaplan Meier estimate of 10-year survival and 95% CI
    m2 <- data.table(estimate=1-kms$surv, L95=1-kms$upper, U95 =1-kms$lower)               # Expected 10-year risk and 95% CI
    m3 <- m2[, .(estimate=-log(1 - estimate)/10, L95=-log(1-L95)/10, U95=-log(1-U95)/10)]  # Annual incidence and 95% CI
    m4 <- m3[, .(estimate=estimate*1000, L95=L95*1000, U95=U95*1000)]                      # Rates per 1000-person years and 95% CI

    metrics <- rbind(idcol="metric",
      "10-year survival"=m1,
      "Expected 10-year risk"=m2,
      "Annual incidence"=m3,
      "Rates per 1000 person years"=m4
    )

    info <- data.table(sex=this_sex, age_group=this_age_group)
    cbind(info, metrics)
}
fwrite(ukb_inci, sep="\t", quote=FALSE, file="analyses/risk_recalibration/UKB_KaplanMeier_incidence_and_risk.txt")

####### Now plot sanity checks

# Compute risk quantiles in each age group for each model and risk calculation
risk_comp <- melt(pred_risk, measure.vars=c("ukb_absrisk", "score2_absrisk", "uk_risk", "recalibrated_risk"), variable.name="risk_model", value.name="risk")
risk_comp <- risk_comp[,.(mean_risk=mean(risk), L95=quantile(risk, 0.025), U95=quantile(risk, 0.975)),
  by=.(sex, age_group, model, risk_model)]
risk_comp[, risk_model := fcase(
  risk_model == "ukb_absrisk", "Risk predicted from Cox models fit in UK Biobank (Aki's code adapted to cross-validation)",
  risk_model == "score2_absrisk", "Risk predicted using uncalibrated SCORE2 formula (baseline hazards estimated across 44 cohorts)",
  risk_model == "uk_risk", "Risk predicted using SCORE2 European CVD low-risk region formula",
  risk_model == "recalibrated_risk", "Predicted risk recalibrated to CPRD-derived estimates of CVD incidence rates"
)]

# Add in CPRD-dervied estimates of CVD incidence published in Sun et al. 2021
toadd <- foreach(this_model = c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %do% {
  CPRD[, .(sex, age_group, model=this_model, risk_model="CVD incidence rates in CPRD (Sun et al. 2021)", mean_risk=expected_risk)]
}
risk_comp <- rbind(risk_comp, toadd, fill=TRUE)

# Add in Kaplan-Meier estimates of observed risk in UK Biobank
toadd <- foreach(this_model = c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %do% {
  ukb_inci[metric == "Expected 10-year risk", .(sex, age_group, model=this_model, risk_model="Kaplan-Meier estimate of 10-year risk in UK Biobank", mean_risk=estimate, L95, U95)]
}
risk_comp <- rbind(risk_comp, toadd)

# Write out risk tables
fwrite(risk_comp, sep="\t", quote=FALSE, file="analyses/risk_recalibration/risk_table_by_sex_age_model.txt")

# Prepare for plotting
risk_comp[, sex := factor(paste0(sex, "s"), levels=c("Males", "Females"))]
risk_comp[, age_group := factor(age_group, levels=c("40-44", "45-49", "50-54", "55-59", "60-64", "65-69"))]
risk_comp[, model := factor(model, levels=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
risk_comp[, risk_model := factor(risk_model, levels=c(
  "Kaplan-Meier estimate of 10-year risk in UK Biobank",
  "Risk predicted from Cox models fit in UK Biobank (Aki's code adapted to cross-validation)",
  "Risk predicted using uncalibrated SCORE2 formula (baseline hazards estimated across 44 cohorts)",
  "Risk predicted using SCORE2 European CVD low-risk region formula",
  "Predicted risk recalibrated to CPRD-derived estimates of CVD incidence rates",
  "CVD incidence rates in CPRD (Sun et al. 2021)"
))]

# Sanity check SCORE2
ggdt <- risk_comp[model == "SCORE2"]
g <- ggplot(ggdt) +
  aes(x=age_group, y=mean_risk, ymin=L95, ymax=U95, color=risk_model) +
  facet_grid(model ~ sex) +
  geom_errorbar(width=0, position=position_dodge(width=0.75)) +
  geom_point(shape=23, fill="white", size=1.2, position=position_dodge(width=0.75)) +
  geom_errorbar(data=ggdt[risk_model == "Kaplan-Meier estimate of 10-year risk in UK Biobank"], width=0, position=position_nudge(x=-0.315), show.legend=FALSE) +
  scale_color_manual(values=c(
    "Kaplan-Meier estimate of 10-year risk in UK Biobank"="#1f78b4",
    "Risk predicted from Cox models fit in UK Biobank (Aki's code adapted to cross-validation)"="#a6cee3",
    "Risk predicted using uncalibrated SCORE2 formula (baseline hazards estimated across 44 cohorts)"="#b2df8a",
    "Risk predicted using SCORE2 European CVD low-risk region formula"="#33a02c",
    "Predicted risk recalibrated to CPRD-derived estimates of CVD incidence rates"="#fb9a99",
    "CVD incidence rates in CPRD (Sun et al. 2021)"="#e31a1c"
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
    legend.position="bottom", legend.justification="left", legend.box.margin=margin(-0.3, 0, 0, -1, unit="cm"),
    legend.title=element_blank(), legend.text=element_text(size=6),
    legend.spacing.y=unit(-0.1, 'cm')
  ) +
  guides(color=guide_legend(order=1, nrow=6, byrow=TRUE))
ggsave(g, width=7.2, height=5, file="analyses/risk_recalibration/SCORE2_risk_computation_check.pdf")
  
# Sanity check all models
g <- ggplot(risk_comp) +
  aes(x=age_group, y=mean_risk, ymin=L95, ymax=U95, color=risk_model) +
  facet_grid(sex ~ model) +
  geom_errorbar(width=0, position=position_dodge(width=0.75)) +
  geom_point(shape=23, fill="white", size=1.2, position=position_dodge(width=0.75)) +
  geom_errorbar(data=risk_comp[risk_model == "Kaplan-Meier estimate of 10-year risk in UK Biobank"], width=0, position=position_nudge(x=-0.315), show.legend=FALSE) +
  scale_color_manual(values=c(
    "Kaplan-Meier estimate of 10-year risk in UK Biobank"="#1f78b4",
    "Risk predicted from Cox models fit in UK Biobank (Aki's code adapted to cross-validation)"="#a6cee3",
    "Risk predicted using uncalibrated SCORE2 formula (baseline hazards estimated across 44 cohorts)"="#b2df8a",
    "Risk predicted using SCORE2 European CVD low-risk region formula"="#33a02c",
    "Predicted risk recalibrated to CPRD-derived estimates of CVD incidence rates"="#fb9a99",
    "CVD incidence rates in CPRD (Sun et al. 2021)"="#e31a1c"
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
    legend.position="bottom", legend.justification="left", legend.box.margin=margin(-0.3, 0, 0, -1, unit="cm"),
    legend.title=element_blank(), legend.text=element_text(size=6),
    legend.spacing.y=unit(-0.1, 'cm')
  ) +
  guides(color=guide_legend(order=1, nrow=6, byrow=TRUE))
ggsave(g, width=7.2, height=5, file="analyses/risk_recalibration/risk_computation_check.pdf")

# Show selected predictions of interest
ggdt <- risk_comp[!(risk_model %in% c(
  "Risk predicted from Cox models fit in UK Biobank (Aki's code adapted to cross-validation)",
  "Risk predicted using uncalibrated SCORE2 formula (baseline hazards estimated across 44 cohorts)"
))]
g <- ggplot(ggdt) +
  aes(x=age_group, y=mean_risk, ymin=L95, ymax=U95, color=risk_model) +
  facet_grid(sex ~ model) +
  geom_errorbar(width=0, position=position_dodge(width=0.75)) +
  geom_point(shape=23, fill="white", size=1.2, position=position_dodge(width=0.75)) +
  geom_errorbar(data=risk_comp[risk_model == "Kaplan-Meier estimate of 10-year risk in UK Biobank"], width=0, position=position_nudge(x=-0.285), show.legend=FALSE) +
  scale_color_manual(values=c(
    "Kaplan-Meier estimate of 10-year risk in UK Biobank"="#377eb8",
    "Risk predicted using SCORE2 European CVD low-risk region formula"="#33a02c",
    "Predicted risk recalibrated to CPRD-derived estimates of CVD incidence rates"="#ff7f00",
    "CVD incidence rates in CPRD (Sun et al. 2021)"="#e41a1c"
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
    legend.position="bottom", legend.justification="left", legend.box.margin=margin(-0.3, 0, 0, -1, unit="cm"),
    legend.title=element_blank(), legend.text=element_text(size=6),
    legend.spacing.y=unit(-0.1, 'cm')
  ) +
  guides(color=guide_legend(order=1, nrow=6, byrow=TRUE))
ggsave(g, width=7.2, height=5, file="analyses/risk_recalibration/risk_and_recalibration.pdf")

# And for SCORE2
ggdt <- risk_comp[model == "SCORE2" & risk_model != "Risk predicted from Cox models fit in UK Biobank (Aki's code adapted to cross-validation)"]
g <- ggplot(ggdt) +
  aes(x=age_group, y=mean_risk, ymin=L95, ymax=U95, color=risk_model) +
  facet_grid(model ~ sex) +
  geom_errorbar(width=0, position=position_dodge(width=0.75)) +
  geom_point(shape=23, fill="white", size=1.2, position=position_dodge(width=0.75)) +
  geom_errorbar(data=ggdt[risk_model == "Kaplan-Meier estimate of 10-year risk in UK Biobank"], width=0, position=position_nudge(x=-0.3), show.legend=FALSE) +
  scale_color_manual(values=c(
    "Kaplan-Meier estimate of 10-year risk in UK Biobank"="#377eb8",
    "Risk predicted using uncalibrated SCORE2 formula (baseline hazards estimated across 44 cohorts)"="#b2df8a",
    "Risk predicted using SCORE2 European CVD low-risk region formula"="#33a02c",
    "Predicted risk recalibrated to CPRD-derived estimates of CVD incidence rates"="#ff7f00",
    "CVD incidence rates in CPRD (Sun et al. 2021)"="#e41a1c"
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
    legend.position="bottom", legend.justification="left", legend.box.margin=margin(-0.3, 0, 0, -1, unit="cm"),
    legend.title=element_blank(), legend.text=element_text(size=6),
    legend.spacing.y=unit(-0.1, 'cm')
  ) +
  guides(color=guide_legend(order=1, nrow=6, byrow=TRUE))
ggsave(g, width=7.2, height=5, file="analyses/risk_recalibration/SCORE2_risk_computation_check.pdf")

