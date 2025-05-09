library(data.table)
library(ggplot2)
library(scales)

# Load predicted risks
pred_scores <- rbind(idcol="cohort",
  discovery=fread("analyses/CVD_weight_training/CVD_linear_predictors_and_risk.txt"),
  replication=fread("analyses/CVD_weight_training/phase3_CVD_linear_predictors_and_risk.txt")
)
pred_scores <- pred_scores[endpoint == "cvd"] # linear predictors identical regardless of endpoint

# Extract absrisks for SCORE2
absrisk <- pred_scores[score == "SCORE2", .(eid, age, age_group, sex, incident_cvd, model, cohort, absrisk=uk_calibrated_risk)]

# Reformat age group for display
absrisk[, age_group := fcase(
  age_group %in% c("40-44", "45-49"), "40-49 years",
  age_group %in% c("50-54", "55-59"), "50-59 years",
  age_group %in% c("60-64", "65-69"), "60-69 years"
)]

# Compute emperical cumulative distribution functions in:
#  (1) Males and females separately
#  (2) For all participants, CVD cases, and non-cases
#  (3) For each decade of age
eCDFs <- rbind(idcol="incident_cvd",
  "All participants"=absrisk[,.(eid, age_group, sex, model, cohort, absrisk)], 
  "Incident CVD"=absrisk[(incident_cvd),.(eid, age_group, sex, model, cohort, absrisk)], 
  "Non-cases"=absrisk[!(incident_cvd),.(eid, age_group, sex, model, cohort, absrisk)]
)

eCDFs[, eCDF := ecdf(absrisk)(absrisk), by=.(incident_cvd, age_group, sex, model)]

eCDFs[, sex := factor(sex, levels=c("Male", "Female"))]
eCDFs[, age_group := factor(age_group, levels=c("40-49 years", "50-59 years", "60-69 years"))]
eCDFs[, incident_cvd := factor(incident_cvd)]
eCDFs[, model := factor(model, levels=c("SCORE2 + NMR scores", "SCORE2 + Biochemistry", "SCORE2 + PRS",
  "SCORE2 + NMR scores + PRS", "SCORE2 + Biochemistry + PRS", "SCORE2"))]

# Plot probability of exceeding X% risk
risk_thresholds <- rbind(
  data.table(age_group="40-49 years", threshold=0.025, threshold_type="low/medium risk"),
  data.table(age_group="40-49 years", threshold=0.075, threshold_type="medium/high risk"),
  data.table(age_group="50-59 years", threshold=0.05, threshold_type="low/medium risk"),
  data.table(age_group="50-59 years", threshold=0.10, threshold_type="medium/high risk"),
  data.table(age_group="60-69 years", threshold=0.05, threshold_type="low/medium risk"),
  data.table(age_group="60-69 years", threshold=0.10, threshold_type="medium/high risk")
)

g <- ggplot(eCDFs) +
  aes(x=absrisk, y=1-eCDF, color=model, linetype=cohort) +
  facet_grid(age_group ~ incident_cvd + sex) + 
  geom_line(linewidth=0.3) +
  geom_vline(data=risk_thresholds[threshold_type == "low/medium risk"], aes(xintercept=threshold), color="#feb24c", linetype=2, linewidth=0.6) + 
  geom_vline(data=risk_thresholds[threshold_type == "medium/high risk"], aes(xintercept=threshold), color="#fc4e2a", linetype=2, linewidth=0.6) + 
  scale_x_continuous("Predicted 10-year CVD risk (calibrated to low-risk region)", labels=scales::percent, limits=c(0,0.3), oob=squish, expand=expansion(mult=0, add=c(0.01, 0.03))) +
  scale_y_continuous("Probability of 10-year CVD risk exceeding X%") +
  scale_linetype_manual(values=c("discovery"="solid", "replication"="dashed")) +
	scale_color_manual(values=c(
		"SCORE2 + NMR scores"="#e41a1c", "SCORE2 + PRS"="#377eb8", "SCORE2 + Biochemistry"="#ff7f00",
		"SCORE2 + NMR scores + PRS"="#4daf4a", "SCORE2 + Biochemistry + PRS"="#984ea3", "SCORE2"="#000000"
	)) +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8), strip.background=element_blank(),
    strip.text=element_text(size=7, face="bold"), legend.position="bottom",
    legend.text=element_text(size=7), legend.title=element_blank()
  )
ggsave(g, width=6.8, height=4.44, file="analyses/test/SCORE2_risk_distribution_comparison.pdf")

# Repeat for QRISK3
absrisk <- pred_scores[score == "QRISK3", .(eid, age, age_group, sex, incident_cvd, model, cohort, absrisk=uk_calibrated_risk)]

absrisk[, age_group := fcase(
  age_group %in% c("40-44", "45-49"), "40-49 years",
  age_group %in% c("50-54", "55-59"), "50-59 years",
  age_group %in% c("60-64", "65-69"), "60-69 years"
)]

eCDFs <- rbind(idcol="incident_cvd",
  "All participants"=absrisk[,.(eid, age_group, sex, model, cohort, absrisk)],
  "Incident CVD"=absrisk[(incident_cvd),.(eid, age_group, sex, model, cohort, absrisk)],
  "Non-cases"=absrisk[!(incident_cvd),.(eid, age_group, sex, model, cohort, absrisk)]
)

eCDFs[, eCDF := ecdf(absrisk)(absrisk), by=.(incident_cvd, age_group, sex, model)]

eCDFs[, sex := factor(sex, levels=c("Male", "Female"))]
eCDFs[, age_group := factor(age_group, levels=c("40-49 years", "50-59 years", "60-69 years"))]
eCDFs[, incident_cvd := factor(incident_cvd)]
eCDFs[, model := factor(model, levels=c("QRISK3 + NMR scores", "QRISK3 + Biochemistry", "QRISK3 + PRS",
  "QRISK3 + NMR scores + PRS", "QRISK3 + Biochemistry + PRS", "QRISK3"))]

g <- ggplot(eCDFs) +
  aes(x=absrisk, y=1-eCDF, color=model, linetype=cohort) +
  facet_grid(age_group ~ incident_cvd + sex) +
  geom_line(linewidth=0.3) +
  geom_vline(xintercept=0.1, color="#fc4e2a", linetype=2, linewidth=0.6) +  # QRISK3 uses a single risk threshold of 10%
  scale_x_continuous("Predicted 10-year CVD risk (calibrated to low-risk region)", labels=scales::percent, limits=c(0,0.3), oob=oob_censor, expand=expansion(mult=0, add=c(0.01, 0.03))) +
  scale_y_continuous("Probability of 10-year CVD risk exceeding X%") +
  scale_linetype_manual(values=c("discovery"="solid", "replication"="dashed")) +
  scale_color_manual(values=c(
    "QRISK3 + NMR scores"="#e41a1c", "QRISK3 + PRS"="#377eb8", "QRISK3 + Biochemistry"="#ff7f00",
    "QRISK3 + NMR scores + PRS"="#4daf4a", "QRISK3 + Biochemistry + PRS"="#984ea3", "QRISK3"="#000000"
  )) +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8), strip.background=element_blank(),
    strip.text=element_text(size=7, face="bold"), legend.position="bottom",
    legend.text=element_text(size=7), legend.title=element_blank()
  )
ggsave(g, width=6.8, height=4.44, file="analyses/test/QRISK3_risk_distribution_comparison.pdf")

# What % of people were classified as high risk in each age group by SCORE2 vs. QRISK3?
risk_strata <- pred_scores[model_type == ""]
risk_strata[, high_risk := fcase(
  score == "QRISK3" & uk_calibrated_risk >= 0.1, TRUE,
  score != "QRISK3" & age >= 50 & uk_calibrated_risk >= 0.1, TRUE,
  score != "QRISK3" & age < 50 & uk_calibrated_risk >= 0.075, TRUE,
  default=FALSE)]
prop_high_risk <- risk_strata[,.(pct_high_risk=sum(high_risk)/.N),by=.(score, sex, age_group)]
prop_high_risk <- dcast(prop_high_risk, age_group ~ score + sex, value.var="pct_high_risk")
prop_high_risk <- prop_high_risk[order(age_group)]
fwrite(prop_high_risk, sep="\t", quote=FALSE, file="analyses/test/prop_high_risk_by_age_sex_score.txt")











