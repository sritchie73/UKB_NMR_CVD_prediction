library(data.table)
library(ggplot2)
library(scales)

# Load predicted risks
pred_scores <- fread("analyses/CVD_weight_training/phase3_CVD_linear_predictors_and_risk.txt")
pred_scores <- pred_scores[score_type == "non-derived"]

# Extract absrisks
absrisk <- pred_scores[,.(eid, age, age_group, sex, incident_cvd, score_type, model, absrisk=uk_calibrated_risk)]

# Compute emperical cumulative distribution functions in:
#  (1) Males and females separately
#  (2) For all participants, CVD cases, and non-cases
#  (3) For each decade of age
absrisk[, age_group := fcase(
  age_group %in% c("40-44", "45-49"), "40-49 years",
  age_group %in% c("50-54", "55-59"), "50-59 years",
  age_group %in% c("60-64", "65-69"), "60-69 years"
)]

eCDFs <- rbind(idcol="incident_cvd",
  "All participants"=absrisk[,.(eid, age_group, sex, model, absrisk)], 
  "Incident CVD"=absrisk[(incident_cvd),.(eid, age_group, sex, model, absrisk)], 
  "Non-cases"=absrisk[!(incident_cvd),.(eid, age_group, sex, model, absrisk)]
)

eCDFs[, eCDF := ecdf(absrisk)(absrisk), by=.(incident_cvd, age_group, sex, model)]

eCDFs[, sex := factor(sex, levels=c("Male", "Female"))]
eCDFs[, age_group := factor(age_group, levels=c("40-49 years", "50-59 years", "60-69 years"))]
eCDFs[, incident_cvd := factor(incident_cvd)]
eCDFs[, model := factor(model, levels=c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs", "SCORE2"))]

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
  aes(x=absrisk, y=1-eCDF, color=model) +
  facet_grid(age_group ~ incident_cvd + sex) + 
  geom_line(linewidth=0.6) +
  scale_color_manual(values=c("SCORE2"="black", "SCORE2 + NMR scores"="#e41a1c", "SCORE2 + PRSs"="#377eb8", "SCORE2 + NMR scores + PRSs"="#4daf4a")) +
  geom_vline(data=risk_thresholds[threshold_type == "low/medium risk"], aes(xintercept=threshold), color="#feb24c", linetype=2, linewidth=0.6) + 
  geom_vline(data=risk_thresholds[threshold_type == "medium/high risk"], aes(xintercept=threshold), color="#fc4e2a", linetype=2, linewidth=0.6) + 
  scale_x_continuous("Predicted 10-year CVD risk (calibrated to low-risk region)", labels=scales::percent, limits=c(0,0.3), oob=squish, expand=expansion(mult=0, add=c(0.01, 0.03))) +
  scale_y_continuous("Probability of 10-year CVD risk exceeding X%") +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8), strip.background=element_blank(),
    strip.text=element_text(size=7, face="bold"), legend.position="bottom",
    legend.text=element_text(size=7), legend.title=element_blank()
  )
ggsave(g, width=6.8, height=4.2, file="analyses/test/phase3_risk_distribution_comparison.pdf")








