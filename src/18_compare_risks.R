library(data.table)
library(ggplot2)
library(scales)
library(forcats)

# Load predicted risk
pred_risk <- fread("analyses/CVD_weight_training/CVD_linear_predictors_and_risk.txt")
pred_risk <- pred_risk[score_type == "non-derived"]

# Compute number and percentage of people at different risk thresholds
pred_risk <- pred_risk[,.(sex, model, age_group, uk_calibrated_risk)]
pred_risk[, risk_group := fcase(
  uk_calibrated_risk >= 0.25, "≥25%",
  uk_calibrated_risk > 0.15, "15%–<25%",
  uk_calibrated_risk > 0.1, "10%–<15%",
  uk_calibrated_risk > 0.075, "7.5%–<10%",
  uk_calibrated_risk > 0.05, "5%–<7.5%",
  uk_calibrated_risk > 0.025, "2.5%–<5%",
  default = "<2.5%"
)]

risk_strata <- pred_risk[, .N, by=.(sex, model, age_group, risk_group)]
risk_strata[pred_risk[,.N,by=.(sex, model, age_group)], on = .(sex, model, age_group), pct := N/i.N]

# Write out
fwrite(risk_strata, sep="\t", quote=FALSE, file="analyses/test/stratification_by_SCORE2_risk_groups.txt")

# Set up labels for plotting
risk_strata[, risk_group := factor(risk_group, levels=c("<2.5%", "2.5%–<5%", "5%–<7.5%", "7.5%–<10%", "10%–<15%", "15%–<25%", "≥25%"))]
risk_strata[, sex := factor(paste0(sex, "s"), levels=c("Males", "Females"))]
risk_strata[, age_group := factor(age_group, levels=c("40-44", "45-49", "50-54", "55-59", "60-64", "65-69"))]
risk_strata[, model := factor(model, levels=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]

# Plot
g <- ggplot(risk_strata) +
  aes(x=age_group, y=N, fill=risk_group) +
  facet_grid(sex ~ model) +
  geom_col(position="fill") +
  geom_hline(yintercept=seq(0,1,by=0.1), color="white", alpha=0.6, linewidth=0.4) +
  scale_y_continuous("UK Biobank participants", labels=percent, expand=c(0,0.04)) +
  scale_x_discrete("Age-group", expand=c(0,0.6)) +
  scale_fill_manual("Predicted 10-year risk of CVD recalibrated to low-risk European region", values=c(
    "<2.5%"="#34b340", "2.5%–<5%"="#71ff50", "5%–<7.5%"="#efff00", "7.5%–<10%"="#eacc00", 
    "10%–<15%"="#e41f1b", "15%–<25%"="#b01815", "≥25%"="#680e0c"
  )) +
  guides(fill=guide_legend(title.position="top", nrow=1, reverse=TRUE)) +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=6), axis.title.y=element_text(size=8),
    axis.text.x=element_text(size=6, angle=45, hjust=1), axis.title.x=element_text(size=8),
    strip.text.x=element_text(size=7, face="bold"), strip.text.y=element_text(size=8, face="bold"),
    strip.background=element_blank(), panel.grid=element_blank(),
    legend.title=element_text(size=8), legend.text=element_text(size=6), 
    legend.box.margin=margin(-0.5, 0, 0, -1, unit="cm"),
    legend.position="bottom", legend.justification="left", 
    legend.box.background=element_blank(), legend.background=element_blank()
  )
ggsave(g, width=7.2, height=4, file="analyses/test/stratification_by_SCORE2_risk_groups_by_age_group.pdf", device=cairo_pdf)

# Plot aggregate across age groups
risk_strata <- risk_strata[,.(N=sum(N)), by=.(sex, model, risk_group)]

g <- ggplot(risk_strata) +
  aes(y=fct_rev(model), x=N, fill=risk_group) +
  facet_wrap(~ sex) +
  geom_col(position="fill") +
  geom_vline(xintercept=seq(0,1,by=0.1), color="white", alpha=0.6, linewidth=0.6) +
  scale_x_continuous("UK Biobank participants", labels=percent, expand=c(0,0.01)) +
  scale_y_discrete("", expand=c(0,0.6)) +
  scale_fill_manual("Predicted 10-year risk of CVD recalibrated to low-risk European region", values=c(
    "<2.5%"="#34b340", "2.5%–<5%"="#71ff50", "5%–<7.5%"="#efff00", "7.5%–<10%"="#eacc00", 
    "10%–<15%"="#e41f1b", "15%–<25%"="#b01815", "≥25%"="#680e0c"
  )) +
  guides(fill=guide_legend(title.position="top", nrow=1, reverse=TRUE)) +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"),
    strip.background=element_blank(), panel.grid=element_blank(),
    legend.title=element_text(size=8), legend.text=element_text(size=6), 
    legend.box.margin=margin(-0.5, 0, 0, -4.2, unit="cm"),
    legend.position="bottom", legend.justification="left", 
    legend.box.background=element_blank(), legend.background=element_blank()
  )
ggsave(g, width=7.2, height=2, file="analyses/test/stratification_by_SCORE2_risk_groups.pdf", device=cairo_pdf)

