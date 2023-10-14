library(data.table)
library(ggplot2)
library(scales)
library(ggstance)
library(forcats)
library(cowplot)

###############################
# Sex-stratified main figure
################################

# Load required data
cinds <- fread("analyses/test/cindices.txt")
risk <- fread("analyses/CVD_weight_training/CVD_linear_predictors_and_risk.txt")
nri <- fread("analyses/test/categorical_nri_estimates.txt")

# Prepare for plotting
cinds <- cinds[sex == "Sex-stratified" & score_type == "non-derived"]
cinds[, model := factor(model, levels=c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]

risk <- risk[score_type == "non-derived"]
risk[, model := factor(model, levels=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
risk[, risk_group := fcase(
  age < 50 & uk_calibrated_risk < 0.025, "Low-to-moderate risk",
  age < 50 & uk_calibrated_risk >= 0.025 & uk_calibrated_risk < 0.075, "High risk",
  age < 50 & uk_calibrated_risk >= 0.075, "Very high risk",
  age >= 50 & uk_calibrated_risk < 0.05, "Low-to-moderate risk",
  age >= 50 & uk_calibrated_risk >= 0.05 & uk_calibrated_risk < 0.10, "High risk",
  age >= 50 & uk_calibrated_risk >= 0.10, "Very high risk"
)]
risk[, risk_group := factor(risk_group, levels=c("Low-to-moderate risk", "High risk", "Very high risk"))]
risk[, cvd_group := ifelse(incident_cvd, "Cases", "Non-cases")]
risk <- risk[, .N, by=.(cvd_group, model, risk_group)]

nri <- nri[sex == "Sex-stratified" & metric %in% c("NRI+", "NRI-")]
nri[, model := gsub("SCORE2 vs. ", "", model_comparison)]
nri[, model := factor(model, levels=c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
nri[, status:= ifelse(metric == "NRI+", "Incident CVD cases", "Non-cases")]

# Make plot
g1 <- ggplot(cinds) +
  aes(x=deltaC, xmin=deltaC.L95, xmax=deltaC.U95, y=fct_rev(model)) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, size=2, fill="white") +
  ylab("") + xlab("Change in C-index (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )

g2 <- ggplot(risk) +
  aes(x=N, y=fct_rev(model), fill=risk_group) +
  facet_grid(. ~ cvd_group) +
  geom_col(position="fill") +
  geom_vline(xintercept=seq(0,1,by=0.1), color="white", alpha=0.6, linewidth=0.4) +
  scale_x_continuous("UK Biobank participants", labels=percent, expand=c(0,0.01)) +
  scale_y_discrete(expand=c(0,0.4)) +
  scale_fill_manual("Predicted ESC 2021 10-year CVD risk category", values=c(
    "Low-to-moderate risk"="#34b340", "High risk"="#eacc00", "Very high risk"="#b01815"
  )) +
  guides(fill=guide_legend(title.position="top", nrow=1, byrow=TRUE, reverse=TRUE)) +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8), panel.grid=element_blank(),
    legend.title=element_text(size=8), legend.text=element_text(size=6),
    legend.box.margin=margin(-0.5, 0, 0, -3, unit="cm"),
    legend.position="bottom", legend.justification="left",
    legend.box.background=element_blank(), legend.background=element_blank()
  )

g3 <- ggplot(nri) +
  aes(x=Estimate, xmin=Lower, xmax=Upper, y=fct_rev(model), color=status) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0, position=position_dodgev(height=0.3)) +
  geom_point(shape=23, fill="white", size=1.2, position=position_dodgev(height=0.3)) +
  scale_color_manual(values=c("Incident CVD cases"="#c51b7d", "Non-cases"="#4d9221")) +
  scale_x_continuous("Categorical NRI, % reclassified (95% CI)", labels=percent) +
  ylab("") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.title=element_blank(), legend.text=element_text(size=8, color="black"),
    legend.box.margin=margin(-0.5, 0, 0, -6, unit="cm"),
    legend.position="bottom", legend.justification="left",
    legend.box.background=element_blank(), legend.background=element_blank()
  )

g <- plot_grid(g1, g2, g3, nrow=1, align="h", axis="l")
ggsave(g, width=10, height=3, file="analyses/test/incremental_improvement.pdf")

###############################
# Sex-specific supp figure
################################

# Load required data
cinds <- fread("analyses/test/cindices.txt")
risk <- fread("analyses/CVD_weight_training/CVD_linear_predictors_and_risk.txt")
nri <- fread("analyses/test/categorical_nri_estimates.txt")

# Prepare for plotting
cinds <- cinds[sex != "Sex-stratified" & score_type == "non-derived"]
cinds[, model := factor(model, levels=c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
cinds[, sex := factor(sex, levels=c("Males", "Females"))]

risk <- risk[score_type == "non-derived"]
risk[, model := factor(model, levels=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
risk[, risk_group := fcase(
  age < 50 & uk_calibrated_risk < 0.025, "Low-to-moderate risk",
  age < 50 & uk_calibrated_risk >= 0.025 & uk_calibrated_risk < 0.075, "High risk",
  age < 50 & uk_calibrated_risk >= 0.075, "Very high risk",
  age >= 50 & uk_calibrated_risk < 0.05, "Low-to-moderate risk",
  age >= 50 & uk_calibrated_risk >= 0.05 & uk_calibrated_risk < 0.10, "High risk",
  age >= 50 & uk_calibrated_risk >= 0.10, "Very high risk"
)]
risk[, risk_group := factor(risk_group, levels=c("Low-to-moderate risk", "High risk", "Very high risk"))]
risk[, cvd_group := ifelse(incident_cvd, "Cases", "Non-cases")]
risk <- risk[, .N, by=.(sex, cvd_group, model, risk_group)]
risk[, sex := factor(sex, levels=c("Male", "Female"))]

nri <- nri[sex != "Sex-stratified" & metric %in% c("NRI+", "NRI-")]
nri[, model := gsub("SCORE2 vs. ", "", model_comparison)]
nri[, model := factor(model, levels=c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
nri[, sex := factor(sex, levels=c("Males", "Females"))]
nri[, status:= ifelse(metric == "NRI+", "Incident CVD cases", "Non-cases")]

# Make plot
g1 <- ggplot(cinds) +
  aes(x=deltaC, xmin=deltaC.L95, xmax=deltaC.U95, y=fct_rev(model)) +
  facet_grid(sex ~ .) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, size=2, fill="white") +
  ylab("") + xlab("Change in C-index (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.background=element_blank(), strip.text=element_text(size=8, face="bold"),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )

g2 <- ggplot(risk) +
  aes(x=N, y=fct_rev(model), fill=risk_group) +
  facet_grid(sex ~ cvd_group) +
  geom_col(position="fill") +
  geom_vline(xintercept=seq(0,1,by=0.1), color="white", alpha=0.6, linewidth=0.4) +
  scale_x_continuous("UK Biobank participants", labels=percent, expand=c(0,0.01)) +
  scale_y_discrete(expand=c(0,0.4)) +
  scale_fill_manual("Predicted ESC 2021 10-year CVD risk category", values=c(
    "Low-to-moderate risk"="#34b340", "High risk"="#eacc00", "Very high risk"="#b01815"
  )) +
  guides(fill=guide_legend(title.position="top", nrow=1, byrow=TRUE, reverse=TRUE)) +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8), panel.grid=element_blank(),
    strip.background=element_blank(), strip.text=element_text(size=8, face="bold"),
    legend.title=element_text(size=8), legend.text=element_text(size=6),
    legend.box.margin=margin(-0.5, 0, 0, -3, unit="cm"),
    legend.position="bottom", legend.justification="left",
    legend.box.background=element_blank(), legend.background=element_blank()
  )

g3 <- ggplot(nri) +
  aes(x=Estimate, xmin=Lower, xmax=Upper, y=fct_rev(model), color=status) +
  facet_grid(sex ~ .) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0, position=position_dodgev(height=0.3)) +
  geom_point(shape=23, fill="white", size=1.2, position=position_dodgev(height=0.3)) +
  scale_color_manual(values=c("Incident CVD cases"="#c51b7d", "Non-cases"="#4d9221")) +
  scale_x_continuous("Categorical NRI, % reclassified (95% CI)", labels=percent) +
  ylab("") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.background=element_blank(), strip.text=element_text(size=8, face="bold"),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.title=element_blank(), legend.text=element_text(size=8, color="black"),
    legend.box.margin=margin(-0.5, 0, 0, -6, unit="cm"),
    legend.position="bottom", legend.justification="left",
    legend.box.background=element_blank(), legend.background=element_blank()
  )

g <- plot_grid(g1, g2, g3, nrow=1, align="h", axis="l")
ggsave(g, width=10, height=5, file="analyses/test/incremental_improvement_sex_specific.pdf")

 

