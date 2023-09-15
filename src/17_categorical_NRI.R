# Takes 15-20 minutes to run on a compute node
library(data.table)
library(foreach)
library(nricens)
library(ggplot2)
library(forcats)
library(ggstance)
library(scales)
library(cowplot)

# Load required data
dat <- fread('analyses/CVD_weight_training/CVD_linear_predictors_and_risk.txt')
dat <- dat[score_type == "non-derived"]

# Wrapper function for categorical NRI test
nri.test <- function(data, base_model, new_model, absrisk_column="uk_calibrated_risk") {
  # Extract predicted 10 year risk
  comp_risk <- data[model == base_model | model == new_model]
  comp_risk[, model := fcase(model == base_model, "base", model == new_model, "new")]
  comp_risk <- dcast(comp_risk, eid + age + incident_cvd_followup + incident_cvd ~ model, value.var=absrisk_column)

  # The ESC 2021 guidelines use different risk thresolds for different age groups:
  #   <2.5%, 2.5-7.5%, and >7.5% in people under 50
  #   <5%, 5-10%, and >10% in people over 50
  # To compute a single NRI across all ages, we can therefor simply add 2.5% risk to everyone under 50 so
  # we can then use the same thresholds in everyone
  comp_risk[age < 50, base := base + 0.025]
  comp_risk[age < 50, new := new + 0.025]

  # Run NRI analysis
  NRI <- nricens(event = comp_risk$incident_cvd, time = comp_risk$incident_cvd_followup,
                 p.std = comp_risk$base, p.new = comp_risk$new,
                 updown = "category", cut = c(0.05, 0.1), t0 = 10, niter = 1000)
  system("rm -f Rplots.pdf")

  # Add in sample size and case numbers
  NRI$n <- comp_risk[,.N]
  NRI$nevent <- comp_risk[, sum(incident_cvd)]
  
  return(NRI)
}

nri_lists <- list(
  "Males"=list(
    "SCORE2 vs. SCORE2 + NMR scores"=nri.test(dat[sex == "Male"], "SCORE2", "SCORE2 + NMR scores"),
    "SCORE2 vs. SCORE2 + PRSs"=nri.test(dat[sex == "Male"], "SCORE2", "SCORE2 + PRSs"),
    "SCORE2 vs. SCORE2 + NMR scores + PRSs"=nri.test(dat[sex == "Male"], "SCORE2", "SCORE2 + NMR scores + PRSs")
  ),
  "Females"=list(
    "SCORE2 vs. SCORE2 + NMR scores"=nri.test(dat[sex == "Female"], "SCORE2", "SCORE2 + NMR scores"),
    "SCORE2 vs. SCORE2 + PRSs"=nri.test(dat[sex == "Female"], "SCORE2", "SCORE2 + PRSs"),
    "SCORE2 vs. SCORE2 + NMR scores + PRSs"=nri.test(dat[sex == "Female"], "SCORE2", "SCORE2 + NMR scores + PRSs")
  ),
  "Sex-stratified"=list(
    "SCORE2 vs. SCORE2 + NMR scores"=nri.test(dat, "SCORE2", "SCORE2 + NMR scores"),
    "SCORE2 vs. SCORE2 + PRSs"=nri.test(dat, "SCORE2", "SCORE2 + PRSs"),
    "SCORE2 vs. SCORE2 + NMR scores + PRSs"=nri.test(dat, "SCORE2", "SCORE2 + NMR scores + PRSs")
  )
)
saveRDS(nri_lists, file="analyses/test/categorical_nri.rds")

# Extract tables of estimates
nri_estimates <- rbindlist(idcol="sex", fill=TRUE, lapply(nri_lists, function(l1) {
  rbindlist(idcol="model_comparison", fill=TRUE, lapply(l1, function(l2) {
    cbind(samples=l2$n, cases=l2$nevent, as.data.table(l2$nri, keep.rownames="metric"))
  }))
}))


# Extract tables of reclassifications for categorical nris
reclassified_cases <- rbindlist(idcol="sex", fill=TRUE, lapply(nri_lists, function(l1) {
  rbindlist(idcol="model_comparison", fill=TRUE, lapply(l1, function(l2) {
    as.data.table(l2[["rtab.case"]])
  }))
}))

reclassified <- rbindlist(idcol="sex", fill=TRUE, lapply(nri_lists, function(l1) {
  rbindlist(idcol="model_comparison", fill=TRUE, lapply(l1, function(l2) {
    as.data.table(l2[["rtab"]])
  }))
}))

# Merge
setnames(reclassified, "N", "All")
reclassified[reclassified_cases, on = .(sex, model_comparison, Standard, New), Cases := N]
setnames(reclassified, "Standard", "Old")

# Add in total sample size and total cases to reclassified table
reclassified[nri_estimates, on = .(sex, model_comparison), Total_Samples := i.samples]
reclassified[nri_estimates, on = .(sex, model_comparison), Total_Cases := i.cases]

fwrite(nri_estimates, sep="\t", quote=FALSE, file="analyses/test/categorical_nri_estimates.txt")
fwrite(reclassified, sep="\t", quote=FALSE, file="analyses/test/categorical_nri_reclassified.txt")

# Plot
ggdt <- nri_estimates[metric %in% c("NRI+", "NRI-")]
ggdt[, type := ifelse(metric == "NRI+", "CVD cases", "Non-cases")]
ggdt[, type := factor(type, levels=c("CVD cases", "Non-cases"))]
ggdt[, sex := factor(sex, levels=c("Males", "Females", "Sex-stratified"))]
ggdt[, model_comparison := factor(model_comparison, levels=c("SCORE2 vs. SCORE2 + NMR scores", "SCORE2 vs. SCORE2 + PRSs", "SCORE2 vs. SCORE2 + NMR scores + PRSs"))]

g <- ggplot(ggdt) +
  aes(x=Estimate, xmin=Lower, xmax=Upper, y=fct_rev(model_comparison), color=type) +
  facet_grid(~ sex) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0, position=position_dodgev(height=0.3)) +
  geom_point(shape=23, fill="white", size=1.2, position=position_dodgev(height=0.3)) +
  scale_color_manual("", values=c("CVD cases"="#c51b7d", "Non-cases"="#4d9221")) +
  scale_x_continuous("Categorical NRI, % reclassified (95% CI)", labels=percent) +
  ylab("") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.title=element_blank(), legend.text=element_text(size=8, color="black"),
    legend.box.margin=margin(-0.5, 0, 0, -6, unit="cm"),
    legend.position="bottom", legend.justification="left",
    legend.box.background=element_blank(), legend.background=element_blank()
  )
ggsave(g, width=7.2, height=3.5, file="analyses/test/categorical_NRI.pdf", device=cairo_pdf)

# Create table of risk strata
risk <- dat[score_type == "non-derived"]
risk[, model := factor(model, levels=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
risk[, age_group := ifelse(age < 50, "40-<50 years of age", "50-<70 years of age")]
risk[, age_group := factor(age_group, levels=c("40-<50 years of age", "50-<70 years of age"))]
risk[, sex := factor(sex, levels=c("Male", "Female"))]
risk[, risk_group := fcase(
  age < 50 & uk_calibrated_risk < 0.025, "Low-to-moderate risk",
  age < 50 & uk_calibrated_risk >= 0.025 & uk_calibrated_risk < 0.075, "High risk",
  age < 50 & uk_calibrated_risk >= 0.075, "Very high risk",
  age >= 50 & uk_calibrated_risk < 0.05, "Low-to-moderate risk",
  age >= 50 & uk_calibrated_risk >= 0.05 & uk_calibrated_risk < 0.10, "High risk",
  age >= 50 & uk_calibrated_risk >= 0.10, "Very high risk"
)]
risk[, risk_group := factor(risk_group, levels=c("Low-to-moderate risk", "High risk", "Very high risk"))]
risk <- risk[, .N, by=.(sex, age_group, model, risk_group)]

group_totals <- risk[, .(total=sum(N)), by=.(sex, age_group, model)]
risk[group_totals, on = .(sex, age_group, model), pct := N/total]
risk[, text := sprintf("%s (%.2f%%)", format(N, big.mark=",", trim=TRUE), round(pct*100, digits=2))]

risk <- dcast(risk, age_group + risk_group ~ sex + model, value.var="text", fill="0 (0.00%)")
fwrite(risk, sep="\t", quote=FALSE, file="analyses/test/ESC_2021_risk_strata.txt")

# Create NRI table for supp
dt <- nri_estimates[,.(sex, model_comparison)]
dt <- unique(dt)
dt <- rbind(unique(dt[,.(sex, model_comparison="SCORE2")]), dt)
dt[, model := gsub("SCORE2 vs. ", "", model_comparison)]
dt[, sex := factor(sex, levels=c("Sex-stratified", "Males", "Females"))]
dt[, model := factor(model, levels=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
dt <- dt[order(model)][order(sex)]

# How many cases are at high risk?
score2_high_risk <- reclassified[Old == "< 0.1" & model_comparison == "SCORE2 vs. SCORE2 + PRSs", .(model_comparison="SCORE2", Cases=sum(Cases)), by=.(sex)]
new_high_risk <- reclassified[New == "< 0.1", .(Cases=sum(Cases)), by=.(sex, model_comparison)]
dt[score2_high_risk, on = .(sex, model_comparison), Cases_high_risk := i.Cases]
dt[new_high_risk, on = .(sex, model_comparison), Cases_high_risk := i.Cases]

# How many low-to-moderate risk are reclassified as high risk?
reclassified_high_risk <- reclassified[Old == "< 0.05" & New == "< 0.1", .(Cases=sum(Cases)), by=.(sex, model_comparison)]
dt[reclassified_high_risk, on = .(sex, model_comparison), Cases_high_risk_reclassified := i.Cases]

# How many cases are at very high risk?
score2_very_high_risk <- reclassified[Old == ">= 0.1" & model_comparison == "SCORE2 vs. SCORE2 + PRSs", .(model_comparison="SCORE2", Cases=sum(Cases)), by=.(sex)]
new_very_high_risk <- reclassified[New == ">= 0.1", .(Cases=sum(Cases)), by=.(sex, model_comparison)]
dt[score2_very_high_risk, on = .(sex, model_comparison), Cases_very_high_risk := i.Cases]
dt[new_very_high_risk, on = .(sex, model_comparison), Cases_very_high_risk := i.Cases]

# How many are reclassified as very high risk
reclassified_very_high_risk <- reclassified[Old != ">= 0.1" & New == ">= 0.1", .(Cases=sum(Cases)), by=.(sex, model_comparison)]
dt[reclassified_very_high_risk, on = .(sex, model_comparison), Cases_very_high_risk_reclassified := i.Cases]

# Add in case NRI
dt[nri_estimates[metric == "NRI+"], on = .(sex, model_comparison), c("Case_NRI", "Case_NRI_L95", "Case_NRI_U95") := .(Estimate, Lower, Upper)]

# Compute bootstrap P-value
case_NRI_pval <- rbindlist(idcol="sex", fill=TRUE, lapply(nri_lists, function(l1) {
  rbindlist(idcol="model_comparison", fill=TRUE, lapply(l1, function(l2) {
    nulls <- sum(l2$bootstrapsample[, "NRI+"] <= 0)
    n <- nrow(l2$bootstrapsample)
    data.table("pval"=(nulls + 1)/(n + 1))
  }))
}))
dt[case_NRI_pval, on = .(sex, model_comparison), Case_NRI_pval := i.pval]

# How many controls are at high risk?
score2_high_risk <- reclassified[Old == "< 0.1" & model_comparison == "SCORE2 vs. SCORE2 + PRSs", .(model_comparison="SCORE2", Controls=sum(All) - sum(Cases)), by=.(sex)]
new_high_risk <- reclassified[New == "< 0.1", .(Controls=sum(All) - sum(Cases)), by=.(sex, model_comparison)]
dt[score2_high_risk, on = .(sex, model_comparison), Controls_high_risk := i.Controls]
dt[new_high_risk, on = .(sex, model_comparison), Controls_high_risk := i.Controls]

# How many very high risk are reclassified as high risk?
reclassified_high_risk <- reclassified[Old == ">= 0.1" & New == "< 0.1", .(Controls=sum(All) - sum(Cases)), by=.(sex, model_comparison)]
dt[reclassified_high_risk, on = .(sex, model_comparison), Controls_high_risk_reclassified := i.Controls]

# How many controls are at low-to-moderate risk
score2_low_risk <- reclassified[Old == "< 0.05" & model_comparison == "SCORE2 vs. SCORE2 + PRSs", .(model_comparison="SCORE2", Controls=sum(All) - sum(Cases)), by=.(sex)]
new_low_risk <- reclassified[New == "< 0.05", .(Controls=sum(All) - sum(Cases)), by=.(sex, model_comparison)]
dt[score2_low_risk, on = .(sex, model_comparison), Controls_low_risk := i.Controls]
dt[new_low_risk, on = .(sex, model_comparison), Controls_low_risk := i.Controls]

# How many very low risk are reclassified as low risk?
reclassified_low_risk <- reclassified[Old != "< 0.05" & New == "< 0.05", .(Controls=sum(All) - sum(Cases)), by=.(sex, model_comparison)]
dt[reclassified_low_risk, on = .(sex, model_comparison), Controls_low_risk_reclassified := i.Controls]

# Add in control NRI
dt[nri_estimates[metric == "NRI-"], on = .(sex, model_comparison), c("Control_NRI", "Control_NRI_L95", "Control_NRI_U95") := .(Estimate, Lower, Upper)]

# Compute bootstrap P-value
case_NRI_pval <- rbindlist(idcol="sex", fill=TRUE, lapply(nri_lists, function(l1) {
  rbindlist(idcol="model_comparison", fill=TRUE, lapply(l1, function(l2) {
    nulls <- sum(l2$bootstrapsample[, "NRI-"] <= 0)
    n <- nrow(l2$bootstrapsample)
    data.table("pval"=(nulls + 1)/(n + 1))
  }))
}))
dt[case_NRI_pval, on = .(sex, model_comparison), Control_NRI_pval := i.pval]

# Write out
fwrite(dt, sep="\t", quote=FALSE, file="analyses/test/categorical_NRI_for_supp.txt")




















