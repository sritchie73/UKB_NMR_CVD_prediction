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
