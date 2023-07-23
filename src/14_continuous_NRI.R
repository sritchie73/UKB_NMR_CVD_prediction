# Needs to be run on compute nodes, takes at least 1:0:0
library(data.table)
library(foreach)
library(nricens)
library(ggplot2)
library(forcats)
library(ggstance)
library(scales)

# Load required data
dat <- fread('analyses/risk_recalibration/CVD_linear_predictors_and_risk.txt')

# Wrapper function for continuous NRI test
nri.test <- function(data, base_model, new_model, absrisk_column) {
  # Extract predicted 10 year risk
  comp_risk <- data[model == base_model | model == new_model]
  comp_risk[, model := fcase(model == base_model, "base", model == new_model, "new")]
  comp_risk <- dcast(comp_risk, eid + incident_cvd_followup + incident_cvd ~ model, value.var=absrisk_column)

  # Run NRI analysis
  contNRI <- nricens(event = comp_risk$incident_cvd, time = comp_risk$incident_cvd_followup,
                     p.std = comp_risk$base, p.new = comp_risk$new,
                     updown = "diff", cut = 0, t0 = 10, niter = 1000)
  system("rm -f Rplots.pdf")

  # Add in sample size and case numbers
  contNRI$n <- comp_risk[,.N]
  contNRI$nevent <- comp_risk[, sum(incident_cvd)]
 
  return(contNRI)
}

nri_lists <- list(
  "uk_risk"=list(
    "Males"=list(
      "SCORE2 vs. SCORE2 + NMR scores"=nri.test(dat[sex == "Male"], "SCORE2", "SCORE2 + NMR scores", "uk_risk"),
      "SCORE2 vs. SCORE2 + PRSs"=nri.test(dat[sex == "Male"], "SCORE2", "SCORE2 + PRSs", "uk_risk"),
      "SCORE2 vs. SCORE2 + NMR scores + PRSs"=nri.test(dat[sex == "Male"], "SCORE2", "SCORE2 + NMR scores + PRSs", "uk_risk")
    ),
    "Females"=list(
      "SCORE2 vs. SCORE2 + NMR scores"=nri.test(dat[sex == "Female"], "SCORE2", "SCORE2 + NMR scores", "uk_risk"),
      "SCORE2 vs. SCORE2 + PRSs"=nri.test(dat[sex == "Female"], "SCORE2", "SCORE2 + PRSs", "uk_risk"),
      "SCORE2 vs. SCORE2 + NMR scores + PRSs"=nri.test(dat[sex == "Female"], "SCORE2", "SCORE2 + NMR scores + PRSs", "uk_risk")
    )
  ),
  "recalibrated_risk"=list(
    "Males"=list(
      "SCORE2 vs. SCORE2 + NMR scores"=nri.test(dat[sex == "Male"], "SCORE2", "SCORE2 + NMR scores", "recalibrated_risk"),
      "SCORE2 vs. SCORE2 + PRSs"=nri.test(dat[sex == "Male"], "SCORE2", "SCORE2 + PRSs", "recalibrated_risk"),
      "SCORE2 vs. SCORE2 + NMR scores + PRSs"=nri.test(dat[sex == "Male"], "SCORE2", "SCORE2 + NMR scores + PRSs", "recalibrated_risk")
    ),
    "Females"=list(
      "SCORE2 vs. SCORE2 + NMR scores"=nri.test(dat[sex == "Female"], "SCORE2", "SCORE2 + NMR scores", "recalibrated_risk"),
      "SCORE2 vs. SCORE2 + PRSs"=nri.test(dat[sex == "Female"], "SCORE2", "SCORE2 + PRSs", "recalibrated_risk"),
      "SCORE2 vs. SCORE2 + NMR scores + PRSs"=nri.test(dat[sex == "Female"], "SCORE2", "SCORE2 + NMR scores + PRSs", "recalibrated_risk")
    )
  )
)
saveRDS(nri_lists, file="analyses/test/nri.rds")

# Extract tables of estimates
nri_estimates <- rbindlist(idcol="risk_col", fill=TRUE, lapply(nri_lists, function(l1) {
  rbindlist(idcol="sex", fill=TRUE, lapply(l1, function(l2) {
    rbindlist(idcol="model_comparison", fill=TRUE, lapply(l2, function(l3) {
      cbind(samples=l3$n, cases=l3$nevent, as.data.table(l3$nri, keep.rownames="metric"))
    }))
  }))
}))
fwrite(nri_estimates, sep="\t", quote=FALSE, file="analyses/test/nri_estimates.txt")

# Plot
ggdt <- nri_estimates[risk_col == "uk_risk" & metric %in% c("NRI+", "NRI-")]
ggdt[, type := ifelse(metric == "NRI+", "CVD cases", "Non-cases")]
ggdt[, type := factor(type, levels=c("CVD cases", "Non-cases"))]
ggdt[, sex := factor(sex, levels=c("Males", "Females"))]
ggdt[, model_comparison := factor(model_comparison, levels=c("SCORE2 vs. SCORE2 + NMR scores", "SCORE2 vs. SCORE2 + PRSs", "SCORE2 vs. SCORE2 + NMR scores + PRSs"))]

g <- ggplot(ggdt) +
  aes(x=Estimate, xmin=Lower, xmax=Upper, y=fct_rev(model_comparison), color=type) + 
  facet_wrap(~ sex) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0, position=position_dodgev(height=0.3)) +
  geom_point(shape=23, fill="white", size=1.2, position=position_dodgev(height=0.3)) +
  scale_color_manual("", values=c("CVD cases"="#c51b7d", "Non-cases"="#4d9221")) +
  scale_x_continuous("Continuous NRI", labels=percent) +
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
ggsave(g, width=7.2, height=2, file="analyses/test/continuous_NRI.pdf")

