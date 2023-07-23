# Takes 15-20 minutes to run on a compute node
library(data.table)
library(foreach)
library(nricens)
library(ggplot2)
library(forcats)
library(ggstance)
library(scales)

# Load required data
dat <- fread("analyses/risk_recalibration/CVD_linear_predictors_and_risk.txt")

# Wrapper function for categorical NRI test
nri.test <- function(data, base_model, new_model, absrisk_column="uk_risk") {
  # Extract predicted 10 year risk
  comp_risk <- data[model == base_model | model == new_model]
  comp_risk[, model := fcase(model == base_model, "base", model == new_model, "new")]
  comp_risk <- dcast(comp_risk, eid + age + incident_cvd_followup + incident_cvd ~ model, value.var=absrisk_column)

  # Split into age groups
  comp_risk_lt50 <- comp_risk[age < 50]
  comp_risk_lt70 <- comp_risk[age >= 50 & age < 70] 

  # Run NRI analysis
  lt50_NRI <- nricens(event = comp_risk_lt50$incident_cvd, time = comp_risk_lt50$incident_cvd_followup,
                      p.std = comp_risk_lt50$base, p.new = comp_risk_lt50$new,
                      updown = "category", cut = 0.075, t0 = 10, niter = 1000)

  lt70_NRI <- nricens(event = comp_risk_lt70$incident_cvd, time = comp_risk_lt70$incident_cvd_followup,
                      p.std = comp_risk_lt70$base, p.new = comp_risk_lt70$new,
                      updown = "category", cut = 0.1, t0 = 10, niter = 1000)
  system("rm -f Rplots.pdf")

  # Add in sample size and case numbers
  lt50_NRI$n <- comp_risk_lt50[,.N]
  lt50_NRI$nevent <- comp_risk_lt50[, sum(incident_cvd)]

  lt70_NRI$n <- comp_risk_lt70[,.N]
  lt70_NRI$nevent <- comp_risk_lt70[, sum(incident_cvd)]

  list("<50"=lt50_NRI, "50–<70"=lt70_NRI)
}

nri_lists <- list(
  "uk_risk"=list(
    "Males"=list(
      "SCORE2 vs. SCORE2 + NMR scores"=nri.test(dat[sex == "Male"], "SCORE2", "SCORE2 + NMR scores"),
      "SCORE2 vs. SCORE2 + PRSs"=nri.test(dat[sex == "Male"], "SCORE2", "SCORE2 + PRSs"),
      "SCORE2 vs. SCORE2 + NMR scores + PRSs"=nri.test(dat[sex == "Male"], "SCORE2", "SCORE2 + NMR scores + PRSs")
    ),
    "Females"=list(
      "SCORE2 vs. SCORE2 + NMR scores"=nri.test(dat[sex == "Female"], "SCORE2", "SCORE2 + NMR scores"),
      "SCORE2 vs. SCORE2 + PRSs"=nri.test(dat[sex == "Female"], "SCORE2", "SCORE2 + PRSs"),
      "SCORE2 vs. SCORE2 + NMR scores + PRSs"=nri.test(dat[sex == "Female"], "SCORE2", "SCORE2 + NMR scores + PRSs")
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
saveRDS(nri_lists, file="analyses/test/categorical_nri.rds")

# Extract tables of estimates
nri_estimates <- rbindlist(idcol="risk_col", fill=TRUE, lapply(nri_lists, function(l1) {
  rbindlist(idcol="sex", fill=TRUE, lapply(l1, function(l2) {
    rbindlist(idcol="model_comparison", fill=TRUE, lapply(l2, function(l3) {
      rbindlist(idcol="age_group", fill=TRUE, lapply(l3, function(l4) {
        cbind(samples=l4$n, cases=l4$nevent, as.data.table(l4$nri, keep.rownames="metric"))
      }))
    }))
  }))
}))

# Extract tables of reclassifications for categorical nris
reclassified_cases <- rbindlist(idcol="risk_col", fill=TRUE, lapply(nri_lists, function(l1) {
  rbindlist(idcol="sex", fill=TRUE, lapply(l1, function(l2) {
    rbindlist(idcol="model_comparison", fill=TRUE, lapply(l2, function(l3) {
      rbindlist(idcol="age_group", fill=TRUE, lapply(l3, function(l4) {
        as.data.table(l4[["rtab.case"]])
      }))
    }))
  }))
}))

reclassified <- rbindlist(idcol="risk_col", fill=TRUE, lapply(nri_lists, function(l1) {
  rbindlist(idcol="sex", fill=TRUE, lapply(l1, function(l2) {
    rbindlist(idcol="model_comparison", fill=TRUE, lapply(l2, function(l3) {
      rbindlist(idcol="age_group", fill=TRUE, lapply(l3, function(l4) {
        as.data.table(l4[["rtab"]])
      }))
    }))
  }))
}))

# Merge
setnames(reclassified, "N", "All")
reclassified[reclassified_cases, on = .(risk_col, sex, model_comparison, age_group, Standard, New), Cases := N]
setnames(reclassified, "Standard", "Old")

# Add in total sample size and total cases to reclassified table
reclassified[nri_estimates, on = .(risk_col, sex, model_comparison, age_group), Total_Samples := i.samples]
reclassified[nri_estimates, on = .(risk_col, sex, model_comparison, age_group), Total_Cases := i.cases]

fwrite(nri_estimates, sep="\t", quote=FALSE, file="analyses/test/categorical_nri_estimates.txt")
fwrite(reclassified, sep="\t", quote=FALSE, file="analyses/test/categorical_nri_reclassified.txt")

# Plot
ggdt <- nri_estimates[metric %in% c("NRI+", "NRI-")]
ggdt[, type := ifelse(metric == "NRI+", "CVD cases", "Non-cases")]
ggdt[, type := factor(type, levels=c("CVD cases", "Non-cases"))]
ggdt[, sex := factor(sex, levels=c("Males", "Females"))]
ggdt[, model_comparison := factor(model_comparison, levels=c("SCORE2 vs. SCORE2 + NMR scores", "SCORE2 vs. SCORE2 + PRSs", "SCORE2 vs. SCORE2 + NMR scores + PRSs"))]
ggdt[, age_group := factor(paste(age_group, "years"), levels=paste(c("<50", "50–<70"), "years"))]

g <- ggplot(ggdt[risk_col == "uk_risk"]) +
  aes(x=Estimate, xmin=Lower, xmax=Upper, y=fct_rev(model_comparison), color=type) +
  facet_grid(age_group ~ sex) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0, position=position_dodgev(height=0.3)) +
  geom_point(shape=23, fill="white", size=1.2, position=position_dodgev(height=0.3)) +
  geom_errorbarh(data=ggdt[risk_col == "uk_risk" & metric == "NRI-"], height=0, position=position_nudge(y=0.075), show.legend=FALSE) +
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

g <- ggplot(ggdt[risk_col == "recalibrated_risk"]) +
  aes(x=Estimate, xmin=Lower, xmax=Upper, y=fct_rev(model_comparison), color=type) +
  facet_grid(age_group ~ sex) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0, position=position_dodgev(height=0.3)) +
  geom_point(shape=23, fill="white", size=1.2, position=position_dodgev(height=0.3)) +
  geom_errorbarh(data=ggdt[risk_col == "recalibrated_risk" & metric == "NRI-"], height=0, position=position_nudge(y=0.075), show.legend=FALSE) +
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
ggsave(g, width=7.2, height=3.5, file="analyses/test/categorical_NRI_after_recalibration.pdf", device=cairo_pdf)


