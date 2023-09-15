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
stop()


#####################
# Build main figure
#####################

# Compute number and percentage of people at different risk thresholds
dat[, risk_group := fcase(
  uk_calibrated_risk >= 0.25, "≥25%",
  uk_calibrated_risk > 0.15, "15%–<25%",
  uk_calibrated_risk > 0.1, "10%–<15%",
  uk_calibrated_risk > 0.075, "7.5%–<10%",
  uk_calibrated_risk > 0.05, "5%–<7.5%",
  uk_calibrated_risk > 0.025, "2.5%–<5%",
  default = "<2.5%"
)]
dat[, age_group := ifelse(age < 50, "40–<50 years", "50–<70 years")]

risk_strata <- dat[, .N, by=.(sex, model, age_group, risk_group)]

# Set up labels for plotting
risk_strata[, risk_group := factor(risk_group, levels=c("<2.5%", "2.5%–<5%", "5%–<7.5%", "7.5%–<10%", "10%–<15%", "15%–<25%", "≥25%"))]
risk_strata[, sex := factor(paste0(sex, "s"), levels=c("Males", "Females"))]
risk_strata[, age_group := factor(age_group, levels=c("40–<50 years", "50–<70 years"))]
risk_strata[, model := factor(model, levels=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]

# Build risk strata plot
g1 <- ggplot(risk_strata) +
  aes(x=N, y=fct_rev(model), fill=risk_group) +
  facet_grid(age_group ~ sex) +
  geom_col(position="fill") +
  geom_vline(xintercept=seq(0,1,by=0.1), color="white", alpha=0.6, linewidth=0.4) +
  scale_x_continuous("UK Biobank participants", labels=percent, expand=c(0,0.01)) +
  scale_y_discrete(expand=c(0,0.4)) +
  scale_fill_manual("Predicted absolute 10-year CVD risk", values=c(
    "<2.5%"="#34b340", "2.5%–<5%"="#71ff50", "5%–<7.5%"="#efff00", "7.5%–<10%"="#eacc00",
    "10%–<15%"="#e41f1b", "15%–<25%"="#b01815", "≥25%"="#680e0c"
  )) +
  guides(fill=guide_legend(title.position="top", nrow=1, byrow=TRUE, reverse=TRUE)) +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text.x=element_text(size=8, face="bold"), strip.text.y=element_text(size=8, face="bold"),
    strip.background=element_blank(), panel.grid=element_blank(),
    legend.title=element_text(size=8), legend.text=element_text(size=6),
    legend.box.margin=margin(-0.5, 0, 0, -3, unit="cm"),
    legend.position="bottom", legend.justification="left",
    legend.box.background=element_blank(), legend.background=element_blank()
  )

# Build NRI plot
nri_ggdt <- ggdt[risk_col == "uk_calibrated_risk" & sex != "Sex-stratified"]
nri_ggdt[, age_group := fct_recode(age_group, "40–<50 years"="<50 years")]
g2 <- ggplot(nri_ggdt) +
  aes(x=Estimate, xmin=Lower, xmax=Upper, y=fct_rev(model_comparison), color=type) +
  facet_grid(age_group ~ sex) +
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

# Combine both and render
g <- plot_grid(g1, g2, nrow=2, labels=c("A", "B"), label_size=10, rel_heights=c(0.53, 0.47))
ggsave(g, width=7.2, height=5.5, file="analyses/test/risk_strata_and_NRI.pdf", device=cairo_pdf)

# output risk strata numbers for supp table
dt1 <- copy(risk_strata)
dt1[, sex := gsub("s", "", sex)]
dt1[dat[,.N,by=.(sex, age_group, model)], on = .(sex, age_group, model), pct := N/i.N]
dt1[, text := sprintf("%s (%.2f%%)", format(N, big.mark=",", trim=TRUE), round(pct*100, digits=2))]
dt1 <- dcast(dt1, age_group + risk_group ~ sex + model, value.var="text", fill="0 (0.00%)")
fwrite(dt1, sep="\t", quote=FALSE, file="analyses/test/risk_strata_for_supp.txt")

# Collate and output NRI details for supp
dt <- nri_estimates[risk_col == "uk_calibrated_risk" & sex != "Sex-stratified"]

dt2 <- reclassified[risk_col == "uk_calibrated_risk" & sex != "Sex-stratified" & Old != New]
dt2[, type := ifelse(Old %like% "<", "Cases", "Non-cases")]
dt2[type == "Cases", Reclassified := Cases]
dt2[type == "Non-cases", Reclassified := All - Cases]
dt2 <- dt2[, .(sex, model_comparison, age_group, samples=Total_Samples, cases=Total_Cases, metric=ifelse(type == "Cases", "Reclassified+", "Reclassified-"), Estimate=Reclassified)]

dt3 <- reclassified[risk_col == "uk_calibrated_risk" & sex != "Sex-stratified"]
dt3 <- dt3[, .(All=sum(All), Cases=sum(Cases)), by=.(sex, model_comparison, age_group, samples=Total_Samples, cases=Total_Cases, New)]
dt3[, type := ifelse(New %like% "<", "Non-cases", "Cases")]
dt3[type == "Cases", Stratified := Cases]
dt3[type == "Non-cases", Stratified := All - Cases]
dt3 <- dt3[, .(sex, model_comparison, age_group, samples, cases, metric=ifelse(type == "Cases", "Stratified+", "Stratified-"), Estimate=Stratified)]

dt4 <- reclassified[risk_col == "uk_calibrated_risk" & sex != "Sex-stratified"]
dt4 <- dt4[, .(All=sum(All), Cases=sum(Cases)), by=.(sex, model_comparison, age_group, samples=Total_Samples, cases=Total_Cases, Old)]
dt4[, type := ifelse(Old %like% "<", "Non-cases", "Cases")]
dt4[type == "Cases", SCORE2 := Cases]
dt4[type == "Non-cases", SCORE2 := All - Cases]
dt4 <- dt4[, .(sex, model_comparison, age_group, samples, cases, metric=ifelse(type == "Cases", "SCORE2+", "SCORE2-"), Estimate=SCORE2)]

dt <- rbind(dt, dt2, dt3, dt4, fill=TRUE)

dt[, age_group := fct_recode(age_group, "40–<50"="<50")]
dt[, model := gsub(".* vs. ", "", model_comparison)] 
dt[, model := factor(model, levels=c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
dt[, sex := factor(sex, levels=c("Males", "Females"))]

dt <- dcast(dt, sex + age_group + samples + cases + model ~ metric, value.var=c("Estimate", "Lower", "Upper"))
dt <- dt[,.(sex, age_group, samples, cases, model,
  `Estimate_SCORE2+`, `Estimate_Stratified+`, `Estimate_Reclassified+`, `Estimate_NRI+`, `Lower_NRI+`, `Upper_NRI+`,
  `Estimate_SCORE2-`, `Estimate_Stratified-`, `Estimate_Reclassified-`, `Estimate_NRI-`, `Lower_NRI-`, `Upper_NRI-`
)]

fwrite(dt, sep="\t", quote=FALSE, file="analyses/test/collated_categorical_NRI_for_supp.txt")







