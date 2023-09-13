library(data.table)
library(foreach)
library(survival)
library(ggplot2)
library(cowplot)

# Make output directory
system("mkdir -p analyses/test", wait=TRUE)

# Load required data
dat <- fread("data/cleaned/analysis_cohort.txt", select=c("eid", "sex", "age", "incident_cvd_followup", "incident_cvd",  "smoking", "sbp", "hdl", "tchol", "SCORE2_excl_UKB", "CAD_metaGRS", "Stroke_metaGRS"))
setnames(dat, "SCORE2_excl_UKB", "SCORE2")

# Add in test scores
test_scores <- fread("analyses/nmr_score_training/aggregate_test_non_derived_NMR_scores.txt")
setnames(test_scores, gsub("NMR", "NMR", names(test_scores)))
dat <- dat[test_scores, on = .(eid)]

# Code smoking status factor
dat[, smoking := fcase(
  is.na(smoking), "other",
  smoking == FALSE, "other",
  smoking == TRUE, "current"
)]
dat[, smoking := factor(smoking, levels=c("other", "current"))]

# Fit multivariable model with SCORE2 as an offset
offset_hrs <- foreach(this_sex = c("Sex-stratified", "Males", "Females"), .combine=rbind) %do% {
  if (this_sex == "Sex-stratified") {
    cx <- coxph(data=dat,
      Surv(incident_cvd_followup, incident_cvd) ~ strata(sex) + offset(SCORE2) +
        scale(CAD_NMR_score) + scale(Stroke_NMR_score) + scale(CAD_metaGRS) + scale(Stroke_metaGRS)
    )
  } else {
    if (this_sex == "Males") {
      this_dat <- dat[sex == "Male"]
    } else {
      this_dat <- dat[sex == "Female"]
    }
    cx <- coxph(data=this_dat,
      Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) +
        scale(CAD_NMR_score) + scale(Stroke_NMR_score) + scale(CAD_metaGRS) + scale(Stroke_metaGRS)
    )
  }

  this_hrs <- as.data.table(coef(summary(cx)), keep.rownames=TRUE)
  setnames(this_hrs, c("variable", "logHR", "HR", "SE", "Z", "pval"))
  this_hrs[, variable := c("CAD NMR score", "Stroke NMR score", "CAD PRS", "Stroke PRS")]
  this_hrs[, c("L95", "U95") :=  as.data.table(exp(confint(cx)))]
  this_hrs[, .(sex=this_sex, variable, HR, L95, U95, pval)]
}
fwrite(offset_hrs, sep="\t", quote=FALSE, file="analyses/test/score_HRs_independent_of_score2.txt")

# Fit multivariable models with SCORE2 component risk factors
multi_hrs <- foreach(this_sex = c("Sex-stratified", "Males", "Females"), .combine=rbind) %do% {
  if (this_sex == "Sex-stratified") {
    cx <- coxph(data=dat,
      Surv(incident_cvd_followup, incident_cvd) ~ strata(sex) + 
        scale(age)*smoking + scale(age)*scale(sbp) + scale(age)*scale(tchol) + scale(age)*scale(hdl) +
        scale(CAD_NMR_score) + scale(Stroke_NMR_score) + scale(CAD_metaGRS) + scale(Stroke_metaGRS)
    )
  } else {
    if (this_sex == "Males") {
      this_dat <- dat[sex == "Male"]
    } else {
      this_dat <- dat[sex == "Female"]
    }
    cx <- coxph(data=this_dat,
      Surv(incident_cvd_followup, incident_cvd) ~
        scale(age)*smoking + scale(age)*scale(sbp) + scale(age)*scale(tchol) + scale(age)*scale(hdl) +
        scale(CAD_NMR_score) + scale(Stroke_NMR_score) + scale(CAD_metaGRS) + scale(Stroke_metaGRS)
    )
  }

  this_hrs <- as.data.table(coef(summary(cx)), keep.rownames=TRUE)
  setnames(this_hrs, c("variable", "logHR", "HR", "SE", "Z", "pval"))
  this_hrs[, variable := c(
    "Age", "Smoker", "SBP", "Total cholesterol", "HDL cholesterol",
    "CAD NMR score", "Stroke NMR score", "CAD PRS", "Stroke PRS",
    "Age x smoker", "Age x SBP", "Age x total cholesterol", "Age x HDL cholesterol"
  )]
  this_hrs[, c("L95", "U95") :=  as.data.table(exp(confint(cx)))]
  this_hrs[, .(sex=this_sex, variable, HR, L95, U95, pval)]
}
fwrite(multi_hrs, sep="\t", quote=FALSE, file="analyses/test/score_HRs_independent_of_score2_risk_factors.txt")

# Create sex-stratified plot
ggdt <- rbind(idcol="model", "offset"=offset_hrs, "risk_factors"=multi_hrs)
ggdt <- ggdt[sex == "Sex-stratified"]
ggdt[, model := factor(model, levels=c("offset", "risk_factors"))]
ggdt <- ggdt[order(pval)]
ggdt[, xorder := factor(1:.N)]

g <- ggplot(ggdt) +
  aes(x=xorder, y=HR, ymin=L95, ymax=U95) +
  facet_grid(~ model, space="free_x", scales="free_x") +
  geom_hline(yintercept=1, linetype=2) +
  geom_errorbar(width=0) +
  geom_point(shape=23, fill="white") +
  ylab("HR per SD (95% CI)") +
  scale_x_discrete(labels=structure(ggdt$variable, names=as.character(ggdt$xorder))) +
  theme_bw() + 
  theme(
    axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5), axis.text.y=element_text(size=6),
    axis.title.x=element_blank(), axis.title.y=element_text(size=8),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
    strip.background=element_blank(), strip.text=element_blank()
  )
ggsave(g, width=7, height=3, file="analyses/test/sex_stratified_HRs.pdf")

# Construct table for supplement
dt <- ggdt[order(model), .(model, variable, HR, L95, U95, pval)]
fwrite(dt, sep="\t", quote=FALSE, file="analyses/test/sex_stratified_HRs.txt")

# Create sex-specific plots
ggdt <- rbind(idcol="model", "offset"=offset_hrs, "risk_factors"=multi_hrs)
ggdt[, model := factor(model, levels=c("offset", "risk_factors"))]
ggdt[, sex := factor(sex, levels=c("Sex-stratified", "Males", "Females"))]
ggdt <- ggdt[order(pval)]
ggdt[sex == "Sex-stratified", xorder := factor(1:.N)]
ggdt[ggdt[sex == "Sex-stratified"], on = .(model, variable), xorder := i.xorder]

g <- ggplot(ggdt) +
  aes(x=xorder, y=HR, ymin=L95, ymax=U95, color=sex) +
  facet_grid(~ model, space="free_x", scales="free_x") +
  geom_hline(yintercept=1, linetype=2) +
  geom_errorbar(width=0, position=position_dodge(width=0.6)) +
  geom_point(shape=23, fill="white", position=position_dodge(width=0.6)) +
  ylab("HR per SD (95% CI)") +
  scale_x_discrete(labels=unique(structure(ggdt$variable, names=as.character(ggdt$xorder)))) +
  scale_color_manual("Sex", values=c("Males"="#e41a1c", "Females"="#377eb8", "Sex-stratified"="#006d2c")) +
  theme_bw() + 
  theme(
    axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5), axis.text.y=element_text(size=6),
    axis.title.x=element_blank(), axis.title.y=element_text(size=8),
    panel.grid.minor.x=element_blank(),
    strip.background=element_blank(), strip.text=element_blank(),
    legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=8)
  )

ggsave(g, width=7, height=5, file="analyses/test/sex_specific_HRs.pdf")

# Construct table for supplement
dt <- dcast(ggdt, model + variable  ~ sex, value.var=c("HR", "L95", "U95", "pval"))
dt <- dt[, .(model, variable, HR_Males, L95_Males, U95_Males, pval_Males, HR_Females, L95_Females, U95_Females, pval_Females)]
dt <- dt[ggdt[sex == "Sex-stratified", .(model, variable)], on = .(model, variable)][order(model)]
fwrite(dt, sep="\t", quote=FALSE, file="analyses/test/sex_specific_HRs.txt")

