library(data.table)
library(foreach)
library(caret)
library(survival)
library(ggplot2)
source("src/utils/cv_coxph.R")
source("src/utils/cox_test.R")
source("src/utils/SCORE2.R")

# Load required data
dat <- fread("data/cleaned/analysis_cohort.txt", select=c("eid", "sex", "incident_cvd_followup", "incident_cvd", "CAD_metaGRS", "Stroke_metaGRS", "SCORE2"))
nmr_scores <- fread("analyses/nmr_score_training/aggregate_test_NMR_scores.txt")
dat <- dat[nmr_scores, on = .(eid)]

# Standardise PRSs
prs_scaling <- fread("data/standardised/prs_scaling_factors.txt")
dat[, CAD_metaGRS := (CAD_metaGRS - prs_scaling[PRS == "CAD_metaGRS", mean])/prs_scaling[PRS == "CAD_metaGRS", sd]]
dat[, Stroke_metaGRS := (Stroke_metaGRS - prs_scaling[PRS == "Stroke_metaGRS", mean])/prs_scaling[PRS == "Stroke_metaGRS", sd]]
 
# Fit Cox proportial hazards models in 10-fold cross-validation
dat[, coxph_cv_foldid := createFolds(paste(incident_cvd, sex), k=10, list=FALSE)]
fwrite(dat[,.(eid, coxph_cv_foldid)], sep="\t", quote=FALSE, file="analyses/CVD_score_weighting/cross_validation_fold_allocation.txt")

cv_cx_lists <- list(
  "Male"=list(
    "SCORE2 + NMR_scores"=cv.coxph("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score", dat[sex == "Male"], "coxph_cv_foldid"),
    "SCORE2 + PRS"=cv.coxph("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_metaGRS + Stroke_metaGRS", dat[sex == "Male"], "coxph_cv_foldid"),
    "SCORE2 + NMR_scores + PRS"=cv.coxph("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score + CAD_metaGRS + Stroke_metaGRS", dat[sex == "Male"], "coxph_cv_foldid")
  ),
  "Female"=list(
    "SCORE2 + NMR_scores"=cv.coxph("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score", dat[sex == "Female"], "coxph_cv_foldid"),
    "SCORE2 + PRS"=cv.coxph("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_metaGRS + Stroke_metaGRS", dat[sex == "Female"], "coxph_cv_foldid"),
    "SCORE2 + NMR_scores + PRS"=cv.coxph("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score + CAD_metaGRS + Stroke_metaGRS", dat[sex == "Female"], "coxph_cv_foldid")
  )
)

# Extract linear predictors predicted in each test fold for downstream analyses
pred_scores <- foreach(this_sex = c("Male", "Female"), .combine=rbind) %:% 
  foreach(this_model = c("SCORE2", "SCORE2 + NMR_scores", "SCORE2 + PRS", "SCORE2 + NMR_scores + PRS"), .combine=rbind) %do% {
    if (this_model == "SCORE2") {
      dat[sex == this_sex, .(model=this_model, eid, sex, linear_predictor=SCORE2)]
    } else {
      this_dat <- dat[sex == this_sex]
      this_formula <- fcase(
        this_model == "SCORE2 + NMR_scores", "Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score",
        this_model == "SCORE2 + PRS", "Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_metaGRS + Stroke_metaGRS",
        this_model == "SCORE2 + NMR_scores + PRS", "Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score + CAD_metaGRS + Stroke_metaGRS"
      )
      this_lp <- cv.linear.predictor(this_formula, this_dat, "coxph_cv_foldid", cv_cx_lists[[this_sex]][[this_model]])
      this_dat[as.integer(names(this_lp)), .(model=this_model, eid, sex, age, incident_cvd_followup, incident_cvd, linear_predictor=this_lp)] 
   }
}

# Compute absolute risk using SCORE2 baseline hazards
pred_scores[, absrisk := score2_absrisk(sex, linear_predictor), by=model]

# Compute absolute risk recalibrated to low-risk european region (including UK)
pred_scores[, uk_risk := score2_recalibration(sex, absrisk, "low"), by=model]

# Write out
fwrite(pred_scores, sep="\t", quote=FALSE, file="analyses/CVD_score_weighting/CVD_linear_predictors_and_risk.txt")

# Extract hazard ratios
cv_hrs <- rbindlist(idcol="sex", lapply(cv_cx_lists, function(l) {
  rbindlist(idcol="model", lapply(l, function(cx_list) {
    rbindlist(idcol="coxph_cv_foldid", lapply(cx_list, function(cx) {
      cf <- coef(summary(cx))
      ci <- confint(cx)
      data.table(score=rownames(cf), logHR=cf[,1], SE=cf[,3], L95=ci[,1], U95=ci[,2], Pval=cf[,5])
    }))
  }))
}))

# Write out set of hazard ratios from cross-validation
fwrite(cv_hrs, sep="\t", quote=FALSE, file="analyses/CVD_score_weighting/cross_validation_hazard_ratios.txt")

# Compute mean hazard ratio to use as score weight
weights <- cv_hrs[,.(weight=mean(logHR)), by=.(model, sex, score)]
weights <- dcast(weights, model + score~ sex, value.var="weight")

# Write out weights and transformations
fwrite(weights, sep="\t", quote=FALSE, file="analyses/CVD_score_weighting/score_weights.txt")
fwrite(prs_scaling, sep="\t", quote=FALSE, file="analyses/CVD_score_weighting/prs_scaling_factors.txt")

# Compute hazard ratios in full data
cx_lists <- list(
  "Male"=list(
    "SCORE2 + NMR_scores"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score", "incident_cvd", dat[sex == "Male"]),
    "SCORE2 + PRS"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat[sex == "Male"]),
    "SCORE2 + NMR_scores + PRS"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat[sex == "Male"])
  ), 
  "Female"=list(
    "SCORE2 + NMR_scores"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score", "incident_cvd", dat[sex == "Female"]),
    "SCORE2 + PRS"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat[sex == "Female"]),
    "SCORE2 + NMR_scores + PRS"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat[sex == "Female"])
  )
)

# Extract hazard ratios
full_hrs <- rbindlist(idcol="sex", lapply(cx_lists, function(l) {
  rbindlist(idcol="model", lapply(l, function(cx) {
    cx$coefficients[,.(score=coefficient, logHR, SE, L95=log(L95), U95=log(U95), Pval=Pvalue)]
  }))
}))

# Plot HRs for weight selection
weights_with_ci <- cv_hrs[,.(logHR=mean(logHR), L95=mean(L95), U95=mean(U95)), by=.(model, sex, score)]
ggdt <- rbind(idcol="hr_type", "10-fold cross-validation"=cv_hrs, "full dataset"=full_hrs, "cross-validation mean"=weights_with_ci, fill=TRUE)
ggdt[, score := gsub("_", " ", score)]
ggdt[, score := factor(score, levels=rev(c("CAD NMR score", "Stroke NMR score", "CAD metaGRS", "Stroke metaGRS")))]
ggdt[, model := gsub("_", " ", model)]
ggdt[, model := factor(model, levels=c("SCORE2 + NMR scores", "SCORE2 + PRS", "SCORE2 + NMR scores + PRS"))]
ggdt[, sex := factor(sex, levels=c("Male", "Female"))]

jitterer <- position_jitter(height=0.25, seed=123)
g <- ggplot(ggdt) +
  aes(x=logHR, xmin=L95, xmax=U95, y=score) +
  facet_grid(model ~ sex, scales="free", space="free_y") + 
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(data=ggdt[hr_type == "10-fold cross-validation"], position=jitterer, height=0, aes(color=hr_type, alpha=hr_type)) +
  geom_point(data=ggdt[hr_type == "10-fold cross-validation"], position=jitterer, aes(color=hr_type, shape=hr_type, size=hr_type, fill=hr_type, alpha=hr_type)) +
  geom_errorbarh(data=ggdt[hr_type == "full dataset"], height=0, aes(color=hr_type, alpha=hr_type), position=position_nudge(y=-0.3)) +
  geom_point(data=ggdt[hr_type == "full dataset"], aes(color=hr_type, shape=hr_type, size=hr_type, fill=hr_type, alpha=hr_type), position=position_nudge(y=-0.3)) +
  geom_errorbarh(data=ggdt[hr_type == "cross-validation mean"], height=0, aes(color=hr_type, alpha=hr_type)) +
  geom_point(data=ggdt[hr_type == "cross-validation mean"], aes(color=hr_type, shape=hr_type, size=hr_type, fill=hr_type, alpha=hr_type)) +
  scale_color_manual(name="Source", values=c("10-fold cross-validation"="#525252", "full dataset"="#810f7c", "cross-validation mean"="black")) +
  scale_size_manual(name="Source", values=c("10-fold cross-validation"=1, "full dataset"=1, "cross-validation mean"=1.2)) + 
  scale_shape_manual(name="Source", values=c("10-fold cross-validation"=19, "full dataset"=19, "cross-validation mean"=23)) +
  scale_fill_manual(name="Source", values=c("10-fold cross-validation"="#525252", "full dataset"="#810f7c", "cross-validation mean"="#41ab5d")) +
  scale_alpha_manual(name="Source", values=c("10-fold cross-validation"=0.6, "full dataset"=1, "cross-validation mean"=1)) +
  ylab("") + xlab("Log Hazard Ratio (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="bottom", legend.justification="left", legend.box.margin=margin(0, 0, 0, -2.5, unit="cm"),
    legend.title=element_text(size=8), legend.text=element_text(size=6)
  ) + 
  guides(
    color=guide_legend(title.position="top"), size=guide_legend(title.position="top"),
    shape=guide_legend(title.position="top"), fill=guide_legend(title.position="top"), 
    alpha=guide_legend(title.position="top")
  )
ggsave(g, width=4, height=6.4, file="analyses/CVD_score_weighting/score_weights_in_cross_validation.pdf")


# For interest, also quantify HRs on a consistent 1 SD scale (not to be used for weights)
dat[, SCORE2 := scale(SCORE2)]
dat[, CAD_NMR_score := scale(CAD_NMR_score)]
dat[, Stroke_NMR_score := scale(Stroke_NMR_score)]

cv_cx_lists2 <- list(
  "Male"=list(
    "SCORE2 + NMR_scores"=cv.coxph("Surv(incident_cvd_followup, incident_cvd) ~ SCORE2 + CAD_NMR_score + Stroke_NMR_score", dat[sex == "Male"], "coxph_cv_foldid"),
    "SCORE2 + PRS"=cv.coxph("Surv(incident_cvd_followup, incident_cvd) ~ SCORE2 + CAD_metaGRS + Stroke_metaGRS", dat[sex == "Male"], "coxph_cv_foldid"),
    "SCORE2 + NMR_scores + PRS"=cv.coxph("Surv(incident_cvd_followup, incident_cvd) ~ SCORE2 + CAD_NMR_score + Stroke_NMR_score + CAD_metaGRS + Stroke_metaGRS", dat[sex == "Male"], "coxph_cv_foldid")
  ),
  "Female"=list(
    "SCORE2 + NMR_scores"=cv.coxph("Surv(incident_cvd_followup, incident_cvd) ~ SCORE2 + CAD_NMR_score + Stroke_NMR_score", dat[sex == "Female"], "coxph_cv_foldid"),
    "SCORE2 + PRS"=cv.coxph("Surv(incident_cvd_followup, incident_cvd) ~ SCORE2 + CAD_metaGRS + Stroke_metaGRS", dat[sex == "Female"], "coxph_cv_foldid"),
    "SCORE2 + NMR_scores + PRS"=cv.coxph("Surv(incident_cvd_followup, incident_cvd) ~ SCORE2 + CAD_NMR_score + Stroke_NMR_score + CAD_metaGRS + Stroke_metaGRS", dat[sex == "Female"], "coxph_cv_foldid")
  )
)

# Extract hazard ratios
cv_hrs2 <- rbindlist(idcol="sex", lapply(cv_cx_lists2, function(l) {
  rbindlist(idcol="model", lapply(l, function(cx_list) {
    rbindlist(idcol="coxph_cv_foldid", lapply(cx_list, function(cx) {
      cf <- coef(summary(cx))
      ci <- confint(cx)
      data.table(score=rownames(cf), logHR=cf[,1], SE=cf[,3], L95=ci[,1], U95=ci[,2], Pval=cf[,5])
    }))
  }))
}))


# Compute hazard ratios in full data
cx_lists2 <- list(
  "Male"=list(
    "SCORE2 + NMR_scores"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ SCORE2 + CAD_NMR_score + Stroke_NMR_score", "incident_cvd", dat[sex == "Male"]),
    "SCORE2 + PRS"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ SCORE2 + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat[sex == "Male"]),
    "SCORE2 + NMR_scores + PRS"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ SCORE2 + CAD_NMR_score + Stroke_NMR_score + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat[sex == "Male"])
  ), 
  "Female"=list(
    "SCORE2 + NMR_scores"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ SCORE2 + CAD_NMR_score + Stroke_NMR_score", "incident_cvd", dat[sex == "Female"]),
    "SCORE2 + PRS"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ SCORE2 + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat[sex == "Female"]),
    "SCORE2 + NMR_scores + PRS"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ SCORE2 + CAD_NMR_score + Stroke_NMR_score + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat[sex == "Female"])
  )
)

# Extract hazard ratios
full_hrs2 <- rbindlist(idcol="sex", lapply(cx_lists2, function(l) {
  rbindlist(idcol="model", lapply(l, function(cx) {
    cx$coefficients[,.(score=coefficient, logHR, SE, L95=log(L95), U95=log(U95), Pval=Pvalue)]
  }))
}))

# Plot HRs for weight selection
weights_with_ci2 <- cv_hrs2[,.(logHR=mean(logHR), L95=mean(L95), U95=mean(U95)), by=.(model, sex, score)]
ggdt <- rbind(idcol="hr_type", "10-fold cross-validation"=cv_hrs2, "full dataset"=full_hrs2, "cross-validation mean"=weights_with_ci2, fill=TRUE)
ggdt[, score := gsub("_", " ", score)]
ggdt[, score := factor(score, levels=rev(c("SCORE2", "CAD NMR score", "Stroke NMR score", "CAD metaGRS", "Stroke metaGRS")))]
ggdt[, model := gsub("_", " ", model)]
ggdt[, model := factor(model, levels=c("SCORE2 + NMR scores", "SCORE2 + PRS", "SCORE2 + NMR scores + PRS"))]
ggdt[, sex := factor(sex, levels=c("Male", "Female"))]

g <- ggplot(ggdt) +
  aes(x=exp(logHR), xmin=exp(L95), xmax=exp(U95), y=score) +
  facet_grid(model ~ sex, scales="free", space="free_y") + 
  geom_vline(xintercept=1, linetype=2) +
  geom_errorbarh(data=ggdt[hr_type == "10-fold cross-validation"], position=jitterer, height=0, aes(color=hr_type, alpha=hr_type)) +
  geom_point(data=ggdt[hr_type == "10-fold cross-validation"], position=jitterer, aes(color=hr_type, shape=hr_type, size=hr_type, fill=hr_type, alpha=hr_type)) +
  geom_errorbarh(data=ggdt[hr_type == "full dataset"], height=0, aes(color=hr_type, alpha=hr_type), position=position_nudge(y=-0.3)) +
  geom_point(data=ggdt[hr_type == "full dataset"], aes(color=hr_type, shape=hr_type, size=hr_type, fill=hr_type, alpha=hr_type), position=position_nudge(y=-0.3)) +
  geom_errorbarh(data=ggdt[hr_type == "cross-validation mean"], height=0, aes(color=hr_type, alpha=hr_type)) +
  geom_point(data=ggdt[hr_type == "cross-validation mean"], aes(color=hr_type, shape=hr_type, size=hr_type, fill=hr_type, alpha=hr_type)) +
  scale_color_manual(name="Source", values=c("10-fold cross-validation"="#525252", "full dataset"="#810f7c", "cross-validation mean"="black")) +
  scale_size_manual(name="Source", values=c("10-fold cross-validation"=1, "full dataset"=1, "cross-validation mean"=1.2)) + 
  scale_shape_manual(name="Source", values=c("10-fold cross-validation"=19, "full dataset"=19, "cross-validation mean"=23)) +
  scale_fill_manual(name="Source", values=c("10-fold cross-validation"="#525252", "full dataset"="#810f7c", "cross-validation mean"="#41ab5d")) +
  scale_alpha_manual(name="Source", values=c("10-fold cross-validation"=0.6, "full dataset"=1, "cross-validation mean"=1)) +
  ylab("") + xlab("Hazard Ratio (95% CI) per SD") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="bottom", legend.justification="left", legend.box.margin=margin(0, 0, 0, -2.5, unit="cm"),
    legend.title=element_text(size=8), legend.text=element_text(size=6)
  ) + 
  guides(
    color=guide_legend(title.position="top"), size=guide_legend(title.position="top"),
    shape=guide_legend(title.position="top"), fill=guide_legend(title.position="top"), 
    alpha=guide_legend(title.position="top")
  )
ggsave(g, width=4, height=6.4, file="analyses/CVD_score_weighting/HRs_in_cross_validation.pdf")

# repeat, omitting SCORE2
ggdt <- ggdt[score != "SCORE2"]
g <- ggplot(ggdt) +
  aes(x=exp(logHR), xmin=exp(L95), xmax=exp(U95), y=score) +
  facet_grid(model ~ sex, scales="free", space="free_y") + 
  geom_vline(xintercept=1, linetype=2) +
  geom_errorbarh(data=ggdt[hr_type == "10-fold cross-validation"], position=jitterer, height=0, aes(color=hr_type, alpha=hr_type)) +
  geom_point(data=ggdt[hr_type == "10-fold cross-validation"], position=jitterer, aes(color=hr_type, shape=hr_type, size=hr_type, fill=hr_type, alpha=hr_type)) +
  geom_errorbarh(data=ggdt[hr_type == "full dataset"], height=0, aes(color=hr_type, alpha=hr_type), position=position_nudge(y=-0.3)) +
  geom_point(data=ggdt[hr_type == "full dataset"], aes(color=hr_type, shape=hr_type, size=hr_type, fill=hr_type, alpha=hr_type), position=position_nudge(y=-0.3)) +
  geom_errorbarh(data=ggdt[hr_type == "cross-validation mean"], height=0, aes(color=hr_type, alpha=hr_type)) +
  geom_point(data=ggdt[hr_type == "cross-validation mean"], aes(color=hr_type, shape=hr_type, size=hr_type, fill=hr_type, alpha=hr_type)) +
  scale_color_manual(name="Source", values=c("10-fold cross-validation"="#525252", "full dataset"="#810f7c", "cross-validation mean"="black")) +
  scale_size_manual(name="Source", values=c("10-fold cross-validation"=1, "full dataset"=1, "cross-validation mean"=1.2)) + 
  scale_shape_manual(name="Source", values=c("10-fold cross-validation"=19, "full dataset"=19, "cross-validation mean"=23)) +
  scale_fill_manual(name="Source", values=c("10-fold cross-validation"="#525252", "full dataset"="#810f7c", "cross-validation mean"="#41ab5d")) +
  scale_alpha_manual(name="Source", values=c("10-fold cross-validation"=0.6, "full dataset"=1, "cross-validation mean"=1)) +
  ylab("") + xlab("Hazard Ratio (95% CI) per SD") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="bottom", legend.justification="left", legend.box.margin=margin(0, 0, 0, -2.5, unit="cm"),
    legend.title=element_text(size=8), legend.text=element_text(size=6)
  ) + 
  guides(
    color=guide_legend(title.position="top"), size=guide_legend(title.position="top"),
    shape=guide_legend(title.position="top"), fill=guide_legend(title.position="top"), 
    alpha=guide_legend(title.position="top")
  )
ggsave(g, width=4, height=6.4, file="analyses/CVD_score_weighting/HRs_in_cross_validation_SCORE2_hidden.pdf")


  







  







