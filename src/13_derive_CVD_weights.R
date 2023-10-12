library(data.table)
library(foreach)
library(survival)
library(ggplot2)
library(cowplot)
source("src/utils/SCORE2.R")

# Make output directory
system("mkdir -p analyses/CVD_weight_training", wait=TRUE)

# Load required data
dat <- fread("data/cleaned/analysis_cohort.txt", select=c("eid", "sex", "age", "smoking", "sbp", "tchol", "hdl", "incident_cvd_followup", "incident_cvd", "cvd_prediction_foldid", "CAD_metaGRS", "Stroke_metaGRS", "SCORE2_excl_UKB"))
setnames(dat, "SCORE2_excl_UKB", "SCORE2")

# Load NMR scores from model training
train_scores <- rbind(idcol="type", 
  "non-derived"=fread("analyses/nmr_score_training/non_derived_NMR_scores.txt"),
  "clinical"=fread("analyses/nmr_score_training/clinical_NMR_scores.txt")
)

# Define population matching SCORE2 risk factor reference distributions so that
# we can center NMR scores and PRSs to match appropriate baseline survival used
# by SCORE2 compute absolute risk and its subsequent recalibration
ref_dat <- dat[
  (age >= 55 & age <= 65) &
  (is.na(smoking) | !(smoking)) &
  (sbp >= 100 & sbp <= 140) & 
  (tchol >= 5 & tchol <= 7) &
  (hdl >= 0.8 & hdl <= 1.8)
]

# Explore distributions of SCORE2, PRSs, and NMR scores in the full training and reference populations
comp <- train_scores[type == "non-derived"]
comp[, dataset := ifelse(cvd_prediction_foldid == prediction_cv_testfold, "test", "training")]
comp[, score := paste(sex, endpoint, "NMR_score", prediction_cv_testfold, dataset, sep="_")]
comp <- comp[,.(eid, score, linear_predictor)]
comp <- rbind(comp, melt(dat, id.vars="eid", measure.vars=c("CAD_metaGRS", "Stroke_metaGRS", "SCORE2"), variable.name="score", value.name="linear_predictor"))
comp[dat, on = .(eid), sex := i.sex]
comp[score %in% c("CAD_metaGRS", "Stroke_metaGRS", "SCORE2"), score := paste0(sex, "_", score)]
score_order <- c(paste0("CAD_NMR_score_", 1:5, "_training"), paste0("CAD_NMR_score_", 1:5, "_test"), paste0("Stroke_NMR_score_", 1:5, "_training"), paste0("Stroke_NMR_score_", 1:5, "_test"))
score_order <- c("SCORE2", "CAD_metaGRS", "Stroke_metaGRS", score_order)
score_order <- c(paste0("Male_", score_order), paste0("Female_", score_order))
comp[,score := factor(score, levels=score_order)]
comp <- rbind(idcol="population",
  "All 168,517 study participants"=comp,
  "18,673 study participants within 1 SD of SCORE2 risk factor means"=comp[eid %in% ref_dat$eid]
)
comp[, population := factor(population, levels=c("All 168,517 study participants", "18,673 study participants within 1 SD of SCORE2 risk factor means"))]

g <- ggplot(comp) + 
  aes(x=linear_predictor, color=population) +
  facet_wrap(~ score, nrow=6) +
  geom_density() +
  geom_vline(xintercept=0, linetype=2) +
  xlab("Linear predictor") +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8),
    strip.background=element_blank(), strip.text=element_text(size=5, face="bold"),
    legend.position="bottom", legend.text=element_text(size=7), legend.title=element_blank()
  )
ggsave(g, width=12, height=8, file="analyses/CVD_weight_training/score_distributions.pdf")
  
# Compute sex-specific centering factors in the reference group matched to SCORE2 risk factor reference
# so that baseline hazards are appropriate later when computing absolute risk from combined scores
prs_scaling <- melt(ref_dat, id.vars=c("eid", "sex"), measure.vars=c("CAD_metaGRS", "Stroke_metaGRS"), variable.name="score", value.name="level")
prs_scaling <- prs_scaling[, .(mean=mean(level)), by=.(sex, score)]

nmr_score_scaling <- train_scores[eid %in% ref_dat$eid & cvd_prediction_foldid != prediction_cv_testfold]
nmr_score_scaling <- nmr_score_scaling[, .(mean=mean(linear_predictor)), by=.(type, sex, endpoint, prediction_cv_testfold)]

# Write out scaling factors
fwrite(prs_scaling, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/prs_SCORE2_reference_centering.txt")
fwrite(nmr_score_scaling, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/nmr_score_SCORE2_reference_centering_in_5fold_cross_validation.txt")

# Apply scaling factors to data
dat[prs_scaling[score == "CAD_metaGRS"], on = .(sex), CAD_metaGRS := (CAD_metaGRS - i.mean)]
dat[prs_scaling[score == "Stroke_metaGRS"], on = .(sex), Stroke_metaGRS := (Stroke_metaGRS - i.mean)]
train_scores[nmr_score_scaling, on = .(type, sex, endpoint, prediction_cv_testfold), linear_predictor := (linear_predictor - i.mean)]

# Explore distributions of SCORE2, PRSs, and NMR scores in the full training and reference populations after centering
comp <- train_scores[type == "non-derived"]
comp[, dataset := ifelse(cvd_prediction_foldid == prediction_cv_testfold, "test", "training")]
comp[, score := paste(sex, endpoint, "NMR_score", prediction_cv_testfold, dataset, sep="_")]
comp <- comp[,.(eid, score, linear_predictor)]
comp <- rbind(comp, melt(dat, id.vars="eid", measure.vars=c("CAD_metaGRS", "Stroke_metaGRS", "SCORE2"), variable.name="score", value.name="linear_predictor"))
comp[dat, on = .(eid), sex := i.sex]
comp[score %in% c("CAD_metaGRS", "Stroke_metaGRS", "SCORE2"), score := paste0(sex, "_", score)]
score_order <- c(paste0("CAD_NMR_score_", 1:5, "_training"), paste0("CAD_NMR_score_", 1:5, "_test"), paste0("Stroke_NMR_score_", 1:5, "_training"), paste0("Stroke_NMR_score_", 1:5, "_test"))
score_order <- c("SCORE2", "CAD_metaGRS", "Stroke_metaGRS", score_order)
score_order <- c(paste0("Male_", score_order), paste0("Female_", score_order))
comp[,score := factor(score, levels=score_order)]
comp <- rbind(idcol="population",
  "All 168,517 study participants"=comp,
  "18,673 study participants within 1 SD of SCORE2 risk factor means"=comp[eid %in% ref_dat$eid]
)
comp[, population := factor(population, levels=c("All 168,517 study participants", "18,673 study participants within 1 SD of SCORE2 risk factor means"))]

g <- ggplot(comp) + 
  aes(x=linear_predictor, color=population) +
  facet_wrap(~ score, nrow=6) +
  geom_density() +
  geom_vline(xintercept=0, linetype=2) +
  xlab("Linear predictor") +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8),
    strip.background=element_blank(), strip.text=element_text(size=5, face="bold"),
    legend.position="bottom", legend.text=element_text(size=7), legend.title=element_blank()
  )
ggsave(g, width=12, height=8, file="analyses/CVD_weight_training/score_distributions_after_centering.pdf")

# Estimate per-score weights to use when combining with SCORE2 using Cox proportional hazards models
cvd_weights <- foreach(this_score_type = c("non-derived", "clinical"), .combine=rbind) %:%
  foreach(this_model = c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %:%
    foreach(this_sex = c("Male", "Female"), .combine=rbind) %:%
      foreach(this_test_fold = 1:5, .combine=rbind) %do% {
        # extract training data
        this_dat <- dat[cvd_prediction_foldid != this_test_fold & sex == this_sex]

        # add in respective nmr scores
        this_train_scores <- train_scores[prediction_cv_testfold == this_test_fold & prediction_cv_testfold != cvd_prediction_foldid & sex == this_sex & type == this_score_type]
        this_dat[this_train_scores[endpoint == "CAD"], on = .(eid), CAD_NMR_score := i.linear_predictor]
        this_dat[this_train_scores[endpoint == "Stroke"], on = .(eid), Stroke_NMR_score := i.linear_predictor]
  
        # Build model formula
        mf <- "Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2)"
        if (this_model %in% c("SCORE2 + NMR scores", "SCORE2 + NMR scores + PRSs")) {
          mf <- paste(mf, "+ CAD_NMR_score + Stroke_NMR_score")
        }
        if (this_model %in% c("SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs")) {
          mf <- paste(mf, "+ CAD_metaGRS + Stroke_metaGRS")
        }
   
        # Fit cox model
        cx <- coxph(as.formula(mf), data=this_dat)

        # Extract relevant coefficients
        ci <- confint(cx)
        cf <- as.data.table(coef(summary(cx)), keep.rownames="score")
        cf[, score := gsub("_", " ", score)]
        cf <- cf[, .(score, weight=coef, L95=ci[,1], U95=ci[,2], pval=`Pr(>|z|)`)]

        # add information
        cbind(test_fold = this_test_fold, sex = this_sex, model = this_model, score_type = this_score_type, cf)
}

# Write out
fwrite(cvd_weights, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/cvd_score_weights_in_5fold_cross_validation.txt")

# Build supp table of avg weights and transformations
dt <- nmr_score_scaling[type == "non-derived", .(centering_factor=mean(mean)), by=.(sex, score=paste(endpoint, "NMR score"))]
dt <- rbind(dt, prs_scaling[,.(sex, score=gsub("_", " ", score), centering_factor=mean)])
dt2 <- cvd_weights[score_type == "non-derived", .(weight=mean(weight)), by=.(model, sex, score)]
dt <- dt[dt2, on = .(sex, score)]
fwrite(dt, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/avg_cvd_score_weights_supp_table.txt")

# Factor levels for plot ordering
cvd_weights[, model := factor(model, levels=c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
cvd_weights[, score := factor(score, levels=c("CAD NMR score", "Stroke NMR score", "CAD metaGRS", "Stroke metaGRS"))]

# Plot weight consistencies
plot_weights <- function(this_sex, this_score_type) {
  ggplot(cvd_weights[sex == this_sex & score_type == this_score_type]) +
    aes(x = score, y = weight, ymin = L95, ymax = U95, color = paste("Test fold", test_fold)) +
    facet_grid(. ~ model, space="free_x", scales="free_x") +
    geom_errorbar(position=position_dodge(width=0.75), alpha=0.6, width=0) +
    geom_point(shape=19, alpha=0.6, position=position_dodge(width=0.75)) +
    geom_hline(yintercept=0, linetype=2) +
    xlab("") + ylab("Weight") +
    guides(color=guide_legend(title="")) +
    theme_bw() +
    theme(
      axis.text.y=element_text(size=6), axis.title.y=element_text(size=8),
      axis.text.x=element_text(size=6, angle=90, hjust=1, vjust=0.5),
      strip.background=element_blank(), strip.text=element_text(size=8, face="bold"),
      panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
      legend.position="bottom", legend.text=element_text(size=6),
      legend.box.margin=margin(-7,-3,0,-3), legend.margin=margin(-3,-3,-3,-3)
    )
}

g1 <- plot_weights("Male", "non-derived")
g2 <- plot_weights("Female", "non-derived")
g <- plot_grid(g1, g2, nrow=2, labels=c("Males  ", "Females"), label_size=10)
ggsave(g, width=7.2, height=4.5, file="analyses/CVD_weight_training/non_derived_NMR_scores_weights.pdf")

g1 <- plot_weights("Male", "clinical")
g2 <- plot_weights("Female", "clinical")
g <- plot_grid(g1, g2, nrow=2, labels=c("Males  ", "Females"), label_size=10)
ggsave(g, width=7.2, height=4.5, file="analyses/CVD_weight_training/clinical_NMR_scores_weights.pdf")

# Build combined linear predictor for each model from aggregate predictions using the scaling factors and
# weights from the 5-fold cross-validation ctraining in the withheld test fold.
pred_scores <- foreach(this_score_type = c("non-derived", "clinical"), .combine=rbind) %:%
  foreach(this_model = c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %:%
    foreach(this_sex = c("Male", "Female"), .combine=rbind) %:%
      foreach(this_test_fold = 1:5, .combine=rbind) %do% {
        # extract test data
        this_dat <- dat[cvd_prediction_foldid == this_test_fold & sex == this_sex]

        # Add in relevant NMR scores (previously centered using appropriate means from training data)
        this_nmr_scores <- train_scores[cvd_prediction_foldid == prediction_cv_testfold & cvd_prediction_foldid == this_test_fold & sex == this_sex & type == this_score_type]
        this_nmr_scores[, score := paste0(endpoint, "_NMR_score")]
        this_nmr_scores <- dcast(this_nmr_scores, eid ~ score, value.var="linear_predictor")
        this_dat <- this_nmr_scores[this_dat, on = .(eid)]

        # Melt columns of interest to long format for summing
        lp_cols <- "SCORE2"
        if (this_model %like% "NMR") {
          lp_cols <- c(lp_cols, "CAD_NMR_score", "Stroke_NMR_score")
        } 
        if (this_model %like% "PRS") {
          lp_cols <- c(lp_cols, "CAD_metaGRS", "Stroke_metaGRS")
        }
        this_dat <- melt(this_dat, id.vars=c("eid", "sex", "age", "incident_cvd", "incident_cvd_followup", "cvd_prediction_foldid"), measure.vars=lp_cols, variable.name="score")

        # Multiply scores by weights
        this_weights <- cvd_weights[test_fold == this_test_fold & sex == this_sex & model == this_model & score_type == this_score_type]
        this_weights[, score := gsub(" ", "_", score)]
        this_dat[this_weights, on = .(score), value := value * i.weight]
      
        # Sum to obtain linear predictor
        this_dat <- this_dat[, .(linear_predictor=sum(value)), by=.(eid, sex, age, incident_cvd, incident_cvd_followup, cvd_prediction_foldid)]
   
        # add model info
        cbind(score_type = this_score_type, model = this_model, this_dat)
}

# Add five year age group (downstream analyses uses these)
pred_scores[, age_group := sprintf("%s-%s", age %/% 5 * 5, age %/% 5 * 5 + 4)]

# Compute absolute risk using SCORE2 baseline hazards
pred_scores[, uncalibrated_risk := score2_absrisk(sex, linear_predictor)]

# Compute absolute risk recalibrated to low-risk european region (including UK)
pred_scores[, uk_calibrated_risk := score2_recalibration(sex, uncalibrated_risk, "low")]

# Reorganize columns
pred_scores <- pred_scores[, .(eid, sex, age, age_group, incident_cvd, incident_cvd_followup, cvd_prediction_foldid,
  model, score_type, linear_predictor, uncalibrated_risk, uk_calibrated_risk)]

# Write out
fwrite(pred_scores, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/CVD_linear_predictors_and_risk.txt")

# Compare abs-risk distributions 
comp <- pred_scores[score_type == "non-derived"]
comp <- dcast(comp, eid + sex ~ model, value.var="uk_calibrated_risk")
comp <- melt(comp, id.vars=c("eid", "sex", "SCORE2"), variable.name="model", value.name="new_absrisk")

comp[, sex := factor(sex, levels=c("Male", "Female"))]
comp[, model := factor(model, levels=c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]

g <- ggplot(comp) +
  aes(x = SCORE2, y = new_absrisk) + 
  facet_grid(sex ~ model) +
  geom_hex() + 
  geom_abline(slope=1, intercept=0, linetype=2, color="red") +
  scale_fill_gradient(name="Participants", low="lightblue1", high="darkblue", trans="log10") +
  xlab("10-year CVD risk with SCORE2 alone") +
  ylab("10-year CVD risk with new model") +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    legend.title=element_text(size=7), legend.text=element_text(size=6)
  )
ggsave(g, width=7.2, height=5, file="analyses/CVD_weight_training/absrisk_comparisons.pdf")





