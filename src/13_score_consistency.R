library(data.table)
library(foreach)
library(ggplot2)
library(forcats)
library(ggrastr)
library(ggpp)
library(ggh4x)
library(ggthemes)
library(cowplot)
library(survival)

options("ggrastr.default.dpi" = 300) 

# Load dataset
dat <- fread("data/cleaned/analysis_cohort.txt")

# Load biomarker information sheet
nmr_info <- fread("data/ukb/NMR_metabolomics/biomarker_information.txt")

# Load scaling factors for NMR data
nmr_scaling <- fread("data/standardised/nmr_scaling_factors.txt")

# Load NMR score coefficients
nmr_coef <- fread("analyses/nmr_score_training/best_fits_coef.txt")

# Build long table of NMR coefficients and age interactions standardised using the given scaling factors
nmr_long <- melt(dat, 
  id.vars=c("eid", "sex", "age", "cvd_prediction_foldid"),
  measure.vars=nmr_coef[!(coef %like% "^age"), unique(coef)])
nmr_long[nmr_scaling, on = .(sex, variable=biomarker), value := (value - i.mean)/i.sd]
nmr_long[, age := (age - 60)/5]
nmr_age_interactions <- copy(nmr_long)
nmr_age_interactions[, value := value * age]
nmr_age_interactions[, variable := paste0("age:", variable)]
nmr_long <- rbind(nmr_long, nmr_age_interactions)

# Compute NMR biomarker scores
models <- unique(nmr_coef[,.(prediction_cv_testfold, endpoint, lambda.fit, type)])
scores <- foreach(midx = models[,.I], .combine=rbind) %do% {
  this_minfo <- models[midx]
  this_coef <- nmr_coef[this_minfo, on = .(prediction_cv_testfold, endpoint, lambda.fit, type)]
  this_sinfo <- dat[,.(eid, sex, cvd_prediction_foldid)] 

  this_score <- nmr_long[this_coef, on = .(sex, variable=coef)][, .(
    linear_predictor = sum(beta * value),
    lp_simplified = sum(round(beta, digits=3) * value)
  ), by=.(eid)]

  this_score <- this_sinfo[this_score, on = .(eid), nomatch=0]

  if (nrow(this_score) > 0) {
    return(cbind(this_minfo, this_score))
  }
}

# Compute average score across all 5 training runs
avg_nmr_coef <- fread("analyses/nmr_score_training/avg_best_coefs.txt")
models <- unique(models[,.(endpoint, lambda.fit, type)])
avg_scores <- foreach(midx = models[,.I], .combine=rbind) %do% {
  this_minfo <- models[midx]
  this_coef <- avg_nmr_coef[this_minfo, on = .(endpoint, lambda.fit, type), nomatch=0]
  this_sinfo <- dat[,.(eid, sex, cvd_prediction_foldid)] 

  this_score <- nmr_long[this_coef, on = .(sex, variable=coef)][, .(
    linear_predictor = sum(beta * value),
    lp_simplified = sum(round(beta, digits=3) * value)
  ), by=.(eid)]

  this_score <- this_sinfo[this_score, on = .(eid), nomatch=0]

  if (nrow(this_score) > 0) {
    return(cbind(this_minfo, this_score))
  }
}

# Assess consistency of NMR biomarker scores across training runs
score_comp <- foreach(this_test_fold = c(1:5, "avg"), .combine=rbind) %do% {
  if (this_test_fold == "avg") {
    this_x <- avg_scores[, .(eid, cvd_prediction_foldid, sex, endpoint, lambda.fit, type,
        x_model="avg", x_score=linear_predictor)]
  } else {
    this_x <- scores[prediction_cv_testfold == this_test_fold,  
      .(eid, cvd_prediction_foldid, sex, endpoint, lambda.fit, type,
        x_model=prediction_cv_testfold, x_score=linear_predictor)]
  } 
  this_y <- scores[, .(eid, cvd_prediction_foldid, sex, endpoint, lambda.fit, type,
      y_model=prediction_cv_testfold, y_score=linear_predictor)]
  this_y <- rbind(this_y, avg_scores[,
    .(eid, cvd_prediction_foldid, sex, endpoint, lambda.fit, type,
      y_model="avg", y_score=linear_predictor)])
  merge(this_x, this_y, by=c("eid", "cvd_prediction_foldid", "sex", "endpoint", "lambda.fit", "type"), all=FALSE)
}

# Color annotation for plots
score_comp[, color_anno := fcase(
  x_model != cvd_prediction_foldid & y_model != cvd_prediction_foldid, "Training samples",
  x_model == cvd_prediction_foldid, "X-axis test samples", 
  y_model == cvd_prediction_foldid, "Y-axis test samples" 
)]

# Create scatterplots of scores
score_scatter <- function(score_comp) {
  # Better facet names
  facet_anno <- function(vars) { 
    sapply(vars, function(ii) {
      if (ii == "avg") {
        "Average"
      } else {
        paste("Test fold", ii) 
      }
    })
  }

  # correlation annoations
  cor_anno <- score_comp[,.(label=sprintf("r=%.2f\np=%.2f", 
    cor(x_score, y_score), 
    cor(x_score, y_score, method="spearman")
  )), by=.(x_model, y_model)]

  ggplot(score_comp) +
    aes(x=x_score, y=y_score, color=color_anno) +
    rasterise(geom_point(data=score_comp[x_model != cvd_prediction_foldid & y_model != cvd_prediction_foldid], alpha=0.5, shape=19, size=0.1)) +
    rasterise(geom_point(data=score_comp[x_model == cvd_prediction_foldid & y_model != cvd_prediction_foldid], alpha=0.5, shape=19, size=0.1)) +
    rasterise(geom_point(data=score_comp[x_model != cvd_prediction_foldid & y_model == cvd_prediction_foldid], alpha=0.5, shape=19, size=0.1)) +
    rasterise(geom_point(data=score_comp[x_model == cvd_prediction_foldid & y_model == cvd_prediction_foldid], alpha=0.5, shape=19, size=0.1)) +
    scale_color_manual(name="", values=c("Training samples"="#252525", "X-axis test samples"="#d94801", "Y-axis test samples"="#006d2c")) +
    guides(color = guide_legend(override.aes = list(size = 1))) +
    geom_abline(intercept=0, slope=1, linetype=2, color="red") +
    geom_text_npc(data=cor_anno, npcx="left", npcy="top", aes(label=label), color="red", size=6*0.352777778) +
    facet_grid2(fct_rev(y_model) ~ x_model, scales="free", independent="all", labeller=as_labeller(facet_anno)) +
    xlab("NMR score") + ylab("NMR score") +
    theme_bw() +
    theme(
      axis.text=element_text(size=6), axis.title=element_text(size=8),
      strip.text=element_text(size=8, face="bold"), strip.background=element_blank(), 
      panel.grid=element_blank(), legend.text=element_text(size=6)
    )
}
g1 <- score_scatter(score_comp[endpoint == "CAD" & sex == "Male" & lambda.fit == "lambda.min" & type == "non-derived" & x_model != "avg" & y_model != "avg"])
g2 <- score_scatter(score_comp[endpoint == "CAD" & sex == "Female" & lambda.fit == "lambda.min" & type == "non-derived" & x_model != "avg" & y_model != "avg"])
g3 <- score_scatter(score_comp[endpoint == "Stroke" & sex == "Male" & lambda.fit == "lambda.min" & type == "non-derived" & x_model != "avg" & y_model != "avg"])
g4 <- score_scatter(score_comp[endpoint == "Stroke" & sex == "Female" & lambda.fit == "lambda.min" & type == "non-derived" & x_model != "avg" & y_model != "avg"])

ggsave(g1, width=7, height=6, file="analyses/nmr_score_training/CAD_male_non_derived_score_consistency.pdf")
ggsave(g2, width=7, height=6, file="analyses/nmr_score_training/CAD_female_non_derived_score_consistency.pdf")
ggsave(g3, width=7, height=6, file="analyses/nmr_score_training/Stroke_male_non_derived_score_consistency.pdf")
ggsave(g4, width=7, height=6, file="analyses/nmr_score_training/Stroke_female_non_derived_score_consistency.pdf")

g1 <- score_scatter(score_comp[endpoint == "CAD" & sex == "Male" & lambda.fit == "lambda.min" & type == "clinical" & x_model != "avg" & y_model != "avg"])
g2 <- score_scatter(score_comp[endpoint == "CAD" & sex == "Female" & lambda.fit == "lambda.min" & type == "clinical" & x_model != "avg" & y_model != "avg"])
g3 <- score_scatter(score_comp[endpoint == "Stroke" & sex == "Male" & lambda.fit == "lambda.min" & type == "clinical" & x_model != "avg" & y_model != "avg"])
g4 <- score_scatter(score_comp[endpoint == "Stroke" & sex == "Female" & lambda.fit == "lambda.min" & type == "clinical" & x_model != "avg" & y_model != "avg"])

ggsave(g1, width=7, height=6, file="analyses/nmr_score_training/CAD_male_clinical_score_consistency.pdf")
ggsave(g2, width=7, height=6, file="analyses/nmr_score_training/CAD_female_clinical_score_consistency.pdf")
ggsave(g3, width=7, height=6, file="analyses/nmr_score_training/Stroke_male_clinical_score_consistency.pdf")
ggsave(g4, width=7, height=6, file="analyses/nmr_score_training/Stroke_female_clinical_score_consistency.pdf")

# Make reduce plot comparing predicted scores in each 20% test fold to training scores 
score_scatter2 <- function(lambda, score_type) {
  ggdt <- score_comp[lambda.fit == lambda & type == score_type]
  ggdt <- ggdt[x_model == cvd_prediction_foldid & y_model != x_model]
  ggdt <- ggdt[, color_anno := sprintf("Test fold %s", y_model)]
  ggdt[color_anno == "Test fold avg", color_anno := "Average score"]
  ggdt[, color_annon := factor(color_anno, levels=c(paste("Test fold", 1:5), "Average score"))]
  ggdt[, facet_anno := sprintf("Test fold %s", x_model)]

  ggplot(ggdt) +
    aes(x=x_score, y=y_score, color=color_anno) +
    facet_grid2(endpoint + fct_rev(sex) ~ facet_anno, scales="free", independent="all") +
    rasterise(geom_point(alpha=0.5, shape=19, size=0.1)) +
    geom_abline(intercept=0, slope=1, linetype=2, color="red") +
    scale_color_manual(values=structure(c(ggthemes::calc_pal()(5), "black"), names=c(paste("Test fold", 1:5), "Average score"))) +
    xlab("NMR score predicted in 20% of withheld test-fold data") +
    ylab("NMR scores including training data") + 
    guides(color=guide_legend("Withheld test fold for score training", override.aes=list(size=1, alpha=1), title.position="top", nrow=1)) +
    theme_bw() +
    theme(
      axis.text=element_text(size=6), axis.title=element_text(size=8),
      strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
      panel.grid=element_blank(), legend.text=element_text(size=6), legend.title=element_text(size=8),
      legend.position="bottom"
    )
}

g <- score_scatter2("lambda.min", "non-derived")
ggsave(g, width=7.2, height=6.6, file="analyses/nmr_score_training/non_derived_score_consistency_simplified.pdf")

g <- score_scatter2("lambda.min", "clinical")
ggsave(g, width=7.2, height=6.6, file="analyses/nmr_score_training/clinical_score_consistency_simplified.pdf")

# Compare average score to test score
avg_vs_test_scatter <- function(lambda, score_type) {
  cor_anno <- score_comp[x_model == "avg" & y_model == cvd_prediction_foldid & lambda.fit == lambda & type == score_type,
    .(label=sprintf("r=%.2f\np=%.2f",
      cor(x_score, y_score),
      cor(x_score, y_score, method="spearman")
    )), by=.(endpoint, sex)]

  ggplot(score_comp[x_model == "avg" & y_model == cvd_prediction_foldid & lambda.fit == lambda & type == score_type]) +
    aes(x=x_score, y=y_score, color=paste("Test fold", y_model)) +
    rasterise(geom_point(alpha=0.5, shape=19, size=0.1)) +
    guides(color = guide_legend(title="", override.aes = list(size = 1))) +
    geom_abline(intercept=0, slope=1, linetype=2, color="red") +
    geom_text_npc(data=cor_anno, npcx="left", npcy="top", aes(label=label), color="red", size=6*0.352777778) +
    facet_grid2(fct_rev(sex) ~ endpoint, scales="free", independent="all") +
    xlab("NMR score from coefficient averages") + ylab("NMR scores aggregated across test folds") +
    theme_bw() +
    theme(
      axis.text=element_text(size=6), axis.title=element_text(size=8),
      strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
      panel.grid=element_blank(), legend.text=element_text(size=6), legend.title=element_text(size=8)
    )
}

g <- avg_vs_test_scatter("lambda.min", "non-derived")
ggsave(g, width=7, height=5.5, file="analyses/nmr_score_training/aggregate_test_vs_coef_average_non_derived_scores.pdf")

g <- avg_vs_test_scatter("lambda.min", "clinical")
ggsave(g, width=7, height=5.5, file="analyses/nmr_score_training/aggregate_test_vs_coef_average_clinical_scores.pdf")

# Can we report coefficients to 3 decimal places (and drop any with effect < 0.001)
full_vs_truncated_scatter <- function(lambda, score_type) {
  cor_anno <- avg_scores[lambda.fit == lambda & type == score_type,
    .(label=sprintf("r=%.2f\np=%.2f",
      cor(linear_predictor, lp_simplified),
      cor(linear_predictor, lp_simplified, method="spearman")
    )), by=.(endpoint, sex)]

  ggplot(avg_scores[lambda.fit == lambda & type == score_type]) +
    aes(x=linear_predictor, y=lp_simplified, color=paste("Test fold", cvd_prediction_foldid)) +
    rasterise(geom_point(alpha=0.5, shape=19, size=0.1)) +
    guides(color = guide_legend(title="", override.aes = list(size = 1))) +
    geom_abline(intercept=0, slope=1, linetype=2, color="red") +
    geom_text_npc(data=cor_anno, npcx="left", npcy="top", aes(label=label), color="red", size=6*0.352777778) +
    facet_grid2(fct_rev(sex) ~ endpoint, scales="free", independent="all") +
    xlab("Average NMR score from full precision coefficients") + 
    ylab("Average NMR score from coefficients rounded to 3 decimal places") +
    theme_bw() +
    theme(
      axis.text=element_text(size=6), axis.title=element_text(size=8),
      strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
      panel.grid=element_blank(), legend.text=element_text(size=6), legend.title=element_text(size=8)
    )
}

g <- full_vs_truncated_scatter("lambda.min", "non-derived")
ggsave(g, width=7, height=5.5, file="analyses/nmr_score_training/full_coef_vs_truncated_coef_non_derived_scores.pdf")

g <- full_vs_truncated_scatter("lambda.min", "clinical")
ggsave(g, width=7, height=5.5, file="analyses/nmr_score_training/full_coef_vs_truncated_coef_clinical_scores.pdf")

# Compare NMR scores derived from clinically calibrated biomarkers to NMR scores from all non-derivable components
clinical_vs_non_derived_scatter <- function(lambda, outcome) {
  ggdt <- score_comp[lambda.fit == lambda & x_model == y_model & endpoint == outcome &(x_model == cvd_prediction_foldid | x_model == "avg")]
  ggdt <- dcast(ggdt, eid + cvd_prediction_foldid + sex + endpoint + lambda.fit + x_model ~ type, value.var="x_score")
  setnames(ggdt, "non-derived", "non_derived")
  setnames(ggdt, "x_model", "score")
  ggdt[, score := sprintf("Score %s", score)]

  cor_anno <- ggdt[,
    .(label=sprintf("r=%.2f\np=%.2f",
      cor(non_derived, clinical),
      cor(non_derived, clinical, method="spearman")
    )), by=.(sex, score)]

  ggplot(ggdt) +
    aes(x=clinical, y=non_derived) +
    facet_grid2(fct_rev(sex) ~ score, scales="free", independent="all") +
    rasterise(geom_point(alpha=0.5, shape=19, size=0.1)) +
    geom_abline(intercept=0, slope=1, linetype=2, color="red") +
    geom_text_npc(data=cor_anno, npcx="left", npcy="top", aes(label=label), color="red", size=6*0.352777778) +
    xlab("NMR score from 21 clinically accredited biomarkers") + 
    ylab("NMR score from 106 non-derived biomarkers") +
    theme_bw() +
    theme(
      axis.text=element_text(size=6), axis.title=element_text(size=8),
      strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
      panel.grid=element_blank(), legend.text=element_text(size=6), legend.title=element_text(size=8)
    )
}

g1 <- clinical_vs_non_derived_scatter("lambda.min", "CAD")
g2 <- clinical_vs_non_derived_scatter("lambda.min", "Stroke")
g <- plot_grid(g1, g2, nrow=2, labels=c("CAD   ", "Stroke"), label_size=10)
ggsave(g, width=7, height=5.5, file="analyses/nmr_score_training/clinical_vs_non_derived_biomarker_scores.pdf")

# Write out full scores including training data
full_scores <- scores[lambda.fit == "lambda.min" & type == "non-derived"]
full_scores <- full_scores[,.(eid, sex, cvd_prediction_foldid, prediction_cv_testfold, endpoint, linear_predictor)]
fwrite(full_scores, sep="\t", quote=FALSE, file="analyses/nmr_score_training/non_derived_NMR_scores.txt")

full_scores <- scores[lambda.fit == "lambda.min" & type == "clinical"]
full_scores <- full_scores[,.(eid, sex, cvd_prediction_foldid, prediction_cv_testfold, endpoint, linear_predictor)]
fwrite(full_scores, sep="\t", quote=FALSE, file="analyses/nmr_score_training/clinical_NMR_scores.txt")

# Write out aggregate scores
test_scores <- scores[lambda.fit == "lambda.min" & prediction_cv_testfold == cvd_prediction_foldid & type == "non-derived"]
test_scores <- dcast(test_scores, eid  ~ endpoint, value.var="linear_predictor")
setnames(test_scores, c("CAD", "Stroke"), c("CAD_NMR_score", "Stroke_NMR_score"))
fwrite(test_scores, sep="\t", quote=FALSE, file="analyses/nmr_score_training/aggregate_test_non_derived_NMR_scores.txt")

test_scores <- scores[lambda.fit == "lambda.min" & prediction_cv_testfold == cvd_prediction_foldid & type == "clinical"]
test_scores <- dcast(test_scores, eid  ~ endpoint, value.var="linear_predictor")
setnames(test_scores, c("CAD", "Stroke"), c("CAD_NMR_score", "Stroke_NMR_score"))
fwrite(test_scores, sep="\t", quote=FALSE, file="analyses/nmr_score_training/aggregate_test_clinical_NMR_scores.txt")

# Also average scores for downstream checks
avg_scores_wide <- dcast(avg_scores[lambda.fit == "lambda.min" & type == "non-derived"], eid ~ endpoint, value.var="linear_predictor")
setnames(avg_scores_wide, c("CAD", "Stroke"), c("CAD_NMR_score", "Stroke_NMR_score"))
fwrite(avg_scores_wide, sep="\t", quote=FALSE, file="analyses/nmr_score_training/coef_avg_non_derived_NMR_scores.txt")

avg_scores_wide <- dcast(avg_scores[lambda.fit == "lambda.min" & type == "clinical"], eid ~ endpoint, value.var="linear_predictor")
setnames(avg_scores_wide, c("CAD", "Stroke"), c("CAD_NMR_score", "Stroke_NMR_score"))
fwrite(avg_scores_wide, sep="\t", quote=FALSE, file="analyses/nmr_score_training/coef_avg_clinical_NMR_scores.txt")

# For each score, how many coefficients are included across N training runs
possible <- nmr_info[Type == "Non-derived" & Biomarker != "Clinical_LDL_C", Biomarker]
possible <- c(possible, paste0("age:", possible))
coef_stats <- expand.grid(endpoint=c("CAD", "Stroke"), sex=c("Male", "Female"), lambda.fit=c("lambda.min", "lambda.1se"), coef=possible, N_scores=0, N_scores_sig=0)
setDT(coef_stats)
coef_stats[nmr_coef[type == "non-derived",.N,by=.(endpoint, sex, lambda.fit, coef)], on = .(endpoint, sex, lambda.fit, coef), N_scores := i.N]
coef_stats[nmr_coef[type == "non-derived",.(N=sum(round(beta, digits=3) != 0)), by=.(endpoint, sex, lambda.fit, coef)], on = .(endpoint, sex, lambda.fit, coef), N_scores_sig := i.N]
coef_stats[, N_score_sig_incl_avg := N_scores_sig]
coef_stats[!avg_nmr_coef[type == "non-derived" & round(beta, digits=3) != 0], on = .(endpoint, sex, lambda.fit, coef), N_score_sig_incl_avg := 0]

ggdt <- coef_stats[lambda.fit=="lambda.min", .(N_coef = .N), by=.(endpoint, sex, N_scores)]
g <- ggplot(ggdt) + 
  aes(x=factor(N_scores), y=N_coef) +
  facet_grid(sex ~ endpoint, drop=FALSE, shrink=FALSE) + 
  geom_col() +
  scale_x_discrete("Number of scores coefficient is included in") +
  scale_y_continuous("Number of coefficients", limits=c(0, 212)) +
  theme_bw() + 
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid=element_blank()
  )
ggsave(g, width=7, height=5.5, file="analyses/nmr_score_training/coefficient_inclusion_matrix_non_derived_scores.pdf")

ggdt <- coef_stats[lambda.fit=="lambda.min", .(N_coef = .N), by=.(endpoint, sex, N_scores_sig)]
g <- ggplot(ggdt) + 
  aes(x=factor(N_scores_sig), y=N_coef) +
  facet_grid(sex ~ endpoint, drop=FALSE, shrink=FALSE) + 
  geom_col() +
  scale_x_discrete("Number of scores coefficient is included in after rounding to 3 decimal places") +
  scale_y_continuous("Number of coefficients", limits=c(0, 212)) +
  theme_bw() + 
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid=element_blank()
  )
ggsave(g, width=7, height=5.5, file="analyses/nmr_score_training/sig_coefficient_inclusion_matrix_non_derived_scores.pdf")

ggdt <- coef_stats[lambda.fit=="lambda.min", .(N_coef = .N), by=.(endpoint, sex, N_score_sig_incl_avg)]
g <- ggplot(ggdt) + 
  aes(x=factor(N_score_sig_incl_avg), y=N_coef) +
  facet_grid(sex ~ endpoint, drop=FALSE, shrink=FALSE) + 
  geom_col() +
  scale_x_discrete("Number of scores coefficient is included in after rounding to 3 decimal places") +
  scale_y_continuous("Number of coefficients", limits=c(0, 212)) +
  theme_bw() + 
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid=element_blank()
  )
ggsave(g, width=7, height=5.5, file="analyses/nmr_score_training/sig_avg_coefficient_inclusion_matrix_non_derived_scores.pdf")

# Do the same for clinical biomarkers
possible <- c("VLDL_C", "Total_TG", "ApoB", "ApoA1", "Omega_3", "Omega_6", "MUFA", "SFA", "DHA", 
  "Ala", "Gly", "His", "Ile", "Leu", "Val", "Phe", "Tyr", "Glucose", "Creatinine", "Albumin", "GlycA")
possible <- c(possible, paste0("age:", possible))
coef_stats <- expand.grid(endpoint=c("CAD", "Stroke"), sex=c("Male", "Female"), lambda.fit=c("lambda.min", "lambda.1se"), coef=possible, N_scores=0, N_scores_sig=0)
setDT(coef_stats)
coef_stats[nmr_coef[type == "clinical",.N,by=.(endpoint, sex, lambda.fit, coef)], on = .(endpoint, sex, lambda.fit, coef), N_scores := i.N]
coef_stats[nmr_coef[type == "clinical",.(N=sum(round(beta, digits=3) != 0)), by=.(endpoint, sex, lambda.fit, coef)], on = .(endpoint, sex, lambda.fit, coef), N_scores_sig := i.N]
coef_stats[, N_score_sig_incl_avg := N_scores_sig]
coef_stats[!avg_nmr_coef[type == "clinical" & round(beta, digits=3) != 0], on = .(endpoint, sex, lambda.fit, coef), N_score_sig_incl_avg := 0]

ggdt <- coef_stats[lambda.fit=="lambda.min", .(N_coef = .N), by=.(endpoint, sex, N_scores)]
g <- ggplot(ggdt) + 
  aes(x=factor(N_scores), y=N_coef) +
  facet_grid(sex ~ endpoint, drop=FALSE, shrink=FALSE) + 
  geom_col() +
  scale_x_discrete("Number of scores coefficient is included in") +
  scale_y_continuous("Number of coefficients", limits=c(0, 212)) +
  theme_bw() + 
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid=element_blank()
  )
ggsave(g, width=7, height=5.5, file="analyses/nmr_score_training/coefficient_inclusion_matrix_clinical_scores.pdf")

ggdt <- coef_stats[lambda.fit=="lambda.min", .(N_coef = .N), by=.(endpoint, sex, N_scores_sig)]
g <- ggplot(ggdt) + 
  aes(x=factor(N_scores_sig), y=N_coef) +
  facet_grid(sex ~ endpoint, drop=FALSE, shrink=FALSE) + 
  geom_col() +
  scale_x_discrete("Number of scores coefficient is included in after rounding to 3 decimal places") +
  scale_y_continuous("Number of coefficients", limits=c(0, 212)) +
  theme_bw() + 
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid=element_blank()
  )
ggsave(g, width=7, height=5.5, file="analyses/nmr_score_training/sig_coefficient_inclusion_matrix_clinical_scores.pdf")

ggdt <- coef_stats[lambda.fit=="lambda.min", .(N_coef = .N), by=.(endpoint, sex, N_score_sig_incl_avg)]
g <- ggplot(ggdt) + 
  aes(x=factor(N_score_sig_incl_avg), y=N_coef) +
  facet_grid(sex ~ endpoint, drop=FALSE, shrink=FALSE) + 
  geom_col() +
  scale_x_discrete("Number of scores coefficient is included in after rounding to 3 decimal places") +
  scale_y_continuous("Number of coefficients", limits=c(0, 212)) +
  theme_bw() + 
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid=element_blank()
  )
ggsave(g, width=7, height=5.5, file="analyses/nmr_score_training/sig_avg_coefficient_inclusion_matrix_clinical_scores.pdf")

# Plot consistency of coefficient estimates
coef_plot <- function(this_endpoint, this_sex, this_type, this_lambda="lambda.min", round=FALSE) {
  ggdt_ribbon <- avg_nmr_coef[endpoint == this_endpoint & sex == this_sex & lambda.fit == this_lambda & type == this_type]
  if (round) {
    ggdt_ribbon[, beta := round(beta, digits=3)]
    ggdt_ribbon <- ggdt_ribbon[beta != 0]
  } 

  ggdt_points <- nmr_coef[endpoint == this_endpoint & sex == this_sex & lambda.fit == this_lambda & type == this_type]
  if (round) {
    ggdt_points[, beta := round(beta, digits=3)]
    ggdt_points <- ggdt_points[beta != 0]
  } 

  if (round) {
    shared <- fintersect(
      ggdt_ribbon[,.(endpoint, sex, coef)],
      unique(ggdt_points[,.(endpoint, sex, coef)])
    )
    ggdt_ribbon <- ggdt_ribbon[shared, on = .(endpoint, sex, coef)]
    ggdt_points <- ggdt_points[shared, on = .(endpoint, sex, coef)]
  }
  
  ggdt_ribbon <- ggdt_ribbon[order(-beta)]
  ggdt_ribbon[, coef := factor(coef, levels=unique(coef))]
  ggdt_ribbon[, xorder := .I]

  ggdt_points[ggdt_ribbon, on = .(coef), xorder := i.xorder]
  ggdt_points[ggdt_ribbon, on = .(coef), discordant := sign(beta) != sign(i.beta)]

  ggplot(ggdt_points) +
    aes(x=xorder, y=exp(beta)) +
    geom_point(data=ggdt_points[!(discordant)], shape=19, fill="#525252", size=0.3) +
    geom_line(data=ggdt_ribbon, color="#6a51a3") +
    geom_hline(yintercept=1, linetype=2) +
    geom_point(data=ggdt_points[(discordant)], shape=1, color="#b30000", size=0.8) +
    scale_y_continuous("Hazard Ratio") +
    scale_x_continuous("", breaks=ggdt_ribbon$xorder, labels=ggdt_ribbon$coef, expand=expansion(mult=0, add=0.5)) +
    theme_bw() +
    theme(
      axis.text.y=element_text(size=6), axis.title.y=element_text(size=8), 
      axis.text.x=element_text(size=6, angle=90, hjust=1, vjust=0.5),
      axis.ticks.x=element_blank(),
      panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
      plot.margin=margin(10, 5.5, 5.5, 5.5)
    )
} 

g1 <- coef_plot("CAD", "Male", "non-derived")
g2 <- coef_plot("CAD", "Female", "non-derived")
g3 <- coef_plot("Stroke", "Male", "non-derived")
g4 <- coef_plot("Stroke", "Female", "non-derived")
g <- plot_grid(g1, g2, g3, g4, ncol=1, vjust=1.5, hjust=-0.2, label_size=8, align="hv",
  labels=c(
    "NMR score for CAD in males               ", 
    "NMR score for CAD in females             ", 
    "NMR score for ischaemic stroke in males  ", 
    "NMR score for ischaemic stroke in females"
  ))

ggsave(g, width=14, height=8, file="analyses/nmr_score_training/coefficient_consistency_non_derived.pdf")

g1 <- coef_plot("CAD", "Male", "non-derived", round=TRUE)
g2 <- coef_plot("CAD", "Female", "non-derived", round=TRUE)
g3 <- coef_plot("Stroke", "Male", "non-derived", round=TRUE)
g4 <- coef_plot("Stroke", "Female", "non-derived", round=TRUE)
g <- plot_grid(g1, g2, g3, g4, ncol=1, vjust=1.5, hjust=-0.2, label_size=8, align="hv",
  labels=c(
    "NMR score for CAD in males               ", 
    "NMR score for CAD in females             ", 
    "NMR score for ischaemic stroke in males  ", 
    "NMR score for ischaemic stroke in females"
  ))

ggsave(g, width=14, height=8, file="analyses/nmr_score_training/sig_coefficient_consistency_non_derived.pdf")

g1 <- coef_plot("CAD", "Male", "clinical")
g2 <- coef_plot("CAD", "Female", "clinical")
g3 <- coef_plot("Stroke", "Male", "clinical")
g4 <- coef_plot("Stroke", "Female", "clinical")
g <- plot_grid(g1, g2, g3, g4, ncol=1, vjust=1.5, hjust=-0.2, label_size=8, align="hv",
  labels=c(
    "NMR score for CAD in males               ", 
    "NMR score for CAD in females             ", 
    "NMR score for ischaemic stroke in males  ", 
    "NMR score for ischaemic stroke in females"
  ))

ggsave(g, width=14, height=8, file="analyses/nmr_score_training/coefficient_consistency_clinical.pdf")

g1 <- coef_plot("CAD", "Male", "clinical", round=TRUE)
g2 <- coef_plot("CAD", "Female", "clinical", round=TRUE)
g3 <- coef_plot("Stroke", "Male", "clinical", round=TRUE)
g4 <- coef_plot("Stroke", "Female", "clinical", round=TRUE)
g <- plot_grid(g1, g2, g3, g4, ncol=1, vjust=1.5, hjust=-0.2, label_size=8, align="hv",
  labels=c(
    "NMR score for CAD in males               ", 
    "NMR score for CAD in females             ", 
    "NMR score for ischaemic stroke in males  ", 
    "NMR score for ischaemic stroke in females"
  ))

ggsave(g, width=14, height=8, file="analyses/nmr_score_training/sig_coefficient_consistency_clinical.pdf")

# Make density plot of coefficients
ggdt <- nmr_coef[lambda.fit == "lambda.min" & type == "non-derived"]
ggdt[, sex := factor(paste0(sex, "s"), levels=c("Males", "Females"))]
ggdt[, endpoint := factor(endpoint, levels=c("CAD", "Stroke"))]
g <- ggplot(ggdt) +
  aes(x=beta, fill=paste("Score", prediction_cv_testfold), color=paste("Score", prediction_cv_testfold)) +
  facet_grid2(sex ~ endpoint, scales="free_x", independent="x") +
  geom_density() +
  geom_vline(xintercept=0, linetype=2) +
  scale_color_manual(values=c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")) +
  scale_fill_manual(values=paste0(c("#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4", "#fed9a6"), "dd")) +
  xlab("Log Hazard Ratio (elasticnet)") + 
  ylab("Density") +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8),
    legend.title=element_blank(), legend.text=element_text(size=6),
    strip.background=element_blank(), strip.text=element_text(size=8, face="bold")
  )
ggsave(g, width=7.2, height=6, file="analyses/nmr_score_training/coefficient_denisities_non_derived.pdf")

ggdt <- nmr_coef[lambda.fit == "lambda.min" & type == "clinical"]
ggdt[, sex := factor(paste0(sex, "s"), levels=c("Males", "Females"))]
ggdt[, endpoint := factor(endpoint, levels=c("CAD", "Stroke"))]
g <- ggplot(ggdt) +
  aes(x=beta, fill=paste("Score", prediction_cv_testfold), color=paste("Score", prediction_cv_testfold)) +
  facet_grid2(sex ~ endpoint, scales="free_x", independent="x") +
  geom_density() +
  geom_vline(xintercept=0, linetype=2) +
  scale_color_manual(values=c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")) +
  scale_fill_manual(values=paste0(c("#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4", "#fed9a6"), "dd")) +
  xlab("Log Hazard Ratio (elasticnet)") + 
  ylab("Density") +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8),
    legend.title=element_blank(), legend.text=element_text(size=6),
    strip.background=element_blank(), strip.text=element_text(size=8, face="bold")
  )
ggsave(g, width=7.2, height=6, file="analyses/nmr_score_training/coefficient_denisities_clinical.pdf")

# Write out table of coefficients for supp
dt <- avg_nmr_coef[lambda.fit == "lambda.min" & type == "non-derived"]
dt[nmr_scaling, on = .(coef=biomarker, sex), c("scaling_mean", "scaling_sd") := .(i.mean, i.sd)]
dt[, sex := factor(sex, levels=c("Male", "Female"))]
dt <- dcast(dt, coef ~ endpoint + sex, value.var=c("scaling_mean", "scaling_sd", "beta", "prop_r2"), fill=0)
dt <- dt[,.(coef, scaling_mean_Male=scaling_mean_CAD_Male, scaling_sd_Male=scaling_sd_CAD_Male, CAD_Male=beta_CAD_Male, Stroke_Male=beta_Stroke_Male, prop_r2_CAD_Male, prop_r2_Stroke_Male,
  scaling_mean_Female=scaling_mean_CAD_Female, scaling_sd_Female=scaling_sd_CAD_Female, CAD_Female=beta_CAD_Female, Stroke_Female=beta_Stroke_Female, prop_r2_CAD_Female, prop_r2_Stroke_Female)]
dt <- dt[order(tolower(coef))]
fwrite(dt, sep="\t", quote=FALSE, file="analyses/nmr_score_training/supp_table_coefficients.txt")
  
