library(data.table)
library(foreach)
library(ukbnmr)
library(ggplot2)
library(forcats)
library(ggrastr)
library(ggpp)
library(ggh4x)
library(hexbin)
library(cowplot)
source("src/utils/SCORE2.R")

options("ggrastr.default.dpi" = 300) 

# Load dataset with NMR biomarkers imputed
dat <- fread("data/imputed/analysis_cohort.txt")

# Load scaling factors for NMR data
nmr_scaling <- fread("data/standardised/nmr_scaling_factors.txt")

# Load NMR score coefficients
nmr_coef <- fread("analyses/nmr_score_training/best_fits_coef.txt")

# Build long table of NMR coefficients and age interactions standardised using the given scaling factors
nmr_long <- melt(dat, id.vars=c("eid", "sex", "age", "prediction_cv_foldid"), measure.vars=ukbnmr::nmr_info[Type == "Non-derived", Biomarker])
nmr_long[nmr_scaling, on = .(variable=biomarker), value := (value - i.mean)/i.sd]
nmr_long[, age := (age - 60)/5]
nmr_age_interactions <- copy(nmr_long)
nmr_age_interactions[, value := value * age]
nmr_age_interactions[, variable := paste0("age:", variable)]
nmr_long <- rbind(nmr_long, nmr_age_interactions)

# Compute NMR biomarker scores
models <- unique(nmr_coef[,.(prediction_cv_testfold, endpoint, lambda.fit)])
nmr <- foreach(midx = models[,.I], .combine=rbind) %do% {
  this_info <- cbind(score="NMR", models[midx])
  this_coef <- nmr_coef[this_info, on = .(prediction_cv_testfold, endpoint, lambda.fit)]
  this_score <- nmr_long[this_coef, on = .(sex, variable=coef)][, .(linear_predictor=sum(beta * value)), by=.(eid, sex, prediction_cv_foldid)]
  cbind(this_info, this_score)
}

# Compute average score across all 5 training runs
avg_nmr_coef <- fread("analyses/nmr_score_training/avg_best_coefs.txt")
models <- unique(models[,.(endpoint, lambda.fit)])
avg_nmr <- foreach(midx = models[,.I], .combine=rbind) %do% {
  this_info <- cbind(score="NMR", models[midx])
  this_coef <- avg_nmr_coef[this_info, on = .(endpoint, lambda.fit)]
  this_score <- nmr_long[this_coef, on = .(sex, variable=coef)][, .(linear_predictor=sum(beta * value)), by=.(eid, sex, prediction_cv_foldid)]
  cbind(this_info, this_score)
}

# Assess consistency of NMR biomarker scores across training runs
score_comp <- foreach(this_test_fold = c(1:5, "avg"), .combine=rbind) %do% {
  if (this_test_fold == "avg") {
    this_x <- avg_nmr[, .(eid, prediction_cv_foldid, sex, endpoint, lambda.fit,
        x_model="avg", x_score=linear_predictor)]
  } else {
    this_x <- nmr[prediction_cv_testfold == this_test_fold,  
      .(eid, prediction_cv_foldid, sex, endpoint, lambda.fit, 
        x_model=prediction_cv_testfold, x_score=linear_predictor)]
  } 
  this_y <- nmr[, .(eid, prediction_cv_foldid, sex, endpoint, lambda.fit,
      y_model=prediction_cv_testfold, y_score=linear_predictor)]
  this_y <- rbind(this_y, avg_nmr[,
    .(eid, prediction_cv_foldid, sex, endpoint, lambda.fit,
      y_model="avg", y_score=linear_predictor)])
  merge(this_x, this_y, by=c("eid", "prediction_cv_foldid", "sex", "endpoint", "lambda.fit"), all=FALSE)
}

# Color annotation for plots
score_comp[, color_anno := fcase(
  x_model != prediction_cv_foldid & y_model != prediction_cv_foldid, "Training samples",
  x_model == prediction_cv_foldid, "X-axis test samples", 
  y_model == prediction_cv_foldid, "Y-axis test samples" 
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
    rasterise(geom_point(data=score_comp[x_model != prediction_cv_foldid & y_model != prediction_cv_foldid], alpha=0.5, shape=19, size=0.1)) +
    rasterise(geom_point(data=score_comp[x_model == prediction_cv_foldid & y_model != prediction_cv_foldid], alpha=0.5, shape=19, size=0.1)) +
    rasterise(geom_point(data=score_comp[x_model != prediction_cv_foldid & y_model == prediction_cv_foldid], alpha=0.5, shape=19, size=0.1)) +
    rasterise(geom_point(data=score_comp[x_model == prediction_cv_foldid & y_model == prediction_cv_foldid], alpha=0.5, shape=19, size=0.1)) +
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
g1 <- score_scatter(score_comp[endpoint == "CHD" & sex == "Male" & lambda.fit == "lambda.min"])
g2 <- score_scatter(score_comp[endpoint == "CHD" & sex == "Female" & lambda.fit == "lambda.min"])
g3 <- score_scatter(score_comp[endpoint == "Stroke" & sex == "Male" & lambda.fit == "lambda.min"])
g4 <- score_scatter(score_comp[endpoint == "Stroke" & sex == "Female" & lambda.fit == "lambda.min"])

ggsave(g1, width=7, height=6, file="analyses/nmr_score_training/CHD_male_score_consistency.pdf")
ggsave(g2, width=7, height=6, file="analyses/nmr_score_training/CHD_female_score_consistency.pdf")
ggsave(g3, width=7, height=6, file="analyses/nmr_score_training/Stroke_male_score_consistency.pdf")
ggsave(g4, width=7, height=6, file="analyses/nmr_score_training/Stroke_female_score_consistency.pdf")

# Compare average score to test score
avg_vs_test_scatter <- function(lambda) {
  cor_anno <- score_comp[x_model == "avg" & y_model == prediction_cv_foldid & lambda.fit == lambda,
    .(label=sprintf("r=%.2f\np=%.2f",
      cor(x_score, y_score),
      cor(x_score, y_score, method="spearman")
    )), by=.(endpoint, sex)]

  ggplot(score_comp[x_model == "avg" & y_model == prediction_cv_foldid & lambda.fit == lambda]) +
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

g <- avg_vs_test_scatter("lambda.min")
ggsave(g, width=7, height=5.5, file="analyses/nmr_score_training/test_aggregate_vs_coef_average.pdf")

# Create dataset with optimal scores
system("mkdir -p data/scores/")
dat[, SCORE2 := score2(sex, age, smoking, sbp, tchol, hdl, type="linear predictor")]
test_scores <- nmr[lambda.fit == "lambda.min" & prediction_cv_testfold == prediction_cv_foldid]
dat[test_scores[endpoint == "CHD"], on = .(eid), CHD_NMR_score := i.linear_predictor]
dat[test_scores[endpoint == "Stroke"], on = .(eid), Stroke_NMR_score := i.linear_predictor]
fwrite(dat, sep="\t", quote=FALSE, file="data/scores/analysis_cohort.txt")

