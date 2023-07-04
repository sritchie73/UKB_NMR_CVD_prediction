library(data.table)
library(foreach)
library(ukbnmr)
library(ggplot2)
library(forcats)
library(ggrastr)
library(ggpp)
library(ggh4x)
library(cowplot)
library(survival)

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
g1 <- score_scatter(score_comp[endpoint == "CAD" & sex == "Male" & lambda.fit == "lambda.min"])
g2 <- score_scatter(score_comp[endpoint == "CAD" & sex == "Female" & lambda.fit == "lambda.min"])
g3 <- score_scatter(score_comp[endpoint == "Stroke" & sex == "Male" & lambda.fit == "lambda.min"])
g4 <- score_scatter(score_comp[endpoint == "Stroke" & sex == "Female" & lambda.fit == "lambda.min"])

ggsave(g1, width=7, height=6, file="analyses/nmr_score_training/CAD_male_score_consistency.pdf")
ggsave(g2, width=7, height=6, file="analyses/nmr_score_training/CAD_female_score_consistency.pdf")
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
ggsave(g, width=7, height=5.5, file="analyses/nmr_score_training/aggregate_test_vs_coef_average.pdf")

# Write out aggregate scores
test_scores <- nmr[lambda.fit == "lambda.min" & prediction_cv_testfold == prediction_cv_foldid]
test_scores <- dcast(test_scores, eid  ~ endpoint, value.var="linear_predictor")
setnames(test_scores, c("CAD", "Stroke"), c("CAD_NMR_score", "Stroke_NMR_score"))
fwrite(test_scores, sep="\t", quote=FALSE, file="analyses/nmr_score_training/aggregate_test_NMR_scores.txt")

# Also average scores for downstream checks
avg_scores <- dcast(avg_nmr[lambda.fit == "lambda.min"], eid ~ endpoint, value.var="linear_predictor")
setnames(avg_scores, c("CAD", "Stroke"), c("CAD_NMR_score", "Stroke_NMR_score"))
fwrite(avg_scores, sep="\t", quote=FALSE, file="analyses/nmr_score_training/coef_avg_NMR_scores.txt")

# For each score, how many coefficients are included across N training runs?
empty <- expand.grid(endpoint=c("CAD", "Stroke"), sex=c("Male", "Female"), scores=0:5)
setDT(empty)
empty[, coef := 0]
coef_stats <- nmr_coef[lambda.fit == "lambda.min", .(scores=.N), by=.(endpoint, sex, coef)]
coef_stats <- coef_stats[,.(coef=.N), by=.(endpoint, sex, scores)]
empty[coef_stats, on = .(endpoint, sex, scores), coef := i.coef]
coef_stats <- empty
coef_stats[, scores := factor(scores, levels=0:5)]

g <- ggplot(coef_stats) + 
  aes(x=scores, y=coef) +
  facet_grid(sex ~ endpoint, drop=FALSE, shrink=FALSE) + 
  geom_col() +
  scale_x_discrete("Number of scores coefficient is included in") +
  scale_y_continuous("Number of coefficients", limits=c(0, 214)) +
  theme_bw() + 
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid=element_blank()
  )
ggsave(g, width=7, height=5.5, file="analyses/nmr_score_training/coefficient_inclusion_matrix.pdf")

# Plot consistency of coefficient estimates
coef_plot <- function(this_endpoint, this_sex, this_lambda="lambda.min") {
  ggdt_ribbon <- avg_nmr_coef[endpoint == this_endpoint & sex == this_sex & lambda.fit == this_lambda]
  
  ggdt_ribbon <- ggdt_ribbon[order(-beta)]
  ggdt_ribbon[, coef := factor(coef, levels=unique(coef))]
  ggdt_ribbon[, xorder := .I]

  ggdt_points <- nmr_coef[endpoint == this_endpoint & sex == this_sex & lambda.fit == this_lambda]
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

g1 <- coef_plot("CAD", "Male")
g2 <- coef_plot("CAD", "Female")
g3 <- coef_plot("Stroke", "Male")
g4 <- coef_plot("Stroke", "Female")
g <- plot_grid(g1, g2, g3, g4, ncol=1, vjust=1.5, hjust=-0.2, label_size=8, align="hv",
  labels=c(
    "NMR score for CAD in males               ", 
    "NMR score for CAD in females             ", 
    "NMR score for ischaemic stroke in males  ", 
    "NMR score for ischaemic stroke in females"
  ))

ggsave(g, width=14, height=8, file="analyses/nmr_score_training/coefficient_consistency.pdf")
  
