library(data.table)
library(foreach)
library(forcats)
library(ggplot2)
library(scales)
library(ggstance)
library(patchwork)

# Load and collate results
res <- rbind(idcol="cohort",
  "discovery"=fread("analyses/public_health_modelling/discovery_screening.txt"),
  "replication"=fread("analyses/public_health_modelling/replication_screening.txt"),
  "pooled"=fread("analyses/public_health_modelling/pooled_screening.txt")
)

# Build Figure 4
ggdt <- res[cohort == "pooled" & endpoint == "cvd" & score == "SCORE2_excl_UKB" & model_sex == "Sex-stratified" & model != "Risk score"]
ggdt <- ggdt[metric %in% c("delta_high_risk", "delta_cvd_high_risk", "delta_cvd_prevented", "delta_NNS", "delta_NNT")]
ggdt[, metric := factor(metric, levels=c("delta_high_risk", "delta_cvd_high_risk", "delta_cvd_prevented", "delta_NNS", "delta_NNT"))]

ggdt[, model_name := fcase(
  model == "Risk score + NMR scores", "SCORE2 + NMR scores for CHD and IS",
  model == "Risk score + clinical biomarkers", "SCORE2 + 11 clinical chemistry biomarkers",
  model == "Risk score + PRSs", "SCORE2 + PRSs for CHD and IS",
  model == "Risk score + NMR scores + PRSs", "SCORE2 + NMR scores for CHD and IS + PRSs for CHD and IS",
  model == "Risk score + clinical biomarkers + PRSs", "SCORE2 + 11 clinical chemistry biomarkers + PRSs for CHD and IS"
)]
ggdt[,model_name := factor(model_name, levels=c("SCORE2 + NMR scores for CHD and IS", "SCORE2 + 11 clinical chemistry biomarkers",
  "SCORE2 + PRSs for CHD and IS", "SCORE2 + NMR scores for CHD and IS + PRSs for CHD and IS",
  "SCORE2 + 11 clinical chemistry biomarkers + PRSs for CHD and IS"))]

# Requires some assembly so we can control axes
screening_plot <- function(ggdt) {
	ggplot(ggdt) +
		aes(x=estimate, xmin=L95, xmax=U95, y=fct_rev(model_name), color=metric) +
		geom_vline(xintercept=0, linetype=2) +
		geom_errorbarh(height=0) +
		geom_point(shape=23, fill="white", size=1.2) +
		scale_color_manual(values=c(
			"delta_high_risk"="#f8766d",
			"delta_cvd_high_risk"="#a3a500",
			"delta_cvd_prevented"="#00bf7d",
			"delta_NNS"="#00b0f6",
			"delta_NNT"="#e76bf3"
		)) +    
		xlab("Change relative to population-wide screening with SCORE2 alone (95% CI)") +
		theme_bw() +
		theme(
			axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
			axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
			strip.background=element_blank(), strip.text=element_blank(),
			panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
			legend.position="none", plot.margin=unit(c(1, 1, 1, 1), units="pt")
		)
}

g1 <- screening_plot(ggdt[metric == "delta_high_risk" & strategy == "blanket"]) + scale_x_continuous(limits=c(0, 4700), breaks=c(0, 2000, 4000))
g2 <- screening_plot(ggdt[metric == "delta_cvd_high_risk" & strategy == "blanket"]) + scale_x_continuous(limits=c(0, 1200), breaks=c(0, 500, 1000))
g3 <- screening_plot(ggdt[metric == "delta_cvd_prevented" & strategy == "blanket"]) + scale_x_continuous(limits=c(0, 110), breaks=c(0, 50, 100))
g4 <- screening_plot(ggdt[metric == "delta_NNS" & strategy == "blanket"]) + scale_x_continuous(limits=c(NA, 0), breaks=c(-400, -200, 0))
g5 <- screening_plot(ggdt[metric == "delta_NNT" & strategy == "blanket"]) + scale_x_continuous(limits=c(-3, 3), breaks=c(-2, 0, 2))
g6 <- screening_plot(ggdt[metric == "delta_high_risk" & strategy == "targeted"]) + scale_x_continuous(limits=c(0, 4700), breaks=c(0, 2000, 4000))
g7 <- screening_plot(ggdt[metric == "delta_cvd_high_risk" & strategy == "targeted"]) + scale_x_continuous(limits=c(0, 1200), breaks=c(0, 500, 1000))
g8 <- screening_plot(ggdt[metric == "delta_cvd_prevented" & strategy == "targeted"]) + scale_x_continuous(limits=c(0, 110), breaks=c(0, 50, 100))
g9 <- screening_plot(ggdt[metric == "delta_NNS" & strategy == "targeted"]) + scale_x_continuous(limits=c(NA, 0), breaks=c(-400, -200, 0))
g10 <- screening_plot(ggdt[metric == "delta_NNT" & strategy == "targeted"]) + scale_x_continuous(limits=c(-3, 3), breaks=c(-2, 0, 2))
g <- g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8 + g9 + g10 + plot_layout(axes="collect", nrow=2, byrow=TRUE)

ggsave(g, width=7.2, height=3.5, file="analyses/public_health_modelling/screening_comparison.pdf")

# Build supp table for the main results
supp1 <- res[endpoint == "cvd" & score == "SCORE2_excl_UKB" & model_sex == "Sex-stratified" & model != "Risk score"]
supp1[, model := gsub("Risk score", "SCORE2", model)]
supp1[, metric_type := fcase(
  metric %like% "^ref_", "ref",
  metric %like% "^alt_", "alt",
  metric %like% "^delta_", "delta"
)]
supp1 <- supp1[!is.na(metric_type)]
supp1[, metric := gsub("(^ref_)|(^alt_)|(^delta_)", "", metric)]
supp1 <- dcast(supp1, strategy + model + cohort + metric ~ metric_type, value.var=c("estimate", "se", "L95", "U95", "pval"))
supp1 <- supp1[, .(strategy, model, cohort, metric, 
  estimate_ref, se_ref, L95_ref, U95_ref,
  estimate_alt, se_alt, L95_alt, U95_alt,
  estimate_delta, se_delta, L95_delta, U95_delta, pval_delta)]
supp1[, strategy := ifelse(strategy == "blanket", "Population-wide", "Targeted")]
supp1[, cohort := fcase(
  cohort == "discovery", "Discovery",
  cohort == "replication", "Replication",
  cohort == "pooled", "Discovery + Replication")]
supp1[, cohort := factor(cohort, levels=c("Discovery", "Replication", "Discovery + Replication"))]
supp1[, model := factor(model, levels=c(
  "SCORE2 + NMR scores",
  "SCORE2 + clinical biomarkers",
  "SCORE2 + PRSs",
  "SCORE2 + NMR scores + PRSs",
  "SCORE2 + clinical biomarkers + PRSs"))]
supp1[, metric := fcase(
  metric == "high_risk", "People classified as high risk",
  metric == "cvd_high_risk", "Incident CVD events in the high-risk group",
  metric == "cvd_prevented", "Expected CVD events prevented by statin initiation",
  metric == "NNS", "Number needed to screen to prevent one CVD event (NNS)",
  metric == "NNT", "Number of statins prescribed per CVD event prevented (NNT)"
)]
supp1[, metric := factor(metric, levels=c(
  "People classified as high risk",
  "Incident CVD events in the high-risk group",
  "Expected CVD events prevented by statin initiation",
  "Number needed to screen to prevent one CVD event (NNS)",
  "Number of statins prescribed per CVD event prevented (NNT)"))]
supp1 <- supp1[order(cohort)][order(metric)][order(model)][order(strategy)]
fwrite(supp1, sep="\t", quote=FALSE, file="analyses/public_health_modelling/screening_comparison_main_analysis.txt")

# Build supp table for sensitivity analyses
supp2 <- res[score != "SCORE2" & model != "Risk score"]
supp2[, metric_type := fcase(
  metric %like% "^ref_", "ref",
  metric %like% "^alt_", "alt",
  metric %like% "^delta_", "delta"
)]
supp2 <- supp2[!is.na(metric_type)]
supp2[, metric := gsub("(^ref_)|(^alt_)|(^delta_)", "", metric)]
supp2 <- dcast(supp2, strategy + model + cohort + model_sex + score + endpoint + metric ~ metric_type, value.var=c("estimate", "se", "L95", "U95", "pval"))
supp2 <- supp2[, .(strategy, model, cohort, model_sex, score, endpoint, metric,
  estimate_ref, se_ref, L95_ref, U95_ref,
  estimate_alt, se_alt, L95_alt, U95_alt,
  estimate_delta, se_delta, L95_delta, U95_delta, pval_delta)]
supp2[, strategy := ifelse(strategy == "blanket", "Population-wide", "Targeted")]
supp2[, cohort := fcase(
  cohort == "discovery", "Discovery",
  cohort == "replication", "Replication",
  cohort == "pooled", "Discovery + Replication")]
supp2[, cohort := factor(cohort, levels=c("Discovery", "Replication", "Discovery + Replication"))]
supp2[, model := factor(model, levels=c(
  "Risk score + NMR scores",
  "Risk score + clinical biomarkers",
  "Risk score + PRSs",
  "Risk score + NMR scores + PRSs",
  "Risk score + clinical biomarkers + PRSs"))]
supp2[, metric := fcase(
  metric == "high_risk", "People classified as high risk",
  metric == "cvd_high_risk", "Incident CVD events in the high-risk group",
  metric == "cvd_prevented", "Expected CVD events prevented by statin initiation",
  metric == "NNS", "Number needed to screen to prevent one CVD event (NNS)",
  metric == "NNT", "Number of statins prescribed per CVD event prevented (NNT)"
)]
supp2[, metric := factor(metric, levels=c(
  "People classified as high risk",
  "Incident CVD events in the high-risk group",
  "Expected CVD events prevented by statin initiation",
  "Number needed to screen to prevent one CVD event (NNS)",
  "Number of statins prescribed per CVD event prevented (NNT)"))]
supp2[, model_sex := factor(model_sex, levels=c("Sex-stratified", "Males", "Females"))]
supp2[score == "SCORE2_excl_UKB", score := "SCORE2"]
supp2[, score := factor(score, levels=c("SCORE2", "QRISK3"))]
supp2[, endpoint := fcase(
  endpoint == "cvd", "Broad",
  endpoint == "cvd_narrow", "Narrow")]
supp2 <- supp2[order(score)][order(model_sex)][order(cohort)][order(metric)][order(model)][order(strategy)]
fwrite(supp2, sep="\t", quote=FALSE, file="analyses/public_health_modelling/screening_comparison_sensitivity_analysis.txt")

# Run pairwise comparisons to assess similarity of results across different settings
##model_info <- unique(res[metric %like% "^delta_",.(strategy, cohort, endpoint, model_sex, score, metric)])
##comp_stats <- foreach(modelIdx1=model_info[,.I], .combine=rbind) %do% {
##  foreach(modelIdx2=model_info[,.I], .combine=rbind) %do% {
##   if (modelIdx1 >= modelIdx2) return(NULL)
##
##    model_info1 <- model_info[modelIdx1]
##    model_info2 <- model_info[modelIdx2]
##
##    res1 <- res[model_info1, on=.(strategy, cohort, endpoint, model_sex, score, metric)]
##    res2 <- res[model_info2, on=.(strategy, cohort, endpoint, model_sex, score, metric)]
##
##    res1 <- res1[, .(model, estimate)]
##    res2 <- res2[, .(model, estimate)]
##
##    comp <- merge(res1, res2, by="model", suffixes=c(".1", ".2"))
##
##    corr <- comp[,cor(estimate.1, estimate.2)]
##    beta <- comp[,coef(lm(estimate.2 ~ 0 + estimate.1))[1]]
##
##    this_comp_stats <- data.table(pearson=corr, lm_beta=beta)
##
##    setnames(model_info1, paste0(names(model_info1), ".1"))
##    setnames(model_info2, paste0(names(model_info2), ".2"))
##    return(cbind(model_info1, model_info2, this_comp_stats))
##  }
##}
##fwrite(comp_stats, sep="\t", quote=FALSE, file="analyses/public_health_modelling/comparison_sensitivity.txt")

# Compare the five metrics across the discovery, replication, and pooled cohorts
cohort_comp <- res[endpoint == "cvd" & model_sex == "Sex-stratified" & score == "SCORE2_excl_UKB" & model != "Risk score"]
cohort_comp <- dcast(cohort_comp, strategy + model + metric ~ cohort, value.var=c("estimate", "L95", "U95"))

common_gg_parts <- function(g, limits, breaks) {
  g + aes(shape=strategy, fill=model, color=model) +
    geom_abline(intercept=0, slope=1, linetype=2) +
    geom_hline(yintercept=0, linetype=3) +
    geom_vline(xintercept=0, linetype=3) +
    geom_errorbarh(height=0, alpha=0.7) + geom_errorbar(width=0, alpha=0.7) +
    geom_point(fill="white") +
    scale_shape_manual(values=c("blanket"=23, "targeted"=22)) +
    scale_color_manual(values=c(
      "Risk score + NMR scores"="#e41a1c", "Risk score + PRSs"="#377eb8", "Risk score + clinical biomarkers"="#ff7f00",
      "Risk score + NMR scores + PRSs"="#4daf4a", "Risk score + clinical biomarkers + PRSs"="#984ea3"
    )) +
    scale_x_continuous(limits=limits, breaks=breaks, labels=comma) +
    scale_y_continuous(limits=limits, breaks=breaks, labels=comma) +
    theme_bw() +
    theme(
      axis.text=element_text(size=6), axis.title=element_text(size=7),
      axis.text.x=element_text(size=6), plot.title=element_text(size=8, face="bold"), 
      legend.position="none"
    )
}

g1 <- ggplot(cohort_comp[metric == "delta_high_risk"]) + aes(
    x=estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=estimate_pooled, ymin=L95_pooled, ymax=U95_pooled, color=model
  ) + xlab("Discovery cohort") + ylab("Pooled cohort") + ggtitle("High Risk")
g1 <- common_gg_parts(g1, c(0, 4500), c(0, 2000, 4000))

g2 <- ggplot(cohort_comp[metric == "delta_high_risk"]) + aes(
    x=estimate_replication, xmin=L95_replication, xmax=U95_replication,
    y=estimate_pooled, ymin=L95_pooled, ymax=U95_pooled, color=model
  ) + xlab("Replication cohort") + ylab("Pooled cohort")
g2 <- common_gg_parts(g2, c(0, 4500), c(0, 2000, 4000))

g3 <- ggplot(cohort_comp[metric == "delta_high_risk"]) + aes(
    x=estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=estimate_replication, ymin=L95_replication, ymax=U95_replication, color=model
  ) + xlab("Discovery cohort") + ylab("Replication cohort") 
g3 <- common_gg_parts(g3, c(0, 4500), c(0, 2000, 4000))


g4 <- ggplot(cohort_comp[metric == "delta_cvd_high_risk"]) + aes(
    x=estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=estimate_pooled, ymin=L95_pooled, ymax=U95_pooled, color=model
  ) + xlab("Discovery cohort") + ylab("Pooled cohort") + ggtitle("CVD high")
g4 <- common_gg_parts(g4, c(0, 1100), c(0, 500, 1000))

g5 <- ggplot(cohort_comp[metric == "delta_cvd_high_risk"]) + aes(
    x=estimate_replication, xmin=L95_replication, xmax=U95_replication,
    y=estimate_pooled, ymin=L95_pooled, ymax=U95_pooled, color=model
  ) + xlab("Replication cohort") + ylab("Pooled cohort")
g5 <- common_gg_parts(g5, c(0, 1100), c(0, 500, 1000))

g6 <- ggplot(cohort_comp[metric == "delta_cvd_high_risk"]) + aes(
    x=estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=estimate_replication, ymin=L95_replication, ymax=U95_replication, color=model
  ) + xlab("Discovery cohort") + ylab("Replication cohort") 
g6 <- common_gg_parts(g6, c(0, 1100), c(0, 500, 1000))


g7 <- ggplot(cohort_comp[metric == "delta_cvd_prevented"]) + aes(
    x=estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=estimate_pooled, ymin=L95_pooled, ymax=U95_pooled, color=model
  ) + xlab("Discovery cohort") + ylab("Pooled cohort") + ggtitle("Prevented")
g7 <- common_gg_parts(g7, c(0, 110), c(0, 50, 100))

g8 <- ggplot(cohort_comp[metric == "delta_cvd_prevented"]) + aes(
    x=estimate_replication, xmin=L95_replication, xmax=U95_replication,
    y=estimate_pooled, ymin=L95_pooled, ymax=U95_pooled, color=model
  ) + xlab("Replication cohort") + ylab("Pooled cohort")
g8 <- common_gg_parts(g8, c(0, 110), c(0, 50, 100))

g9 <- ggplot(cohort_comp[metric == "delta_cvd_prevented"]) + aes(
    x=estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=estimate_replication, ymin=L95_replication, ymax=U95_replication, color=model
  ) + xlab("Discovery cohort") + ylab("Replication cohort") 
g9 <- common_gg_parts(g9, c(0, 110), c(0, 50, 100))


g10 <- ggplot(cohort_comp[metric == "delta_NNS"]) + aes(
    x=estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=estimate_pooled, ymin=L95_pooled, ymax=U95_pooled, color=model
  ) + xlab("Discovery cohort") + ylab("Pooled cohort") + ggtitle("NNS")
g10 <- common_gg_parts(g10, c(-550, 0), c(-400, -200, 0))

g11 <- ggplot(cohort_comp[metric == "delta_NNS"]) + aes(
    x=estimate_replication, xmin=L95_replication, xmax=U95_replication,
    y=estimate_pooled, ymin=L95_pooled, ymax=U95_pooled, color=model
  ) + xlab("Replication cohort") + ylab("Pooled cohort")
g11 <- common_gg_parts(g11, c(-550, 0), c(-400, -200, 0))

g12 <- ggplot(cohort_comp[metric == "delta_NNS"]) + aes(
    x=estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=estimate_replication, ymin=L95_replication, ymax=U95_replication, color=model
  ) + xlab("Discovery cohort") + ylab("Replication cohort") 
g12 <- common_gg_parts(g12, c(-550, 0), c(-400, -200, 0))


g13 <- ggplot(cohort_comp[metric == "delta_NNT"]) + aes(
    x=estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=estimate_pooled, ymin=L95_pooled, ymax=U95_pooled, color=model
  ) + xlab("Discovery cohort") + ylab("Pooled cohort") + ggtitle("NNT")
g13 <- common_gg_parts(g13, c(-4, 4), c(-2, 0, 2))

g14 <- ggplot(cohort_comp[metric == "delta_NNT"]) + aes(
    x=estimate_replication, xmin=L95_replication, xmax=U95_replication,
    y=estimate_pooled, ymin=L95_pooled, ymax=U95_pooled, color=model
  ) + xlab("Replication cohort") + ylab("Pooled cohort")
g14 <- common_gg_parts(g14, c(-4, 4), c(-2, 0, 2))

g15 <- ggplot(cohort_comp[metric == "delta_NNT"]) + aes(
    x=estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=estimate_replication, ymin=L95_replication, ymax=U95_replication, color=model
  ) + xlab("Discovery cohort") + ylab("Replication cohort") 
g15 <- common_gg_parts(g15, c(-4, 4), c(-2, 0, 2))


g <- g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8 + g9 + g10 + g11 + g12 + g13 + g14 + g15 + plot_layout(nrow=3, byrow=FALSE)
ggsave(g, width=7.2, height=4.2, file="analyses/public_health_modelling/cohort_sensitivity.pdf")


