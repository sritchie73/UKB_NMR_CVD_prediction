library(data.table)
library(foreach)
library(forcats)
library(ggplot2)
library(scales)
library(ggstance)
library(patchwork)

# Load and collate results
blanket <- rbind(idcol="cohort",
  "discovery"=fread("analyses/public_health_modelling/discovery_blanket_screening.txt"),
  "replication"=fread("analyses/public_health_modelling/replication_blanket_screening.txt"),
  "pooled"=fread("analyses/public_health_modelling/pooled_blanket_screening.txt")
)

targeted <- rbind(idcol="cohort",
  "discovery"=fread("analyses/public_health_modelling/discovery_targeted_screening.txt"),
  "replication"=fread("analyses/public_health_modelling/replication_targeted_screening.txt"),
  "pooled"=fread("analyses/public_health_modelling/pooled_targeted_screening.txt")
)

res <- rbind(blanket, targeted)

# Build Figure 5
ggdt <- res[cohort == "pooled" & endpoint == "cvd" & score == "SCORE2_excl_UKB" & model_sex == "Sex-stratified" & model != "Risk score"]

ggdt[, model_name := fcase(
  model == "Risk score + NMR scores", "SCORE2 + NMR scores for CHD and IS",
  model == "Risk score + clinical biomarkers", "SCORE2 + 11 clinical chemistry biomarkers",
  model == "Risk score + PRSs", "SCORE2 + PRSs for CHD and IS",
  model == "Risk score + NMR scores + PRSs", "SCORE2 + NMR scores + PRSs",
  model == "Risk score + clinical biomarkers + PRSs", "SCORE2 + 11 clinical biomarkers + PRSs",
  model == "Risk score + NMR scores + clinical biomarkers", "SCORE2 + NMR scores + 11 clinical biomarkers",
  model == "Risk score + NMR scores + clinical biomarkers + PRSs", "SCORE2 + NMR scores + 11 clinical biomarkers + PRSs"
)]
ggdt[,model_name := factor(model_name, levels=c("SCORE2 + PRSs for CHD and IS", "SCORE2 + NMR scores for CHD and IS",
  "SCORE2 + 11 clinical chemistry biomarkers", "SCORE2 + NMR scores + 11 clinical biomarkers", "SCORE2 + NMR scores + PRSs",
  "SCORE2 + 11 clinical biomarkers + PRSs", "SCORE2 + NMR scores + 11 clinical biomarkers + PRSs"))]

# Requires some assembly so we can control axes
screening_plot <- function(ggdt) {
	ggplot(ggdt) +
		aes(x=estimate, xmin=L95, xmax=U95, y=fct_rev(model_name), color=metric) +
		geom_vline(xintercept=0, linetype=2) +
		geom_errorbarh(height=0) +
		geom_point(shape=23, fill="white", size=1.2) +
		scale_color_manual(values=c(
			"delta_high_risk"="#f8766d", "alt_delta_ref_high_risk"="#f8766d", "alt_delta_null_high_risk"="#f8766d", 
			"delta_cvd_high_risk"="#a3a500", "alt_delta_ref_cvd_high_risk"="#a3a500", "alt_delta_null_cvd_high_risk"="#a3a500",
			"delta_cvd_prevented"="#00bf7d", "alt_delta_ref_cvd_prevented"="#00bf7d", "alt_delta_null_cvd_prevented"="#00bf7d",
			"delta_NNS"="#00b0f6", "alt_delta_ref_NNS"="#00b0f6", "alt_delta_null_NNS"="#00b0f6",
			"delta_NNT"="#e76bf3", "alt_delta_ref_NNT"="#e76bf3", "alt_delta_null_NNT"="#e76bf3"
		)) +    
		xlab("Change relative to SCORE2 alone (95% CI)") +
		theme_bw() +
		theme(
			axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
			axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
			strip.background=element_blank(), strip.text=element_blank(),
			panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
			legend.position="none", plot.margin=unit(c(1, 1, 1, 1), units="pt")
		)
}

g1 <- screening_plot(ggdt[metric == "delta_high_risk" & strategy == "blanket"]) + scale_x_continuous(limits=c(-200, 4700), breaks=c(0, 2000, 4000))
g2 <- screening_plot(ggdt[metric == "delta_cvd_high_risk" & strategy == "blanket"]) + scale_x_continuous(limits=c(0, 1200), breaks=c(0, 500, 1000))
g3 <- screening_plot(ggdt[metric == "delta_cvd_prevented" & strategy == "blanket"]) + scale_x_continuous(limits=c(0, 220), breaks=c(0, 100, 200))
g4 <- screening_plot(ggdt[metric == "delta_NNS" & strategy == "blanket"]) + scale_x_continuous(limits=c(-275, 0), breaks=c(-200, -100, 0))
g5 <- screening_plot(ggdt[metric == "delta_NNT" & strategy == "blanket"]) + scale_x_continuous(limits=c(-4.25, 1.5), breaks=c(-4, -2, 0))

g6 <- screening_plot(ggdt[metric == "alt_delta_ref_high_risk" & strategy == "targeted"]) + scale_x_continuous(limits=c(-200, 4700), breaks=c(0, 2000, 4000))
g7 <- screening_plot(ggdt[metric == "alt_delta_ref_cvd_high_risk" & strategy == "targeted"]) + scale_x_continuous(limits=c(0, 1200), breaks=c(0, 500, 1000))
g8 <- screening_plot(ggdt[metric == "alt_delta_ref_cvd_prevented" & strategy == "targeted"]) + scale_x_continuous(limits=c(0, 220), breaks=c(0, 100, 200))
g9 <- screening_plot(ggdt[metric == "alt_delta_ref_NNS" & strategy == "targeted"]) + scale_x_continuous(limits=c(-275, 0), breaks=c(-200, -100, 0))
g10 <- screening_plot(ggdt[metric == "alt_delta_ref_NNT" & strategy == "targeted"]) + scale_x_continuous(limits=c(-4.25, 1.5), breaks=c(-4, -2, 0))

g11 <- screening_plot(ggdt[metric == "alt_delta_null_high_risk" & strategy == "targeted"]) + scale_x_continuous(limits=c(-200, 4700), breaks=c(0, 2000, 4000))
g12 <- screening_plot(ggdt[metric == "alt_delta_null_cvd_high_risk" & strategy == "targeted"]) + scale_x_continuous(limits=c(0, 1200), breaks=c(0, 500, 1000))
g13 <- screening_plot(ggdt[metric == "alt_delta_null_cvd_prevented" & strategy == "targeted"]) + scale_x_continuous(limits=c(0, 220), breaks=c(0, 100, 200))
g14 <- screening_plot(ggdt[metric == "alt_delta_null_NNS" & strategy == "targeted"]) + scale_x_continuous(limits=c(-275, 0), breaks=c(-200, -100, 0))
g15 <- screening_plot(ggdt[metric == "alt_delta_null_NNT" & strategy == "targeted"]) + scale_x_continuous(limits=c(-4.25, 1.5), breaks=c(-4, -2, 0))

g <- g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8 + g9 + g10 + g11 + g12 + g13 + g14 + g15 + plot_layout(axes="collect", nrow=3, byrow=TRUE)
ggsave(g, width=7.2, height=4.425, file="analyses/public_health_modelling/screening_comparison.pdf")

# Build supp for main analyses - population wide screening
supp1 <- res[endpoint == "cvd" & score == "SCORE2_excl_UKB" & model_sex == "Sex-stratified" & model != "Risk score" & strategy == "blanket"]
supp1[, model := gsub("Risk score", "SCORE2", model)]
supp1[metric %like% "^delta", metric := gsub("delta", "alt_delta_ref", metric)]
supp1[, metric_type := fcase(
  metric %like% "^alt_delta_ref_", "alt_delta_ref",
  metric %like% "^ref_", "ref", 
  metric %like% "^alt_", "alt"
)]
supp1 <- supp1[!is.na(metric_type)]
supp1[, metric := gsub("(^alt_delta_ref_)|(^ref_)|(^alt_)", "", metric)]
supp1 <- dcast(supp1, strategy + model + cohort + metric ~ metric_type, value.var=c("estimate", "se", "L95", "U95", "pval"))
supp1 <- supp1[, .(strategy, model, cohort, metric, 
  estimate_ref, se_ref, L95_ref, U95_ref,
  estimate_alt, se_alt, L95_alt, U95_alt,
  estimate_alt_delta_ref, se_alt_delta_ref, L95_alt_delta_ref, U95_alt_delta_ref, pval_alt_delta_ref
)]
supp1[, strategy := ifelse(strategy == "blanket", "Population-wide", "Targeted")]
supp1[, cohort := fcase(
  cohort == "discovery", "Discovery",
  cohort == "replication", "Replication",
  cohort == "pooled", "Discovery + Replication")]
supp1[, cohort := factor(cohort, levels=c("Discovery", "Replication", "Discovery + Replication"))]
supp1[, model := factor(model, levels=c(
  "SCORE2 + PRSs",
  "SCORE2 + NMR scores",
  "SCORE2 + clinical biomarkers",
  "SCORE2 + NMR scores + clinical biomarkers",
  "SCORE2 + NMR scores + PRSs",
  "SCORE2 + clinical biomarkers + PRSs",
  "SCORE2 + NMR scores + clinical biomarkers + PRSs"
))]
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
fwrite(supp1, sep="\t", quote=FALSE, file="analyses/public_health_modelling/population_screening_comparison_main_analysis.txt")

# Build supp for main analyses - targeted screening
supp2 <- res[endpoint == "cvd" & score == "SCORE2_excl_UKB" & model_sex == "Sex-stratified" & model != "Risk score" & strategy == "targeted"]
supp2[, model := gsub("Risk score", "SCORE2", model)]
supp2[metric %like% "^delta", metric := gsub("delta", "alt_delta_ref", metric)]
supp2[, metric_type := fcase(
  metric %like% "^alt_delta_ref_", "alt_delta_ref",
  metric %like% "^alt_delta_null_", "alt_delta_null",
  metric %like% "^null_delta_ref_", "null_delta_ref",
  metric %like% "^ref_", "ref", 
  metric %like% "^alt_", "alt",
  metric %like% "^null_", "null"
)]
supp2 <- supp2[!is.na(metric_type)]
supp2[, metric := gsub("(^alt_delta_ref_)|(^alt_delta_null_)|(^null_delta_ref_)|(^ref_)|(^alt_)|(^null_)", "", metric)]
supp2 <- dcast(supp2, strategy + model + cohort + metric ~ metric_type, value.var=c("estimate", "se", "L95", "U95", "pval"))
supp2 <- supp2[, .(strategy, model, cohort, metric, 
  estimate_ref, se_ref, L95_ref, U95_ref,
  estimate_alt, se_alt, L95_alt, U95_alt,
  estimate_null, se_null, L95_null, U95_null,
  estimate_alt_delta_ref, se_alt_delta_ref, L95_alt_delta_ref, U95_alt_delta_ref, pval_alt_delta_ref,
  estimate_null_delta_ref, se_null_delta_ref, L95_null_delta_ref, U95_null_delta_ref, pval_null_delta_ref,
  estimate_alt_delta_null, se_alt_delta_null, L95_alt_delta_null, U95_alt_delta_null, pval_alt_delta_null
)]
supp2[, strategy := ifelse(strategy == "blanket", "Population-wide", "Targeted")]
supp2[, cohort := fcase(
  cohort == "discovery", "Discovery",
  cohort == "replication", "Replication",
  cohort == "pooled", "Discovery + Replication")]
supp2[, cohort := factor(cohort, levels=c("Discovery", "Replication", "Discovery + Replication"))]
supp2[, model := factor(model, levels=c(
  "SCORE2 + PRSs",
  "SCORE2 + NMR scores",
  "SCORE2 + clinical biomarkers",
  "SCORE2 + NMR scores + clinical biomarkers",
  "SCORE2 + NMR scores + PRSs",
  "SCORE2 + clinical biomarkers + PRSs",
  "SCORE2 + NMR scores + clinical biomarkers + PRSs"
))]
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
supp2 <- supp2[order(cohort)][order(metric)][order(model)][order(strategy)]
fwrite(supp2, sep="\t", quote=FALSE, file="analyses/public_health_modelling/targeted_screening_comparison_main_analysis.txt")

# Build supp for sensitivity analyses - population wide screening
supp3 <- res[score != "SCORE2" & model != "Risk score" & strategy == "blanket"]
supp3[metric %like% "^delta", metric := gsub("delta", "alt_delta_ref", metric)]
supp3[, metric_type := fcase(
  metric %like% "^alt_delta_ref_", "alt_delta_ref",
  metric %like% "^ref_", "ref", 
  metric %like% "^alt_", "alt"
)]
supp3 <- supp3[!is.na(metric_type)]
supp3[, metric := gsub("(^alt_delta_ref_)|(^ref_)|(^alt_)", "", metric)]
supp3 <- dcast(supp3, strategy + model + cohort + model_sex + score + endpoint + metric ~ metric_type, value.var=c("estimate", "se", "L95", "U95", "pval"))
supp3 <- supp3[, .(strategy, model, cohort, model_sex, score, endpoint, metric,
  estimate_ref, se_ref, L95_ref, U95_ref,
  estimate_alt, se_alt, L95_alt, U95_alt,
  estimate_alt_delta_ref, se_alt_delta_ref, L95_alt_delta_ref, U95_alt_delta_ref, pval_alt_delta_ref
)]
supp3[, strategy := ifelse(strategy == "blanket", "Population-wide", "Targeted")]
supp3[, cohort := fcase(
  cohort == "discovery", "Discovery",
  cohort == "replication", "Replication",
  cohort == "pooled", "Discovery + Replication")]
supp3[, cohort := factor(cohort, levels=c("Discovery", "Replication", "Discovery + Replication"))]
supp3[, model := factor(model, levels=c(
  "Risk score + PRSs",
  "Risk score + NMR scores",
  "Risk score + clinical biomarkers",
  "Risk score + NMR scores + clinical biomarkers",
  "Risk score + NMR scores + PRSs",
  "Risk score + clinical biomarkers + PRSs",
  "Risk score + NMR scores + clinical biomarkers + PRSs"
))]
supp3[, metric := fcase(
  metric == "high_risk", "People classified as high risk",
  metric == "cvd_high_risk", "Incident CVD events in the high-risk group",
  metric == "cvd_prevented", "Expected CVD events prevented by statin initiation",
  metric == "NNS", "Number needed to screen to prevent one CVD event (NNS)",
  metric == "NNT", "Number of statins prescribed per CVD event prevented (NNT)"
)]
supp3[, metric := factor(metric, levels=c(
  "People classified as high risk",
  "Incident CVD events in the high-risk group",
  "Expected CVD events prevented by statin initiation",
  "Number needed to screen to prevent one CVD event (NNS)",
  "Number of statins prescribed per CVD event prevented (NNT)"))]
supp3[, model_sex := factor(model_sex, levels=c("Sex-stratified", "Males", "Females"))]
supp3[score == "SCORE2_excl_UKB", score := "SCORE2"]
supp3[, score := factor(score, levels=c("SCORE2", "QRISK3"))]
supp3[, endpoint := fcase(
  endpoint == "cvd", "Broad",
  endpoint == "cvd_narrow", "Narrow")]
supp3 <- supp3[order(score)][order(model_sex)][order(cohort)][order(metric)][order(model)][order(strategy)]
fwrite(supp3, sep="\t", quote=FALSE, file="analyses/public_health_modelling/population_screening_comparison_sensitivity_analysis.txt")

# Build supp for sensitivity analyses - targeted screening
supp4 <- res[score != "SCORE2" & model != "Risk score" & strategy == "targeted"]
supp4[metric %like% "^delta", metric := gsub("delta", "alt_delta_ref", metric)]
supp4[, metric_type := fcase(
  metric %like% "^alt_delta_ref_", "alt_delta_ref",
  metric %like% "^alt_delta_null_", "alt_delta_null",
  metric %like% "^null_delta_ref_", "null_delta_ref",
  metric %like% "^ref_", "ref", 
  metric %like% "^alt_", "alt",
  metric %like% "^null_", "null"
)]
supp4 <- supp4[!is.na(metric_type)]
supp4[, metric := gsub("(^alt_delta_ref_)|(^alt_delta_null_)|(^null_delta_ref_)|(^ref_)|(^alt_)|(^null_)", "", metric)]
supp4 <- dcast(supp4, strategy + model + cohort + model_sex + score + endpoint + metric ~ metric_type, value.var=c("estimate", "se", "L95", "U95", "pval"))
supp4 <- supp4[, .(strategy, model, cohort, model_sex, score, endpoint, metric,
  estimate_ref, se_ref, L95_ref, U95_ref,
  estimate_alt, se_alt, L95_alt, U95_alt,
  estimate_null, se_null, L95_null, U95_null,
  estimate_alt_delta_ref, se_alt_delta_ref, L95_alt_delta_ref, U95_alt_delta_ref, pval_alt_delta_ref,
  estimate_null_delta_ref, se_null_delta_ref, L95_null_delta_ref, U95_null_delta_ref, pval_null_delta_ref,
  estimate_alt_delta_null, se_alt_delta_null, L95_alt_delta_null, U95_alt_delta_null, pval_alt_delta_null
)]
supp4[, strategy := ifelse(strategy == "blanket", "Population-wide", "Targeted")]
supp4[, cohort := fcase(
  cohort == "discovery", "Discovery",
  cohort == "replication", "Replication",
  cohort == "pooled", "Discovery + Replication")]
supp4[, cohort := factor(cohort, levels=c("Discovery", "Replication", "Discovery + Replication"))]
supp4[, model := factor(model, levels=c(
  "Risk score + PRSs",
  "Risk score + NMR scores",
  "Risk score + clinical biomarkers",
  "Risk score + NMR scores + clinical biomarkers",
  "Risk score + NMR scores + PRSs",
  "Risk score + clinical biomarkers + PRSs",
  "Risk score + NMR scores + clinical biomarkers + PRSs"
))]
supp4[, metric := fcase(
  metric == "high_risk", "People classified as high risk",
  metric == "cvd_high_risk", "Incident CVD events in the high-risk group",
  metric == "cvd_prevented", "Expected CVD events prevented by statin initiation",
  metric == "NNS", "Number needed to screen to prevent one CVD event (NNS)",
  metric == "NNT", "Number of statins prescribed per CVD event prevented (NNT)"
)]
supp4[, metric := factor(metric, levels=c(
  "People classified as high risk",
  "Incident CVD events in the high-risk group",
  "Expected CVD events prevented by statin initiation",
  "Number needed to screen to prevent one CVD event (NNS)",
  "Number of statins prescribed per CVD event prevented (NNT)"))]
supp4[, model_sex := factor(model_sex, levels=c("Sex-stratified", "Males", "Females"))]
supp4[score == "SCORE2_excl_UKB", score := "SCORE2"]
supp4[, score := factor(score, levels=c("SCORE2", "QRISK3"))]
supp4[, endpoint := fcase(
  endpoint == "cvd", "Broad",
  endpoint == "cvd_narrow", "Narrow")]
supp4 <- supp4[order(score)][order(model_sex)][order(cohort)][order(metric)][order(model)][order(strategy)]
fwrite(supp4, sep="\t", quote=FALSE, file="analyses/public_health_modelling/targeted_screening_comparison_sensitivity_analysis.txt")

# Run pairwise comparisons to assess similarity of results across different settings
model_info <- unique(res[metric %like% "delta",.(strategy, cohort, endpoint, model_sex, score, metric)])
comp_stats <- foreach(modelIdx1=model_info[,.I], .combine=rbind) %do% {
  foreach(modelIdx2=model_info[,.I], .combine=rbind) %do% {
   if (modelIdx1 >= modelIdx2) return(NULL)

    model_info1 <- model_info[modelIdx1]
    model_info2 <- model_info[modelIdx2]

    res1 <- res[model_info1, on=.(strategy, cohort, endpoint, model_sex, score, metric)]
    res2 <- res[model_info2, on=.(strategy, cohort, endpoint, model_sex, score, metric)]

    res1 <- res1[, .(model, estimate)]
    res2 <- res2[, .(model, estimate)]

    comp <- merge(res1, res2, by="model", suffixes=c(".1", ".2"))

    corr <- comp[,cor(estimate.1, estimate.2)]
    beta <- comp[,coef(lm(estimate.2 ~ 0 + estimate.1))[1]]

    this_comp_stats <- data.table(pearson=corr, lm_beta=beta)

    setnames(model_info1, paste0(names(model_info1), ".1"))
    setnames(model_info2, paste0(names(model_info2), ".2"))
    return(cbind(model_info1, model_info2, this_comp_stats))
  }
}
fwrite(comp_stats, sep="\t", quote=FALSE, file="analyses/public_health_modelling/comparison_sensitivity.txt")

# Compare the five metrics across the discovery, replication, and pooled cohorts
cohort_comp <- res[endpoint == "cvd" & model_sex == "Sex-stratified" & score == "SCORE2_excl_UKB" & model != "Risk score"]
cohort_comp[, group := fcase(
  strategy == "blanket" & metric %like% "delta", "blanket",
  strategy == "targeted" & metric %like% "alt_delta_ref", "targeted",
  strategy == "targeted" & metric %like% "alt_delta_null", "targeted2"
)]
cohort_comp <- cohort_comp[!is.na(group)]
cohort_comp[, metric := gsub("(^alt_delta_ref_)|(^alt_delta_null_)|(^delta_)", "", metric)]
cohort_comp <- dcast(cohort_comp, group + model + metric ~ cohort, value.var=c("estimate", "L95", "U95"))

common_gg_parts <- function(g, limits, breaks) {
  g + aes(shape=group, fill=model, color=model) +
    geom_abline(intercept=0, slope=1, linetype=2) +
    geom_hline(yintercept=0, linetype=3) +
    geom_vline(xintercept=0, linetype=3) +
    geom_errorbarh(height=0, alpha=0.7) + geom_errorbar(width=0, alpha=0.7) +
    geom_point(fill="white") +
    scale_shape_manual(values=c("blanket"=23, "targeted"=22, "targeted2"=21)) +
    scale_color_manual(values=c(
      "Risk score + NMR scores"="#e41a1c", "Risk score + PRSs"="#377eb8", "Risk score + clinical biomarkers"="#ff7f00",
      "Risk score + NMR scores + PRSs"="#4daf4a", "Risk score + clinical biomarkers + PRSs"="#984ea3", "Risk score + NMR scores + clinical biomarkers"="#f781bf",
      "Risk score + NMR scores + clinical biomarkers + PRSs"="#a65628"
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

g1 <- ggplot(cohort_comp[metric == "high_risk"]) + aes(
    x=estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=estimate_pooled, ymin=L95_pooled, ymax=U95_pooled
  ) + xlab("Discovery cohort") + ylab("Pooled cohort") + ggtitle("High Risk")
g1 <- common_gg_parts(g1, c(-200, 4800), c(0, 2000, 4000))

g2 <- ggplot(cohort_comp[metric == "high_risk"]) + aes(
    x=estimate_replication, xmin=L95_replication, xmax=U95_replication,
    y=estimate_pooled, ymin=L95_pooled, ymax=U95_pooled
  ) + xlab("Replication cohort") + ylab("Pooled cohort")
g2 <- common_gg_parts(g2, c(-200, 4800), c(0, 2000, 4000))

g3 <- ggplot(cohort_comp[metric == "high_risk"]) + aes(
    x=estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=estimate_replication, ymin=L95_replication, ymax=U95_replication
  ) + xlab("Discovery cohort") + ylab("Replication cohort") 
g3 <- common_gg_parts(g3, c(-200, 4800), c(0, 2000, 4000))


g4 <- ggplot(cohort_comp[metric == "cvd_high_risk"]) + aes(
    x=estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=estimate_pooled, ymin=L95_pooled, ymax=U95_pooled
  ) + xlab("Discovery cohort") + ylab("Pooled cohort") + ggtitle("CVD high")
g4 <- common_gg_parts(g4, c(0, 1100), c(0, 500, 1000))

g5 <- ggplot(cohort_comp[metric == "cvd_high_risk"]) + aes(
    x=estimate_replication, xmin=L95_replication, xmax=U95_replication,
    y=estimate_pooled, ymin=L95_pooled, ymax=U95_pooled
  ) + xlab("Replication cohort") + ylab("Pooled cohort")
g5 <- common_gg_parts(g5, c(0, 1100), c(0, 500, 1000))

g6 <- ggplot(cohort_comp[metric == "cvd_high_risk"]) + aes(
    x=estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=estimate_replication, ymin=L95_replication, ymax=U95_replication
  ) + xlab("Discovery cohort") + ylab("Replication cohort") 
g6 <- common_gg_parts(g6, c(0, 1100), c(0, 500, 1000))


g7 <- ggplot(cohort_comp[metric == "cvd_prevented"]) + aes(
    x=estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=estimate_pooled, ymin=L95_pooled, ymax=U95_pooled
  ) + xlab("Discovery cohort") + ylab("Pooled cohort") + ggtitle("Prevented")
g7 <- common_gg_parts(g7, c(0, 220), c(0, 100, 200))

g8 <- ggplot(cohort_comp[metric == "cvd_prevented"]) + aes(
    x=estimate_replication, xmin=L95_replication, xmax=U95_replication,
    y=estimate_pooled, ymin=L95_pooled, ymax=U95_pooled
  ) + xlab("Replication cohort") + ylab("Pooled cohort")
g8 <- common_gg_parts(g8, c(0, 220), c(0, 100, 200))

g9 <- ggplot(cohort_comp[metric == "cvd_prevented"]) + aes(
    x=estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=estimate_replication, ymin=L95_replication, ymax=U95_replication
  ) + xlab("Discovery cohort") + ylab("Replication cohort") 
g9 <- common_gg_parts(g9, c(0, 220), c(0, 100, 200))


g10 <- ggplot(cohort_comp[metric == "NNS"]) + aes(
    x=estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=estimate_pooled, ymin=L95_pooled, ymax=U95_pooled
  ) + xlab("Discovery cohort") + ylab("Pooled cohort") + ggtitle("NNS")
g10 <- common_gg_parts(g10, c(-285, 0), c(-200, -100, 0))

g11 <- ggplot(cohort_comp[metric == "NNS"]) + aes(
    x=estimate_replication, xmin=L95_replication, xmax=U95_replication,
    y=estimate_pooled, ymin=L95_pooled, ymax=U95_pooled
  ) + xlab("Replication cohort") + ylab("Pooled cohort")
g11 <- common_gg_parts(g11, c(-285, 0), c(-200, -100, 0))

g12 <- ggplot(cohort_comp[metric == "NNS"]) + aes(
    x=estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=estimate_replication, ymin=L95_replication, ymax=U95_replication
  ) + xlab("Discovery cohort") + ylab("Replication cohort") 
g12 <- common_gg_parts(g12, c(-285, 0), c(-200, -100, 0))


g13 <- ggplot(cohort_comp[metric == "NNT"]) + aes(
    x=estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=estimate_pooled, ymin=L95_pooled, ymax=U95_pooled
  ) + xlab("Discovery cohort") + ylab("Pooled cohort") + ggtitle("NNT")
g13 <- common_gg_parts(g13, c(-4.5, 2), c(-4, -2, 0, 2))

g14 <- ggplot(cohort_comp[metric == "NNT"]) + aes(
    x=estimate_replication, xmin=L95_replication, xmax=U95_replication,
    y=estimate_pooled, ymin=L95_pooled, ymax=U95_pooled
  ) + xlab("Replication cohort") + ylab("Pooled cohort")
g14 <- common_gg_parts(g14, c(-4.5, 2), c(-4, -2, 0, 2))

g15 <- ggplot(cohort_comp[metric == "NNT"]) + aes(
    x=estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=estimate_replication, ymin=L95_replication, ymax=U95_replication
  ) + xlab("Discovery cohort") + ylab("Replication cohort") 
g15 <- common_gg_parts(g15, c(-4.5, 2), c(-4, -2, 0, 2))


g <- g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8 + g9 + g10 + g11 + g12 + g13 + g14 + g15 + plot_layout(nrow=3, byrow=FALSE)
ggsave(g, width=7.2, height=4.2, file="analyses/public_health_modelling/cohort_sensitivity.pdf")

# Compare broad CVD to narrow CVD
cvd_comp <- res[cohort == "pooled" & score == "SCORE2_excl_UKB" & model_sex == "Sex-stratified" & model != "Risk score"]
cvd_comp[, group := fcase(
  strategy == "blanket" & metric %like% "delta", "blanket",
  strategy == "targeted" & metric %like% "alt_delta_ref", "targeted",
  strategy == "targeted" & metric %like% "alt_delta_null", "targeted2"
)]
cvd_comp <- cvd_comp[!is.na(group)]
cvd_comp[, metric := gsub("(^alt_delta_ref_)|(^alt_delta_null_)|(^delta_)", "", metric)]
cvd_comp <- dcast(cvd_comp, group + model + metric ~ endpoint, value.var=c("estimate", "L95", "U95"))

g1 <- ggplot(cvd_comp[metric == "high_risk"]) + aes(
    x=estimate_cvd_narrow, xmin=L95_cvd_narrow, xmax=U95_cvd_narrow,
    y=estimate_cvd, ymin=L95_cvd, ymax=U95_cvd
  ) + xlab("Narrow CVD") + ylab("Broad CVD") + ggtitle("High Risk")
g1 <- common_gg_parts(g1, c(-200, 4800), c(0, 2000, 4000))

g4 <- ggplot(cvd_comp[metric == "cvd_high_risk"]) + aes(
    x=estimate_cvd_narrow, xmin=L95_cvd_narrow, xmax=U95_cvd_narrow,
    y=estimate_cvd, ymin=L95_cvd, ymax=U95_cvd
  ) + xlab("Narrow CVD") + ylab("Broad CVD") + ggtitle("CVD high")
g4 <- common_gg_parts(g4, c(0, 1100), c(0, 500, 1000))

g7 <- ggplot(cvd_comp[metric == "cvd_prevented"]) + aes(
    x=estimate_cvd_narrow, xmin=L95_cvd_narrow, xmax=U95_cvd_narrow,
    y=estimate_cvd, ymin=L95_cvd, ymax=U95_cvd
  ) + xlab("Narrow CVD") + ylab("Broad CVD") + ggtitle("Prevented")
g7 <- common_gg_parts(g7, c(0, 220), c(0, 100, 200))

g10 <- ggplot(cvd_comp[metric == "NNS"]) + aes(
    x=estimate_cvd_narrow, xmin=L95_cvd_narrow, xmax=U95_cvd_narrow,
    y=estimate_cvd, ymin=L95_cvd, ymax=U95_cvd
  ) + xlab("Narrow CVD") + ylab("Broad CVD") + ggtitle("NNS")
g10 <- common_gg_parts(g10, c(-275, 0), c(-200, -100, 0))

g13 <- ggplot(cvd_comp[metric == "NNT"]) + aes(
    x=estimate_cvd_narrow, xmin=L95_cvd_narrow, xmax=U95_cvd_narrow,
    y=estimate_cvd, ymin=L95_cvd, ymax=U95_cvd
  ) + xlab("Narrow CVD") + ylab("Broad CVD") + ggtitle("NNT")
g13 <- common_gg_parts(g13, c(-4.5, 2), c(-4, -2, 0, 2))

g <- g1 + g4 + g7 + g10 + g13 + plot_layout(nrow=1)
ggsave(g, width=7.2, height=1.62181102, file="analyses/public_health_modelling/endpoint_sensitivity.pdf")

# Build Figure similar to main plot showing results in Males
ggdt <- res[cohort == "pooled" & endpoint == "cvd" & score == "SCORE2_excl_UKB" & model_sex == "Males" & model != "Risk score"]
ggdt[, model_name := fcase(
  model == "Risk score + NMR scores", "SCORE2 + NMR scores for CHD and IS",
  model == "Risk score + clinical biomarkers", "SCORE2 + 11 clinical chemistry biomarkers",
  model == "Risk score + PRSs", "SCORE2 + PRSs for CHD and IS",
  model == "Risk score + NMR scores + PRSs", "SCORE2 + NMR scores + PRSs",
  model == "Risk score + clinical biomarkers + PRSs", "SCORE2 + 11 clinical biomarkers + PRSs",
  model == "Risk score + NMR scores + clinical biomarkers", "SCORE2 + NMR scores + 11 clinical biomarkers",
  model == "Risk score + NMR scores + clinical biomarkers + PRSs", "SCORE2 + NMR scores + 11 clinical biomarkers + PRSs"
)]
ggdt[,model_name := factor(model_name, levels=c("SCORE2 + PRSs for CHD and IS", "SCORE2 + NMR scores for CHD and IS",
  "SCORE2 + 11 clinical chemistry biomarkers", "SCORE2 + NMR scores + 11 clinical biomarkers", "SCORE2 + NMR scores + PRSs",
  "SCORE2 + 11 clinical biomarkers + PRSs", "SCORE2 + NMR scores + 11 clinical biomarkers + PRSs"))]

g1 <- screening_plot(ggdt[metric == "delta_high_risk" & strategy == "blanket"]) + scale_x_continuous(limits=c(-600, 7800), breaks=c(0, 3000, 6000))
g2 <- screening_plot(ggdt[metric == "delta_cvd_high_risk" & strategy == "blanket"]) + scale_x_continuous(limits=c(0, 1800), breaks=c(0, 1000))
g3 <- screening_plot(ggdt[metric == "delta_cvd_prevented" & strategy == "blanket"]) + scale_x_continuous(limits=c(0, 360), breaks=c(0, 100, 200, 300))
g4 <- screening_plot(ggdt[metric == "delta_NNS" & strategy == "blanket"]) + scale_x_continuous(limits=c(-140, 0), breaks=c(-100, -50, 0))
g5 <- screening_plot(ggdt[metric == "delta_NNT" & strategy == "blanket"]) + scale_x_continuous(limits=c(-4.25, 4.25), breaks=c(-4, -2, 0, 2, 4))

g6 <- screening_plot(ggdt[metric == "alt_delta_ref_high_risk" & strategy == "targeted"]) + scale_x_continuous(limits=c(-600, 7800), breaks=c(0, 3000, 6000))
g7 <- screening_plot(ggdt[metric == "alt_delta_ref_cvd_high_risk" & strategy == "targeted"]) + scale_x_continuous(limits=c(0, 1800), breaks=c(0, 1000))
g8 <- screening_plot(ggdt[metric == "alt_delta_ref_cvd_prevented" & strategy == "targeted"]) + scale_x_continuous(limits=c(0, 360), breaks=c(0, 100, 200, 300))
g9 <- screening_plot(ggdt[metric == "alt_delta_ref_NNS" & strategy == "targeted"]) + scale_x_continuous(limits=c(-140, 0), breaks=c(-100, -50, 0))
g10 <- screening_plot(ggdt[metric == "alt_delta_ref_NNT" & strategy == "targeted"]) + scale_x_continuous(limits=c(-4.25, 4.25), breaks=c(-4, -2, 0, 2, 4))

g11 <- screening_plot(ggdt[metric == "alt_delta_null_high_risk" & strategy == "targeted"]) + scale_x_continuous(limits=c(-600, 7800), breaks=c(0, 3000, 6000))
g12 <- screening_plot(ggdt[metric == "alt_delta_null_cvd_high_risk" & strategy == "targeted"]) + scale_x_continuous(limits=c(0, 1800), breaks=c(0, 1000))
g13 <- screening_plot(ggdt[metric == "alt_delta_null_cvd_prevented" & strategy == "targeted"]) + scale_x_continuous(limits=c(0, 360), breaks=c(0, 100, 200, 300))
g14 <- screening_plot(ggdt[metric == "alt_delta_null_NNS" & strategy == "targeted"]) + scale_x_continuous(limits=c(-140, 0), breaks=c(-100, -50, 0))
g15 <- screening_plot(ggdt[metric == "alt_delta_null_NNT" & strategy == "targeted"]) + scale_x_continuous(limits=c(-4.25, 4.25), breaks=c(-4, -2, 0, 2, 4))

g <- g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8 + g9 + g10 + g11 + g12 + g13 + g14 + g15 + plot_layout(axes="collect", nrow=3, byrow=TRUE)
ggsave(g, width=7.2, height=4.425, file="analyses/public_health_modelling/screening_comparison_males.pdf")

# Build Figure similar to main plot showing results in Females
ggdt <- res[cohort == "pooled" & endpoint == "cvd" & score == "SCORE2_excl_UKB" & model_sex == "Females" & model != "Risk score"]
ggdt[, model_name := fcase(
  model == "Risk score + NMR scores", "SCORE2 + NMR scores for CHD and IS",
  model == "Risk score + clinical biomarkers", "SCORE2 + 11 clinical chemistry biomarkers",
  model == "Risk score + PRSs", "SCORE2 + PRSs for CHD and IS",
  model == "Risk score + NMR scores + PRSs", "SCORE2 + NMR scores + PRSs",
  model == "Risk score + clinical biomarkers + PRSs", "SCORE2 + 11 clinical biomarkers + PRSs",
  model == "Risk score + NMR scores + clinical biomarkers", "SCORE2 + NMR scores + 11 clinical biomarkers",
  model == "Risk score + NMR scores + clinical biomarkers + PRSs", "SCORE2 + NMR scores + 11 clinical biomarkers + PRSs"
)]
ggdt[,model_name := factor(model_name, levels=c("SCORE2 + PRSs for CHD and IS", "SCORE2 + NMR scores for CHD and IS",
  "SCORE2 + 11 clinical chemistry biomarkers", "SCORE2 + NMR scores + 11 clinical biomarkers", "SCORE2 + NMR scores + PRSs",
  "SCORE2 + 11 clinical biomarkers + PRSs", "SCORE2 + NMR scores + 11 clinical biomarkers + PRSs"))]

g1 <- screening_plot(ggdt[metric == "delta_high_risk" & strategy == "blanket"]) + scale_x_continuous(limits=c(-200, 2800), breaks=c(0, 1000, 2000))
g2 <- screening_plot(ggdt[metric == "delta_cvd_high_risk" & strategy == "blanket"]) + scale_x_continuous(limits=c(-50, 500), breaks=c(0, 200, 400))
g3 <- screening_plot(ggdt[metric == "delta_cvd_prevented" & strategy == "blanket"]) + scale_x_continuous(limits=c(-10, 100), breaks=c(0, 40, 80))
g4 <- screening_plot(ggdt[metric == "delta_NNS" & strategy == "blanket"]) + scale_x_continuous(limits=c(-3250, 250), breaks=c(-2000, 0))
g5 <- screening_plot(ggdt[metric == "delta_NNT" & strategy == "blanket"]) + scale_x_continuous(limits=c(-9, 9), breaks=c(-8, -4, 0, 4, 8))

g6 <- screening_plot(ggdt[metric == "alt_delta_ref_high_risk" & strategy == "targeted"]) + scale_x_continuous(limits=c(-200, 2800), breaks=c(0, 1000, 2000))
g7 <- screening_plot(ggdt[metric == "alt_delta_ref_cvd_high_risk" & strategy == "targeted"]) + scale_x_continuous(limits=c(-50, 500), breaks=c(0, 200, 400))
g8 <- screening_plot(ggdt[metric == "alt_delta_ref_cvd_prevented" & strategy == "targeted"]) + scale_x_continuous(limits=c(-10, 100), breaks=c(0, 40, 80))
g9 <- screening_plot(ggdt[metric == "alt_delta_ref_NNS" & strategy == "targeted"]) + scale_x_continuous(limits=c(-3250, 250), breaks=c(-2000, 0))
g10 <- screening_plot(ggdt[metric == "alt_delta_ref_NNT" & strategy == "targeted"]) + scale_x_continuous(limits=c(-9, 9), breaks=c(-8, -4, 0, 4, 8))

g11 <- screening_plot(ggdt[metric == "alt_delta_null_high_risk" & strategy == "targeted"]) + scale_x_continuous(limits=c(-200, 2800), breaks=c(0, 1000, 2000))
g12 <- screening_plot(ggdt[metric == "alt_delta_null_cvd_high_risk" & strategy == "targeted"]) + scale_x_continuous(limits=c(-50, 500), breaks=c(0, 200, 400))
g13 <- screening_plot(ggdt[metric == "alt_delta_null_cvd_prevented" & strategy == "targeted"]) + scale_x_continuous(limits=c(-10, 100), breaks=c(0, 40, 80))
g14 <- screening_plot(ggdt[metric == "alt_delta_null_NNS" & strategy == "targeted"]) + scale_x_continuous(limits=c(-3250, 250), breaks=c(-2000, 0))
g15 <- screening_plot(ggdt[metric == "alt_delta_null_NNT" & strategy == "targeted"]) + scale_x_continuous(limits=c(-9, 9), breaks=c(-8, -4, 0, 4, 8))

g <- g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8 + g9 + g10 + g11 + g12 + g13 + g14 + g15 + plot_layout(axes="collect", nrow=3, byrow=TRUE)
ggsave(g, width=7.2, height=4.425, file="analyses/public_health_modelling/screening_comparison_females.pdf")

# Build Figure similar to main plot showing results with QRISK3
ggdt <- res[cohort == "pooled" & endpoint == "cvd" & score == "QRISK3" & model_sex == "Sex-stratified" & model != "Risk score"]
ggdt[, model_name := fcase(
  model == "Risk score + NMR scores", "QRISK3 + NMR scores for CHD and IS",
  model == "Risk score + clinical biomarkers", "QRISK3 + 11 clinical chemistry biomarkers",
  model == "Risk score + PRSs", "QRISK3 + PRSs for CHD and IS",
  model == "Risk score + NMR scores + PRSs", "QRISK3 + NMR scores + PRSs",
  model == "Risk score + clinical biomarkers + PRSs", "QRISK3 + 11 clinical biomarkers + PRSs",
  model == "Risk score + NMR scores + clinical biomarkers", "QRISK3 + NMR scores + 11 clinical biomarkers",
  model == "Risk score + NMR scores + clinical biomarkers + PRSs", "QRISK3 + NMR scores + 11 clinical biomarkers + PRSs"
)]
ggdt[,model_name := factor(model_name, levels=c("QRISK3 + PRSs for CHD and IS", "QRISK3 + NMR scores for CHD and IS",
  "QRISK3 + 11 clinical chemistry biomarkers", "QRISK3 + NMR scores + 11 clinical biomarkers", "QRISK3 + NMR scores + PRSs",
  "QRISK3 + 11 clinical biomarkers + PRSs", "QRISK3 + NMR scores + 11 clinical biomarkers + PRSs"))]

g1 <- screening_plot(ggdt[metric == "delta_high_risk" & strategy == "blanket"]) + scale_x_continuous(limits=c(-150, 7000), breaks=c(0, 3000, 6000))
g2 <- screening_plot(ggdt[metric == "delta_cvd_high_risk" & strategy == "blanket"]) + scale_x_continuous(limits=c(0, 750), breaks=c(0, 300, 600))
g3 <- screening_plot(ggdt[metric == "delta_cvd_prevented" & strategy == "blanket"]) + scale_x_continuous(limits=c(0, 142), breaks=c(0, 50, 100))
g4 <- screening_plot(ggdt[metric == "delta_NNS" & strategy == "blanket"]) + scale_x_continuous(limits=c(-11, 0), breaks=c(-10, -5, 0))
g5 <- screening_plot(ggdt[metric == "delta_NNT" & strategy == "blanket"]) + scale_x_continuous(limits=c(-1.2, 2.2), breaks=c(-1, 0, 1, 2))

g6 <- screening_plot(ggdt[metric == "alt_delta_ref_high_risk" & strategy == "targeted"]) + scale_x_continuous(limits=c(-150, 7000), breaks=c(0, 3000, 6000))
g7 <- screening_plot(ggdt[metric == "alt_delta_ref_cvd_high_risk" & strategy == "targeted"]) + scale_x_continuous(limits=c(0, 750), breaks=c(0, 300, 600))
g8 <- screening_plot(ggdt[metric == "alt_delta_ref_cvd_prevented" & strategy == "targeted"]) + scale_x_continuous(limits=c(0, 142), breaks=c(0, 50, 100))
g9 <- screening_plot(ggdt[metric == "alt_delta_ref_NNS" & strategy == "targeted"]) + scale_x_continuous(limits=c(-11, 0), breaks=c(-10, -5, 0))
g10 <- screening_plot(ggdt[metric == "alt_delta_ref_NNT" & strategy == "targeted"]) + scale_x_continuous(limits=c(-1.2, 2.2), breaks=c(-1, 0, 1, 2))

g11 <- screening_plot(ggdt[metric == "alt_delta_null_high_risk" & strategy == "targeted"]) + scale_x_continuous(limits=c(-150, 7000), breaks=c(0, 3000, 6000))
g12 <- screening_plot(ggdt[metric == "alt_delta_null_cvd_high_risk" & strategy == "targeted"]) + scale_x_continuous(limits=c(0, 750), breaks=c(0, 300, 600))
g13 <- screening_plot(ggdt[metric == "alt_delta_null_cvd_prevented" & strategy == "targeted"]) + scale_x_continuous(limits=c(0, 142), breaks=c(0, 50, 100))
g14 <- screening_plot(ggdt[metric == "alt_delta_null_NNS" & strategy == "targeted"]) + scale_x_continuous(limits=c(-11, 0), breaks=c(-10, -5, 0))
g15 <- screening_plot(ggdt[metric == "alt_delta_null_NNT" & strategy == "targeted"]) + scale_x_continuous(limits=c(-1.2, 2.2), breaks=c(-1, 0, 1, 2))

g <- g1 + g2 + g3 + g4 + g5 + g6 + g7 + g8 + g9 + g10 + g11 + g12 + g13 + g14 + g15 + plot_layout(axes="collect", nrow=3, byrow=TRUE)
ggsave(g, width=7.150787, height=4.425, file="analyses/public_health_modelling/screening_comparison_QRISK3.pdf")



