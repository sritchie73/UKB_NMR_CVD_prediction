library(foreach)
library(forcats)
library(ggplot2)
library(ggstance)
library(scales)
library(patchwork)

# Load and collate results
res <- rbind(idcol="cohort",
  "discovery"=fread("analyses/test/discovery_categorical_nri_estimates.txt"),
  "replication"=fread("analyses/test/replication_categorical_nri_estimates.txt"),
  "pooled"=fread("analyses/test/pooled_categorical_nri_estimates.txt")
)

res2 <- rbind(idcol="cohort",
  "discovery"=fread("analyses/test/discovery_categorical_nri_reclassified.txt"),
  "replication"=fread("analyses/test/replication_categorical_nri_reclassified.txt"),
  "pooled"=fread("analyses/test/pooled_categorical_nri_reclassified.txt")
)

# Build main figure 3B
ggdt <- res[cohort == "pooled" & model_sex == "Sex-stratified" & endpoint == "cvd" & score == "SCORE2_excl_UKB"]
ggdt <- ggdt[metric %in% c("NRI+", "NRI-")]
ggdt[, model_name := fcase(
  model == "SCORE2 + NMR scores", "SCORE2 + NMR scores for CHD and IS",
  model == "SCORE2 + Biochemistry", "SCORE2 + 11 clinical chemistry biomarkers",
  model == "SCORE2 + PRS", "SCORE2 + PRSs for CHD and IS",
  model == "SCORE2 + NMR scores + PRS", "SCORE2 + NMR scores for CHD and IS + PRSs for CHD and IS",
  model == "SCORE2 + Biochemistry + PRS", "SCORE2 + 11 clinical chemistry biomarkers + PRSs for CHD and IS"
)]
ggdt[,model_name := factor(model_name, levels=unique(model_name))]

g <- ggplot(ggdt) +
  aes(x=Estimate, xmin=L95, xmax=U95, y=fct_rev(model_name), color=metric) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0, position=position_dodgev(height=0.3)) +
  geom_point(shape=23, size=2, fill="white", position=position_dodgev(height=0.3)) +
  scale_color_manual(values=c("NRI+"="#c51b7d", "NRI-"="#4d9221")) +
  scale_x_continuous("Categorical NRI, % reclassified (95% CI)", labels=percent, limits=c(-0.05, 0.17)) +
  ylab("") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )
ggsave(g, width=5.3, height=2, file="analyses/test/main_categorical_nri.pdf")

# Build supp tables for main results
supp1 <- res[endpoint == "cvd" & model_sex == "Sex-stratified" & score == "SCORE2_excl_UKB"]
supp1 <- supp1[metric %in% c("NRI+", "NRI-")]
supp1[, metric := ifelse(metric == "NRI+", "Net % CVD cases correctly reclassified", "Net % non-cases correctly reclassified")]
supp1[, model := factor(model, levels=unique(model))]
supp1 <- supp1[order(model)][order(metric)]
supp1[, c("endpoint", "model_sex", "score", "model_type", "Lower", "Upper") := NULL]
fwrite(supp1, sep="\t", quote=FALSE, file="analyses/test/categorical_nri_main_analysis.txt")

supp2 <- res2[endpoint == "cvd" & model_sex == "Sex-stratified" & score == "SCORE2_excl_UKB"]
supp2[, Non_cases := All - Cases]
supp2[, Old := fcase( 
  Old == "< 0.05", "Low risk",
  Old == "< 0.1", "Medium risk",
  Old == ">= 0.1", "High risk"
)]
supp2[, New := fcase(
  New == "< 0.05", "Low risk",
  New == "< 0.1", "Medium risk",
  New == ">= 0.1", "High risk"
)]
supp2[,Old := factor(Old, levels=c("Low risk", "Medium risk", "High risk"))]
supp2[,New := factor(New, levels=c("Low risk", "Medium risk", "High risk"))]
supp2 <- supp2[,.(model, cohort, SCORE2=Old, New, Cases, Non_cases)]
supp2[, model := factor(model, levels=unique(model))]
supp2 <- supp2[order(SCORE2)][order(New)][order(model)]
fwrite(supp2, sep="\t", quote=FALSE, file="analyses/test/categorical_nri_main_analysis_reclassified_numbers.txt")

# Build supp tables for sensitivity analyses
supp3 <- res[score != "SCORE2"]
supp3[score == "SCORE2_excl_UKB", score := "SCORE2"]
supp3 <- supp3[metric %in% c("NRI+", "NRI-")]
supp3[, metric := ifelse(metric == "NRI+", "Net % CVD cases correctly reclassified", "Net % non-cases correctly reclassified")]
supp3[, model := gsub("SCORE2", "Risk score", model)]
supp3[, model := gsub("QRISK3", "Risk score", model)]
supp3[, model := factor(model, levels=unique(model))]
supp3[, endpoint := fcase(
  endpoint == "cvd", "Broad",
  endpoint == "cvd_narrow", "Narrow")]
supp3[, model_sex := factor(model_sex, levels=c("Sex-stratified", "Males", "Females"))]
supp3[, cohort := factor(cohort, levels=c("discovery", "replication", "pooled"))]
supp3[, score := factor(score, levels=c("SCORE2", "QRISK3"))]
supp3 <- supp3[order(model_sex)][order(cohort)][order(model)][order(metric)]
fwrite(supp3, sep="\t", quote=FALSE, file="analyses/test/categorical_nri_sensitivity_analysis.txt")

supp4 <- res2[score != "SCORE2"]
supp4[score == "SCORE2_excl_UKB", score := "SCORE2"]
supp4[, Non_cases := All - Cases]
supp4[, Old := fcase( 
  Old == "< 0.05", "Low risk",
  Old == "< 0.1" & score == "QRISK3", "Low risk",
  Old == "< 0.1" & score == "SCORE2", "Medium risk",
  Old == ">= 0.1", "High risk"
)]
supp4[, New := fcase(
  New == "< 0.05", "Low risk",
  New == "< 0.1" & score == "QRISK3", "Low risk",
  New == "< 0.1" & score == "SCORE2", "Medium risk",
  New == ">= 0.1", "High risk"
)]
supp4[,Old := factor(Old, levels=c("Low risk", "Medium risk", "High risk"))]
supp4[,New := factor(New, levels=c("Low risk", "Medium risk", "High risk"))]
supp4 <- supp4[,.(model, cohort, model_sex, score, endpoint, Old, New, Cases, Non_cases)]
supp4[, model := factor(model, levels=unique(model))]
supp4 <- supp4[order(Old)][order(New)][order(model_sex)][order(cohort)][order(model)]
fwrite(supp4, sep="\t", quote=FALSE, file="analyses/test/categorical_nri_sensitivity_analysis_reclassified_numbers.txt")

# Run pairwise comparisons to assess similarity of results across different settings
model_info <- unique(res[,.(cohort, endpoint, model_sex, score)])
comp_stats <- foreach(modelIdx1=model_info[,.I], .combine=rbind) %do% {
  foreach(modelIdx2=model_info[,.I], .combine=rbind) %do% {
   if (modelIdx1 == modelIdx2) return(NULL)

    model_info1 <- model_info[modelIdx1]
    model_info2 <- model_info[modelIdx2]

    res1 <- res[model_info1, on=.(cohort, endpoint, model_sex, score)]
    res2 <- res[model_info2, on=.(cohort, endpoint, model_sex, score)]

    res1 <- res1[, .(metric, model_type, Estimate)]
    res2 <- res2[, .(metric, model_type, Estimate)]

    comp <- merge(res1, res2, by=c("model_type", "metric"), suffixes=c(".1", ".2"))

    corr1 <- comp[metric == "NRI+", cor(Estimate.1, Estimate.2)]
    beta1 <- comp[metric == "NRI+", coef(lm(Estimate.2 ~ 0 + Estimate.1))[1]]

    corr2 <- comp[metric == "NRI-", cor(Estimate.1, Estimate.2)]
    beta2 <- comp[metric == "NRI-", coef(lm(Estimate.2 ~ 0 + Estimate.1))[1]]

    this_comp_stats <- data.table(metric=c("NRI+", "NRI-"), pearson=c(corr1, corr2), lm_beta=c(beta1, beta2))

    setnames(model_info1, paste0(names(model_info1), ".1"))
    setnames(model_info2, paste0(names(model_info2), ".2"))
    return(cbind(model_info1, model_info2, this_comp_stats))
  }
}
fwrite(comp_stats, sep="\t", quote=FALSE, file="analyses/test/categorical_nri_comparison_sensitivity.txt")

# Create plot of key comparisons of interest
common_gg_parts <- function(g) {
  g +
    geom_abline(intercept=0, slope=1, linetype=2) +
    geom_hline(yintercept=0, linetype=2) +
    geom_vline(xintercept=0, linetype=2) +
    geom_errorbarh(height=0, alpha=0.7) + geom_errorbar(width=0, alpha=0.7) +
    geom_point(shape=23) +
    scale_x_continuous(limits=c(-0.05, 0.25), oob=scales::squish, breaks=c(0, 0.1, 0.2), labels=percent) +
    scale_y_continuous(limits=c(-0.05, 0.25), oob=scales::squish, breaks=c(0, 0.1, 0.2), labels=percent) +
    scale_color_manual(values=c(
      "Risk score + NMR scores"="#e41a1c", "Risk score + PRS"="#377eb8", "Risk score + Biochemistry"="#ff7f00",
      "Risk score + NMR scores + PRS"="#4daf4a", "Risk score + Biochemistry + PRS"="#984ea3"
    )) +
    scale_fill_manual(values=c(
      "Risk score + NMR scores"="#e41a1c", "Risk score + PRS"="#377eb8", "Risk score + Biochemistry"="#ff7f00",
      "Risk score + NMR scores + PRS"="#4daf4a", "Risk score + Biochemistry + PRS"="#984ea3"
    )) +
    theme_bw() +
    theme(
      axis.text=element_text(size=6), axis.title=element_text(size=8),
      axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
      legend.position="bottom", legend.text=element_text(size=7), legend.title=element_blank()
    )
}

res[, model := paste("Risk score +", model_type)]
res <- res[metric %in% c("NRI+", "NRI-")]

# Compare Pooled cohort to discovery and replication
cohort_comp <- res[endpoint == "cvd" & model_sex == "Sex-stratified" & score == "SCORE2_excl_UKB"]
cohort_comp <- dcast(cohort_comp, metric + model ~ cohort, value.var=c("Estimate", "L95", "U95"))

g1 <- ggplot(cohort_comp) + aes(
    x=Estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=Estimate_pooled, ymin=L95_pooled, ymax=U95_pooled, color=model, fill=model
  ) + xlab("Discovery cohort net % reclassified") + ylab("Pooled cohort net % reclassified")
g1 <- common_gg_parts(g1)

g2 <- ggplot(cohort_comp) + aes(
    x=Estimate_replication, xmin=L95_replication, xmax=U95_replication,
    y=Estimate_pooled, ymin=L95_pooled, ymax=U95_pooled, color=model, fill=model
  ) + xlab("Replication cohort net % reclassified") + ylab("Pooled cohort net % reclassified")
g2 <- common_gg_parts(g2)

g3 <- ggplot(cohort_comp) + aes(
    x=Estimate_discovery, xmin=L95_discovery, xmax=U95_discovery,
    y=Estimate_replication, ymin=L95_replication, ymax=U95_replication, color=model, fill=model
  ) + xlab("Discovery cohort net % reclassified") + ylab("Replication cohort net % reclassified")
g3 <- common_gg_parts(g3)

# Compare Sex-stratified to sex-specific
sex_comp <- res[cohort == "pooled" & endpoint == "cvd" & score == "SCORE2_excl_UKB"]
sex_comp <- dcast(sex_comp, metric + model ~ model_sex, value.var=c("Estimate", "L95", "U95"))

g4 <- ggplot(sex_comp) + aes(
  x=Estimate_Males, xmin=L95_Males, xmax=U95_Males,
  y=`Estimate_Sex-stratified`, ymin=`L95_Sex-stratified`, ymax=`U95_Sex-stratified`,
  color=model, fill=model) + xlab("Net % reclassified in males") + ylab("Sex-stratified net % reclassified")
g4 <- common_gg_parts(g4)

g5 <- ggplot(sex_comp) + aes(
  x=Estimate_Females, xmin=L95_Females, xmax=U95_Females,
  y=`Estimate_Sex-stratified`, ymin=`L95_Sex-stratified`, ymax=`U95_Sex-stratified`,
  color=model, fill=model) + xlab("Net % reclassified in females") + ylab("Sex-stratified net % reclassified")
g5 <- common_gg_parts(g5)

g6 <- ggplot(sex_comp) + aes(
  x=Estimate_Males, xmin=L95_Males, xmax=U95_Males,
  y=Estimate_Females, ymin=L95_Females, ymax=U95_Females,
  color=model, fill=model) + xlab("Net % reclassified in males") + ylab("Net % reclassified in females")
g6 <- common_gg_parts(g6)

# Compare SCORE2 to QRISK3
score_comp <- res[cohort == "pooled" & endpoint == "cvd" & model_sex == "Sex-stratified" & score != "SCORE2"]
score_comp <- dcast(score_comp, metric + model ~ score, value.var=c("Estimate", "L95", "U95"))

g7 <- ggplot(score_comp) + aes(
  x=Estimate_QRISK3, xmin=L95_QRISK3, xmax=U95_QRISK3,
  y=Estimate_SCORE2_excl_UKB, ymin=L95_SCORE2_excl_UKB, ymax=U95_SCORE2_excl_UKB,
  color=model, fill=model) + xlab("Net % reclassified relative to QRISK3") + ylab("Net % reclassified relative to SCORE2")
g7 <- common_gg_parts(g7)

# Compare broad CVD to narrow CVD
cvd_comp <- res[cohort == "pooled" & score == "SCORE2_excl_UKB" & model_sex == "Sex-stratified"]
cvd_comp <- dcast(cvd_comp, metric + model ~ endpoint, value.var=c("Estimate", "L95", "U95"))

g8 <- ggplot(cvd_comp) + aes(
  x=Estimate_cvd_narrow, xmin=L95_cvd_narrow, xmax=U95_cvd_narrow,
  y=Estimate_cvd, ymin=L95_cvd, ymax=U95_cvd, color=model, fill=model) +
  xlab("Narrow CVD definition net % reclassified") + ylab("Broad CVD definition net % reclassified")
g8 <- common_gg_parts(g8)

# Add all the parts together into a single plot
g <- (g1 | g2 | g3) / (g4 | g5 | g6) / (g7 | g8 | g8)
ggsave(g, width=7.2, height=8.1, file="analyses/test/categorical_nri_comparison_sensitivity.pdf")

