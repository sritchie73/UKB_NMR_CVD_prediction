library(data.table)
library(foreach)
library(forcats)
library(ggplot2)
library(ggstance)
library(patchwork)

# Load and collate results
res <- rbind(idcol="cohort",
  "discovery"=fread("analyses/test/discovery_cindices.txt"),
  "replication"=fread("analyses/test/replication_cindices.txt"),
  "pooled"=fread("analyses/test/pooled_cindices.txt")
)

# Build main figure 3A
ggdt <- res[cohort == "pooled" & model_sex == "Sex-stratified" & endpoint == "cvd" & score == "SCORE2_excl_UKB"]
ggdt[, model_name := fcase(
  model == "SCORE2 + NMR scores", "SCORE2 + NMR scores for CHD and IS",
  model == "SCORE2 + Biochemistry", "SCORE2 + 11 clinical chemistry biomarkers",
  model == "SCORE2 + PRS", "SCORE2 + PRSs for CHD and IS",
  model == "SCORE2 + NMR scores + PRS", "SCORE2 + NMR scores for CHD and IS + PRSs for CHD and IS",
  model == "SCORE2 + Biochemistry + PRS", "SCORE2 + 11 clinical chemistry biomarkers + PRSs for CHD and IS"
)]
ggdt[,model_name := factor(model_name, levels=model_name)]

g <- ggplot(ggdt) +
  aes(x=deltaC, xmin=deltaC.L95, xmax=deltaC.U95, y=fct_rev(model_name), color=model) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, size=2, fill="white") +
  scale_x_continuous("Change in C-index (95% CI)", limits=c(0, 0.03), breaks=c(0, 0.01, 0.02, 0.03)) +
	scale_color_manual(values=c(
		"SCORE2 + NMR scores"="#e41a1c", "SCORE2 + PRS"="#377eb8", "SCORE2 + Biochemistry"="#ff7f00",
		"SCORE2 + NMR scores + PRS"="#4daf4a", "SCORE2 + Biochemistry + PRS"="#984ea3"
	)) +
  ylab("") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )
ggsave(g, width=5.3, height=2, file="analyses/test/main_delta_cindices.pdf")

# Build supp table of main results
supp1 <- res[endpoint == "cvd" & model_sex == "Sex-stratified" & score == "SCORE2_excl_UKB"]
supp1[, model := factor(model, levels=unique(model))]
supp1 <- supp1[order(model)]
supp1[, c("endpoint", "model_sex", "score", "model_type") := NULL]
fwrite(supp1, sep="\t", quote=FALSE, file="analyses/test/cindices_main_analysis.txt")

# Build supp table of sensitivity analysis results
supp2 <- res[score != "SCORE2"]
supp2[score == "SCORE2_excl_UKB", score := "SCORE2"]
supp2[, model := gsub("SCORE2", "Risk score", model)]
supp2[, model := gsub("QRISK3", "Risk score", model)]
supp2[, model := factor(model, levels=unique(model))]
supp2[, endpoint := fcase(
  endpoint == "cvd", "Broad",
  endpoint == "cvd_narrow", "Narrow")]
supp2[, model_sex := factor(model_sex, levels=c("Sex-stratified", "Males", "Females"))]
supp2[, cohort := factor(cohort, levels=c("discovery", "replication", "pooled"))]
supp2[, score := factor(score, levels=c("SCORE2", "QRISK3"))]
supp2 <- supp2[order(model_sex)][order(cohort)][order(model)]
fwrite(supp2, sep="\t", quote=FALSE, file="analyses/test/cindices_sensitivity_analysis.txt")

# Run pairwise comparisons to assess similarity of results across different settings
model_info <- unique(res[,.(cohort, endpoint, model_sex, score)])
comp_stats <- foreach(modelIdx1=model_info[,.I], .combine=rbind) %do% {
  foreach(modelIdx2=model_info[,.I], .combine=rbind) %do% {
   if (modelIdx1 == modelIdx2) return(NULL)

    model_info1 <- model_info[modelIdx1]
    model_info2 <- model_info[modelIdx2]

    res1 <- res[model_info1, on=.(cohort, endpoint, model_sex, score)]
    res2 <- res[model_info2, on=.(cohort, endpoint, model_sex, score)]

    res1 <- res1[, .(model_type, deltaC)]
    res2 <- res2[, .(model_type, deltaC)]

    comp <- merge(res1, res2, by="model_type", suffixes=c(".1", ".2"))

    corr <- comp[,cor(deltaC.1, deltaC.2)]
    beta <- comp[,coef(lm(deltaC.2 ~ 0 + deltaC.1))[1]]

    this_comp_stats <- data.table(pearson=corr, lm_beta=beta)

    setnames(model_info1, paste0(names(model_info1), ".1"))
    setnames(model_info2, paste0(names(model_info2), ".2"))
    return(cbind(model_info1, model_info2, this_comp_stats))
  }
}
fwrite(comp_stats, sep="\t", quote=FALSE, file="analyses/test/deltaC_comparison_sensitivity.txt")

# Create plot of key comparisons of interest
common_gg_parts <- function(g) {
  g +
    geom_abline(intercept=0, slope=1, linetype=2) +
    geom_errorbarh(height=0, alpha=0.7) + geom_errorbar(width=0, alpha=0.7) +
    geom_point(shape=23, fill="white") +
    scale_x_continuous(limits=c(0, 0.03), oob=scales::squish, breaks=c(0, 0.01, 0.02, 0.03)) +
    scale_y_continuous(limits=c(0, 0.03), oob=scales::squish, breaks=c(0, 0.01, 0.02, 0.03)) +
		scale_color_manual(values=c(
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

# Compare Pooled cohort to discovery and replication
cohort_comp <- res[endpoint == "cvd" & model_sex == "Sex-stratified" & score == "SCORE2_excl_UKB"]
cohort_comp <- dcast(cohort_comp, model ~ cohort, value.var=c("deltaC", "deltaC.L95", "deltaC.U95"))

g1 <- ggplot(cohort_comp) + aes(
    x=deltaC_discovery, xmin=deltaC.L95_discovery, xmax=deltaC.U95_discovery,
    y=deltaC_pooled, ymin=deltaC.L95_pooled, ymax=deltaC.U95_pooled, color=model
  ) + xlab("Discovery cohort Delta C-index") + ylab("Pooled cohort Delta C-index")
g1 <- common_gg_parts(g1)

g2 <- ggplot(cohort_comp) + aes(
    x=deltaC_replication, xmin=deltaC.L95_replication, xmax=deltaC.U95_replication,
    y=deltaC_pooled, ymin=deltaC.L95_pooled, ymax=deltaC.U95_pooled, color=model
  ) + xlab("Replication cohort Delta C-index") + ylab("Pooled cohort Delta C-index")
g2 <- common_gg_parts(g2)

g3 <- ggplot(cohort_comp) + aes(
    x=deltaC_discovery, xmin=deltaC.L95_discovery, xmax=deltaC.U95_discovery,
    y=deltaC_replication, ymin=deltaC.L95_replication, ymax=deltaC.U95_replication, color=model
  ) + xlab("Discovery cohort Delta C-index") + ylab("Replication cohort Delta C-index")
g3 <- common_gg_parts(g3)

# Compare Sex-stratified to sex-specific
sex_comp <- res[cohort == "pooled" & endpoint == "cvd" & score == "SCORE2_excl_UKB"]
sex_comp <- dcast(sex_comp, model ~ model_sex, value.var=c("deltaC", "deltaC.L95", "deltaC.U95"))

g4 <- ggplot(sex_comp) + aes(
  x=deltaC_Males, xmin=deltaC.L95_Males, xmax=deltaC.U95_Males,
  y=`deltaC_Sex-stratified`, ymin=`deltaC.L95_Sex-stratified`, ymax=`deltaC.U95_Sex-stratified`,
  color=model) + xlab("Delta C-index in males") + ylab("Sex-stratified Delta C-index")
g4 <- common_gg_parts(g4)

g5 <- ggplot(sex_comp) + aes(
  x=deltaC_Females, xmin=deltaC.L95_Females, xmax=deltaC.U95_Females,
  y=`deltaC_Sex-stratified`, ymin=`deltaC.L95_Sex-stratified`, ymax=`deltaC.U95_Sex-stratified`,
  color=model) + xlab("Delta C-index in females") + ylab("Sex-stratified Delta C-index")
g5 <- common_gg_parts(g5)

g6 <- ggplot(sex_comp) + aes(
  x=deltaC_Males, xmin=deltaC.L95_Males, xmax=deltaC.U95_Males,
  y=deltaC_Females, ymin=deltaC.L95_Females, ymax=deltaC.U95_Females,
  color=model) + xlab("Delta C-index in males") + ylab("Delta C-index in females")
g6 <- common_gg_parts(g6)

# Compare SCORE2 to QRISK3
score_comp <- res[cohort == "pooled" & endpoint == "cvd" & model_sex == "Sex-stratified" & score != "SCORE2"]
score_comp <- dcast(score_comp, model ~ score, value.var=c("deltaC", "deltaC.L95", "deltaC.U95"))

g7 <- ggplot(score_comp) + aes(
  x=deltaC_QRISK3, xmin=deltaC.L95_QRISK3, xmax=deltaC.U95_QRISK3,
  y=deltaC_SCORE2_excl_UKB, ymin=deltaC.L95_SCORE2_excl_UKB, ymax=deltaC.U95_SCORE2_excl_UKB,
  color=model) + xlab("Delta C-index relative to QRISK3") + ylab("Delta C-index relative to SCORE2")
g7 <- common_gg_parts(g7)

# Compare broad CVD to narrow CVD
cvd_comp <- res[cohort == "pooled" & score == "SCORE2_excl_UKB" & model_sex == "Sex-stratified"]
cvd_comp <- dcast(cvd_comp, model ~ endpoint, value.var=c("deltaC", "deltaC.L95", "deltaC.U95"))

g8 <- ggplot(cvd_comp) + aes(
  x=deltaC_cvd_narrow, xmin=deltaC.L95_cvd_narrow, xmax=deltaC.U95_cvd_narrow,
  y=deltaC_cvd, ymin=deltaC.L95_cvd, ymax=deltaC.U95_cvd, color=model) +
  xlab("Narrow CVD definition Delta C-index") + ylab("Broad CVD definition Delta C-index")
g8 <- common_gg_parts(g8)

# Add all the parts together into a single plot
g <- (g1 | g2 | g3) / (g4 | g5 | g6) / (g7 | g8 | g8)
ggsave(g, width=7.2, height=8.1, file="analyses/test/deltaC_comparison_sensitivity.pdf")


