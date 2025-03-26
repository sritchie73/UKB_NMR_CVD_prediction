library(data.table)
library(foreach)
library(ggplot2)
library(ggstance)
library(patchwork)

# Load and collate results
res <- rbind(idcol="cohort",
  "discovery"=fread("analyses/univariate/discovery_analysis.txt"),
  "replication"=fread("analyses/univariate/replication_analysis.txt"),
  "pooled"=fread("analyses/univariate/pooled_analysis.txt")
)

# Plot C-index/deltaC-index for top 10 biomarkers in pooled analysis + HRs from discovery analysis
ggdt <- res[cohort == "pooled" & sex == "Sex-stratified" & endpoint == "cvd" & score == "SCORE2_excl_UKB"]
ggdt <- ggdt[order(-deltaC)]
ggdt <- ggdt[model %in% ggdt[, unique(model)][1:10]]
ggdt[res[cohort == "discovery"], on = .(sex, endpoint, score, biomarker), 
  c("biomarker.HR", "biomarker.HR.L95", "biomarker.HR.U95") := 
  .(i.biomarker.HR, i.biomarker.HR.L95, i.biomarker.HR.U95)]
ggdt <- rbind(
  ggdt[,.(model, model_type, metric="HR", estimate=biomarker.HR, L95=biomarker.HR.L95, U95=biomarker.HR.U95)],
  ggdt[,.(model, model_type, metric="deltaC", estimate=deltaC, L95=deltaC.L95, U95=deltaC.U95)]
)
ggdt[, model := factor(model, levels=rev(unique(model)))]
ggdt[, metric := factor(metric, levels=c("HR", "deltaC"))]
ggref <- data.table(metric=c("HR", "HR", "HR", "deltaC", "deltaC"), null=c(0.8, 1, 1.2, 0, 0.01)) 

g <- ggplot(ggdt) +
  aes(x=estimate, xmin=L95, xmax=U95, y=model, color=model_type) +
  facet_grid(~ metric, scales="free_x") +
  geom_vline(data=ggref, aes(xintercept=null), linetype=2) +
  geom_errorbarh(height=0, position=position_dodgev(height=0.6)) +
  geom_point(shape=23, size=2, fill="white", position=position_dodgev(height=0.6)) +
  scale_x_continuous("estimate (95% CI)") +
  scale_color_manual(values=c("assays"="#aa4400", "NMR"="#8800aa")) +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="bottom", legend.text=element_text(size=7), legend.title=element_blank()
  )
ggsave(g, width=7.2, height=3, file="analyses/univariate/top10_biomarkers.pdf")

# Build supp table of main results
supp1 <- res[endpoint == "cvd" & sex == "Sex-stratified" & score == "SCORE2_excl_UKB" & biomarker != ""]
supp1o <- supp1[cohort == "pooled"][order(deltaC.fdr, na.last=FALSE)][,.(biomarker)]
supp1 <- supp1[supp1o, on = .(biomarker)]
supp1[, model_type := fcase(
  model_type == "assays", "Clinical biochemistry assay", 
  model_type == "NMR", "NMR spectroscopy",
  default="-")]
supp1[, c("endpoint", "sex", "score", "biomarker") := NULL]
fwrite(supp1, sep="\t", quote=FALSE, file="analyses/univariate/cindices_main_analysis.txt")

# Build supp table of sensitivity analysis results
supp2 <- res[score != "SCORE2" & biomarker != ""]
supp2[score == "SCORE2_excl_UKB", score := "SCORE2"]
supp2[, model := gsub("SCORE2", "Risk score", model)]
supp2[, model := gsub("QRISK3", "Risk score", model)]
supp2[, model_type := fcase(
  model_type == "assays", "Clinical biochemistry assay", 
  model_type == "NMR", "NMR spectroscopy",
  default="-")]
supp2[, endpoint := fcase(
  endpoint == "cvd", "ASCVD", 
  endpoint == "cvd_narrow", "CHD + Stroke")]
supp2[, sex := factor(sex, levels=c("Sex-stratified", "Males", "Females"))]
supp2[, cohort := factor(cohort, levels=c("discovery", "replication", "pooled"))]
supp2[, score := factor(score, levels=c("SCORE2", "QRISK3"))]
supp2 <- supp2[order(sex)][order(cohort)]
supp2 <- supp2[supp1o, on = .(biomarker)]
fwrite(supp2, sep="\t", quote=FALSE, file="analyses/univariate/cindices_sensitivity_analysis.txt")

# Run pairwise comparisons to assess similarity of results across different settings
model_info <- unique(res[,.(cohort, endpoint, sex, score)])
comp_stats <- foreach(modelIdx1=model_info[,.I], .combine=rbind) %do% {
  foreach(modelIdx2=model_info[,.I], .combine=rbind) %do% {
   if (modelIdx1 == modelIdx2) return(NULL)

    model_info1 <- model_info[modelIdx1]
    model_info2 <- model_info[modelIdx2]

    res1 <- res[model_info1, on=.(cohort, endpoint, sex, score)]
    res2 <- res[model_info2, on=.(cohort, endpoint, sex, score)]
 
    res1 <- res1[biomarker != "", .(model_type, biomarker, deltaC)]
    res2 <- res2[biomarker != "", .(model_type, biomarker, deltaC)]

    comp <- merge(res1, res2, by=c("model_type", "biomarker"), suffixes=c(".1", ".2")) 

    corr1 <- comp[,cor(deltaC.1, deltaC.2)]
    beta1 <- comp[,coef(lm(deltaC.2 ~ 0 + deltaC.1))[1]]

    corr2 <- comp[model_type == "NMR", cor(deltaC.1, deltaC.2)]
    beta2 <- comp[model_type == "NMR", coef(lm(deltaC.2 ~ 0 + deltaC.1))[1]]
    
    corr3 <- comp[model_type == "assays", cor(deltaC.1, deltaC.2)]
    beta3 <- comp[model_type == "assays", coef(lm(deltaC.2 ~ 0 + deltaC.1))[1]]

    this_comp_stats <- rbind(idcol="model_type", 
      "all"=data.table(pearson=corr1, lm_beta=beta1),
      "NMR"=data.table(pearson=corr2, lm_beta=beta2),
      "assays"=data.table(pearson=corr3, lm_beta=beta3))

    setnames(model_info1, paste0(names(model_info1), ".1"))
    setnames(model_info2, paste0(names(model_info2), ".2"))
    return(cbind(model_info1, model_info2, this_comp_stats))
  }
}
fwrite(comp_stats, sep="\t", quote=FALSE, file="analyses/univariate/deltaC_comparison_sensitivity.txt")

# Create plot of key comparisons of interest
common_gg_parts <- function(g) { 
  g + 
    geom_abline(intercept=0, slope=1, linetype=2) +
		geom_errorbarh(height=0, alpha=0.7) + geom_errorbar(width=0, alpha=0.7) +
		geom_point(shape=23, fill="white") +
		scale_x_continuous(limits=c(0, 0.01), oob=scales::squish, breaks=c(0, 0.005, 0.01)) +
		scale_y_continuous(limits=c(0, 0.01), oob=scales::squish, breaks=c(0, 0.005, 0.01)) +
		scale_color_manual(values=c("assays"="#aa4400", "NMR"="#8800aa")) +
		theme_bw() +
		theme(
			axis.text=element_text(size=6), axis.title=element_text(size=8),
			axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
			legend.position="bottom", legend.text=element_text(size=7), legend.title=element_blank()
		)
}

# Compare Pooled cohort to discovery and replication
cohort_comp <- res[endpoint == "cvd" & sex == "Sex-stratified" & score == "SCORE2_excl_UKB" & biomarker != ""]
cohort_comp <- dcast(cohort_comp, model + model_type ~ cohort, value.var=c("deltaC", "deltaC.L95", "deltaC.U95"))

g1 <- ggplot(cohort_comp) + aes(
		x=deltaC_discovery, xmin=deltaC.L95_discovery, xmax=deltaC.U95_discovery,
		y=deltaC_pooled, ymin=deltaC.L95_pooled, ymax=deltaC.U95_pooled, color=model_type
  ) + xlab("Discovery cohort Delta C-index") + ylab("Pooled cohort Delta C-index")
g1 <- common_gg_parts(g1)

g2 <- ggplot(cohort_comp) + aes(
		x=deltaC_replication, xmin=deltaC.L95_replication, xmax=deltaC.U95_replication,
		y=deltaC_pooled, ymin=deltaC.L95_pooled, ymax=deltaC.U95_pooled, color=model_type
  ) + xlab("Replication cohort Delta C-index") + ylab("Pooled cohort Delta C-index")
g2 <- common_gg_parts(g2)

g3 <- ggplot(cohort_comp) + aes(
    x=deltaC_discovery, xmin=deltaC.L95_discovery, xmax=deltaC.U95_discovery,
		y=deltaC_replication, ymin=deltaC.L95_replication, ymax=deltaC.U95_replication, color=model_type
  ) + xlab("Discovery cohort Delta C-index") + ylab("Replication cohort Delta C-index")
g3 <- common_gg_parts(g3)

# Compare Sex-stratified to sex-specific
sex_comp <- res[cohort == "pooled" & endpoint == "cvd" & score == "SCORE2_excl_UKB" & biomarker != ""]
sex_comp <- dcast(sex_comp, model + model_type ~ sex, value.var=c("deltaC", "deltaC.L95", "deltaC.U95"))

g4 <- ggplot(sex_comp) + aes(
  x=deltaC_Males, xmin=deltaC.L95_Males, xmax=deltaC.U95_Males,
  y=`deltaC_Sex-stratified`, ymin=`deltaC.L95_Sex-stratified`, ymax=`deltaC.U95_Sex-stratified`,
  color=model_type) + xlab("Delta C-index in males") + ylab("Sex-stratified Delta C-index")
g4 <- common_gg_parts(g4)

g5 <- ggplot(sex_comp) + aes(
  x=deltaC_Females, xmin=deltaC.L95_Females, xmax=deltaC.U95_Females,
  y=`deltaC_Sex-stratified`, ymin=`deltaC.L95_Sex-stratified`, ymax=`deltaC.U95_Sex-stratified`,
  color=model_type) + xlab("Delta C-index in females") + ylab("Sex-stratified Delta C-index")
g5 <- common_gg_parts(g5)

g6 <- ggplot(sex_comp) + aes(
  x=deltaC_Males, xmin=deltaC.L95_Males, xmax=deltaC.U95_Males,
  y=deltaC_Females, ymin=deltaC.L95_Females, ymax=deltaC.U95_Females,
  color=model_type) + xlab("Delta C-index in males") + ylab("Delta C-index in females")
g6 <- common_gg_parts(g6)

# Compare SCORE2 to QRISK3
score_comp <- res[cohort == "pooled" & endpoint == "cvd" & sex == "Sex-stratified" & biomarker != "" & score != "SCORE2"]
score_comp <- dcast(score_comp, biomarker + model_type ~ score, value.var=c("deltaC", "deltaC.L95", "deltaC.U95"))

g7 <- ggplot(score_comp) + aes(
  x=deltaC_QRISK3, xmin=deltaC.L95_QRISK3, xmax=deltaC.U95_QRISK3,
  y=deltaC_SCORE2_excl_UKB, ymin=deltaC.L95_SCORE2_excl_UKB, ymax=deltaC.U95_SCORE2_excl_UKB,
  color=model_type) + xlab("Delta C-index relative to QRISK3") + ylab("Delta C-index relative to SCORE2")
g7 <- common_gg_parts(g7)

# Compare broad CVD to narrow CVD
cvd_comp <- res[cohort == "pooled" & score == "SCORE2_excl_UKB" & sex == "Sex-stratified" & biomarker != ""]
cvd_comp <- dcast(cvd_comp, biomarker + model_type ~ endpoint, value.var=c("deltaC", "deltaC.L95", "deltaC.U95"))

g8 <- ggplot(cvd_comp) + aes(
  x=deltaC_cvd_narrow, xmin=deltaC.L95_cvd_narrow, xmax=deltaC.U95_cvd_narrow,
  y=deltaC_cvd, ymin=deltaC.L95_cvd, ymax=deltaC.U95_cvd, color=model_type) +
  xlab("Narrow CVD definition Delta C-index") + ylab("Broad CVD definition Delta C-index")
g8 <- common_gg_parts(g8)
  
# Add all the parts together into a single plot
g <- (g1 | g2 | g3) / (g4 | g5 | g6) / (g7 | g8 | g8)
ggsave(g, width=7.2, height=8, file="analyses/univariate/deltaC_comparison_sensitivity.pdf")

