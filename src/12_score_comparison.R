library(data.table)
library(foreach)
library(ggplot2)
library(forcats)
library(ggrastr)
library(ggpp)
library(ggh4x)

options("ggrastr.default.dpi" = 300)

# Make output directoru
system("mkdir -p analyses/test", wait=TRUE)

# Compare non-derived NMR scores to each other
dat <- fread("data/cleaned/analysis_cohort.txt", select=c("eid", "sex", "CAD_metaGRS", "Stroke_metaGRS", "SCORE2_excl_UKB"))
setnames(dat, "SCORE2_excl_UKB", "SCORE2")

nmr_scores <- fread("analyses/nmr_score_training/aggregate_test_non_derived_NMR_scores.txt")
dat <- dat[nmr_scores, on = .(eid)]
dat <- melt(dat, id.vars=c("eid", "sex"), variable.name="score")

score_comp <- foreach(this_rn = unique(dat$score), .combine=rbind) %:%
  foreach(this_cn = unique(dat$score), .combine=rbind) %do% {
    this_x <- dat[score == this_rn, .(eid, sex, x_score=score, x_value=value)]
    this_y <- dat[score == this_cn, .(eid, sex, y_score=score, y_value=value)]
    merge(this_x, this_y, by=c("eid", "sex"))
}
score_comp[, x_score := factor(x_score, levels=c("SCORE2", "CAD_NMR_score", "Stroke_NMR_score", "CAD_metaGRS", "Stroke_metaGRS"))]
score_comp[, y_score := factor(y_score, levels=c("SCORE2", "CAD_NMR_score", "Stroke_NMR_score", "CAD_metaGRS", "Stroke_metaGRS"))]

# Get pairwise correlation coefficient statistics:
cor_stats <- score_comp[x_score != y_score, {
  cc <- cor.test(x_value, y_value)
  list(estimate=cc$estimate, L95=cc$conf.int[1], U95=cc$conf.int[2], pval=cc$p.value)
}, by=.(sex, x_score, y_score)]
fwrite(cor_stats, sep="\t", quote=FALSE, file="analyses/test/pairwise_correlations_sex_specific_scores.txt")

# Make pairwise scatterplots
score_scatter <- function(score_comp) {
  cor_anno <- score_comp[,.(label=sprintf("r=%.2f\np=%.2f",
    cor(x_value, y_value),
    cor(x_value, y_value, method="spearman")
  )), by=.(x_score, y_score)]

  ggplot(score_comp) +
    aes(x=x_value, y=y_value) +
    rasterise(geom_point(alpha=0.5, shape=19, size=0.1)) +
    geom_smooth(method="lm", color="red", linetype=2) + 
    geom_text_npc(data=cor_anno, npcx="left", npcy="top", aes(label=label), color="red", size=6*0.352777778) +
    facet_grid2(fct_rev(y_score) ~ x_score, scales="free", independent="all") +
    xlab("") + ylab("") +
    theme_bw() +
    theme(
      axis.text=element_text(size=6), axis.title=element_text(size=8),
      strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
      panel.grid=element_blank(), legend.text=element_text(size=6)
    )
}

g1 <- score_scatter(score_comp[sex == "Male"])
g2 <- score_scatter(score_comp[sex == "Female"])

ggsave(g1, width=7, height=7, file="analyses/test/non_derived_score_compare_males.pdf")
ggsave(g2, width=7, height=7, file="analyses/test/non_derived_score_compare_females.pdf")

# Compare clinical NMR scores to each other
dat <- fread("data/cleaned/analysis_cohort.txt", select=c("eid", "sex", "CAD_metaGRS", "Stroke_metaGRS", "SCORE2_excl_UKB"))
setnames(dat, "SCORE2_excl_UKB", "SCORE2")

nmr_scores <- fread("analyses/nmr_score_training/aggregate_test_clinical_NMR_scores.txt")
dat <- dat[nmr_scores, on = .(eid)]
dat <- melt(dat, id.vars=c("eid", "sex"), variable.name="score")

score_comp <- foreach(this_rn = unique(dat$score), .combine=rbind) %:%
  foreach(this_cn = unique(dat$score), .combine=rbind) %do% {
    this_x <- dat[score == this_rn, .(eid, sex, x_score=score, x_value=value)]
    this_y <- dat[score == this_cn, .(eid, sex, y_score=score, y_value=value)]
    merge(this_x, this_y, by=c("eid", "sex"))
}
score_comp[, x_score := factor(x_score, levels=c("SCORE2", "CAD_NMR_score", "Stroke_NMR_score", "CAD_metaGRS", "Stroke_metaGRS"))]
score_comp[, y_score := factor(y_score, levels=c("SCORE2", "CAD_NMR_score", "Stroke_NMR_score", "CAD_metaGRS", "Stroke_metaGRS"))]

g1 <- score_scatter(score_comp[sex == "Male"])
g2 <- score_scatter(score_comp[sex == "Female"])

ggsave(g1, width=7, height=7, file="analyses/test/clinical_score_compare_males.pdf")
ggsave(g2, width=7, height=7, file="analyses/test/clinical_score_compare_females.pdf")

