library(data.table)
library(ggplot2)
library(ggh4x)
library(cowplot)

# Load results stratified by sex
blanket <- fread("analyses/public_health_modelling/blanket_screening/population_screening_by_sex.txt")
targeted <- fread("analyses/public_health_modelling/targeted_screening/population_screening_by_sex.txt")

# Combine
dat <- rbind(idcol="strategy", "blanket"=blanket, "targeted"=targeted)

# Format columns of interest
shownum <- function(num) { format(num, big.mark=",", trim=TRUE) }
dat[, text := sprintf("%s (%s-%s)", shownum(estimate), shownum(L95), shownum(U95))]
dat[, people := shownum(people)]
dat[, cases := shownum(cases)]
dat <- dcast(dat, sex + people + cases + strategy + model ~ number, value.var="text")
dat <- dat[, .(sex, people, cases, strategy, model, high_risk, high_risk_cases, events_prevented, NNS, NNT, delta_high_risk_cases, delta_events_prevented)]

# Impose ordering
dat[, model := factor(model, levels=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
dat[, sex := factor(sex, levels=c("Males", "Females"))]
dat[, strategy := factor(strategy, levels=c("blanket", "targeted"))]
dat <- dat[order(model)][order(strategy)][order(sex)]

# Write out
fwrite(dat, sep="\t", quote=FALSE, file="analyses/public_health_modelling/screening_summary_by_sex.txt")

# Also plot
dat <- rbind(idcol="strategy", "blanket"=blanket, "targeted"=targeted)
dat[, metric := fcase(
  number == "high_risk", "People classified as high risk",
  number == "high_risk_cases", "CVD cases classified as high risk",
  number == "events_prevented", "CVD events prevented by statin prescription",
  number == "NNS", "Number needed to screen to prevent 1 CVD event",
  number == "NNT", "Number of statins prescribed per CVD event prevented",
  default = NA
)]
dat <- dat[!is.na(metric)]
dat[, metric := factor(metric, levels=unique(metric))]
dat[, strategy := fcase(
  strategy == "blanket", "Population screening",
  strategy == "targeted", "Targeted screening"
)]
dat[, model := factor(model, levels=rev(c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs")))]

ref <- dat[model == "SCORE2"]
ref[, strategy := NULL]

g1 <- ggplot(dat[sex == "Males"]) + 
  aes(x=estimate, xmin=L95, xmax=U95, y=model, color=metric) +
  facet_grid(strategy ~ metric, scales="free", labeller=labeller(strategy=label_wrap_gen(width=10), metric=label_wrap_gen(width=20))) +
  geom_vline(data=ref[sex == "Males"], aes(xintercept=estimate), linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, fill="white", size=1.2) +
  scale_x_continuous("Number (95% CI)", n.breaks=4) +
  theme_bw() + 
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.background=element_blank(), strip.text=element_text(size=6, face="bold"),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(), 
    legend.position="none"
  )

g2 <- ggplot(dat[sex == "Females"]) + 
  aes(x=estimate, xmin=L95, xmax=U95, y=model, color=metric) +
  facet_grid(strategy ~ metric, scales="free", labeller=labeller(strategy=label_wrap_gen(width=10), metric=label_wrap_gen(width=20))) +
  geom_vline(data=ref[sex == "Females"], aes(xintercept=estimate), linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, fill="white", size=1.2) +
  scale_x_continuous("Number (95% CI)", n.breaks=4) +
  theme_bw() + 
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.background=element_blank(), strip.text=element_text(size=6, face="bold"),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(), 
    legend.position="none"
  )

g <- plot_grid(g1, g2, labels=c("Males  ", "Females"), label_size=10, ncol=1)
ggsave(g, width=7.2, height=4, file="analyses/public_health_modelling/screening_comparison_by_sex.pdf")

# Load results not stratified by sex
blanket <- fread("analyses/public_health_modelling/blanket_screening/population_screening.txt")
targeted <- fread("analyses/public_health_modelling/targeted_screening/population_screening.txt")

# Combine
dat <- rbind(idcol="strategy", "blanket"=blanket, "targeted"=targeted)

# Format columns of interest
dat[, text := sprintf("%s (%s-%s)", shownum(estimate), shownum(L95), shownum(U95))]
dat[, people := shownum(people)]
dat[, cases := shownum(cases)]
dat <- dcast(dat, people + cases + strategy + model ~ number, value.var="text")
dat <- dat[, .(people, cases, strategy, model, high_risk, high_risk_cases, events_prevented, NNS, NNT, delta_high_risk_cases, delta_events_prevented)]

# Impose ordering
dat[, model := factor(model, levels=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
dat[, strategy := factor(strategy, levels=c("blanket", "targeted"))]
dat <- dat[order(model)][order(strategy)]

# Write out
fwrite(dat, sep="\t", quote=FALSE, file="analyses/public_health_modelling/screening_summary.txt")

# Build supp tables
blanket <- fread("analyses/public_health_modelling/blanket_screening/population_screening_by_sex.txt")
targeted <- fread("analyses/public_health_modelling/targeted_screening/population_screening_by_sex.txt")

dt1 <- dcast(blanket, sex + people + cases + model ~ number, value.var=c("estimate", "L95", "U95"))
dt1[, sex := factor(sex, levels=c("Males", "Females"))]
dt1[, model := factor(model, levels=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
dt1 <- dt1[order(model)][order(sex)]
dt1 <- dt1[,.(sex, people, cases, model, 
  estimate_high_risk, L95_high_risk, U95_high_risk,
  estimate_high_risk_cases, L95_high_risk_cases, U95_high_risk_cases,
  estimate_events_prevented, L95_events_prevented, U95_events_prevented,
  estimate_NNS, L95_NNS, U95_NNS,
  estimate_NNT, L95_NNT, U95_NNT
)]
fwrite(dt1, sep="\t", quote=FALSE, file="analyses/public_health_modelling/blanket_screening_supp.txt")

dt2 <- dcast(targeted, sex + people + cases + model ~ number, value.var=c("estimate", "L95", "U95"))
dt2[, sex := factor(sex, levels=c("Males", "Females"))]
dt2[, model := factor(model, levels=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
dt2 <- dt2[order(model)][order(sex)]
dt2 <- dt2[,.(sex, people, cases, model, 
  estimate_high_risk, L95_high_risk, U95_high_risk,
  estimate_high_risk_cases, L95_high_risk_cases, U95_high_risk_cases,
  estimate_events_prevented, L95_events_prevented, U95_events_prevented,
  estimate_NNS, L95_NNS, U95_NNS,
  estimate_NNT, L95_NNT, U95_NNT
)]
fwrite(dt2, sep="\t", quote=FALSE, file="analyses/public_health_modelling/targeted_screening_supp.txt")

dt3 <- dcast(targeted, sex + people + cases + model ~ number, value.var=c("estimate", "L95", "U95"))
dt3[, sex := factor(sex, levels=c("Males", "Females"))]
dt3[, model := factor(model, levels=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
dt3 <- dt3[order(model)][order(sex)]
dt3 <- dt3[,.(sex,
  estimate_medium_risk_score2, L95_medium_risk_score2, U95_medium_risk_score2,
  estimate_medium_risk_score2_cases, L95_medium_risk_score2_cases, U95_medium_risk_score2_cases,
  model, 
  estimate_reclassified, L95_reclassified, U95_reclassified,
  estimate_reclassified_cases, L95_reclassified_cases, U95_reclassified_cases,
  estimate_delta_events_prevented, L95_delta_events_prevented, U95_delta_events_prevented
)]
fwrite(dt3, sep="\t", quote=FALSE, file="analyses/public_health_modelling/targeted_screening2_supp.txt")  

# Build supp tables split into each age group
dt1 <- fread("analyses/public_health_modelling/blanket_screening/population_screening_by_sex_and_age_group.txt")
dt2 <- fread("analyses/public_health_modelling/targeted_screening/population_screening_by_sex_and_age_group.txt")
dt2[, model := paste("SCORE2 then", model)]
dt <- rbind(idcol="strategy", "population"=dt1, "targeted"=dt2)
dt[, sex := factor(sex, levels=c("Males", "Females"))]

dt1 <- dcast(dt[number == "high_risk_cases"], sex + age_group + cases ~ model, value.var=c("estimate", "L95", "U95"))
dt1 <- dt1[,.(sex, age_group, cases, 
  estimate_SCORE2, L95_SCORE2, U95_SCORE2,
  `estimate_SCORE2 + NMR scores`, `L95_SCORE2 + NMR scores`, `U95_SCORE2 + NMR scores`,
  `estimate_SCORE2 + PRSs`, `L95_SCORE2 + PRSs`, `U95_SCORE2 + PRSs`,
  `estimate_SCORE2 + NMR scores + PRSs`, `L95_SCORE2 + NMR scores + PRSs`, `U95_SCORE2 + NMR scores + PRSs`
)]
fwrite(dt1, sep="\t", quote=FALSE, file="analyses/public_health_modelling/high_risk_cases_by_sex_and_age_for_supp.txt")

dt2 <- dcast(dt[number == "high_risk_non_cases"], sex + age_group + non_cases ~ model, value.var=c("estimate", "L95", "U95"))
dt2 <- dt2[,.(sex, age_group, non_cases, 
  estimate_SCORE2, L95_SCORE2, U95_SCORE2,
  `estimate_SCORE2 + NMR scores`, `L95_SCORE2 + NMR scores`, `U95_SCORE2 + NMR scores`,
  `estimate_SCORE2 + PRSs`, `L95_SCORE2 + PRSs`, `U95_SCORE2 + PRSs`,
  `estimate_SCORE2 + NMR scores + PRSs`, `L95_SCORE2 + NMR scores + PRSs`, `U95_SCORE2 + NMR scores + PRSs`
)]
fwrite(dt2, sep="\t", quote=FALSE, file="analyses/public_health_modelling/high_risk_non_cases_by_sex_and_age_for_supp.txt")

dt3 <- dcast(dt[number == "medium_risk_score2_cases" & model == "SCORE2 then SCORE2 + NMR scores"], age_group ~ sex, value.var=c("estimate", "L95", "U95")) 
dt3 <- dt3[, .(age_group, estimate_Males, L95_Males, U95_Males, estimate_Females, L95_Females, U95_Females)]
fwrite(dt3, sep="\t", quote=FALSE, file="analyses/public_health_modelling/intermediate_risk_score2_cases_by_sex_and_age_for_supp.txt")

dt4 <- dcast(dt[number == "medium_risk_score2_non_cases" & model == "SCORE2 then SCORE2 + NMR scores"], age_group ~ sex, value.var=c("estimate", "L95", "U95")) 
dt4 <- dt4[, .(age_group, estimate_Males, L95_Males, U95_Males, estimate_Females, L95_Females, U95_Females)]
fwrite(dt4, sep="\t", quote=FALSE, file="analyses/public_health_modelling/intermediate_risk_score2_non_cases_by_sex_and_age_for_supp.txt")

dt5 <- dcast(dt[number == "reclassified_cases"], sex + age_group ~ model, value.var=c("estimate", "L95", "U95"))
dt5 <- dt5[,.(sex, age_group, 
  `estimate_SCORE2 then SCORE2 + NMR scores`, `L95_SCORE2 then SCORE2 + NMR scores`, `U95_SCORE2 then SCORE2 + NMR scores`,
  `estimate_SCORE2 then SCORE2 + PRSs`, `L95_SCORE2 then SCORE2 + PRSs`, `U95_SCORE2 then SCORE2 + PRSs`,
  `estimate_SCORE2 then SCORE2 + NMR scores + PRSs`, `L95_SCORE2 then SCORE2 + NMR scores + PRSs`, `U95_SCORE2 then SCORE2 + NMR scores + PRSs`
)]
fwrite(dt5, sep="\t", quote=FALSE, file="analyses/public_health_modelling/reclassified_cases_by_sex_and_age_for_supp.txt")

dt6 <- dcast(dt[number == "reclassified_non_cases"], sex + age_group ~ model, value.var=c("estimate", "L95", "U95"))
dt6 <- dt6[,.(sex, age_group, 
  `estimate_SCORE2 then SCORE2 + NMR scores`, `L95_SCORE2 then SCORE2 + NMR scores`, `U95_SCORE2 then SCORE2 + NMR scores`,
  `estimate_SCORE2 then SCORE2 + PRSs`, `L95_SCORE2 then SCORE2 + PRSs`, `U95_SCORE2 then SCORE2 + PRSs`,
  `estimate_SCORE2 then SCORE2 + NMR scores + PRSs`, `L95_SCORE2 then SCORE2 + NMR scores + PRSs`, `U95_SCORE2 then SCORE2 + NMR scores + PRSs`
)]
fwrite(dt6, sep="\t", quote=FALSE, file="analyses/public_health_modelling/reclassified_non_cases_by_sex_and_age_for_supp.txt")







