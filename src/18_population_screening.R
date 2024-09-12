library(data.table)
library(foreach)
library(survival)
library(boot)
options(boot.parallel="multicore")
options(boot.ncpus=10)
library(ggplot2)
library(forcats)
library(cowplot)
library(doMC)

# create output directory
system("mkdir -p analyses/public_health_modelling")

# Load in predicted risks
pred_risk <- fread("analyses/CVD_weight_training/CVD_linear_predictors_and_risk.txt")
pred_risk <- pred_risk[score_type == "non-derived"]

# Allocate to ESC 2021 risk groups
pred_risk[, risk_group := fcase(
  age < 50 & uk_calibrated_risk < 0.025, "Low-to-moderate risk",
  age < 50 & uk_calibrated_risk >= 0.025 & uk_calibrated_risk < 0.075, "High risk",
  age < 50 & uk_calibrated_risk >= 0.075, "Very high risk",
  age >= 50 & uk_calibrated_risk < 0.05, "Low-to-moderate risk",
  age >= 50 & uk_calibrated_risk >= 0.05 & uk_calibrated_risk < 0.10, "High risk",
  age >= 50 & uk_calibrated_risk >= 0.10, "Very high risk"
)]

# Reformat to wide table so all models have same bootstraps applied
models <- unique(pred_risk[,.(model)])
models <- models[model != "SCORE2"]
models[, colname := paste0("model", .I)]
pred_risk[, colname := model]
pred_risk[models, on = .(model), colname := i.colname]
pred_risk <- dcast(pred_risk, eid + sex + age + age_group + incident_cvd_followup + incident_cvd ~ colname, value.var="risk_group")

# The function called by the bootstrap procedure must return a flat vector, so we need to
# keep track of the statistics separately
boot_stats_info <- as.data.table(expand.grid(
  sex=c("Male", "Female"),
  age_group=c("40-44", "45-49", "50-54", "55-59", "60-64", "65-69"),
  incident_cvd=c(TRUE, FALSE),
  model=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"),
  stats=c("pct_v.high_risk", "pct_high_risk", "pct_reclassified_to_v.high_risk")
))

# Sanity check loop ordering:
check_loop <- foreach(this_stat=c("pct_v.high_risk", "pct_high_risk", "pct_reclassified_to_v.high_risk"), .combine=rbind) %:% 
  foreach(this_model=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %:% 
    foreach(this_cvd=c(TRUE, FALSE), .combine=rbind) %:% 
      foreach(this_age_group=c("40-44", "45-49", "50-54", "55-59", "60-64", "65-69"), .combine=rbind) %:%
        foreach(this_sex=c("Male", "Female"), .combine=rbind) %do% {
   data.table(sex=this_sex, age_group=this_age_group, incident_cvd=this_cvd, model=this_model, stats=this_stat)
}
check_loop[, loop_row_order := .I]
boot_stats_info[, info_row_order := .I]
check_loop[boot_stats_info, on=.(sex, age_group, incident_cvd, model, stats), info_row_order := i.info_row_order]
stopifnot(all(check_loop$loop_row_order == check_loop$info_row_order))
boot_stats_info[, info_row_order := NULL]

# Create bootstrap function which returns a flat vector
boot_func <- function(dt) {
  foreach(this_stat=c("pct_v.high_risk", "pct_high_risk", "pct_reclassified_to_v.high_risk"), .combine=c) %:% 
    foreach(this_model=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=c) %:% 
      foreach(this_cvd=c(TRUE, FALSE), .combine=c) %:% 
        foreach(this_age_group=c("40-44", "45-49", "50-54", "55-59", "60-64", "65-69"), .combine=c) %:%
          foreach(this_sex=c("Male", "Female"), .combine=c) %do% {

    # Filter input dataset:
    this_dt <- dt[sex == this_sex & age_group == this_age_group & !xor(incident_cvd, this_cvd)]

    # Filter to subset classified as high risk by SCORE2 if relevant
    if (this_stat == "pct_reclassified_to_v.high_risk") {
      this_dt <- this_dt[SCORE2 == "High risk"]
    }

    # Return requested statistic, doing look up by column:
    if (this_stat == "pct_reclassified_to_v.high_risk" && this_model == "SCORE2") {
      return(1) # For comparison, this is where we re-allocate all high risk to very high risk
    } else if (this_stat %in% c("pct_reclassified_to_v.high_risk", "pct_v.high_risk")) {
      if (this_model == "SCORE2") {
        return(this_dt[, sum(SCORE2 == "Very high risk")/.N])
      } else if (this_model == "SCORE2 + NMR scores") {
        return(this_dt[, sum(model1 == "Very high risk")/.N])
      } else if (this_model == "SCORE2 + PRSs") {
        return(this_dt[, sum(model2 == "Very high risk")/.N])
      } else if (this_model == "SCORE2 + NMR scores + PRSs") {
        return(this_dt[, sum(model3 == "Very high risk")/.N])
      }
    } else if (this_stat == "pct_high_risk") {
      if (this_model == "SCORE2") {
        return(this_dt[, sum(SCORE2 == "High risk")/.N])
      } else if (this_model == "SCORE2 + NMR scores") {
        return(this_dt[, sum(model1 == "High risk")/.N])
      } else if (this_model == "SCORE2 + PRSs") {
        return(this_dt[, sum(model2 == "High risk")/.N])
      } else if (this_model == "SCORE2 + NMR scores + PRSs") {
        return(this_dt[, sum(model3 == "High risk")/.N])
      }
    }

    return(NA) # indicates error in code above
  }
}

# Run bootstrap analysis
surv_cols_idx <- match(c("incident_cvd_followup", "incident_cvd"), names(pred_risk))
boot_res <- censboot(pred_risk, boot_func, 1000, index=surv_cols_idx) 
saveRDS(boot_res, "analyses/public_health_modelling/bootstraps.rds")

# Extract bootstrap statistics
boot_stats <- foreach(this_bootstrap = 0:1000, .combine=rbind) %do% {
  this_res <- cbind(bootstrap=this_bootstrap, boot_stats_info)
  if(this_bootstrap == 0) {
    this_res[, value := boot_res$t0]
  } else {
    this_res[, value := boot_res$t[this_bootstrap,]]
  }
  return(this_res)
}

# Load hypothetical population 
ons_pop <- fread("analyses/public_health_modelling/ONS_hypothetical_100k_pop_by_age_sex.txt")

# For each model, stratify people into risk groups
pop_boot <- foreach(this_bootstrap = 0:1000, .combine=rbind) %:% 
  foreach(this_model = c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %:%
    foreach(this_strategy = c("population-wide", "targeted"), .combine=rbind) %do% {

  this_boot <- cbind(bootstrap=this_bootstrap, model=this_model, strategy=this_strategy, ons_pop) 
  
  if (this_strategy == "population-wide") {
    this_boot_stats <- boot_stats[this_bootstrap == bootstrap & model == this_model]
  } else {
    this_boot_stats <- boot_stats[this_bootstrap == bootstrap & model == "SCORE2"] 
  }

  this_boot[this_boot_stats[stats == "pct_v.high_risk" & incident_cvd], on = .(sex, age_group), v.high_risk_cases := cases * i.value]
  this_boot[this_boot_stats[stats == "pct_high_risk" & incident_cvd], on = .(sex, age_group), high_risk_cases := cases * i.value]

  this_boot[this_boot_stats[stats == "pct_v.high_risk" & !(incident_cvd)], on = .(sex, age_group), v.high_risk_controls := controls * i.value]
  this_boot[this_boot_stats[stats == "pct_high_risk" & !(incident_cvd)], on = .(sex, age_group), high_risk_controls := controls * i.value]

  if (this_strategy == "targeted") {
    this_boot_stats <- boot_stats[this_bootstrap == bootstrap & model == this_model & stats == "pct_reclassified_to_v.high_risk"]
    this_boot[this_boot_stats[(incident_cvd)], on = .(sex, age_group), reclassified_cases := high_risk_cases * i.value]
    this_boot[this_boot_stats[!(incident_cvd)], on = .(sex, age_group), reclassified_controls := high_risk_controls * i.value]
    this_boot[, reclassified := reclassified_cases + reclassified_controls]
  } else {
    this_boot[, reclassified_cases := 0]
    this_boot[, reclassified_controls := 0]
    this_boot[, reclassified := 0]
  }
  this_boot[, not_reclassified_cases := high_risk_cases - reclassified_cases]
  this_boot[, not_reclassified_controls := high_risk_controls - reclassified_controls]
  this_boot[, not_reclassified := not_reclassified_cases + not_reclassified_controls]

  return(this_boot)
}

# Compute additional statistics
pop_boot[, cases_treated := v.high_risk_cases + reclassified_cases]
pop_boot[, pct_cases_treated := cases_treated / cases]
pop_boot[, N_treated := cases_treated + v.high_risk_controls + reclassified_controls]
pop_boot[, events_prevented := (v.high_risk_cases + reclassified_cases) * 0.2] # assuming 20% reduction due to statins
pop_boot[, additional_prevented := reclassified_cases * 0.2] # Just set of additional events prevented by targeted screening
pop_boot[, NNS := N / events_prevented] # Number needed to screen per event prevented
pop_boot[, NNT := N_treated / events_prevented] # Number of statins prescribed per event prevented

# Write out
fwrite(pop_boot, sep="\t", quote=FALSE, file="analyses/public_health_modelling/simulated_pop_stratified_with_bootstraps.txt")

# Compute aggregate statistics, bootstrap SE, 95% CI and P-values in the whole simulated population
agg_estimates <- melt(pop_boot, id.vars=c("bootstrap", "model", "strategy", "sex", "age_group"), variable.name="metric")
agg_estimates <- agg_estimates[!(metric %in% c("NNS", "NNT", "pct_cases_treated"))] # Need recomputing after summing
agg_estimates <- agg_estimates[, .(value=sum(value)), by=.(bootstrap, model, strategy, metric)]
agg_estimates <- dcast(agg_estimates, bootstrap + model + strategy ~ metric, value.var="value")
agg_estimates[, pct_cases_treated := cases_treated / cases]
agg_estimates[, NNS := N / events_prevented]
agg_estimates[, NNT := N_treated / events_prevented]

agg_estimates <- melt(agg_estimates, id.vars=c("bootstrap", "model", "strategy", "N", "cases", "controls"), variable.name="metric")
agg_ref <- agg_estimates[model == "SCORE2" & strategy == "population-wide"]
agg_deltas <- copy(agg_estimates)

agg_estimates <- agg_estimates[, .(estimate=value[1], SE=sd(value[-1])), by=.(model, strategy, metric)]
agg_estimates[, L95 := estimate - qnorm(1-(0.05/2))*SE]
agg_estimates[, U95 := estimate + qnorm(1-(0.05/2))*SE]

agg_deltas[agg_ref, on = .(bootstrap, metric), ref := i.value]
agg_deltas[, delta := value - ref]
agg_deltas <- agg_deltas[, .(delta=delta[1], delta.SE=sd(delta[-1])), by=.(model, strategy, metric)]
agg_deltas[, delta.L95 := delta - qnorm(1-(0.05/2))*delta.SE]
agg_deltas[, delta.U95 := delta + qnorm(1-(0.05/2))*delta.SE]
agg_deltas[, delta.pval := pmin(1, pnorm(abs(delta/delta.SE), lower.tail=FALSE)*2)]

agg_estimates[agg_deltas, on = .(model, strategy, metric), c("delta", "delta.SE", "delta.L95", "delta.U95", "delta.pval") := .(delta, delta.SE, delta.L95, delta.U95, delta.pval)]
agg_estimates[model == "SCORE2" & strategy == "population-wide", c("delta", "delta.SE", "delta.L95", "delta.U95", "delta.pval") := NA]

fwrite(agg_estimates, sep="\t", quote=FALSE, file="analyses/public_health_modelling/simulated_pop_stratified.txt")

# Compute sex-specific aggregate statistics, bootstrap SE, 95% CI and P-values
sex_estimates <- melt(pop_boot, id.vars=c("bootstrap", "model", "strategy", "sex", "age_group"), variable.name="metric")
sex_estimates <- sex_estimates[!(metric %in% c("NNS", "NNT", "pct_cases_treated"))] # Need recomputing after summing
sex_estimates <- sex_estimates[, .(value=sum(value)), by=.(bootstrap, model, strategy, sex, metric)]
sex_estimates <- dcast(sex_estimates, bootstrap + model + strategy + sex ~ metric, value.var="value")
sex_estimates[, pct_cases_treated := cases_treated / cases]
sex_estimates[, NNS := N / events_prevented]
sex_estimates[, NNT := N_treated / events_prevented]

sex_estimates <- melt(sex_estimates, id.vars=c("bootstrap", "model", "strategy", "sex", "N", "cases", "controls"), variable.name="metric")
sex_ref <- sex_estimates[model == "SCORE2" & strategy == "population-wide"]
sex_deltas <- copy(sex_estimates)

sex_estimates <- sex_estimates[, .(estimate=value[1], SE=sd(value[-1])), by=.(model, strategy, sex, metric)]
sex_estimates[, L95 := estimate - qnorm(1-(0.05/2))*SE]
sex_estimates[, U95 := estimate + qnorm(1-(0.05/2))*SE]

sex_deltas[sex_ref, on = .(bootstrap, sex, metric), ref := i.value]
sex_deltas[, delta := value - ref]
sex_deltas <- sex_deltas[, .(delta=delta[1], delta.SE=sd(delta[-1])), by=.(model, strategy, sex, metric)]
sex_deltas[, delta.L95 := delta - qnorm(1-(0.05/2))*delta.SE]
sex_deltas[, delta.U95 := delta + qnorm(1-(0.05/2))*delta.SE]
sex_deltas[, delta.pval := pmin(1, pnorm(abs(delta/delta.SE), lower.tail=FALSE)*2)]

sex_estimates[sex_deltas, on = .(model, strategy, sex, metric), c("delta", "delta.SE", "delta.L95", "delta.U95", "delta.pval") := .(delta, delta.SE, delta.L95, delta.U95, delta.pval)]
sex_estimates[model == "SCORE2" & strategy == "population-wide", c("delta", "delta.SE", "delta.L95", "delta.U95", "delta.pval") := NA]

fwrite(sex_estimates, sep="\t", quote=FALSE, file="analyses/public_health_modelling/simulated_pop_stratified_sex_specific.txt")

# Build supp tables
long <- rbind(agg_estimates, sex_estimates, fill=TRUE)
long[is.na(sex), sex := "Combined"]

dt <- foreach(this_strategy = c("population-wide", "targeted"), .combine=rbind) %:%
  foreach(this_sex = c("Combined", "Male", "Female"), .combine=rbind) %:%
    foreach(this_model = c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %do% {
      data.table(strategy=this_strategy, sex=this_sex, model=this_model)
}

dt[long[metric == "N_treated"], on = .(strategy, sex, model),
  c("N_treated", "N_treated.SE", "N_treated.L95", "N_treated.U95", "delta.N_treated", "delta.N_treated.SE", "delta.N_treated.L95", "delta.N_treated.U95", "delta.N_treated.pval") :=
  .(i.estimate, i.SE, i.L95, i.U95, i.delta, i.delta.SE, i.delta.L95, i.delta.U95, i.delta.pval)
]

dt[long[metric == "cases_treated"], on = .(strategy, sex, model),
  c("cases_treated", "cases_treated.SE", "cases_treated.L95", "cases_treated.U95", "delta.cases_treated", "delta.cases_treated.SE", "delta.cases_treated.L95", "delta.cases_treated.U95", "delta.cases_treated.pval") :=
  .(i.estimate, i.SE, i.L95, i.U95, i.delta, i.delta.SE, i.delta.L95, i.delta.U95, i.delta.pval)
]

dt[long[metric == "events_prevented"], on = .(strategy, sex, model),
  c("events_prevented", "events_prevented.SE", "events_prevented.L95", "events_prevented.U95", "delta.events_prevented", "delta.events_prevented.SE", "delta.events_prevented.L95", "delta.events_prevented.U95", "delta.events_prevented.pval") :=
  .(i.estimate, i.SE, i.L95, i.U95, i.delta, i.delta.SE, i.delta.L95, i.delta.U95, i.delta.pval)
]

dt[long[metric == "NNS"], on = .(strategy, sex, model),
  c("NNS", "NNS.SE", "NNS.L95", "NNS.U95", "delta.NNS", "delta.NNS.SE", "delta.NNS.L95", "delta.NNS.U95", "delta.NNS.pval") :=
  .(i.estimate, i.SE, i.L95, i.U95, i.delta, i.delta.SE, i.delta.L95, i.delta.U95, i.delta.pval)
]

dt[long[metric == "NNT"], on = .(strategy, sex, model),
  c("NNT", "NNT.SE", "NNT.L95", "NNT.U95", "delta.NNT", "delta.NNT.SE", "delta.NNT.L95", "delta.NNT.U95", "delta.NNT.pval") :=
  .(i.estimate, i.SE, i.L95, i.U95, i.delta, i.delta.SE, i.delta.L95, i.delta.U95, i.delta.pval)
]

fwrite(dt, sep="\t", quote=FALSE, file="analyses/public_health_modelling/screening_table_for_supp.txt")

# Create plot for main figure
ggdt <- agg_estimates[metric %in% c("N_treated", "cases_treated", "events_prevented", "NNS", "NNT")]
ggdt[, metric := factor(metric, levels=c("N_treated", "cases_treated", "events_prevented", "NNS", "NNT"))]
ggdt[, model := factor(model, levels=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]

g <- ggplot(ggdt[model != "SCORE2"]) + 
  aes(x=delta, xmin=delta.L95, xmax=delta.U95, y=fct_rev(model), color=metric) +
  facet_grid(strategy ~ metric, scales="free") + 
  geom_vline(xintercept=0, linetype=2) + 
  geom_errorbarh(height=0) + 
  geom_point(shape=23, fill="white", size=1.2) +
  scale_color_manual(values=c("N_treated"="#f8766d", "cases_treated"="#a3a500", "events_prevented"="#00bf7d", "NNS"="#00b0f6", "NNT"="#e76bf3")) + 
  xlab("Change relative to population-wide screening with SCORE2 alone (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.background=element_blank(), strip.text=element_text(size=6, face="bold"),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )

ggsave(g, width=7.2, height=2, file="analyses/public_health_modelling/screening_comparison.pdf")

# Create supp plot by sex
ggdt <- sex_estimates[metric %in% c("N_treated", "cases_treated", "events_prevented", "NNS", "NNT")]
ggdt[, metric := factor(metric, levels=c("N_treated", "cases_treated", "events_prevented", "NNS", "NNT"))]
ggdt[, model := factor(model, levels=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]

g1 <- ggplot(ggdt[model != "SCORE2" & sex == "Male"]) + 
  aes(x=delta, xmin=delta.L95, xmax=delta.U95, y=fct_rev(model), color=metric) +
  facet_grid(strategy ~ metric, scales="free") + 
  geom_vline(xintercept=0, linetype=2) + 
  geom_errorbarh(height=0) + 
  geom_point(shape=23, fill="white", size=1.2) +
  scale_color_manual(values=c("N_treated"="#f8766d", "cases_treated"="#a3a500", "events_prevented"="#00bf7d", "NNS"="#00b0f6", "NNT"="#e76bf3")) + 
  xlab("Change relative to population-wide screening with SCORE2 alone (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.background=element_blank(), strip.text=element_text(size=6, face="bold"),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )

g2 <- ggplot(ggdt[model != "SCORE2" & sex == "Female"]) + 
  aes(x=delta, xmin=delta.L95, xmax=delta.U95, y=fct_rev(model), color=metric) +
  facet_grid(strategy ~ metric, scales="free") + 
  geom_vline(xintercept=0, linetype=2) + 
  geom_errorbarh(height=0) + 
  geom_point(shape=23, fill="white", size=1.2) +
  scale_color_manual(values=c("N_treated"="#f8766d", "cases_treated"="#a3a500", "events_prevented"="#00bf7d", "NNS"="#00b0f6", "NNT"="#e76bf3")) + 
  xlab("Change relative to population-wide screening with SCORE2 alone (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.background=element_blank(), strip.text=element_text(size=6, face="bold"),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )

g <- plot_grid(g1, g2, nrow=2, align="hv") 
ggsave(g, width=7.2, height=4, file="analyses/public_health_modelling/screening_comparison_sex_specific.pdf")

