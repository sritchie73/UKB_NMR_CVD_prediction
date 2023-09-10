library(data.table)
library(foreach)
library(boot)
library(statmod)
library(ggplot2)
library(ggh4x)
library(ggstance)

# Make output directory
system("mkdir -p analyses/test")

# Load required data
dat <- fread("analyses/CVD_weight_training/CVD_linear_predictors_and_risk.txt")

# Build data.table of models to test
model_info <- foreach(this_sex=c("Males", "Females", "Sex-stratified"), .combine=rbind) %:%
  foreach(this_model=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %:%
    foreach(this_score_type=c("non-derived", "clinical"), .combine=rbind) %do% {
      if (this_model == "SCORE2") {
        if (this_score_type == "clinical") {
          return(NULL)
        } else {
          this_score_type <- NA
        }
      }
      data.table(sex=this_sex,  model=this_model, score_type=this_score_type)
}
model_info[, model_var := paste0("model_", 1:.N), by=sex]

# Convert sex to integer for strata
strata_num <- dat[,.GRP, by=sex]
dat[strata_num, on = .(sex), sex_int := i.GRP]

# Simplify dat to include only one instance of SCORE2
dat <- dat[!(model == "SCORE2" & score_type == "clinical")]
dat[model == "SCORE2", score_type := NA]

# Cast to wide for bootstraps (i.e. we sample with replacement the same sample set for each comparison)
dat[model_info, on = .(model, score_type), model_var := i.model_var]
dat <- dcast(dat, eid + sex + sex_int + incident_cvd + incident_cvd_followup ~ model_var, value.var="linear_predictor")

# Create bootstrap function that: 
# (1) calculates the absolute C-index for each model
# (2) calculates the C-index standard error for each model
# (3) calculates the difference in C-index from SCORE2 for each model
# (4) calculates the % change in C-index from SCORE2 for each model
# For each entry in the model_info table
boot.fun <- function(dt) {
  model_cinds <- foreach(midx = model_info[,.I], .combine=rbind) %dopar% {
    this_minfo <- model_info[midx]
    if (this_minfo$sex == "Sex-stratified") {
      this_y <- Surv(dt[["incident_cvd_followup"]], dt[["incident_cvd"]])
      this_x <- dt[[this_minfo$model_var]]
      cf <- survival::concordancefit(y = this_y, x = this_x, strata = dt[["sex_int"]], reverse = TRUE, timefix = TRUE)
    } else {
      this_dt <- dt[sex == gsub("s$", "", this_minfo$sex)]
      this_y <- Surv(this_dt[["incident_cvd_followup"]], this_dt[["incident_cvd"]])
      this_x <- this_dt[[this_minfo$model_var]]
      cf <- survival::concordancefit(y = this_y, x = this_x, reverse = TRUE, timefix = TRUE)
    }       
    
    cbind(this_minfo, C.index=cf$concordance, SE=sqrt(cf$var)) # jackknife SE computed same way as in coxph()
  }

  # Compute delta C-index from SCORE2
  model_cinds[model_cinds[model == "SCORE2"], on = .(sex), delta.C := C.index - i.C.index]
  
  # Compute % change in C-index from SCORE2 
  model_cinds[model_cinds[model == "SCORE2"], on = .(sex), pct_change := delta.C/(i.C.index - 0.5)*100]
  
  # Need to return as flat vector
  res <- melt(model_cinds, measure.vars=c("C.index", "SE", "delta.C", "pct_change")) 
  res$value
}

# Run bootstrap analysis (should take ~30 minutes across 7 cores)
surv_cols_idx <- match(c("incident_cvd_followup", "incident_cvd"), names(dat))
boot_res <- censboot(dat, boot.fun, 1000, index=surv_cols_idx)
saveRDS(boot_res, fil="analyses/test/cindices_bootstraps.rds")

# Now we need to extract and collate the results
cinds <- foreach(this_metric = c("C.index", "SE", "delta.C", "pct_change"), .combine=rbind) %do% {
  cbind(model_info, metric = this_metric)
}
cinds[, estimate := boot_res$t0]
cinds[, L95 := apply(boot_res$t, 2, function(v) {  sort(v)[25] })] 
cinds[, U95 := apply(boot_res$t, 2, function(v) {  sort(v)[975] })] 
cinds[, pval := apply(boot_res$t, 2, function(v) {  (1 + sum(v <= 0))/(1000 + 1) })] # one-sided

# Cast to wide
cinds <- dcast(cinds, sex + model + score_type ~ metric, value.var=c("estimate", "L95", "U95", "pval"))
cinds <- cinds[,.(sex, model, score_type, 
  C.index = estimate_C.index, C.L95 = L95_C.index, C.U95 = U95_C.index,
  SE = estimate_SE, SE.L95 = L95_SE, SE.U95 = U95_SE, 
  deltaC = estimate_delta.C, deltaC.L95 = L95_delta.C, deltaC.U95 = U95_delta.C, deltaC.pval = pval_delta.C,
  pct_change = estimate_pct_change, pct.L95 = L95_pct_change, pct.U95 = U95_pct_change, pct.pval = pval_pct_change
)]

# Write out
fwrite(cinds, sep="\t", quote=FALSE, file="analyses/test/cindices.txt")

# Extract data for plotting
ggdt <- rbind(
  cinds[,.(sex, model, score_type, metric="C.index", estimate=C.index, L95=C.L95, U95=C.U95)],
  cinds[model != "SCORE2",.(sex, model, score_type, metric="deltaC", estimate=deltaC, L95=deltaC.L95, U95=deltaC.U95)]
)

ggdt <- ggdt[!(model == "SCORE2 + PRSs" & score_type == "clinical")]
ggdt[model == "SCORE2 + PRSs", score_type := NA]

# Encode factors for plotting
ggdt[, model := factor(model, levels=rev(c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs")))]
ggdt[, sex := factor(sex, levels=c("Males", "Females", "Sex-stratified"))]

# Extract reference for SCORE2
ref <- rbind(
  ggdt[model == "SCORE2", .(sex, metric, estimate)],
  ggdt[model == "SCORE2", .(sex, metric="deltaC", estimate=0)]
)

# Plot sex-stratified analysis for full biomarker scores
g <- ggplot(ggdt[(is.na(score_type) | score_type == "non-derived") & sex == "Sex-stratified"]) +
  aes(x=estimate, xmin=L95, xmax=U95, y=model) +
  facet_wrap(~ metric, nrow=1, scales="free_x") +
  geom_vline(data=ref[sex == "Sex-stratified"], aes(xintercept=estimate), linetype=2) + 
  geom_errorbarh(height=0) +
  geom_point(shape=23, size=2, fill="white") + 
  ylab("") + xlab("Estimate (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )
ggsave(g, width=7.2, height=2, file="analyses/test/cindex_sex_stratified.pdf")

# Plot analysis in males and females separately
g <- ggplot(ggdt[(is.na(score_type) | score_type == "non-derived") & sex != "Sex-stratified"]) +
  aes(x=estimate, xmin=L95, xmax=U95, y=model, color=sex) +
  facet_grid2(sex ~ metric, independent="x", scales="free_x") +
  geom_vline(data=ref[sex != "Sex-stratified"], aes(xintercept=estimate), linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, size=2, fill="white") +
  scale_color_manual("Sex", values=c("Males"="#e41a1c", "Females"="#377eb8")) +
  ylab("") + xlab("Estimate (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )
ggsave(g, width=7.2, height=3.5, file="analyses/test/cindex_sex_specific.pdf")

# Generate comparison plot showing clinical scores have worse performance:
g <- ggplot(ggdt) +
  aes(x=estimate, xmin=L95, xmax=U95, y=model, color=score_type) +
  facet_grid2(sex ~ metric, independent="x", scales="free_x") +
  geom_vline(data=ref, aes(xintercept=estimate), linetype=2) +
  geom_errorbarh(height=0, position=position_dodgev(height=0.6)) +
  geom_point(shape=23, size=2, fill="white", position=position_dodgev(height=0.6)) +
  scale_color_manual(values=c("non-derived"="black", "clinical"="red"), labels=c("non-derived"="Full NMR scores", "clinical"="Clinical biomarkers")) +
  ylab("") + xlab("Estimate (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=8)
  )
ggsave(g, width=7.2, height=5, file="analyses/test/cindex_with_clinical_scores.pdf")


# Create formatted table for manuscript
dt <- cinds[is.na(score_type) | score_type == "non-derived"]
dt[, score_type := NULL]
dt[, sex := factor(sex, levels=c("Sex-stratified", "Males", "Females"))]
dt <- dt[order(sex)]
dt[model == "SCORE2", c("deltaC", "deltaC.L95", "deltaC.U95", "deltaC.pval", "pct_change", "pct.L95", "pct.U95") := NA]
dt[, pct.pval := NULL]
dt[, pct_change := pct_change / 100]
dt[, pct.L95 := pct.L95 / 100]
dt[, pct.U95 := pct.U95 / 100]
fwrite(dt, sep="\t", quote=FALSE, file="analyses/test/cindex_for_supp.txt")

