library(data.table)
library(foreach)
library(doMC)
library(survival)
library(boot)
library(ggplot2)
registerDoMC(6) # bootstraps take ~ 30 mins

system("mkdir -p analyses/test/")

# Load data
dat <- fread("data/cleaned/analysis_cohort.txt")

# Build data.table of models to test
model_info <- foreach(this_sex=c("Males", "Females", "Sex-stratified"), .combine=rbind) %:%
  foreach(this_model=c("SCORE2", "SCORE2_excl_UKB"), .combine=rbind) %do% {
    data.table(sex=this_sex,  model=this_model)
}

# Convert sex to integer for strata
strata_num <- dat[,.GRP, by=sex]
dat[strata_num, on = .(sex), sex_int := i.GRP]

# Filter to relevant columns for bootstrap tests
dat <- dat[,.(eid, sex, sex_int, incident_cvd, incident_cvd_followup, SCORE2, SCORE2_excl_UKB)]

# Create bootstrap function that: 
# (1) calculates the absolute C-index for each model
# (2) calculates the C-index standard error for each model
boot.fun <- function(dt) {
  model_cinds <- foreach(midx = model_info[,.I], .combine=rbind) %dopar% {
    this_minfo <- model_info[midx]
    if (this_minfo$sex == "Sex-stratified") {
      this_y <- Surv(dt[["incident_cvd_followup"]], dt[["incident_cvd"]])
      this_x <- dt[[this_minfo$model]]
      cf <- survival::concordancefit(y = this_y, x = this_x, strata = dt[["sex_int"]], reverse = TRUE, timefix = TRUE)
    } else {
      this_dt <- dt[sex == gsub("s$", "", this_minfo$sex)]
      this_y <- Surv(this_dt[["incident_cvd_followup"]], this_dt[["incident_cvd"]])
      this_x <- this_dt[[this_minfo$model]]
      cf <- survival::concordancefit(y = this_y, x = this_x, reverse = TRUE, timefix = TRUE)
    }
    cbind(this_minfo, C.index=cf$concordance, SE=sqrt(cf$var)) # jackknife SE computed same way as in coxph()
  }
  # Need to return as flat vector
  res <- melt(model_cinds, measure.vars=c("C.index", "SE"))
  res$value
}

surv_cols_idx <- match(c("incident_cvd_followup", "incident_cvd"), names(dat))
boot_res <- censboot(dat, boot.fun, 1000, index=surv_cols_idx)
saveRDS(boot_res, fil="analyses/test/SCORE2_overfitting_bootstraps.rds")

# Now we need to extract and collate the results
cinds <- foreach(this_metric = c("C.index", "SE"), .combine=rbind) %do% {
  cbind(model_info, metric = this_metric)
}
cinds[, estimate := boot_res$t0]
cinds[, L95 := apply(boot_res$t, 2, function(v) {  sort(v)[25] })]
cinds[, U95 := apply(boot_res$t, 2, function(v) {  sort(v)[975] })]

# Cast to wide
cinds <- dcast(cinds, sex + model ~ metric, value.var=c("estimate", "L95", "U95"))
cinds <- cinds[,.(sex, SCORE2_method = ifelse(model == "SCORE2", "Weights derived from all datasets", "Weights derived excluding UK Biobank"),
  C.index = estimate_C.index, C.L95 = L95_C.index, C.U95 = U95_C.index,
  SE = estimate_SE, SE.L95 = L95_SE, SE.U95 = U95_SE
)]
fwrite(cinds, sep="\t", quote=FALSE, file="analyses/test/cindex_by_SCORE2_method.txt")

cinds[, sex := factor(sex, levels=c("Sex-stratified", "Males", "Females"))]
cinds[, SCORE2_method := factor(SCORE2_method, levels=c("Weights derived excluding UK Biobank", "Weights derived from all datasets"))]
g <- ggplot(cinds) + 
  aes(x=C.index, xmin=C.L95, xmax=C.U95, y=SCORE2_method, color=SCORE2_method) +
  facet_wrap(~ sex, scales="free_x") +
  geom_errorbarh(height=0) +
  geom_point(shape=18) +
  scale_color_manual(values=c("Weights derived excluding UK Biobank"="#2166ac", "Weights derived from all datasets"="#b2182b")) +
  xlab("C-index (95% CI)") +
  guides(color=guide_legend(title="SCORE2 model", reverse=TRUE)) +
  theme_bw() +
  theme(
    axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    strip.background=element_blank(), strip.text=element_text(size=8, face="bold"), 
    legend.text=element_text(size=6), legend.title=element_text(size=8)
  )
ggsave(g, width=7.2, height=1.5, file="analyses/test/cindex_by_SCORE2_method.pdf")


