library(data.table)
library(ggplot2)

# Load NRI estimates
nri_estimates <- fread("analyses/public_health_modelling/net_reclassification/nri_estimates.txt")

# Load model information
model_info <- fread("analyses/test/model_fit_information.txt")

# NRI estimates compare how case/control are reclassified when changing from one model 
# (e.g. conventional risk factors alone) to another (e.g. conventional risk factors +
#  selected biomarkers).
#
# We've compute various comparisons, using risk thresholds of 0/5/10% (UK guidelines),
# 0/5/7.5% (American guidelines), and continuous NRI.
#
# First, we want to see whether, like the C-indices, the lasso models with lambda =
# lambda.min outperform those with lambda.1se
gg_dt <- nri_estimates[metric %in% c("NRI+", "NRI-")]
gg_dt[, compared_to := paste("Compared to:", ifelse(grepl("^Base \\+ PGS", model), "Conventional RF + PGS", "Conventional RF"))]
gg_dt[, biomarkers := gsub(".*to ", "", model)]
gg_dt[, biomarkers := gsub(" \\(.*\\)", "", biomarkers)]
gg_dt[model_info[lambda == "lambda.min"], on = .(biomarkers=name), long_name := ifelse(lambda == "Best model", i.long_name, NA_character_)]
gg_dt[model_info[lambda == "lambda.1se"], on = .(biomarkers=name), long_name := ifelse(lambda == "Best model with fewest features", i.long_name, x.long_name)]
gg_dt[biomarkers == "PGS", long_name := "Conventional RF + PGS"]
gg_dt[biomarkers == "PGS + Blood" & lambda == "Best model", 
      long_name := model_info[lambda == "lambda.min" & (PGS) & name == "Blood", long_name]]
gg_dt[biomarkers == "PGS + Blood" & lambda == "Best model with fewest features",
      long_name := model_info[lambda == "lambda.1se" & (PGS) & name == "Blood", long_name]]
gg_dt[biomarkers == "PGS + Nightingale" & lambda == "Best model", 
      long_name := model_info[lambda == "lambda.min" & (PGS) & name == "Nightingale", long_name]]
gg_dt[biomarkers == "PGS + Nightingale" & lambda == "Best model with fewest features",
      long_name := model_info[lambda == "lambda.1se" & (PGS) & name == "Nightingale", long_name]]
gg_dt[biomarkers == "PGS + Blood + Nightingale" & lambda == "Best model", 
      long_name := model_info[lambda == "lambda.min" & (PGS) & name == "Blood + Nightingale", long_name]]
gg_dt[biomarkers == "PGS + Blood + Nightingale" & lambda == "Best model with fewest features",
      long_name := model_info[lambda == "lambda.1se" & (PGS) & name == "Blood + Nightingale", long_name]]
gg_dt[biomarkers %like% "CRP", long_name := paste("Conventional RF +", biomarkers)]
gg_dt[biomarkers %like% "GlycA", long_name := paste("Conventional RF +", biomarkers)]
gg_dt[, risk_categories := gsub(" Categorical NRI", "", nri_type)]
gg_dt[, risk_categories := factor(risk_categories, levels=c("Continuous NRI", "ACC/AHA 2019", "NICE 2014"))]
gg_dt[, nri_group := ifelse(metric == "NRI+", "cases", "controls")]
gg_dt[, biomarkers := factor(biomarkers, levels=rev(c("CRP", "GlycA", "Nightingale", "Blood", "Blood + Nightingale", 
  "PGS", "PGS + CRP", "PGS + GlycA", "PGS + Nightingale", "PGS + Blood", "PGS + Blood + Nightingale")))]
gg_dt[, lambda := factor(lambda, levels=rev(c("No feature selection", "Best model with fewest features", "Best model")))]
gg_dt <- gg_dt[order(lambda)][order(biomarkers)]
gg_dt[, long_name := factor(long_name, levels=unique(long_name))]

# curate dataset for showing Conventional RF + PGS intercept line
gg_dt2 <- gg_dt[biomarkers == "PGS"]
gg_dt2 <- rbind(gg_dt2, gg_dt2)
gg_dt2[1:(.N/2), compared_to := "Compared to: Conventional RF + PGS"]

g <- ggplot(gg_dt) +
  aes(x = Estimate*100, xmin=Lower*100, xmax=Upper*100, y=biomarkers, color=lambda, shape=biomarkers) +
  geom_vline(xintercept=0, linetype=2) +
  geom_vline(data=gg_dt2, aes(xintercept=Estimate*100), color="red", linetype=2) +
  geom_errorbarh(height=0, position=position_dodge(width=0.6)) +
  geom_point(size=2, fill="white", position=position_dodge(width=0.6)) +
  facet_grid(compared_to ~ nri_group + risk_categories, scales="free") +
  scale_shape_manual(name="Biomarkers", values=c(
    "CRP"=23, "GlycA"=24, "Blood"=23, "Nightingale"=24, "Blood + Nightingale"=22, "PGS"=17,
    "PGS + CRP"=18, "PGS + GlycA"=17, "PGS + Blood"=18, "PGS + Nightingale"=17, "PGS + Blood + Nightingale"=15
  )) +
  scale_color_manual(name="Feature selection", values=c("No feature selection"="red", "Best model"="#00B050",
                     "Best model with fewest features"="#0070C0")) +
  ylab("") + xlab("Categorical NRI, % reclassified (95% CI)") +
  theme_bw() + theme(legend.position="bottom", legend.box="vertical")
ggsave(g, width=13, height=7, units="in", file="analyses/public_health_modelling/net_reclassification/nri_compare_lambda_min_1se.pdf")

# Now just show the lambda.min models
g <- ggplot(gg_dt[lambda != "Best model with fewest features"]) +
  aes(x = Estimate*100, xmin=Lower*100, xmax=Upper*100, y=long_name, color=lambda, shape=biomarkers) +
  geom_vline(xintercept=0, linetype=2) +
  geom_vline(data=gg_dt2, aes(xintercept=Estimate*100), color="#d95f02", linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(size=2, fill="white") +
  facet_grid(compared_to ~ nri_group + risk_categories, scales="free") +
  scale_shape_manual(name="Biomarkers", values=c(
    "CRP"=23, "GlycA"=24, "Blood"=23, "Nightingale"=24, "Blood + Nightingale"=22, "PGS"=17,
    "PGS + CRP"=18, "PGS + GlycA"=17, "PGS + Blood"=18, "PGS + Nightingale"=17, "PGS + Blood + Nightingale"=15
  )) +
  scale_color_manual(name="Feature selection", values=c("No feature selection"="#d95f02", "Best model"="#7570b3")) +
  ylab("") + xlab("Categorical NRI, % reclassified (95% CI)") +
  theme_bw() + theme(legend.position="bottom")
ggsave(g, width=13, height=7, units="in", file="analyses/public_health_modelling/net_reclassification/nri_compare.pdf")

# Write out table for supp
wide <- dcast(gg_dt[lambda != "Best model with fewest features"], long_name + biomarkers + compared_to ~ nri_type + metric, value.var=c("Estimate", "Lower", "Upper"))
wide <- wide[order(-biomarkers)][order(compared_to)]
wide <- wide[, .(long_name, biomarkers, compared_to, 
  `ContNRI+`=`Estimate_Continuous NRI_NRI+`, `ContNRI+_L95`=`Lower_Continuous NRI_NRI+`, `ContNRI+_U95`=`Upper_Continuous NRI_NRI+`,
  `ContNRI-`=`Estimate_Continuous NRI_NRI-`, `ContNRI-_L95`=`Lower_Continuous NRI_NRI-`, `ContNRI-_U95`=`Upper_Continuous NRI_NRI-`,
  `ACCAHA_NRI+`=`Estimate_ACC/AHA 2019 Categorical NRI_NRI+`, `ACCAHA_NRI+_L95`=`Lower_ACC/AHA 2019 Categorical NRI_NRI+`, `ACCAHA_NRI+_U95`=`Upper_ACC/AHA 2019 Categorical NRI_NRI+`,
  `ACCAHA_NRI-`=`Estimate_ACC/AHA 2019 Categorical NRI_NRI-`, `ACCAHA_NRI-_L95`=`Lower_ACC/AHA 2019 Categorical NRI_NRI-`, `ACCAHA_NRI-_U95`=`Upper_ACC/AHA 2019 Categorical NRI_NRI-`,
  `NICE_NRI+`=`Estimate_NICE 2014 Categorical NRI_NRI+`, `NICE_NRI+_L95`=`Lower_NICE 2014 Categorical NRI_NRI+`, `NICE_NRI+_U95`=`Upper_NICE 2014 Categorical NRI_NRI+`,
  `NICE_NRI-`=`Estimate_NICE 2014 Categorical NRI_NRI-`, `NICE_NRI-_L95`=`Lower_NICE 2014 Categorical NRI_NRI-`, `NICE_NRI-_U95`=`Upper_NICE 2014 Categorical NRI_NRI-`
)]
fwrite(wide, sep="\t", quote=FALSE, file="analyses/public_health_modelling/net_reclassification/nri_compare.txt")
