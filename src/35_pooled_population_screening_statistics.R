library(data.table)
library(foreach)

statin_modifier <- 1/10  # Assume statins prevent 1 in 10 cases, i.e. a 10% risk reduction over 10 years

# Load bootstrapped demographic standardised risk stratification
boot_res <- fread("analyses/public_health_modelling/pooled_standardised_risk_stratification_bootstraps.txt")

# Aggregate into population modelling statistics
res <- foreach(this_strategy = c("blanket", "targeted"), .combine=rbind) %:%
  foreach(this_endpoint = c("cvd", "cvd_narrow"), .combine=rbind) %:%
    foreach(this_score = c("SCORE2", "SCORE2_excl_UKB", "QRISK3"), .combine=rbind) %:%
      foreach(this_model_sex = c("Sex-stratified", "Males", "Females"), .combine=rbind) %:%
        foreach(this_model = unique(boot_res$model_name), .combine=rbind) %do% {
					# Extract relevant bootstrap statistics
					this_res <- boot_res[strategy == this_strategy & endpoint == this_endpoint & score == this_score & model_name == this_model]

					# Extract relevant sex
					if (this_model_sex == "Males") {
						this_res <- this_res[sex == "Male"]
					} else if (this_model_sex == "Females") {
						this_res <- this_res[sex == "Female"]
					}

					# Compute relevant statistics
					total_n <- this_res[bootstrap == 0 & comparitor == "ref", sum(N)]
					total_n_boot <- this_res[bootstrap != 0 & comparitor == "ref", .(N=sum(N)), by=bootstrap][,N] # by definition, all should be identical to total_n
					total_cases <- this_res[bootstrap == 0 & comparitor == "ref" & status == "case", sum(N)]
					total_cases_boot <- this_res[bootstrap != 0 & comparitor == "ref" & status == "case", .(N=sum(N)), by=bootstrap][,N] # by definition, all should be identical to total_cases
					ref_high_risk <- this_res[bootstrap == 0 & comparitor == "ref" & risk_group == "high", sum(N)]
					ref_high_risk_boot <- this_res[bootstrap != 0 & comparitor == "ref" & risk_group == "high", .(N=sum(N)), by=bootstrap][,N]
					ref_cvd_high_risk <- this_res[bootstrap == 0 & comparitor == "ref" & risk_group == "high" & status == "case", sum(N)]
					ref_cvd_high_risk_boot <- this_res[bootstrap != 0 & comparitor == "ref" & risk_group == "high" & status == "case", .(N=sum(N)), by=bootstrap][,N]
					ref_cvd_prevented <- ref_cvd_high_risk * statin_modifier
					ref_cvd_prevented_boot <- ref_cvd_high_risk_boot * statin_modifier
					ref_NNS <- total_n / ref_cvd_prevented
					ref_NNS_boot <- total_n_boot / ref_cvd_prevented_boot
					ref_NNT <- total_cases / ref_cvd_prevented
					ref_NNT_boot <- total_cases_boot / ref_cvd_prevented_boot

					if (this_model != "Risk score") {
						alt_high_risk <- this_res[bootstrap == 0 & comparitor == "alt" & risk_group == "high", sum(N)]
						alt_high_risk_boot <- this_res[bootstrap != 0 & comparitor == "alt" & risk_group == "high", .(N=sum(N)), by=bootstrap][,N]
						alt_cvd_high_risk <- this_res[bootstrap == 0 & comparitor == "alt" & risk_group == "high" & status == "case", sum(N)]
						alt_cvd_high_risk_boot <- this_res[bootstrap != 0 & comparitor == "alt" & risk_group == "high" & status == "case", .(N=sum(N)), by=bootstrap][,N]
						alt_cvd_prevented <- alt_cvd_high_risk * statin_modifier
						alt_cvd_prevented_boot <- alt_cvd_high_risk_boot * statin_modifier
						alt_NNS <- total_n / alt_cvd_prevented
						alt_NNS_boot <- total_n_boot / alt_cvd_prevented_boot
						alt_NNT <- total_cases / alt_cvd_prevented
						alt_NNT_boot <- total_cases_boot / alt_cvd_prevented_boot

						delta_high_risk <- alt_high_risk - ref_high_risk
						delta_high_risk_boot <- alt_high_risk_boot - ref_high_risk_boot
						delta_cvd_high_risk <- alt_cvd_high_risk - ref_cvd_high_risk
						delta_cvd_high_risk_boot <- alt_cvd_high_risk_boot - ref_cvd_high_risk_boot
						delta_cvd_prevented <- alt_cvd_prevented - ref_cvd_prevented
						delta_cvd_prevented_boot <- alt_cvd_prevented_boot - ref_cvd_prevented_boot
						delta_NNS <- alt_NNS - ref_NNS
						delta_NNS_boot <- alt_NNS_boot - ref_NNS_boot
						delta_NNT <- alt_NNT - ref_NNT
						delta_NNT_boot <- alt_NNT_boot - ref_NNT_boot
					}

					# Build output table of statistics
					this_stats <- rbind(use.names=FALSE,
						data.table(metric="total_n", estimate=total_n, se=sd(total_n_boot)),
						data.table("total_cases", total_cases, sd(total_cases_boot)),
						data.table("ref_high_risk", ref_high_risk, sd(ref_high_risk_boot)),
						data.table("ref_cvd_high_risk", ref_cvd_high_risk, sd(ref_cvd_high_risk_boot)),
						data.table("ref_cvd_prevented", ref_cvd_prevented, sd(ref_cvd_prevented_boot)),
						data.table("ref_NNS", ref_NNS, sd(ref_NNS_boot)),
						data.table("ref_NNT", ref_NNT, sd(ref_NNT_boot))
					)

					if (this_model != "Risk score") {
						this_stats <- rbind(this_stats, use.names=FALSE,
							data.table("alt_high_risk", alt_high_risk, sd(alt_high_risk_boot)),
							data.table("alt_cvd_high_risk", alt_cvd_high_risk, sd(alt_cvd_high_risk_boot)),
							data.table("alt_cvd_prevented", alt_cvd_prevented, sd(alt_cvd_prevented_boot)),
							data.table("alt_NNS", alt_NNS, sd(alt_NNS_boot)),
							data.table("alt_NNT", alt_NNT, sd(alt_NNT_boot)),
							data.table("delta_high_risk", delta_high_risk, sd(delta_high_risk_boot)),
							data.table("delta_cvd_high_risk", delta_cvd_high_risk, sd(delta_cvd_high_risk_boot)),
							data.table("delta_cvd_prevented", delta_cvd_prevented, sd(delta_cvd_prevented_boot)),
							data.table("delta_NNS", delta_NNS, sd(delta_NNS_boot)),
							data.table("delta_NNT", delta_NNT, sd(delta_NNT_boot))
						)
					}

					# Compute 95% confidence intervals
					this_stats[, L95 := estimate - qnorm(1-(0.05/2))*se]
					this_stats[, U95 := estimate + qnorm(1-(0.05/2))*se]

					# Compute P-values for delta statistics
					this_stats[metric %like% "delta", pval := pmin(1, pnorm(abs(estimate/se), lower.tail=FALSE)*2)]
				 
					# Add model information and return
					cbind(strategy=this_strategy, endpoint=this_endpoint, score=this_score, model_sex=this_model_sex, model=this_model, this_stats)
}
fwrite(res, sep="\t", quote=FALSE, file="analyses/public_health_modelling/pooled_screening.txt")

