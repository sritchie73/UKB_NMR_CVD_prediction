library(data.table)
library(foreach)

statin_modifier <- 1/5  # Assume statins prevent 1 in 5 cases, i.e. a 20% risk reduction over 10 years

# Load bootstrapped demographic standardised risk stratification
boot_res <- fread("analyses/public_health_modelling/discovery_standardised_targeted_risk_stratification_bootstraps.txt")

# Aggregate into population modelling statistics
res <- foreach(this_strategy = c("targeted"), .combine=rbind) %:%
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

          # Compute sex-specific scaling factors, so that estimates are per 100,000 males and per 100,000 females respectively,
          # Instead of per 100,000 of the male+female population
          total_n <- this_res[bootstrap == 0 & comparitor == "ref", sum(N)]
          scaling_factor <- 1e5/total_n  # will be 1 by definition for Sex-stratified analyses

          # Compute totals for sanity checking
          total_n <- total_n * scaling_factor
					total_n_boot <- this_res[bootstrap != 0 & comparitor == "ref", .(N=sum(N)), by=bootstrap][,N*scaling_factor] # by definition, all should be identical to total_n (i.e 1e5)
					total_cases <- this_res[bootstrap == 0 & comparitor == "ref" & status == "case", sum(N)*scaling_factor]
					total_cases_boot <- this_res[bootstrap != 0 & comparitor == "ref" & status == "case", .(N=sum(N)), by=bootstrap][,N*scaling_factor] # by definition, all should be identical to total_cases

					# Compute relevant statistics
					ref_high_risk <- this_res[bootstrap == 0 & comparitor == "ref" & risk_group == "high", sum(N)*scaling_factor]
					ref_high_risk_boot <- this_res[bootstrap != 0 & comparitor == "ref" & risk_group == "high", .(N=sum(N)), by=bootstrap][,N*scaling_factor]
					ref_cvd_high_risk <- this_res[bootstrap == 0 & comparitor == "ref" & risk_group == "high" & status == "case", sum(N)*scaling_factor]
					ref_cvd_high_risk_boot <- this_res[bootstrap != 0 & comparitor == "ref" & risk_group == "high" & status == "case", .(N=sum(N)), by=bootstrap][,N*scaling_factor]
					ref_cvd_prevented <- ref_cvd_high_risk * statin_modifier
					ref_cvd_prevented_boot <- ref_cvd_high_risk_boot * statin_modifier
					ref_NNS <- total_n / ref_cvd_prevented
					ref_NNS_boot <- total_n_boot / ref_cvd_prevented_boot
					ref_NNT <- ref_high_risk / ref_cvd_prevented
					ref_NNT_boot <- ref_high_risk_boot / ref_cvd_prevented_boot

					if (this_model != "Risk score") {
						alt_high_risk <- this_res[bootstrap == 0 & comparitor == "alt" & risk_group == "high", sum(N)*scaling_factor]
						alt_high_risk_boot <- this_res[bootstrap != 0 & comparitor == "alt" & risk_group == "high", .(N=sum(N)), by=bootstrap][,N*scaling_factor]
						alt_cvd_high_risk <- this_res[bootstrap == 0 & comparitor == "alt" & risk_group == "high" & status == "case", sum(N)*scaling_factor]
						alt_cvd_high_risk_boot <- this_res[bootstrap != 0 & comparitor == "alt" & risk_group == "high" & status == "case", .(N=sum(N)), by=bootstrap][,N*scaling_factor]
						alt_cvd_prevented <- alt_cvd_high_risk * statin_modifier
						alt_cvd_prevented_boot <- alt_cvd_high_risk_boot * statin_modifier
						alt_NNS <- total_n / alt_cvd_prevented
						alt_NNS_boot <- total_n_boot / alt_cvd_prevented_boot
						alt_NNT <- alt_high_risk / alt_cvd_prevented
						alt_NNT_boot <- alt_high_risk_boot / alt_cvd_prevented_boot

						null_high_risk <- this_res[bootstrap == 0 & comparitor == "null" & risk_group == "high", sum(N)*scaling_factor]
						null_high_risk_boot <- this_res[bootstrap != 0 & comparitor == "null" & risk_group == "high", .(N=sum(N)), by=bootstrap][,N*scaling_factor]
						null_cvd_high_risk <- this_res[bootstrap == 0 & comparitor == "null" & risk_group == "high" & status == "case", sum(N)*scaling_factor]
						null_cvd_high_risk_boot <- this_res[bootstrap != 0 & comparitor == "null" & risk_group == "high" & status == "case", .(N=sum(N)), by=bootstrap][,N*scaling_factor]
						null_cvd_prevented <- null_cvd_high_risk * statin_modifier
						null_cvd_prevented_boot <- null_cvd_high_risk_boot * statin_modifier
						null_NNS <- total_n / null_cvd_prevented
						null_NNS_boot <- total_n_boot / null_cvd_prevented_boot
						null_NNT <- null_high_risk / null_cvd_prevented
						null_NNT_boot <- null_high_risk_boot / null_cvd_prevented_boot

						alt_delta_ref_high_risk <- alt_high_risk - ref_high_risk
						alt_delta_ref_high_risk_boot <- alt_high_risk_boot - ref_high_risk_boot
						alt_delta_ref_cvd_high_risk <- alt_cvd_high_risk - ref_cvd_high_risk
						alt_delta_ref_cvd_high_risk_boot <- alt_cvd_high_risk_boot - ref_cvd_high_risk_boot
						alt_delta_ref_cvd_prevented <- alt_cvd_prevented - ref_cvd_prevented
						alt_delta_ref_cvd_prevented_boot <- alt_cvd_prevented_boot - ref_cvd_prevented_boot
						alt_delta_ref_NNS <- alt_NNS - ref_NNS
						alt_delta_ref_NNS_boot <- alt_NNS_boot - ref_NNS_boot
						alt_delta_ref_NNT <- alt_NNT - ref_NNT
						alt_delta_ref_NNT_boot <- alt_NNT_boot - ref_NNT_boot

						alt_delta_null_high_risk <- alt_high_risk - null_high_risk
						alt_delta_null_high_risk_boot <- alt_high_risk_boot - null_high_risk_boot
						alt_delta_null_cvd_high_risk <- alt_cvd_high_risk - null_cvd_high_risk
						alt_delta_null_cvd_high_risk_boot <- alt_cvd_high_risk_boot - null_cvd_high_risk_boot
						alt_delta_null_cvd_prevented <- alt_cvd_prevented - null_cvd_prevented
						alt_delta_null_cvd_prevented_boot <- alt_cvd_prevented_boot - null_cvd_prevented_boot
						alt_delta_null_NNS <- alt_NNS - null_NNS
						alt_delta_null_NNS_boot <- alt_NNS_boot - null_NNS_boot
						alt_delta_null_NNT <- alt_NNT - null_NNT
						alt_delta_null_NNT_boot <- alt_NNT_boot - null_NNT_boot

						null_delta_ref_high_risk <- null_high_risk - ref_high_risk
						null_delta_ref_high_risk_boot <- null_high_risk_boot - ref_high_risk_boot
						null_delta_ref_cvd_high_risk <- null_cvd_high_risk - ref_cvd_high_risk
						null_delta_ref_cvd_high_risk_boot <- null_cvd_high_risk_boot - ref_cvd_high_risk_boot
						null_delta_ref_cvd_prevented <- null_cvd_prevented - ref_cvd_prevented
						null_delta_ref_cvd_prevented_boot <- null_cvd_prevented_boot - ref_cvd_prevented_boot
						null_delta_ref_NNS <- null_NNS - ref_NNS
						null_delta_ref_NNS_boot <- null_NNS_boot - ref_NNS_boot
						null_delta_ref_NNT <- null_NNT - ref_NNT
						null_delta_ref_NNT_boot <- null_NNT_boot - ref_NNT_boot
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
							data.table("null_high_risk", null_high_risk, sd(null_high_risk_boot)),
							data.table("null_cvd_high_risk", null_cvd_high_risk, sd(null_cvd_high_risk_boot)),
							data.table("null_cvd_prevented", null_cvd_prevented, sd(null_cvd_prevented_boot)),
							data.table("null_NNS", null_NNS, sd(null_NNS_boot)),
							data.table("null_NNT", null_NNT, sd(null_NNT_boot)),
							data.table("alt_delta_ref_high_risk", alt_delta_ref_high_risk, sd(alt_delta_ref_high_risk_boot)),
							data.table("alt_delta_ref_cvd_high_risk", alt_delta_ref_cvd_high_risk, sd(alt_delta_ref_cvd_high_risk_boot)),
							data.table("alt_delta_ref_cvd_prevented", alt_delta_ref_cvd_prevented, sd(alt_delta_ref_cvd_prevented_boot)),
							data.table("alt_delta_ref_NNS", alt_delta_ref_NNS, sd(alt_delta_ref_NNS_boot)),
							data.table("alt_delta_ref_NNT", alt_delta_ref_NNT, sd(alt_delta_ref_NNT_boot)),
							data.table("alt_delta_null_high_risk", alt_delta_null_high_risk, sd(alt_delta_null_high_risk_boot)),
							data.table("alt_delta_null_cvd_high_risk", alt_delta_null_cvd_high_risk, sd(alt_delta_null_cvd_high_risk_boot)),
							data.table("alt_delta_null_cvd_prevented", alt_delta_null_cvd_prevented, sd(alt_delta_null_cvd_prevented_boot)),
							data.table("alt_delta_null_NNS", alt_delta_null_NNS, sd(alt_delta_null_NNS_boot)),
							data.table("alt_delta_null_NNT", alt_delta_null_NNT, sd(alt_delta_null_NNT_boot)),
							data.table("null_delta_ref_high_risk", null_delta_ref_high_risk, sd(null_delta_ref_high_risk_boot)),
							data.table("null_delta_ref_cvd_high_risk", null_delta_ref_cvd_high_risk, sd(null_delta_ref_cvd_high_risk_boot)),
							data.table("null_delta_ref_cvd_prevented", null_delta_ref_cvd_prevented, sd(null_delta_ref_cvd_prevented_boot)),
							data.table("null_delta_ref_NNS", null_delta_ref_NNS, sd(null_delta_ref_NNS_boot)),
							data.table("null_delta_ref_NNT", null_delta_ref_NNT, sd(null_delta_ref_NNT_boot))
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
fwrite(res, sep="\t", quote=FALSE, file="analyses/public_health_modelling/discovery_targeted_screening.txt")

