library(data.table)

# Load bootstrap results
boot_res <- fread("analyses/public_health_modelling/pooled_standardised_risk_stratification_bootstraps.txt")

# Filter to representative example for stratifying population by SCORE2
boot_res <- boot_res[strategy == "blanket" & endpoint == "cvd" & score == "SCORE2_excl_UKB" & model_name == "Risk score + NMR scores" & comparitor == "ref"]

# Compute relevant statistics
ref_high_risk <- boot_res[bootstrap == 0 & risk_group == "high", sum(N)]
ref_high_risk_boot <- boot_res[bootstrap != 0 & risk_group == "high", .(N=sum(N)), by=bootstrap][,N]
ref_cvd_high_risk <- boot_res[bootstrap == 0 & risk_group == "high" & status == "case", sum(N)]
ref_cvd_high_risk_boot <- boot_res[bootstrap != 0 & risk_group == "high" & status == "case", .(N=sum(N)), by=bootstrap][,N]
ref_cvd_prevented <- ref_cvd_high_risk/5
ref_cvd_prevented_boot <- ref_cvd_high_risk_boot/5
ref_NNS <- total_n / ref_cvd_prevented
ref_NNS_boot <- total_n_boot / ref_cvd_prevented_boot
ref_NNT <- ref_high_risk / ref_cvd_prevented
ref_NNT_boot <- ref_high_risk_boot / ref_cvd_prevented_boot

alt_high_risk <- boot_res[bootstrap == 0 & risk_group != "low", sum(N)]
alt_high_risk_boot <- boot_res[bootstrap != 0 & risk_group != "low", .(N=sum(N)), by=bootstrap][,N]
alt_cvd_high_risk <- boot_res[bootstrap == 0 & risk_group != "low" & status == "case", sum(N)]
alt_cvd_high_risk_boot <- boot_res[bootstrap != 0 & risk_group != "low" & status == "case", .(N=sum(N)), by=bootstrap][,N]
alt_cvd_prevented <- alt_cvd_high_risk/5
alt_cvd_prevented_boot <- alt_cvd_high_risk_boot/5
alt_NNS <- total_n / alt_cvd_prevented
alt_NNS_boot <- total_n_boot / alt_cvd_prevented_boot
alt_NNT <- alt_high_risk / alt_cvd_prevented
alt_NNT_boot <- alt_high_risk_boot / alt_cvd_prevented_boot

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

# Build output table of statistics
stats <- rbind(use.names=FALSE,
	data.table(metric="total_n", estimate=total_n, se=sd(total_n_boot)),
	data.table("total_cases", total_cases, sd(total_cases_boot)),
	data.table("ref_high_risk", ref_high_risk, sd(ref_high_risk_boot)),
	data.table("ref_cvd_high_risk", ref_cvd_high_risk, sd(ref_cvd_high_risk_boot)),
	data.table("ref_cvd_prevented", ref_cvd_prevented, sd(ref_cvd_prevented_boot)),
	data.table("ref_NNS", ref_NNS, sd(ref_NNS_boot)),
	data.table("ref_NNT", ref_NNT, sd(ref_NNT_boot)),
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

# Compute 95% confidence intervals
stats[, L95 := estimate - qnorm(1-(0.05/2))*se]
stats[, U95 := estimate + qnorm(1-(0.05/2))*se]

# Compute P-values for delta statistics
stats[metric %like% "delta", pval := pmin(1, pnorm(abs(estimate/se), lower.tail=FALSE)*2)]

# Write out 
fwrite(stats, sep="\t", quote=FALSE, file="analyses/public_health_modelling/treat_medium_and_high_risk.txt")





