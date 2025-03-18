library(data.table)
library(ggplot2)
library(ggstance)
library(ggh4x)
library(patchwork)
library(cowplot)

# load in combined scores
comb_scores <- fread("analyses/CVD_weight_training/phase3_CVD_linear_predictors_and_risk.txt")

# Compute densities
comb_densities <- foreach(this_model = c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs")) %do% {
  l2 <- foreach(this_sex = c("Male", "Female")) %do% {
		this_dat <- comb_scores[model == this_model & sex == this_sex]
		this_dat[, density(linear_predictor, na.rm=TRUE)]
  }
  names(l2) <- c("Male", "Female")
  l2
}
names(comb_densities) <- c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs")

# Extract table of just the densities for plotting
dens_check <- foreach(this_model = c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %do% {
  foreach(this_sex = c("Male", "Female"), .combine=rbind) %do% {
		this_dens <- comb_densities[[this_model]][[this_sex]]
		data.table(model=this_model, sex=this_sex, cohort="UKB replication", x=this_dens$x, y=this_dens$y)
  }
}

# Add in distributions from phase 1 and 2
ukb_dens <- readRDS("analyses/replication/UKB_combined_score_densities.rds")

ukb_dens <- foreach(this_model = c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %do% {
	foreach(this_sex = c("Male", "Female"), .combine=rbind) %do% {
		this_dens <- ukb_dens[[this_model]][[this_sex]]
		data.table(model=this_model, sex=this_sex, cohort="UKB discovery", x=this_dens$x, y=this_dens$y)
	}
}
dens_check <- rbind(dens_check, ukb_dens)

# Add in SCORE2 densities
score2_dens <- foreach(this_model = c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %do% {
  foreach(this_sex = c("Male", "Female"), .combine=rbind) %do% {
    this_dat <- comb_scores[model == "SCORE2" & sex == this_sex]
    this_dens <- this_dat[,density(linear_predictor, na.rm=TRUE)]
    data.table(model=this_model, sex=this_sex, cohort="SCORE2 in UKB replication", x=this_dens$x, y=this_dens$y)
  }
}
dens_check <- rbind(dens_check, score2_dens)

disc_comb_scores <- fread("analyses/CVD_weight_training/CVD_linear_predictors_and_risk.txt")
disc_score2_dens <- foreach(this_model = c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %do% {
  foreach(this_sex = c("Male", "Female"), .combine=rbind) %do% {
    this_dat <- disc_comb_scores[model == "SCORE2" & sex == this_sex]
    this_dens <- this_dat[,density(linear_predictor, na.rm=TRUE)]
    data.table(model=this_model, sex=this_sex, cohort="SCORE2 in UKB discovery", x=this_dens$x, y=this_dens$y)
  }
}
dens_check <- rbind(dens_check, disc_score2_dens)

# Create density plot
dens_check[, sex := factor(sex, levels=c("Male", "Female"))]
g <- ggplot(dens_check) +
	aes(x=x, y=y, color=cohort) +
	facet_nested(sex ~ model) +
	geom_line(linewidth=0.5) +
	scale_color_manual("", values=c(
		"SCORE2 in UKB discovery"="black",
		"SCORE2 in UKB replication"="#d73027",
		"UKB discovery"="#542788", 
		"UKB replication"="#1b7837"), 
	guide=guide_legend(ncol=1, reverse=TRUE)) +
	scale_x_continuous("Score value (linear predictor)") +
	scale_y_continuous("Density") +
	theme_bw() +
	theme(
		axis.title=element_text(size=8), axis.text=element_text(size=6),
		strip.text=element_text(size=8, face="bold"),
		legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=8),
		plot.title=element_text(size=10, hjust=0, vjust=0.1, face="bold"),
	)

ggsave(g, width=13, height=5, file="analyses/CVD_weight_training/phase3_combined_score_distribution_checks.pdf")



