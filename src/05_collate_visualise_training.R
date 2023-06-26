library(data.table)
library(foreach)
library(ggplot2)
library(ggh4x)

# Load fits for all parameters
fitdt <- foreach(this_test_fold = 1:5, .combine=rbind) %do% {
  foreach(this_endpoint = c("CHD", "Stroke"), .combine=rbind) %do% {
    foreach(this_sex = c("Female", "Male"), .combine=rbind) %do% {
      model_dir <- sprintf("analyses/nmr_score_training/test_fold_%s/%s/%s", this_test_fold, this_endpoint, this_sex)
      fread(sprintf("%s/cv_coxnet_list_all_fits.txt", model_dir))
    }
  }
}

# Load best fit per alpha
bestfit <- foreach(this_test_fold = 1:5, .combine=rbind) %do% {
  foreach(this_endpoint = c("CHD", "Stroke"), .combine=rbind) %do% {
    foreach(this_sex = c("Female", "Male"), .combine=rbind) %do% {
      model_dir <- sprintf("analyses/nmr_score_training/test_fold_%s/%s/%s", this_test_fold, this_endpoint, this_sex)
      fread(sprintf("%s/cv_coxnet_best_fits.txt", model_dir))
    }
  }
}

# Load optimal fit per trait, sex, and test-fold
bestbestfit <- foreach(this_test_fold = 1:5, .combine=rbind) %do% {
  foreach(this_endpoint = c("CHD", "Stroke"), .combine=rbind) %do% {
    foreach(this_sex = c("Female", "Male"), .combine=rbind) %do% {
      model_dir <- sprintf("analyses/nmr_score_training/test_fold_%s/%s/%s", this_test_fold, this_endpoint, this_sex)
      fread(sprintf("%s/cv_coxnet_best_best_fits.txt", model_dir))
    }
  }
}
bestbestfit <- bestbestfit[,-(7:11)] # some columns duplicated in training code

# Load in non-zero coefficients
bestcoef <- foreach(this_test_fold = 1:5, .combine=rbind) %do% {
  foreach(this_endpoint = c("CHD", "Stroke"), .combine=rbind) %do% {
    foreach(this_sex = c("Female", "Male"), .combine=rbind) %do% {
      model_dir <- sprintf("analyses/nmr_score_training/test_fold_%s/%s/%s", this_test_fold, this_endpoint, this_sex)
      fread(sprintf("%s/best_best_fits_coefficients.txt", model_dir))
    }
  }
}

# Combine coefficients across test folds
avgbestcoef <- bestcoef[,.(beta=sum(beta)/5, sd=sd(beta), min=min(beta), max=max(beta), non_zero=.N),by=.(endpoint, sex, lambda.fit, coef)]
avgbestcoef <- avgbestcoef[order(-abs(beta))][order(lambda.fit)][order(sex)][order(endpoint)]

# Write out collated stats
fwrite(fitdt, sep="\t", quote=FALSE, file="analyses/nmr_score_training/all_elasticnet_fits.txt")
fwrite(bestfit, sep="\t", quote=FALSE, file="analyses/nmr_score_training/best_fits_per_alpha.txt")
fwrite(bestbestfit, sep="\t", quote=FALSE, file="analyses/nmr_score_training/best_fits.txt")
fwrite(bestcoef, sep="\t", quote=FALSE, file="analyses/nmr_score_training/best_fits_coef.txt")
fwrite(avgbestcoef, sep="\t", quote=FALSE, file="analyses/nmr_score_training/avg_best_coefs.txt")

# For each training independent 5-fold training, plot the model fits
fitdt[, sex := factor(sex, levels=c("Male", "Female"))]
for (this_test_fold in 1:5) {
  this_fitdt <- fitdt[prediction_cv_testfold == this_test_fold]
  this_bestfit <- bestfit[prediction_cv_testfold == this_test_fold]
  this_bestbestfit <- bestbestfit[prediction_cv_testfold == this_test_fold & model == "lambda.min"]

  out_dir <- sprintf("analyses/nmr_score_training/test_fold_%s/", this_test_fold)

  # Plot model fits
  g <- ggplot(this_fitdt) +
    aes(x=log(lambda), y=fit_mean, ymin=fit_mean_minus_sd, ymax=fit_mean_plus_sd,
        fill=factor(alpha), colour=factor(alpha)) +
    geom_ribbon(colour="#00000000", alpha=0.3, show.legend=FALSE) +
    geom_line() +
    facet_grid2(endpoint + sex ~ alpha, scales="free", independent="x") +
    geom_vline(data=this_bestfit, aes(xintercept=log(lambda), colour=factor(alpha), linetype=model)) +
    geom_vline(data=this_bestbestfit, aes(xintercept=log(lambda)), linetype="twodash", color="black", show.legend=FALSE) +
    scale_fill_manual(name="Penalty mixing\n(1=lasso, 0=ridge)", values=c("0"="#5e4fa2", "0.1"="#3288bd", "0.25"="#66c2a5",
                      "0.5"="#ffff33", "0.75"="#f46d43", "0.9"="#d53e4f", "1"="#9e0142")) +
    scale_colour_manual(name="Penalty mixing\n(1=lasso, 0=ridge)", values=c("0"="#5e4fa2", "0.1"="#3288bd", "0.25"="#66c2a5",
                      "0.5"="#ffff33", "0.75"="#f46d43", "0.9"="#d53e4f", "1"="#9e0142")) +
    scale_linetype_manual(name="Lambda selection", values=c("lambda.min"="dotted", "lambda.1se"="dashed")) +
    guides(fill=guide_legend(nrow=1), colour=guide_legend(nrow=1)) +
    xlab("Log(lambda)") +
    ylab(sprintf("%s (± SD)", fitdt$fit_metric[1])) +
    theme_bw() + 
    theme(axis.text=element_text(size=6), axis.title=element_text(size=8), 
          legend.title=element_text(size=8), legend.text=element_text(size=6), 
          strip.text.x=element_blank(), strip.background.x=element_blank(), 
          strip.text.y=element_text(size=8), 
          legend.position="bottom", legend.box="vertical", 
          legend.box.background=element_blank(), legend.background=element_blank(),
          legend.box.margin=margin(-7,-3,0,-3), legend.margin=margin(-3,-3,-3,-3))
    ggsave(g, width=7.2, height=6, units="in", file=sprintf("%s/cv_coxnet_fit.pdf", out_dir))
}

# Now do a fixed viewport so we can better see the optimal stroke selection and compare across
# endpoints
for (this_test_fold in 1:5) {
  this_fitdt <- fitdt[prediction_cv_testfold == this_test_fold]
  this_bestfit <- bestfit[prediction_cv_testfold == this_test_fold]
  this_bestbestfit <- bestbestfit[prediction_cv_testfold == this_test_fold & model == "lambda.min"]

  out_dir <- sprintf("analyses/nmr_score_training/test_fold_%s/", this_test_fold)

  # Fixed view
  g <- ggplot(this_fitdt) +
    aes(x=log(lambda), y=fit_mean, ymin=fit_mean_minus_sd, ymax=fit_mean_plus_sd,
        fill=factor(alpha), colour=factor(alpha)) +
    geom_ribbon(colour="#00000000", alpha=0.3, show.legend=FALSE) +
    geom_line() +
    facet_grid2(endpoint + sex ~ alpha, scales="free_x", independent="x") +
    geom_vline(data=this_bestfit, aes(xintercept=log(lambda), colour=factor(alpha), linetype=model)) +
    geom_vline(data=this_bestbestfit, aes(xintercept=log(lambda)), linetype="twodash", color="black", show.legend=FALSE) +
    scale_fill_manual(name="Penalty mixing\n(1=lasso, 0=ridge)", values=c("0"="#5e4fa2", "0.1"="#3288bd", "0.25"="#66c2a5",
                      "0.5"="#ffff33", "0.75"="#f46d43", "0.9"="#d53e4f", "1"="#9e0142")) +
    scale_colour_manual(name="Penalty mixing\n(1=lasso, 0=ridge)", values=c("0"="#5e4fa2", "0.1"="#3288bd", "0.25"="#66c2a5",
                      "0.5"="#ffff33", "0.75"="#f46d43", "0.9"="#d53e4f", "1"="#9e0142")) +
    scale_linetype_manual(name="Lambda selection", values=c("lambda.min"="dotted", "lambda.1se"="dashed")) +
    guides(fill=guide_legend(nrow=1), colour=guide_legend(nrow=1)) +
    xlab("Log(lambda)") +
    scale_y_continuous(name=sprintf("%s (± SD)", fitdt$fit_metric[1]), limits=c(21.8, 23.8), expand=c(0,0), oob=scales::squish) +
    theme_bw() + 
    theme(axis.text=element_text(size=6), axis.title=element_text(size=8), 
          legend.title=element_text(size=8), legend.text=element_text(size=6), 
          strip.text.x=element_blank(), strip.background.x=element_blank(), 
          strip.text.y=element_text(size=8), 
          legend.position="bottom", legend.box="vertical", 
          legend.box.background=element_blank(), legend.background=element_blank(),
          legend.box.margin=margin(-7,-3,0,-3), legend.margin=margin(-3,-3,-3,-3))
    ggsave(g, width=7.2, height=6, units="in", file=sprintf("%s/cv_coxnet_fit_fixed_zoom.pdf", out_dir))
}


