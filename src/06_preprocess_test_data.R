library(data.table)
library(ukbnmr)
library(foreach)
library(ggplot2)
source("src/utils/logit.R")

# Create output directory
system("mkdir -p data/processed/test")

# Load test dataset
test <- fread("data/cleaned/test_data.txt")

# Flag complete data for selected biomarkers
lasso_coef <- fread("analyses/train/lasso_coefficients.txt")
selected_nmr <- lasso_coef[name == "NMR" & lambda == "lambda.min" & coef_type == "NMR Metabolomics", var]
test[, complete_data := complete.cases(test[,.SD,.SDcols=selected_nmr])]

# Convert SBP from integer to numeric
test[, sbp := as.numeric(sbp)]

# Load biomarker information sheets
bio_info <- fread("data/ukb/biomarkers/output/biomarker_info.txt")
nmr_info <- ukbnmr::nmr_info

# Examine distributions of continuous variables:
cont <- c("CAD_metaGRS", "Stroke_metaGRS", "age", "bmi", "sbp", bio_info$var, nmr_info$Biomarker)
cont <- intersect(cont, names(test))
ggdt <- melt(test, id.vars="eid", measure.vars=cont, na.rm=TRUE)

n_var <- length(cont)
n_col <- ceiling(sqrt(n_var))
plot_dim <- n_col * 1.4

g <- ggplot(ggdt, aes(x = value)) +
  geom_density(trim=TRUE, size=0.4) +
  facet_wrap(~ variable, scales="free", ncol=n_col) +
  xlab("Training data raw value") +
  ylab("Density") +
  theme_bw() +
  theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
        strip.background=element_blank(), strip.text=element_text(size=6),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  )

ggsave(g, width=plot_dim, height=plot_dim, units="in", file="data/processed/test/continuous_variable_distributions.png")

# Log or logit transform continuous variables depending on whether they are
# percentages (logit) or not (log)

# Go through each biomarker and apply the correct transformation. We also need to handle cases
# where biomarker concentrations are 0 (in which case we need to apply a small offset so these
# can be logged) and where percentages are 0 or 100 (in which case we need to shift, or shrink
# if both are present, so they can be logit transformed). The returned var_info table contains                                                                                                                                               # information about the range of values present for each biomarker and any shifting/shrinking
# applied.

# This process is the same as what we did for the training data, except here we don't impute 
# missing values to their median
var_info <- foreach(var = cont, .combine=rbind) %do% {
  # This process only applies to biomarkers, we'll handle the other five manually
  if (var %in% c("CAD_metaGRS", "Stroke_metaGRS", "age", "bmi", "sbp")) {
    return(NULL)
  }

  # Extract values
  values <- test[[var]]

  # Percentages are logit transformed whereas other biomarkers are log transformed
  is_percentage <- var %in% c(nmr_info[Units == "%", Biomarker], bio_info[units == "%", var])

  # Get range of values for the biomarker, as well as the min/max values that are not 0 or 100% so we know
  # what offset(s) need to be applied for log/logit transformation
  lim <- data.table(
    variable = var, is_percentage=is_percentage,
    min = min(na.omit(values)),
    non_zero_min = min(na.omit(values[values != 0])),
    max = max(na.omit(values)),
    non_100pct_max = ifelse(is_percentage, max(na.omit(values[values < 100])), NA_real_)
  )

  # Determine offset for transformation
  lim[, log_offset := fcase(
    !(is_percentage) & min == 0, non_zero_min / 2,  # right shift biomarkers with 0 values
    (is_percentage) & max == 100 & min > 0, (100 - non_100pct_max)/2*-1, # left shift % with 100% values
    (is_percentage) & max == 100 & min == 0, non_zero_min / 2, # right shift percentages with 0% and 100%
    default = 0)]

  # When % have both 0% and 100%, we need to shrink the values slightly as well
  lim[, log_mult := ifelse((is_percentage) & max == 100 & min == 0,
    (100 - (100 - non_100pct_max)/2)/(100 + log_offset), # multiplier to bring right shifted % back to below 100%
     1)]

  # Apply any right shift and rescaling required for log or logit transformation
  values <- (values + lim$log_offset)*lim$log_mult

  # Need to convert percentages from 0%-100% to 0.00-1.00 prior to logit transformation
  if (is_percentage) values <- values/100

  # Apply log/logit transformation
  if (is_percentage) {
    values <- logit(values)
  } else {
    values <- log(values)
  }

  # Get a summary of the distribution after transformation
  trans_summary <- as.data.table(as.list(summary(na.omit(values))))
  setnames(trans_summary, c("trans_min", "trans_L25", "trans_median", "trans_mean", "trans_U75", "trans_max"))
  lim <- cbind(lim, trans_summary)

  # Also get the standard deviation
  lim[, trans_sd := sd(values, na.rm=TRUE)]

  # Get missingness
  lim[, missing := sum(is.na(values))]
  lim[, missing_pct := paste0(round(missing/length(values)*100, digits=1), "%")]

  # Winsorize distribution at +/- 5SD
  lim[, lower_5sd := trans_mean - 5*trans_sd]
  lim[, below_5sd := sum(values < lower_5sd, na.rm=TRUE)]
  lim[, below_5sd_pct := paste0(round(below_5sd/length(values)*100, digits=1), "%")]
  values[values < lim$lower_5sd] <- lim$lower_5sd

  lim[, upper_5sd := trans_mean + 5*trans_sd]
  lim[, above_5sd := sum(values > upper_5sd, na.rm=TRUE)]
  lim[, above_5sd_pct := paste0(round(above_5sd/length(values)*100, digits=1), "%")]
  values[values > lim$upper_5sd] <- lim$upper_5sd

  # Standardisation will be applied at the model level - i.e. concentrations preserved in dataset

  # Overwrite column
  test[, c(var) := values]

  # return offset information
  return(lim)
}

# Write out transformation information
fwrite(var_info, sep="\t", quote=FALSE, file="data/processed/test/continuous_variable_transformation.txt")

# BMI and SBP also log transformed
test[, bmi := log(bmi)] # Preserve units, but log transform. Standardisation applied within models
test[, sbp := log(sbp)]

# Check distributions after transformation
ggdt <- melt(test, id.vars="eid", measure.vars=cont, na.rm=TRUE)

g <- ggplot(ggdt, aes(x = value)) +
  geom_density(trim=TRUE, size=0.4) +
  facet_wrap(~ variable, scales="free", ncol=n_col) +
  xlab("Training data transformed value") +
  ylab("Density") +
  theme_bw() +
  theme(axis.text=element_text(size=6), axis.title=element_text(size=10),
        strip.background=element_blank(), strip.text=element_text(size=6),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  )

ggsave(g, width=plot_dim, height=plot_dim, units="in", file="data/processed/test/post_transformation_distributions.png")

# Write out processed test data
fwrite(test, sep="\t", quote=FALSE, file="data/processed/test/processed_test_data.txt")
