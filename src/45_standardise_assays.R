library(data.table)

system("mkdir -p data/standardised/")

# load imputed assay data
assay_dat <- fread("data/imputed/assays.txt")

# load assay metadata
assay_info <- fread('data/ukb/biomarkers/output/biomarker_info.txt')

# Add sex
dat <- fread("data/imputed/analysis_cohort.txt")
assay_dat[dat, on = .(eid), sex := i.sex]

# Compute means and standard deviations of assays
assay_long <- melt(assay_dat, id.vars=c("eid", "sex"), variable.name="biomarker", value.name="concentration")
assay_scaling <- assay_long[,.(mean=mean(concentration), sd=sd(concentration)), by=.(sex, biomarker)]
assay_scaling[assay_info, on = .(biomarker=var), units := i.units]
fwrite(assay_scaling, sep="\t", quote=FALSE, file="data/standardised/assay_scaling_factors.txt")

# Standardise assays
assay_long[assay_scaling, on = .(biomarker, sex), standardised := (concentration - mean)/sd]
assay_scaled <- dcast(assay_long, eid ~ biomarker, value.var="standardised")

# Give same row order as dat
assay_scaled <- assay_scaled[dat[,.(eid)], on = .(eid)]

# Write out
fwrite(assay_scaled, sep="\t", quote=FALSE, file="data/standardised/assay_concentrations.txt")

