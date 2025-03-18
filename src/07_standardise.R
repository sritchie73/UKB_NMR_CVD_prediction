library(data.table)

system("mkdir -p data/standardised/")

# load imputed data
dat <- fread("data/imputed/analysis_cohort.txt")

# Load NMR information
nmr_info <- fread("data/ukb/NMR_metabolomics/biomarker_information.txt")
non_derived <- nmr_info[Type == "Non-derived", Biomarker]
non_derived <- setdiff(non_derived, "Clinical_LDL_C")
composite <- nmr_info[Type == "Composite", Biomarker]

# Compute means and standard deviations of non-derived and composite NMR biomarkers
nmr_long <- melt(dat, id.vars=c("eid", "sex"), measure.vars=c(non_derived, composite),
  variable.name="biomarker", value.name="concentration")

nmr_scaling <- nmr_long[,.(mean=mean(concentration), sd=sd(concentration)), by=.(sex, biomarker)]
nmr_scaling[nmr_info, on = .(biomarker=Biomarker), units := i.Units]
fwrite(nmr_scaling, sep="\t", quote=FALSE, file="data/standardised/nmr_scaling_factors.txt")

# Standardise NMR biomarkers
nmr_long[nmr_scaling, on = .(biomarker, sex), standardised := (concentration - mean)/sd]
nmr_scaled <- dcast(nmr_long, eid ~ biomarker, value.var="standardised")

# Give same row order as dat
nmr_scaled <- nmr_scaled[dat[,.(eid)], on = .(eid)]

# Write out
fwrite(nmr_scaled, sep="\t", quote=FALSE, file="data/standardised/nmr_concentrations.txt")

