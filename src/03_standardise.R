library(data.table)
library(ukbnmr)

system("mkdir -p data/standardised/")

# load imputed data
dat <- fread("data/imputed/analysis_cohort.txt")

# Compute means and standard deviations of non-derived NMR biomarkers
nmr_long <- melt(dat, id.vars="eid", measure.vars=ukbnmr::nmr_info[Type == "Non-derived", Biomarker],
  variable.name="biomarker", value.name="concentration")

nmr_scaling <- nmr_long[,.(mean=mean(concentration), sd=sd(concentration)), by=biomarker]
nmr_scaling[ukbnmr::nmr_info, on = .(biomarker=Biomarker), units := i.Units]
fwrite(nmr_scaling, sep="\t", quote=FALSE, file="data/standardised/nmr_scaling_factors.txt")

# Compute scaling factors for PRS
prs_scaling <- rbind(idcol="PRS", 
  "CAD_metaGRS"=dat[,.(mean=mean(CAD_metaGRS), sd=sd(CAD_metaGRS))],
  "Stroke_metaGRS"=dat[,.(mean=mean(Stroke_metaGRS), sd=sd(Stroke_metaGRS))]
)
fwrite(prs_scaling, sep="\t", quote=FALSE, file="data/standardised/prs_scaling_factors.txt")

# Standardise NMR biomarkers
nmr_long[nmr_scaling, on = .(biomarker), standardised := (concentration - mean)/sd]
nmr_scaled <- dcast(nmr_long, eid ~ biomarker, value.var="standardised")

# Give same row order as dat
nmr_scaled <- nmr_scaled[dat[,.(eid)], on = .(eid)]

# Write out
fwrite(nmr_scaled, sep="\t", quote=FALSE, file="data/standardised/non_derived_nmr.txt")






