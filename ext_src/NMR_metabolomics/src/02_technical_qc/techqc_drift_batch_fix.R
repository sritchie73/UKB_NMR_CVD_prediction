library(data.table)
library(ukbnmr)
library(MASS)

system("mkdir -p output/techqc_drift_batch_fix")

# Load NMR data with basic sample QC and reformat to match expectations of UKBNMR package
bio <- fread("output/basic_sample_QC/biomarker_measurements.txt")
setnames(bio, c("sample_id", "visit"), c("eid", "visit_index"))
bio[, visit_index := fcase(
  visit_index == "Main Phase", 0L,
  visit_index == "Repeat Assessment", 1L)]

# Load biomarker tags and reformat to match expectations of UKBNMR package
bio_qc <- fread("output/basic_sample_QC/biomarker_QC_tags.txt")
setnames(bio_qc, c("sample_id", "visit", "biomarker", "QC_tag"), c("eid", "visit_index", "Biomarker", "value"))
bio_qc[, visit_index := fcase(
  visit_index == "Main Phase", 0L,
  visit_index == "Repeat Assessment", 1L)]
bio_qc[, value := gsub(";", "; ", value)]
bio_qc[, value := gsub("_", " ", value)]

# Load sample information and reformat to match expectations of UKBNMR package
sinfo <- fread("output/basic_sample_QC/sample_information.txt")
sinfo <- sinfo[!(removed), .(eid=sample_id, visit_index=visit, Shipment.Plate=plate_id, Spectrometer=spectrometer, 
  Well.Position.Within.Plate=plate_position, Well.Row=well_row, Well.Column=well_column, Shipping.Batch=batch,
  Sample.Measured.Date=date_MEASURED_AT, Prep.to.Measure.Duration=prep_to_measured, tag)]
sinfo[tag %like% "High_lactate", High.Lactate := "Yes"]
sinfo[tag %like% "High_pyruvate", High.Pyruvate := "Yes"]
sinfo[tag %like% "Low_glucose", Low.Glucose := "Yes"]
sinfo[tag %like% "Low_protein", Low.Protein := "Yes"]
sinfo[, tag := NULL]
sinfo[, visit_index := fcase(
  visit_index == "Main Phase", 0L,
  visit_index == "Repeat Assessment", 1L)]

# Get majority date of sample measurement per shipment plate
plate_measure_dates <- sinfo[,.N,by=list(Shipment.Plate, Sample.Measured.Date)]
plate_measure_dates <- plate_measure_dates[,.SD[which.max(N)],by=Shipment.Plate]
sinfo[plate_measure_dates, on = list(Shipment.Plate), Plate.Measured.Date := i.Sample.Measured.Date]

# Here, we modify the original procedure so that instead of splitting each spectrometer into
# 10 bins, we (1) make sure the bin size is consistent regardless of number of samples measured
# by the spectrometer (which can vary in size dramatically, and increases over data releases)
# and (2) et a hard split for spectrometer 10278626 (number 5 by color) at plate 490000006726 where
# we see visually a spectrometer recalibration event in Alanine concentrations
sinfo <- sinfo[order(Plate.Measured.Date)][order(Spectrometer)]
sinfo[, Spectrometer := as.character(Spectrometer)]
spec_to_split <- sinfo[Shipment.Plate == "490000006726", Spectrometer][1]
sinfo[, row := .I]
idx_to_split <- sinfo[Shipment.Plate == "490000006726", max(row)]
sinfo[Spectrometer == spec_to_split & row > idx_to_split, Spectrometer := paste(Spectrometer, "Group 2")]
sinfo[, row := NULL]

# Split a series of dates into 10 equal size bins, ordered by date
bin_dates <- function(date, n=10) {
  date_as_int <- as.integer(date)
  date_order <- as.integer(factor(date_as_int))
  bins <- cut(unique(date_order), n, labels=FALSE)
  bin_map <- data.table(date=unique(date_order), bin=bins)
  bin_map[data.table(date=date_order), on=list(date), bin]
}

sinfo[, Spectrometer.Date.Bin := bin_dates(Plate.Measured.Date, floor(.N/2000)), by=Spectrometer]
offset <- sinfo[,.(max_bin=max(Spectrometer.Date.Bin)),by=Spectrometer][order(Spectrometer)]
offset[-1, offset := offset[,cumsum(max_bin)[-.N]]]
offset[1, offset := 0]
sinfo[offset, on = .(Spectrometer), Spectrometer.Date.Bin := Spectrometer.Date.Bin + i.offset]

# Melt to long, filtering to non-derived biomarkers and dropping missing values
bio <- melt(bio, id.vars=c("eid", "visit_index"), variable.name="Biomarker", na.rm=TRUE,
  measure.vars=intersect(names(bio), ukbnmr::nmr_info[Type == "Non-derived", Biomarker]))

# add in relevant sample processing information
bio <- sinfo[, .(eid, visit_index, Well.Row, Well.Column, Spectrometer.Date.Bin, Shipping.Batch,
  Spectrometer, Shipment.Plate, Prep.to.Measure.Duration)][bio, on = .(eid, visit_index), nomatch=0]

# Determine offset for log transformation (required for variables with measurements == 0):
# Acetate, Acetoacetate, Albumin, bOHbutyrate, Clinical_LDL_C, Gly, Ile, L_LDL_CE, L_LDL_FC, M_LDL_CE, Phosphatidylc
log_offset <- bio[!is.na(value), .(Minimum=min(value), Minimum.Non.Zero=min(value[value != 0])),by=Biomarker]
log_offset[, Log.Offset := ifelse(Minimum == 0, Minimum.Non.Zero / 2, 0)]

# Get log transformed raw value
bio[log_offset, on = list(Biomarker), log_value := log(value + Log.Offset)]

# Adjust for time between sample prep and sample measurement
bio[, adj := rlm(log_value ~ log(Prep.to.Measure.Duration))$residuals, by=Biomarker]

# Adjust for within plate structure across 96-well plate rows A-H
bio[, adj := rlm(adj ~ ukbnmr:::factor_by_size(Well.Row))$residuals, by=.(Shipping.Batch, Biomarker)]

# Adjust for within plate structure across 96-well plate columns 1-12
bio[, adj := rlm(adj ~ ukbnmr:::factor_by_size(Well.Column))$residuals, by=.(Shipping.Batch, Biomarker)]

# Adjust for drift over time within spectrometer
bio[, adj := rlm(adj ~ ukbnmr:::factor_by_size(Spectrometer.Date.Bin))$residuals, by=list(Biomarker, Spectrometer)]

# Rescale to absolute units. First, shift the residuals to have the same central
# parameter estimate (e.g. mean, estimated via robust linear regression) as the
# log concentrations
bio[, adj := adj + as.vector(coef(rlm(log_value ~ 1)))[1], by=Biomarker]

# Undo the log transform
bio[, adj := exp(adj)]

# Remove the log offset
bio[log_offset, on = list(Biomarker), adj := adj - Log.Offset]

# Some values that were 0 are now < 0, apply small right shift
# for these biomarkers (shift is very small, i.e. the impact on
# the distribution's median is essentially numeric error)
shift <- bio[, list(Right.Shift=-pmin(0, min(adj))), by=Biomarker]
log_offset <- log_offset[shift, on = list(Biomarker)]
bio[log_offset, on = list(Biomarker), adj := adj + Right.Shift]

# Identify and remove outlier plates.
# Model plate medians as a normal distribution (across all 3,190 plates), then
# we can flag outlier plates as those > 3.60 standard deviations from the
# mean where 3.6 is derived from the range of the theoretical normal
# distribution for 3,190 samples
n_plates <- sinfo[, length(unique(Shipment.Plate))]
sdlim <- max(qnorm(ppoints(n_plates))) 
plate_medians <- bio[, list(value = median(adj)), by=list(Biomarker, Shipment.Plate)]

outlier_lim <- plate_medians[, list(
  Lower.Limit = mean(value) - sd(value) * sdlim,
  Mean.Plate.Medians = mean(value),
  Upper.Limit = mean(value) + sd(value) * sdlim
), by=Biomarker]

plate_medians[, outlier := "no"]
plate_medians[outlier_lim, on = list(Biomarker, value < Lower.Limit), outlier := "low"]
plate_medians[outlier_lim, on = list(Biomarker, value > Upper.Limit), outlier := "high"]

# Add outlier plate tags to biomarker qc tags
bio_qc <- bio_qc[Biomarker %in% ukbnmr::nmr_info[Type == "Non-derived", Biomarker]]

outlier_flags <- sinfo[plate_medians[outlier != "no"], on = .(Shipment.Plate),
                     .(eid, visit_index, Biomarker, value=ifelse(
                       outlier == "high", "High outlier plate", "Low outlier plate"))]

bio_qc <- rbind(bio_qc, outlier_flags)
bio_qc <- bio_qc[, .(value = paste(value, collapse="; ")), by=list(eid, visit_index, Biomarker)]

# Create files for intermediate steps for diagnostic plotting
system("mkdir -p output/techqc_drift_batch_fix/intermediate_steps")

# 1. Create version of biomarker measurements keeping outlier plates

# cast to wide
bio_with_outliers <- dcast(bio, eid + visit_index ~ Biomarker, value.var="adj")

# recompute derived biomarkers
bio_with_outliers <- ukbnmr:::nightingale_composite_biomarker_compute(bio_with_outliers)
bio_with_outliers <- ukbnmr:::nightingale_ratio_compute(bio_with_outliers)
bio_with_outliers <- ukbnmr:::extended_ratios_compute(bio_with_outliers)

# reformat to match pre-release data
setnames(bio_with_outliers, c("eid", "visit_index"), c("sample_id", "visit"))
bio_with_outliers[, visit := fcase(
  visit == 0L, "Main Phase",
  visit == 1L, "Repeat Assessment")]

# Write out
fwrite(bio_with_outliers, sep="\t", quote=FALSE, file="output/techqc_drift_batch_fix/intermediate_steps/biomarker_measurements_all_samples.txt")

# 2. Create version with outlier plates removed
bio <- bio[!plate_medians[outlier != "no"], on = .(Biomarker, Shipment.Plate)]
bio <- dcast(bio, eid + visit_index ~ Biomarker, value.var="adj")

bio <- ukbnmr:::nightingale_composite_biomarker_compute(bio)
bio <- ukbnmr:::nightingale_ratio_compute(bio)
bio <- ukbnmr:::extended_ratios_compute(bio)

setnames(bio, c("eid", "visit_index"), c("sample_id", "visit"))
bio[, visit := fcase(
  visit == 0L, "Main Phase",
  visit == 1L, "Repeat Assessment")]

fwrite(bio, sep="\t", quote=FALSE, file="output/techqc_drift_batch_fix/intermediate_steps/biomarker_measurements_all_samples_outlier_plates_removed.txt")

# Reload sinfo and add in extra columns
sinfo_qc <- copy(sinfo)
sinfo_qc[, visit := fcase(
  visit_index == 0L, "Main Phase",
  visit_index == 1L, "Repeat Assessment")]

sinfo <- fread("output/basic_sample_QC/sample_information.txt")
sinfo[sinfo_qc, on = .(sample_id=eid, visit), spectrometer_plate_bin := Spectrometer.Date.Bin]
sinfo[(removed), spectrometer_plate_bin := NA]

fwrite(sinfo, sep="\t", quote=FALSE, file="output/techqc_drift_batch_fix/intermediate_steps/sample_information.txt")

# Curate biomarker tags
bio_qc <- dcast(bio_qc, eid + visit_index ~ Biomarker, value.var="value")

bio_qc <- ukbnmr:::nightingale_composite_biomarker_flags(bio_qc)
bio_qc <- ukbnmr:::nightingale_ratio_flags(bio_qc)
bio_qc <- ukbnmr:::extended_ratios_flags(bio_qc)

setnames(bio_qc, c("eid", "visit_index"), c("sample_id", "visit"))
bio_qc[, visit := fcase(
  visit == 0L, "Main Phase",
  visit == 1L, "Repeat Assessment")]

bio_qc <- melt(bio_qc, id.vars=c("sample_id", "visit"), variable.name="biomarker", value.name="QC_tag", na.rm=TRUE)
fwrite(bio_qc, sep="\t", quote=FALSE, file="output/techqc_drift_batch_fix/intermediate_steps/biomarker_QC_tags.txt")

# 3. Create main dataset with blind duplicates removed
sinfo <- fread("output/basic_participant_QC/sample_information.txt")
sinfo[sinfo_qc, on = .(sample_id=eid, visit), spectrometer_plate_bin := Spectrometer.Date.Bin]
sinfo[(removed), spectrometer_plate_bin := NA]
fwrite(sinfo, sep="\t", quote=FALSE, file="output/techqc_drift_batch_fix/sample_information.txt")

bio_qc <- bio_qc[data.table(biomarker=names(bio)), on = .(biomarker), nomatch=0]
bio_qc <- bio_qc[sinfo[!(removed), .(eid_30418, sample_id, visit)], on = .(sample_id, visit), nomatch=0]
bio_qc <- bio_qc[, .(eid_30418, visit, biomarker, QC_tag)]
fwrite(bio_qc, sep="\t", quote=FALSE, file="output/techqc_drift_batch_fix/biomarker_QC_tags.txt")

bio <- bio[sinfo[!(removed), .(sample_id, visit)], on = .(sample_id, visit)]
bio <- sinfo[!(removed), .(eid_30418, sample_id, visit)][bio, on = .(sample_id, visit)]
bio[, sample_id := NULL]
fwrite(bio, sep="\t", quote=FALSE, file="output/techqc_drift_batch_fix/biomarker_measurements.txt")

fwrite(log_offset, sep="\t", quote=FALSE, file="output/techqc_drift_batch_fix/log_offset.txt")
fwrite(outlier_lim, sep="\t", quote=FALSE, file="output/techqc_drift_batch_fix/outlier_plate_detection_limits.txt")

