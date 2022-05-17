library(data.table)

# Load pre-curated information sheets
biomarker_info <- fread("output/biomarker_info.txt")
sample_info <- fread("output/sample_info.txt")

# Load extracted data
extracted <- fread("data/extracted/biomarkers_extracted.csv", colClasses=c("eid"="character"))

# Split out visit repeat combinations into unique rows
visit_repeats <- setdiff(unique(gsub("^[0-9]+-", "", names(extracted))), "eid")
extracted <- rbindlist(fill=TRUE, use.names=TRUE, lapply(visit_repeats, function(vr) {
  # Find columns matching this visit repeat pair (e.g. ending in -0.0)
  this_cols <- names(extracted)[grepl(pattern=paste0(vr, "$"), names(extracted))]

	# Filter to these columns
	this_vr <- extracted[, .SD, .SDcols=c("eid", this_cols)]

  # Drop repeat visit pair label from column name
  setnames(this_vr, this_cols, gsub(paste0("-", vr, "$"), "", this_cols))

  # Add columns for visit and repeat index
  this_vr[, visit_index := as.integer(gsub("\\..*", "", vr))]
  this_vr[, repeat_index := as.integer(gsub(".*\\.", "", vr))]

  # Move to start of data table
  this_vr <- this_vr[,.SD,.SDcols=c("eid", "visit_index", "repeat_index", gsub("-.*", "", this_cols))]

  # Drop instance and array index combinations with all missing data
  # eid, visit_index, and array_index always non-missing
  this_vr <- this_vr[apply(this_vr, 1, function(row) { sum(!is.na(row)) > 3L })]

  # Return
  this_vr
}))

# Repeat index always 0
extracted[, repeat_index := NULL]

# Extract blood biomarkers
blood <- extracted[, .SD, .SDcols=c("eid", "visit_index",
  biomarker_info[sample_type != "Urine" & !is.na(UKB.Field.ID), UKB.Field.ID])]

setnames(blood, 
  biomarker_info[!is.na(UKB.Field.ID) & sample_type != "Urine", as.character(UKB.Field.ID)], 
  biomarker_info[!is.na(UKB.Field.ID) & sample_type != "Urine", var])

# Extract blood biomarker missingness reason where missing due to being above or below
# detection limits (Reportability field, coding: https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=4917)
miss <- extracted[, .SD, .SDcols=c("eid", "visit_index",
  biomarker_info[sample_type != "Urine" & !is.na(Reportability.Field.ID), Reportability.Field.ID])]

setnames(miss, 
  biomarker_info[!is.na(Reportability.Field.ID) & sample_type != "Urine", as.character(Reportability.Field.ID)], 
  biomarker_info[!is.na(Reportability.Field.ID) & sample_type != "Urine", var])

# Melt to long to determine lower/upper detection limits
blood <- melt(blood, id.vars=c("eid", "visit_index"))
miss <- melt(miss, id.vars=c("eid", "visit_index"), na.rm=TRUE)

blood[, below_detection := FALSE]
blood[miss[value %in% c(2, 4)], on = .(eid, visit_index, variable), below_detection := TRUE]

blood[, above_detection := FALSE]
blood[miss[value %in% c(3, 5)], on = .(eid, visit_index, variable), above_detection := TRUE]

# Set values above/below detection to their respective limits (with small offset of 0.0001)
lims <- blood[!is.na(value), .(min = min(value), max = max(value)), by=variable]
blood <- rbind(
  blood[!(below_detection) & !(above_detection)],
  blood[(below_detection)][lims, on = .(variable), nomatch=0, .(
    eid, visit_index, variable, value = i.min - 0.0001, 
    below_detection, above_detection
  )],
  blood[(above_detection)][lims, on = .(variable), nomatch=0, .(
    eid, visit_index, variable, value = i.max + 0.0001, 
    below_detection, above_detection
  )]
)

# Derive non-hdl
nonhdl <- merge(blood[variable == "tchol"], blood[variable == "hdl"], 
                by=c("eid", "visit_index"), suffixes=c(".tchol", ".hdl"))
nonhdl[, value := value.tchol - value.hdl]
nonhdl[, below_detection := FALSE]
nonhdl[, above_detection := FALSE]

# Upper and lower limits interact in different ways:
#
# Case numbers:
#
#    below_detection.tchol above_detection.tchol below_detection.hdl above_detection.hdl      N
# 1:                 FALSE                 FALSE               FALSE               FALSE 445445
# 2:                  TRUE                 FALSE                TRUE               FALSE      4
# 3:                 FALSE                 FALSE               FALSE                TRUE      1
# 4:                 FALSE                 FALSE                TRUE               FALSE      1
#
# Types:
#  
# 1. Both tchol and hdl are within detection limits, nonhdl is a precise estimate
# 2. Both tchol and hdl are below detection limits, the true value is in range [nonhdl, tchol].
# 3. hdl is above detection, then the true value is in range [0, nonhdl].
# 4. hdl is below detection, then the true value is in range [nonhdl, tchol].
#
nonhdl[below_detection.tchol & below_detection.hdl, above_detection := TRUE]
nonhdl[!(below_detection.tchol) & above_detection.hdl, below_detection := TRUE] 
nonhdl[!(below_detection.tchol) & below_detection.hdl, above_detection := TRUE]
  
# add to main table
nonhdl <- nonhdl[, .(eid, visit_index, variable="nonhdl", value, below_detection, above_detection)]
blood <- rbind(blood, nonhdl)

# Derive ratio of ApoB to ApoA1
apobapoa1 <- merge(blood[variable == "apob"], blood[variable == "apoa1"],
                   by=c("eid", "visit_index"), suffixes=c(".apob", ".apoa1"))
apobapoa1[, value := value.apob / value.apoa1]
apobapoa1[, below_detection := FALSE]
apobapoa1[, above_detection := FALSE]

# Upper and lower limits interact in different ways:
# 
# Case numbers:
#
#    below_detection.apob above_detection.apob below_detection.apoa1 above_detection.apoa1      N
# 1:                FALSE                FALSE                 FALSE                 FALSE 440741
# 2:                FALSE                FALSE                 FALSE                  TRUE   1792
# 3:                FALSE                 TRUE                 FALSE                 FALSE    306
# 4:                 TRUE                FALSE                 FALSE                 FALSE    730
# 5:                FALSE                FALSE                  TRUE                 FALSE     10
# 6:                 TRUE                FALSE                 FALSE                  TRUE     11
# 7:                 TRUE                FALSE                  TRUE                 FALSE     10
# 8:                FALSE                 TRUE                 FALSE                  TRUE      2
# 
# Types:
#
# 1. Both apob and apoa1 are within detection limits, apobapoa1 is a precise ratio
# 2. Only apoa1 is above detection, then the true value is in range [0, apobapoa1]
# 3. Only apob is above detection, then the true value is in range [apobapoa1, inf)
# 4. Only apob is below detection, then the true value is in range [0, apobapoa1]
# 5. Only apoa1 is below detection, then the true value is in range [apobapoa1, inf)
# 6. apob is below detection and apoa1 is above detection, then the true value is in range [0, apobapoa1]
# 7. apob and apoa1 are both below detection, then the ratio is an imprecise estimate, but could truly be anywhere between [0, Inf).
# 8. apob and apoa1 are both above detection, then the ratio is an imprecise estimate, but could truly be anywhere between [0, Inf).
apobapoa1[, below_detection := fcase(
  !is.na(value) & above_detection.apoa1, TRUE,
  !is.na(value) & below_detection.apob, TRUE,
  !is.na(value) & below_detection.apob & above_detection.apoa1, TRUE,
  below_detection.apoa1 & below_detection.apob, FALSE,
  above_detection.apoa1 & above_detection.apob, FALSE,
  default = FALSE)]
apobapoa1[, above_detection := fcase(
  !is.na(value) & below_detection.apoa1, TRUE,
  !is.na(value) & above_detection.apob, TRUE,
  below_detection.apoa1 & below_detection.apob, FALSE,
  above_detection.apoa1 & above_detection.apob, FALSE,
  default = FALSE)]

# add to main table
apobapoa1 <- apobapoa1[, .(eid, visit_index, variable="apobapoa1", value, below_detection, above_detection)]
blood <- rbind(blood, apobapoa1)

# Define and add hba1c_pct
hba1c_pct <- blood[variable == "hba1c"]
hba1c_pct[, value := value/10.929+2.15]
hba1c_pct[, variable := "hba1c_pct"]
blood <- rbind(blood, hba1c_pct)

# Define and add fasting glucose based on:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3186886/
# Basically, we compute the median glucose in males and females separately
# where fasting time is 0 hours, 1 hour, 2 hours, and 3 or more hours, then
# use the sex-specific difference in medians compared to 3 or more hours as
# the adjustment factor to obtain fasting glucose.
fasting_glucose <- blood[variable == "glucose" & !is.na(value)]
fasting_glucose[extracted, on = .(eid, visit_index), fasting_time := i.74]
fasting_glucose <- fasting_glucose[!is.na(fasting_time)]
fasting_glucose[, fast_lt_3hrs := pmin(fasting_time, 3)]
sex <- fread("data/anthropometrics/anthropometrics.txt", select=c("eid", "visit_index", "sex"), colClasses=c("eid"="character"))
fasting_glucose[sex, on = .(eid, visit_index), sex := i.sex]

med_glc_by_fast_sex <- fasting_glucose[, .(med_glc=median(value)), by=.(fast_lt_3hrs, sex)]
med_glc_by_fast_sex[med_glc_by_fast_sex[fast_lt_3hrs == 3], on = .(sex), adj := med_glc - i.med_glc]
fasting_glucose[med_glc_by_fast_sex, on = .(fast_lt_3hrs, sex), value := value - i.adj]

fasting_glucose <- fasting_glucose[, .(eid, visit_index, variable="fasting_glucose", value, below_detection, above_detection)]
blood <- rbind(blood, fasting_glucose)

# Write out information on adjustment
med_glc_by_fast_sex <- med_glc_by_fast_sex[order(-fast_lt_3hrs)][order(sex)]
setnames(med_glc_by_fast_sex, c("fasting_time", "sex", "median_glucose", "adjustment_factor"))
med_glc_by_fast_sex[, fasting_time := fcase( 
  fasting_time == 3, "3 hours or more",
  fasting_time == 2, "2 hours",
  fasting_time == 1, "1 hours",
  fasting_time == 0, "0 hours"
)]
fwrite(med_glc_by_fast_sex, sep="\t", quote=FALSE, file="output/adjusting_glucose_for_fasting.txt")

# Drop missing values
blood <- blood[!is.na(value)]

# Extract urine biomarkers
urine <- extracted[, .SD, .SDcols=c("eid", "visit_index",
  biomarker_info[sample_type == "Urine" & !is.na(UKB.Field.ID), UKB.Field.ID])]

setnames(urine, 
  biomarker_info[!is.na(UKB.Field.ID) & sample_type == "Urine", as.character(UKB.Field.ID)], 
  biomarker_info[!is.na(UKB.Field.ID) & sample_type == "Urine", var])

# Extract blood biomarker missingness reason where missing due to being above or below
# detection limits (Reportability field, coding: https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=4917)
miss <- extracted[, .SD, .SDcols=c("eid", "visit_index",
  biomarker_info[sample_type == "Urine" & !is.na(Reportability.Field.ID), Reportability.Field.ID])]

setnames(miss, 
  biomarker_info[!is.na(Reportability.Field.ID) & sample_type == "Urine", as.character(Reportability.Field.ID)], 
  biomarker_info[!is.na(Reportability.Field.ID) & sample_type == "Urine", var])

# Melt to long to determine lower/upper detection limits
urine[, uriacc := as.numeric(uriacc)] # otherwise integer column throws warning in melt
urine <- melt(urine, id.vars=c("eid", "visit_index"))
miss <- melt(miss, id.vars=c("eid", "visit_index"), na.rm=TRUE)
miss <- miss[value != ""]

urine[, below_detection := FALSE]
urine[miss[value %like% "<"], on = .(eid, visit_index, variable), below_detection := TRUE]

urine[, above_detection := FALSE]
urine[miss[value %like% ">"], on = .(eid, visit_index, variable), above_detection := TRUE]

# Set values above/below detection to their respective limits (with small offset of 0.0001)
# Unlike blood biomarkers, limits are hard coded into the reportability fields
lims <- unique(miss[,.(variable, value)])
lower <- lims[value %like% "<", .(variable, min=as.numeric(gsub("<", "", value)))]
upper <- lims[value %like% ">", .(variable, max=as.numeric(gsub(">", "", value)))]
lims <- merge(lower, upper, by = "variable", all=TRUE)

urine <- rbind(
  urine[!(below_detection) & !(above_detection)],
  urine[(below_detection)][lims, on = .(variable), nomatch=0, .(
    eid, visit_index, variable, value = i.min - 0.0001,
    below_detection, above_detection
  )],
  urine[(above_detection)][lims, on = .(variable), nomatch=0, .(
    eid, visit_index, variable, value = i.max + 0.0001,
    below_detection, above_detection
  )]
)

# Drop missing
urine <- urine[!is.na(value)]

# Build combined biomarker table
biomarkers <- rbind(blood, urine)

imputed_limits <- biomarkers[(below_detection) | (above_detection), 
  .(eid, visit_index, variable, type = ifelse(below_detection, "below detection limit", "above detection limit"))]
fwrite(imputed_limits, sep="\t", quote=FALSE, file="output/samples_outside_detection.txt")

biomarkers <- dcast(biomarkers, eid + visit_index ~ variable, value.var="value")
biomarkers <- biomarkers[extracted[,.(eid, visit_index)], on = .(eid, visit_index), nomatch=0]
fwrite(biomarkers, sep="\t", quote=FALSE, file="output/biomarkers.txt")

# Also get list of samples that were not measured (as opposed to missing data for some other reason)
not_measured <- extracted[, .(eid, visit_index, 
  no_blood_sample = ifelse(!is.na(`20050`), TRUE, FALSE), 
  no_urine_sample = ifelse(!is.na(`20072`), TRUE, FALSE))]
not_measured <- not_measured[(no_blood_sample) | (no_urine_sample)]

fwrite(not_measured, sep="\t", quote=FALSE, file="output/samples_not_measured.txt")

# Write out fasting time
fasting <- extracted[,.(eid, visit_index, fasting_time=`74`)]
fasting <- fasting[biomarkers[,.(eid, visit_index)], on = .(eid, visit_index)]
fwrite(fasting, sep="\t", quote=FALSE, file="output/fasting_time_hours.txt")

