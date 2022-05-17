library(data.table)

system("mkdir -p output")

raw <- fread("data/raw/ukbiobank/extracted/blood_pressure.csv")
info <- fread("data/raw/ukbiobank/extracted/field_info.txt")

# Split out instance (visit) and array index (repeat measure) fields so they
# are rows instead of columns
visit_repeats <- setdiff(unique(gsub("^[0-9]+-", "", names(raw))), "eid")
raw <- rbindlist(fill=TRUE, use.names=TRUE, lapply(visit_repeats, function(vr) {
	# Find columns matching this visit repeat pair (e.g. ending in -0.0)
	this_cols <- names(raw)[grepl(pattern=paste0(vr, "$"), names(raw))]

	# Filter to these columns
	this_raw <- raw[, .SD, .SDcols=c("eid", this_cols)]

	# Drop repeat visit pair label from column name
	setnames(this_raw, this_cols, gsub(paste0("-", vr, "$"), "", this_cols))

	# Add columns for visit and repeat index
	this_raw[, visit_index := as.integer(gsub("\\..*", "", vr))]
	this_raw[, repeat_index := as.integer(gsub(".*\\.", "", vr))]

  # Move to start of data table
  this_raw <- this_raw[,.SD,.SDcols=c("eid", "visit_index", "repeat_index", gsub("-.*", "", this_cols))]

	# Drop instance and array index combinations with all missing data
	# eid, visit_index, and array_index always non-missing
	this_raw <- this_raw[apply(this_raw, 1, function(row) { sum(!is.na(row)) > 3L })]

  # Return
  this_raw
}))

# Set names
setnames(raw, c("93", "94", "95", "102", "4079", "4080"), 
  c("sbp_manual", "dbp_manual", "pulse_rate_manual",
    "pulse_rate_auto", "dbp_automatic", "sbp_automatic"))

# Average repeated measures:
avg <- function(xx) {
  if(all(is.na(xx))) {
    return(NA_real_)
  } else {
    return(mean(na.omit(xx)))
  }
}

raw <- raw[, .(
  sbp_manual = avg(sbp_manual),
  dbp_manual = avg(dbp_manual),
  pulse_rate_manual = avg(pulse_rate_manual),
  sbp_automatic = avg(sbp_automatic),
  dbp_automatic = avg(dbp_automatic),
  pulse_rate_auto = avg(pulse_rate_auto)
), by=.(eid, visit_index)]

# Average automatic and manual measurements where both are present
raw <- raw[, .(
  sbp = avg(c(sbp_manual, sbp_automatic)),
  dbp = avg(c(dbp_manual, dbp_automatic)),
  pulse_rate = avg(c(pulse_rate_manual, pulse_rate_auto))
), by=.(eid, visit_index)]

# Write out
fwrite(raw, sep="\t", quote=FALSE, file="output/blood_pressure.txt")








