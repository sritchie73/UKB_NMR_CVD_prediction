library(data.table)

system("mkdir -p output")

raw <- fread("data/raw/ukbiobank/extracted/smoking.csv")
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

# Convert field ids to variable names
setnames(raw, as.character(info$field.id), info$var)

# Curate data on current smoking and smoking status
# (raw contains the full set of fields relating to smoking if we need them later)
smoking_status <- raw[repeat_index == 0, .(eid, visit_index, smoking_status)]

# Convert field codes to relevant labels (TODO for most variables)
smoking_status[, smoking_status := fcase(
  smoking_status == 0, "Never",
  smoking_status == 1, "Previous",
  smoking_status == 2, "Current",
  smoking_status == -3, "Prefer not to answer",
  default = NA_character_)]

# Define current smoking
smoking_status[, current_smoker := fcase(
  smoking_status == "Current", TRUE,
  smoking_status == "Prefer not to answer", NA,
  default = FALSE)]

# Add number of cigarettes smoked per day
smoking_status[raw[repeat_index == 0], on = .(eid, visit_index), daily_cigarettes := as.numeric(i.daily_cigarettes)]
smoking_status[!(current_smoker), daily_cigarettes := 0]
smoking_status[daily_cigarettes == -1, daily_cigarettes := NA] # Do not know
smoking_status[daily_cigarettes == -3, daily_cigarettes := NA] # Prefer not to answer
smoking_status[daily_cigarettes == -10, daily_cigarettes := 0.5] # Less than one a day

fwrite(smoking_status, sep="\t", quote=FALSE, file="output/smoking_status.txt")
