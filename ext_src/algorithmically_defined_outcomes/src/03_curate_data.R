library(data.table)
library(lubridate)

system("mkdir -p output")

raw <- fread("data/raw/ukbiobank/extracted/algorithmically_derived_outcomes.csv")
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

# No repeat information, so drop column
raw[, repeat_index := NULL]

# Also visit index isn't relevant here either
raw[, visit_index := NULL]

# Convert unknown dates to NAs (note there are huge numbers of NAs already)
for (date_field in info[var %like% "date$", var]) {
  v <- raw[[date_field]]
  v[v == "1900-01-01"] <- NA
  raw[, c(date_field) := v]
  raw[,]
}

# Map source codings
for (source_field in info[var %like% "source$", var]) {
  v <- raw[[source_field]]
  v <- fcase(
    v == 0, "Self-reported only",
    v == 1, "Hospital admission",
    v == 2, "Death only",
    v == 11, "Hospital primary",
    v == 12, "Death primary",
    v == 21, "Hospital secondary",
    v == 22, "Death contributory",
    default = NA_character_
  )
  raw[, c(source_field) := v]
  raw[,]
}


# Write out
fwrite(raw, sep="\t", quote=FALSE, file="output/algorithmically_defined_outcomes.txt")
fwrite(info, sep="\t", quote=FALSE, file="output/algorithmically_defined_outcomes_info.txt")

