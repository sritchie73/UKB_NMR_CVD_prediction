library(data.table)

system("mkdir -p output")

raw <- fread("data/raw/ukbiobank/extracted/medication.csv")
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

# Rename fields - note touchsreen questions for standard medications
# had different fields for males and females, hence the column naming 
# and additional subsequent processing
setnames(raw, c("2492", "137", "20003", "6671", "20199", "6153", "6177"), 
  c("other_prescription_meds", "num_current_meds", "medication_code", 
    "num_antibiotics_last_3mo", "antibiotic_code", "Female", "Male"))


# Extract standard touchscreen survey medication question and convert to
# wide format
tchscrn <- raw[,.(eid, visit_index, Female, Male)]
tchscrn <- tchscrn[!is.na(Female) | !is.na(Male)]

# Convert codes to labels
tchscrn[, Female := fcase(
  Female == 1, "Cholesterol lowering medication",
  Female == 2, "Blood pressure medication",
  Female == 3, "Insulin",
  Female == 4, "Hormone replacement therapy",
  Female == 5, "Oral contraceptive pill or minipill",
  Female == -7, "None of the above",
  Female == -1, "Do not know",
  Female == -3, "Prefer not to answer"
)]

tchscrn[, Male := fcase(
  Male == 1, "Cholesterol lowering medication",
  Male == 2, "Blood pressure medication",
  Male == 3, "Insulin",
  Male == -7, "None of the above",
  Male == -1, "Do not know",
  Male == -3, "Prefer not to answer"
)]

tchscrn[, medications := ifelse(is.na(Male), Female, Male)]

# Split out medication types into a wide table
wide <- unique(tchscrn[,.(eid, visit_index)])

# Curate cholesterol lowering medication
wide[, cholesterol_medication := FALSE]
wide[tchscrn[medications == "Cholesterol lowering medication"], on = .(eid, visit_index), cholesterol_medication := TRUE]
wide[tchscrn[medications %in% c("Do not know", "Prefer not to answer")], on = .(eid, visit_index), cholesterol_medication := NA]

# Curate blood pressure lowering medication
wide[, blood_pressure_medication := FALSE]
wide[tchscrn[medications == "Blood pressure medication"], on = .(eid, visit_index), blood_pressure_medication := TRUE]
wide[tchscrn[medications %in% c("Do not know", "Prefer not to answer")], on = .(eid, visit_index), blood_pressure_medication := NA]

# Curate insulin usage
wide[, insulin_medication := FALSE]
wide[tchscrn[medications == "Insulin"], on = .(eid, visit_index), insulin_medication := TRUE]
wide[tchscrn[medications %in% c("Do not know", "Prefer not to answer")], on = .(eid, visit_index), insulin_medication := NA]

# Hormone replacement therapy (female only question)
wide[, hormone_replacement_therapy := FALSE]
wide[tchscrn[medications == "Hormone replacement therapy"], on = .(eid, visit_index), hormone_replacement_therapy := TRUE]
wide[tchscrn[!is.na(Female) & medications %in% c("Do not know", "Prefer not to answer")], on = .(eid, visit_index), hormone_replacement_therapy := NA]

# Oral contraceptive (female only question)
wide[, oral_contraceptive := FALSE]
wide[tchscrn[medications == "Oral contraceptive pill or minipill"], on = .(eid, visit_index), oral_contraceptive := TRUE]
wide[tchscrn[!is.na(Female) & medications %in% c("Do not know", "Prefer not to answer")], on = .(eid, visit_index), oral_contraceptive := NA]

# Can now rename to tchscrn
tchscrn <- wide

# Add in information about other prescription meds
tchscrn[raw[repeat_index == 0], on = .(eid, visit_index), other_prescription_meds := fcase(
  i.other_prescription_meds == 1, TRUE, # Yes - you will be asked about this later by an interviewer
  i.other_prescription_meds == 0, FALSE, # No
  i.other_prescription_meds == -1, NA, # Do not know
  i.other_prescription_meds == -3, NA, # Prefer not to answer
  is.na(i.other_prescription_meds), NA # No data available
)] 

# Add in number of current medications - note in some instances these are non-zero where
# 'other_prescription_meds' is FALSE and other columns curated above are all FALSE. This 
# can occur when a person is taking over the counter medictions, e.g. painkillers.
tchscrn[raw[repeat_index == 0], on = .(eid, visit_index), num_current_meds := i.num_current_meds]

# Add in number of antibiotics taken in last three months
tchscrn[raw[repeat_index == 0], on = .(eid, visit_index), num_antibiotics_last_3mo := i.num_antibiotics_last_3mo]

# Write out touchscreen survey medication info
fwrite(tchscrn, sep="\t", quote=FALSE, file="output/medications_simple.txt")

# Extract detailed medications from verbal interview
meds <- raw[,.(eid, visit_index, medication_code)]

# Add in labels
med_labels <- fread("data/raw/ukbiobank/docs/field_20003_medication_code_labels.tsv")
meds[med_labels, on = .(medication_code = coding), medication_name := i.meaning]

# Distinguish missing data from cases where the person was not taking any medications,
# use coding -2 for consistency with UKB codings (Not applicable - i.e. person was not
# asked this question due to answering no to early question, other option would be -7,
# none of the above).
#
meds[tchscrn[num_current_meds == 0], on = .(eid, visit_index), medication_code := -2]
meds[medication_code == -2, medication_name := "Not applicable"]

# Drop missing data i.e. where the person was not asked the earlier touchscreen interview
# data, or where 'num_current_meds' is otherwise NA.
meds <- meds[!is.na(medication_code)]

# Write out
fwrite(meds, sep="\t", quote=FALSE, file="output/detailed_medications_field_20003.txt")

# Do the same for antibiotics
anti <- raw[,.(eid, visit_index, antibiotic_code)]
anti_labels <- fread("data/raw/ukbiobank/docs/field_20199_antibiotic_code_labels.tsv")
anti[anti_labels, on = .(antibiotic_code = coding), antibiotic_name := i.meaning]
anti[tchscrn[num_antibiotics_last_3mo == 0], on = .(eid, visit_index), antibiotic_code := -2]
anti[antibiotic_code == -2, antibiotic_name := "Not applicable"]
anti <- anti[!is.na(antibiotic_code)]
fwrite(anti, sep="\t", quote=FALSE, file="output/detailed_antibiotics_field_20199.txt")

# Curate information
info <- rbind(
  data.table(var="cholesterol_medication", description="Answered \"Cholesterol lowering medication\" in field 6153 (females) or 6177 (males)"),
  data.table(var="blood_pressure_medication", description="Answered \"Blood pressure medication\" in field 6153 (females) or 6177 (males)"),
  data.table(var="insulin_medication", description="Answered \"Insulin\" in field 6153 (females) or 6177 (males)"),
  data.table(var="hormone_replacement_therapy", description="Answered \"Hormone replacement therapy\" in field 6153 (females)"),
  data.table(var="oral_contraceptive", description="Answered \"Oral contraceptive pill or minipill\" in field 6153 (females)"),
  data.table(var="other_prescription_meds", description="Answered \"Yes\" to \"Taking other prescription medications\" question (field 2492)"),
  data.table(var="num_current_meds", description="Number of treatments/medications taken (prescription or over the counter) (field 137)"),
  data.table(var="num_antibiotics_last_3mo", description="Number of antibiotics taken in last 3 months (field 6671) (note question asked only at imaging visits; visit_index 2 or 3)")
)

# Write out
fwrite(info, sep="\t", quote=FALSE, file="output/medications_simple_info.txt")

