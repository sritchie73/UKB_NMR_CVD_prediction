library(data.table)

system("mkdir -p output")

raw <- fread("data/raw/ukbiobank/extracted/family_history.csv")
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

# Convert codes to labels
raw[, father := fcase(
  father == 13, "Prostate cancer",
  father == 12, "Severe depression",
  father == 11, "Parkinson's disease",
  father == 10, "Alzheimer's disease/dementia",
  father == 9, "Diabetes",
  father == 8, "High blood pressure",
  father == 6, "Chronic bronchitis/emphysema",
  father == 5, "Breast cancer",
  father == 4, "Bowel cancer",
  father == 3, "Lung cancer",
  father == 2, "Stroke",
  father == 1, "Heart disease",
  father == -11, "Do not know (group 1)",
  father == -13, "Prefer not to answer (group 1)",
  father == -17, "None of the above (group 1)",
  father == -21, "Do not know (group 2)",
  father == -23, "Prefer not to answer (group 2)",
  father == -27, "None of the above (group 2)"
)]

raw[, mother := fcase(
  mother == 13, "Prostate cancer",
  mother == 12, "Severe depression",
  mother == 11, "Parkinson's disease",
  mother == 10, "Alzheimer's disease/dementia",
  mother == 9, "Diabetes",
  mother == 8, "High blood pressure",
  mother == 6, "Chronic bronchitis/emphysema",
  mother == 5, "Breast cancer",
  mother == 4, "Bowel cancer",
  mother == 3, "Lung cancer",
  mother == 2, "Stroke",
  mother == 1, "Heart disease",
  mother == -11, "Do not know (group 1)",
  mother == -13, "Prefer not to answer (group 1)",
  mother == -17, "None of the above (group 1)",
  mother == -21, "Do not know (group 2)",
  mother == -23, "Prefer not to answer (group 2)",
  mother == -27, "None of the above (group 2)"
)]

raw[, siblings := fcase(
  siblings == 13, "Prostate cancer",
  siblings == 12, "Severe depression",
  siblings == 11, "Parkinson's disease",
  siblings == 10, "Alzheimer's disease/dementia",
  siblings == 9, "Diabetes",
  siblings == 8, "High blood pressure",
  siblings == 6, "Chronic bronchitis/emphysema",
  siblings == 5, "Breast cancer",
  siblings == 4, "Bowel cancer",
  siblings == 3, "Lung cancer",
  siblings == 2, "Stroke",
  siblings == 1, "Heart disease",
  siblings == -11, "Do not know (group 1)",
  siblings == -13, "Prefer not to answer (group 1)",
  siblings == -17, "None of the above (group 1)",
  siblings == -21, "Do not know (group 2)",
  siblings == -23, "Prefer not to answer (group 2)",
  siblings == -27, "None of the above (group 2)"
)]

# Build wide datasets for each disease
father <- unique(raw[,.(eid, visit_index)])

father[, no_disease := FALSE]
none <- raw[father == "None of the above (group 1)", .(eid, visit_index)]
none <- none[raw[father == "None of the above (group 2)", .(eid, visit_index)], on = .(eid, visit_index), nomatch=0]
father[none, on = .(eid, visit_index), no_disease := TRUE]
father[raw[father %like% "Do not know" | father %like% "Prefer not to answer"], on = .(eid, visit_index), no_disease := NA]

father[, heart_disease := FALSE]
father[raw[father == "Heart disease"], on = .(eid, visit_index), heart_disease := TRUE]
father[raw[father %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), heart_disease := NA]

father[, stroke := FALSE]
father[raw[father == "Stroke"], on = .(eid, visit_index), stroke := TRUE]
father[raw[father %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), stroke := NA]

father[, cardiovascular_disease := heart_disease | stroke]

father[, copd := FALSE]
father[raw[father == "Chronic bronchitis/emphysema"], on = .(eid, visit_index), copd := TRUE]
father[raw[father %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), copd := NA]

father[, hypertension := FALSE]
father[raw[father == "High blood pressure"], on = .(eid, visit_index), hypertension := TRUE]
father[raw[father %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), hypertension := NA]

father[, diabetes := FALSE]
father[raw[father == "Diabetes"], on = .(eid, visit_index), diabetes := TRUE]
father[raw[father %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), diabetes := NA]

father[, alzheimers := FALSE]
father[raw[father == "Alzheimer's disease/dementia"], on = .(eid, visit_index), alzheimers := TRUE]
father[raw[father %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), alzheimers := NA]

father[, parkinsons := FALSE]
father[raw[father == "Parkinson's disease"], on = .(eid, visit_index), parkinsons := TRUE]
father[raw[father %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), parkinsons := NA]

father[, depression := FALSE]
father[raw[father == "Severe depression"], on = .(eid, visit_index), depression := TRUE]
father[raw[father %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), depression := NA]

father[, lung_cancer := FALSE]
father[raw[father == "Lung cancer"], on = .(eid, visit_index), lung_cancer := TRUE]
father[raw[father %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), lung_cancer := NA]

father[, bowel_cancer := FALSE]
father[raw[father == "Bowel cancer"], on = .(eid, visit_index), bowel_cancer := TRUE]
father[raw[father %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), bowel_cancer := NA]

father[, breast_cancer := FALSE] # Not an option for illness of father field

father[, prostate_cancer := FALSE]
father[raw[father == "Prostate cancer"], on = .(eid, visit_index), prostate_cancer := TRUE]
father[raw[father %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), prostate_cancer := NA]

father[, any_cancer := lung_cancer | bowel_cancer | breast_cancer | prostate_cancer]

# ----------

mother <- unique(raw[,.(eid, visit_index)])

mother[, no_disease := FALSE]
none <- raw[mother == "None of the above (group 1)", .(eid, visit_index)]
none <- none[raw[mother == "None of the above (group 2)", .(eid, visit_index)], on = .(eid, visit_index), nomatch=0]
mother[none, on = .(eid, visit_index), no_disease := TRUE]
mother[raw[mother %like% "Do not know" | mother %like% "Prefer not to answer"], on = .(eid, visit_index), no_disease := NA]

mother[, heart_disease := FALSE]
mother[raw[mother == "Heart disease"], on = .(eid, visit_index), heart_disease := TRUE]
mother[raw[mother %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), heart_disease := NA]

mother[, stroke := FALSE]
mother[raw[mother == "Stroke"], on = .(eid, visit_index), stroke := TRUE]
mother[raw[mother %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), stroke := NA]

mother[, cardiovascular_disease := heart_disease | stroke]

mother[, copd := FALSE]
mother[raw[mother == "Chronic bronchitis/emphysema"], on = .(eid, visit_index), copd := TRUE]
mother[raw[mother %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), copd := NA]

mother[, hypertension := FALSE]
mother[raw[mother == "High blood pressure"], on = .(eid, visit_index), hypertension := TRUE]
mother[raw[mother %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), hypertension := NA]

mother[, diabetes := FALSE]
mother[raw[mother == "Diabetes"], on = .(eid, visit_index), diabetes := TRUE]
mother[raw[mother %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), diabetes := NA]

mother[, alzheimers := FALSE]
mother[raw[mother == "Alzheimer's disease/dementia"], on = .(eid, visit_index), alzheimers := TRUE]
mother[raw[mother %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), alzheimers := NA]

mother[, parkinsons := FALSE]
mother[raw[mother == "Parkinson's disease"], on = .(eid, visit_index), parkinsons := TRUE]
mother[raw[mother %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), parkinsons := NA]

mother[, depression := FALSE]
mother[raw[mother == "Severe depression"], on = .(eid, visit_index), depression := TRUE]
mother[raw[mother %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), depression := NA]

mother[, lung_cancer := FALSE]
mother[raw[mother == "Lung cancer"], on = .(eid, visit_index), lung_cancer := TRUE]
mother[raw[mother %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), lung_cancer := NA]

mother[, bowel_cancer := FALSE]
mother[raw[mother == "Bowel cancer"], on = .(eid, visit_index), bowel_cancer := TRUE]
mother[raw[mother %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), bowel_cancer := NA]

mother[, breast_cancer := FALSE]
mother[raw[mother == "Breast cancer"], on = .(eid, visit_index), breast_cancer := TRUE]
mother[raw[mother %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), breast_cancer := NA]

mother[, prostate_cancer := FALSE] # Not an option for illness of mother field

mother[, any_cancer := lung_cancer | bowel_cancer | breast_cancer | prostate_cancer]

# ----- 

siblings <- unique(raw[,.(eid, visit_index)])

siblings[, no_disease := FALSE]
none <- raw[siblings == "None of the above (group 1)", .(eid, visit_index)]
none <- none[raw[siblings == "None of the above (group 2)", .(eid, visit_index)], on = .(eid, visit_index), nomatch=0]
siblings[none, on = .(eid, visit_index), no_disease := TRUE]
siblings[raw[siblings %like% "Do not know" | siblings %like% "Prefer not to answer"], on = .(eid, visit_index), no_disease := NA]

siblings[, heart_disease := FALSE]
siblings[raw[siblings == "Heart disease"], on = .(eid, visit_index), heart_disease := TRUE]
siblings[raw[siblings %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), heart_disease := NA]

siblings[, stroke := FALSE]
siblings[raw[siblings == "Stroke"], on = .(eid, visit_index), stroke := TRUE]
siblings[raw[siblings %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), stroke := NA]

siblings[, cardiovascular_disease := heart_disease | stroke]

siblings[, copd := FALSE]
siblings[raw[siblings == "Chronic bronchitis/emphysema"], on = .(eid, visit_index), copd := TRUE]
siblings[raw[siblings %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), copd := NA]

siblings[, hypertension := FALSE]
siblings[raw[siblings == "High blood pressure"], on = .(eid, visit_index), hypertension := TRUE]
siblings[raw[siblings %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), hypertension := NA]

siblings[, diabetes := FALSE]
siblings[raw[siblings == "Diabetes"], on = .(eid, visit_index), diabetes := TRUE]
siblings[raw[siblings %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), diabetes := NA]

siblings[, alzheimers := FALSE]
siblings[raw[siblings == "Alzheimer's disease/dementia"], on = .(eid, visit_index), alzheimers := TRUE]
siblings[raw[siblings %in% c("Do not know (group 1)", "Prefer not to answer (group 1)")], on = .(eid, visit_index), alzheimers := NA]

siblings[, parkinsons := FALSE]
siblings[raw[siblings == "Parkinson's disease"], on = .(eid, visit_index), parkinsons := TRUE]
siblings[raw[siblings %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), parkinsons := NA]

siblings[, depression := FALSE]
siblings[raw[siblings == "Severe depression"], on = .(eid, visit_index), depression := TRUE]
siblings[raw[siblings %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), depression := NA]

siblings[, lung_cancer := FALSE]
siblings[raw[siblings == "Lung cancer"], on = .(eid, visit_index), lung_cancer := TRUE]
siblings[raw[siblings %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), lung_cancer := NA]

siblings[, bowel_cancer := FALSE]
siblings[raw[siblings == "Bowel cancer"], on = .(eid, visit_index), bowel_cancer := TRUE]
siblings[raw[siblings %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), bowel_cancer := NA]

siblings[, breast_cancer := FALSE]
siblings[raw[siblings == "Breast cancer"], on = .(eid, visit_index), breast_cancer := TRUE]
siblings[raw[siblings %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), breast_cancer := NA]

siblings[, prostate_cancer := FALSE]
siblings[raw[siblings == "Prostate cancer"], on = .(eid, visit_index), prostate_cancer := TRUE]
siblings[raw[siblings %in% c("Do not know (group 2)", "Prefer not to answer (group 2)")], on = .(eid, visit_index), prostate_cancer := NA]

siblings[, any_cancer := lung_cancer | bowel_cancer | breast_cancer | prostate_cancer]

# ------

p1 <- melt(father, id.vars=c("eid", "visit_index"), value.name="father")
p2 <- melt(mother, id.vars=c("eid", "visit_index"), value.name="mother")
parents <- p1[p2, on = .(eid, visit_index, variable)]
parents[variable != "no_disease", parents := father | mother]
parents[variable == "no_disease", parents := father & mother]
parents <- dcast(parents, eid + visit_index ~ variable, value.var="parents")
parents <- parents[, .SD, .SDcols=names(father)]
parents <- parents[unique(raw[,.(eid, visit_index)]), on = .(eid, visit_index)]

# -----

p1 <- melt(father, id.vars=c("eid", "visit_index"), value.name="father")
p2 <- melt(mother, id.vars=c("eid", "visit_index"), value.name="mother")
p3 <- melt(siblings, id.vars=c("eid", "visit_index"), value.name="siblings")
relatives <- p1[p2, on = .(eid, visit_index, variable)][p3, on = .(eid, visit_index, variable)]
relatives[variable != "no_disease", relatives := father | mother | siblings]
relatives[variable == "no_disease", relatives := father & mother & siblings]
relatives <- dcast(relatives, eid + visit_index ~ variable, value.var="relatives")
relatives <- relatives[, .SD, .SDcols=names(father)]
relatives <- relatives[unique(raw[,.(eid, visit_index)]), on = .(eid, visit_index)]

# Write out
fwrite(father, sep="\t", quote=FALSE, file="output/illness_of_father.txt")
fwrite(mother, sep="\t", quote=FALSE, file="output/illness_of_mother.txt")
fwrite(siblings, sep="\t", quote=FALSE, file="output/illness_of_siblings.txt")
fwrite(parents, sep="\t", quote=FALSE, file="output/illness_of_parents.txt")
fwrite(relatives, sep="\t", quote=FALSE, file="output/illness_of_first_degree_relatives.txt")
