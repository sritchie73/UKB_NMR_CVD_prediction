library(data.table)
library(lubridate)

system("mkdir -p output")

raw <- fread("data/raw/ukbiobank/extracted/anthropometrics.csv")
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

# Some fields (e.g. sex) only available at baseline, copy to repeat visits
raw[raw[visit_index == 0], on = .(eid), sex := i.sex]
raw[raw[visit_index == 0], on = .(eid), birth_year := i.birth_year]
raw[raw[visit_index == 0], on = .(eid), birth_month := i.birth_month]
raw[raw[visit_index == 0], on = .(eid), ethnicity := i.ethnicity]
raw[raw[visit_index == 0], on = .(eid), townsend := i.townsend]

# Convert codes for sex
raw[, sex := ifelse(sex == 1L, "Male", "Female")]

# Convert codes for pregnancy
raw[, pregnant := fcase(
  sex == "Male", "No",
  pregnant == 0L, "No",
  pregnant == 1L, "Yes",
  pregnant == 2L, "Unsure"
)]

# Convert assessment centre codes to locations
raw[, assessment_centre := fcase(
  assessment_centre == "11012", "Barts",
  assessment_centre == "11021", "Birmingham",
  assessment_centre == "11011", "Bristol",
  assessment_centre == "11008", "Bury",
  assessment_centre == "11003", "Cardiff",
  assessment_centre == "11024", "Cheadle (revisit)",
  assessment_centre == "11020", "Croydon",
  assessment_centre == "11005", "Edinburgh",
  assessment_centre == "11004", "Glasgow",
  assessment_centre == "11018", "Hounslow",
  assessment_centre == "11010", "Leeds",
  assessment_centre == "11016", "Liverpool",
  assessment_centre == "11001", "Manchester",
  assessment_centre == "11017", "Middlesborough",
  assessment_centre == "11009", "Newcastle",
  assessment_centre == "11013", "Nottingham",
  assessment_centre == "11002", "Oxford",
  assessment_centre == "11007", "Reading",
  assessment_centre == "11014", "Sheffield",
  assessment_centre == "10003", "Stockport (pilot)",
  assessment_centre == "11006", "Stoke",
  assessment_centre == "11022", "Swansea",
  assessment_centre == "11023", "Wrexham",
  assessment_centre == "11025", "Cheadle (imaging)",
  assessment_centre == "11026", "Reading (imaging)",
  assessment_centre == "11027", "Newcastle (imaging)",
  assessment_centre == "11028", "Bristol (imaging)"
)]

# Code Nation of assessment centre
raw[, assessment_nation := "England"]
raw[assessment_centre %in% c("Cardiff", "Swansea", "Wrexham"), assessment_nation := "Wales"] 
raw[assessment_centre %in% c("Edinburgh", "Glasgow"), assessment_nation := "Scotland"]

# Convert ethnicity codes to labels. Codes and labels are organised
# into a tree, with codes in the 1000s belonging to top level codes
# starting with the same number. E.g. "British" (1001), "Irish" (1002),
# "Any other white background" (1003) are all sub-groups of "White" (1).
# https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=1001
raw[, ethnicity_subgroup := fcase(
  ethnicity == 1, "White",
  ethnicity == 1001, "British",
  ethnicity == 2001, "White and Black Caribbean",
  ethnicity == 3001, "Indian",
  ethnicity == 4001, "Caribbean",
  ethnicity == 2, "Mixed",
  ethnicity == 1002, "Irish",
  ethnicity == 2002, "White and Black African",
  ethnicity == 3002, "Pakistani",
  ethnicity == 4002, "African",
  ethnicity == 3, "Asian or Asian British",
  ethnicity == 1003, "Any other white background",
  ethnicity == 2003, "White and Asian",
  ethnicity == 3003, "Bangladeshi",
  ethnicity == 4003, "Any other Black background",
  ethnicity == 4, "Black or Black British",
  ethnicity == 2004, "Any other mixed background",
  ethnicity == 3004, "Any other Asian background",
  ethnicity == 5, "Chinese",
  ethnicity == 6, "Other ethnic group",
  ethnicity == -1, "Do not know",
  ethnicity == -3, "Prefer not to answer"
)]

raw[, ethnicity_group := fcase(
  ethnicity %in% c(1, 1001, 1002, 1003), "White",
  ethnicity %in% c(2, 2001, 2002, 2003, 2004), "Mixed",
  ethnicity %in% c(3, 3001, 3002, 3003, 3004), "Asian or Asian British",
  ethnicity %in% c(4, 4001, 4002, 4003), "Black or Black British",
  ethnicity == 5, "Chinese",
  ethnicity == 6, "Other ethnic group",
  ethnicity == -1, "Do not know",
  ethnicity == -3, "Prefer not to answer"
)]

# Compute waist to hip ratio
raw[, waist_hip_ratio := waist / hip]

# Compute decimal age which incorporates birth month into age
# Note day of birth not provided for privacy reasons, so we guess
# this as the midpoint of each month
raw[, approx_birth_date := as.IDate(sprintf("%s-%s-%s", birth_year, birth_month, 15))]
raw[, age_decimal := time_length(as.Date(assessment_date) - as.Date(approx_birth_date), unit="years")]
raw[, age_decimal := round(age_decimal, digits=2)] # any further resolution is meaningless due to lack of precision

# Load in genetic information and add
geno <- fread("data/raw/ukbiobank/ukb_sqc_v2.txt")

raw[geno, on = .(eid), genetic_sex := Inferred.Gender]
raw[, genetic_sex := fcase(
  genetic_sex == "M", "Male",
  genetic_sex == "F", "Female",
  default = NA
)]

raw[geno, on = .(eid), genetic_white_british := in.white.British.ancestry.subset]
raw[, genetic_white_british := fcase(
  genetic_white_british == 1L, TRUE,
  genetic_white_british == 0L, FALSE,
  default = NA
)]

# Provide some sort of sensible ordering to columns
raw <- raw[, .(eid, visit_index, assessment_date, assessment_centre, assessment_nation, townsend,
               sex, genetic_sex, pregnant, ethnicity_group, ethnicity_subgroup, genetic_white_british,
               age, age_decimal, birth_year, birth_month, approx_birth_date, height, weight, bmi,
               waist, hip, waist_hip_ratio)]

# Add information on additional fields
info <- rbind(use.names=TRUE, fill=TRUE, info,
  data.table(var="assessment_nation", name="Nation assessment centre is located in"),
  data.table(var="age_decimal", name="Age compute using birth date and birth month"),
  data.table(var="approx_birth_date", name="Date of birth, inferred as the 15th of the 'birth_month' in the 'birth_year'"),
  data.table(var="waist_hip_ratio", name="Waist to hip ratio, computed as waist / hip"),
  data.table(var="ethnicity_group", name="Top-level ethnicity grouping in field #21000"),
  data.table(field.id=21000, var="ethnicity_subgroup", name="Ethnicity"),
  data.table(var="genetic_sex", name="Sex in UK Biobank genotype sample QC file (Resource #531) inferred from sex-specific genotypes"),
  data.table(var="genetic_white_british", name="Sample in White British ancestry subgroup in genotype PCA (see Resource #531)")
)

# Reorder rows and columns
info <- info[names(raw), on = .(var), nomatch=0]
info <- info[, .(var, name, field.id)]

# Write out
fwrite(raw, sep="\t", quote=FALSE, file="output/anthropometrics.txt")
fwrite(info, sep="\t", quote=FALSE, file="output/anthropometrics_info.txt")

