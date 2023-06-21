#!/usr/local/Cluster-Apps/R/R.4.0.3/bin/Rscript --vanilla
########################################################################################
# Load R package dependencies
########################################################################################

# Set the R package location directory to one with the R packages pre-installed.
srcDir = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P30418/curated_phenotypes/endpoints/"
SoftwareDir = sprintf("%s/software_dependencies/", srcDir)
RpackageDir = sprintf("%s/Rpackages/", SoftwareDir)
RpackageDir = sprintf("%s/%s.%s/", RpackageDir, R.version$major, gsub("\\..*", "", R.version$minor))
if (Sys.getenv("SLURM_JOB_PARTITION") %in% c("icelake", "icelake-himem")) RpackageDir = gsub("/$", "-icelake/", RpackageDir)
.libPaths(RpackageDir)

suppressMessages(library("data.table"))
suppressMessages(library("foreach"))
suppressMessages(library("docopt"))
suppressMessages(library("bit64"))
suppressMessages(library("lubridate"))

########################################################################################
# Set up program options
########################################################################################

"Curate incidence and/or prevalence for a custom endpoint

Given a set of ICD codes (hospitalisations and deaths), OPCS codes (operations), and
UK Biobank self reported medical history field codes, along with related options in the
provided endpoint definition file, curates prevalent and/or incident disease events for
each UK Biobank participant at each assessment visit. See README.txt for more details.

Usage:
  curate_endpoint.R --def-file <file> --output <directory>
  curate_endpoint.R --def-file <file>
  curate_endpoint.R --def-file <file> --output <directory> --verbose
  curate_endpoint.R --def-file <file> --verbose
  curate_endpoint.R -h | --help

Options:
  -h --help               Show this screen.
  --def-file <file>       File containing endpoint definition, see template.txt for
                          example, and README.txt for details on options.
  --output <directory>    Directory to save the curated endpoint information. Defaults
                          to the same directory as the definition file.
  --verbose               Prints status updates.
" -> doc

# Parse arguments 
args <- docopt(doc, strict=TRUE)

if (args[["--verbose"]]) {
  message("Input arguments read in by docopt:")
  message(show(args))
}

########################################################################################
# Load and check input def-file
########################################################################################
def_file <- args[["--def-file"]]
if (!file.exists(def_file)) {
  stop(sprintf("--def-file %s: file does not exist", def_file))
}

if (args[["--verbose"]]) {
  message("Parsing endpoint definition file...")
}

opts <- readLines(def_file)
opts <- gsub("#.*$", "", opts) # strip out comments
opts <- gsub("\\s*$", "", opts) # remove trailing whitespace
opts <- opts[opts != ""] # remove empty lines

# Split into fields and codes
fields <- tolower(gsub(":.*", "", opts))
opts <- gsub(".*:\\s*", "", opts)

# Check for valid fields
valid_self_report <- c(
  "self-report field 20002", "self-report field 20004", "self-report field 2443",
	"self-report field 2463", "self-report field 4728", "self-report field 5452",
	"self-report field 5463", "self-report field 5474", "self-report field 5485",
	"self-report field 5496", "self-report field 5507", "self-report field 5518",
	"self-report field 5529", "self-report field 5540", "self-report field 6150",
	"self-report field 6152", "self-report field 20001", "self-report field 3005",
	"self-report field 6151", "self-report field 2986", "self-report field 4041",
  "self-report field 20126", "self-report field 20122", "self-report field 20127", 
  "self-report field 20124", "self-report field 20125", "self-report field 20123", 
  "self-report field 1920", "self-report field 1930", "self-report field 1940", 
  "self-report field 1950", "self-report field 1960", "self-report field 1970", 
  "self-report field 1980", "self-report field 1990", "self-report field 2000", 
  "self-report field 2010", "self-report field 2020", "self-report field 2030", 
  "self-report field 2040", "self-report field 4526", "self-report field 4537", 
  "self-report field 4548", "self-report field 4559", "self-report field 4570", 
  "self-report field 4581", "self-report field 2050", "self-report field 2060", 
  "self-report field 2070", "self-report field 2080", "self-report field 2090", 
  "self-report field 2100", "self-report field 4598", "self-report field 4609", 
  "self-report field 4620", "self-report field 4631", "self-report field 5375", 
  "self-report field 5386", "self-report field 4642", "self-report field 4653", 
  "self-report field 6156", "self-report field 5663", "self-report field 5674", 
  "self-report field 6145"
)

valid_opcs3 <- c("prevalent opcs-3", "prevalent opcs-3 primary cause only")
valid_opcs4 <- c(
  "opcs-4", "prevalent opcs-4", "prevalent opcs-4 primary cause only",
  "incident opcs-4", "incident opcs-4 primary cause only", 
  "non-fatal incident opcs-4", "non-fatal incident opcs-4 primary cause only"
)
valid_icd9 <- c("prevalent icd-9", "prevalent icd-9 primary cause only")
valid_icd10 <- c(
  "icd-10", "prevalent icd-10", "prevalent icd-10 primary cause only",
  "incident icd-10", "incident icd-10 primary cause only", 
  "non-fatal incident icd-10", "non-fatal incident icd-10 primary cause only",
  "fatal incident icd-10", "fatal incident icd-10 primary cause only"
)

valid_options <- c(
  "max follow years", "max follow date", "max follow age",
  "min follow years", "min follow date", "min follow age", "incident only where not prevalent"
)

valid_opcs3 <- c(valid_opcs3, paste("excluding", valid_opcs3))
valid_opcs4 <- c(valid_opcs4, paste("excluding", valid_opcs4))
valid_icd9 <- c(valid_icd9, paste("excluding", valid_icd9))
valid_icd10 <- c(valid_icd10, paste("excluding", valid_icd10))

valid <- c(valid_self_report, valid_opcs3, valid_opcs4, valid_icd9, valid_icd10, valid_options)

invalid <- setdiff(fields, valid)
if (length(invalid) > 0) {
  warning(sprintf("discarding invalid fields from endpoint file: %s", paste(invalid, collapse=", ")))
}

fields <- fields[which(fields %in% valid)]
opts <- opts[which(fields %in% valid)]
if (length(opts) == 0) {
  stop(sprintf("No valid fields found in --def-file %s", def_file))
}

# Flags indicating whether endpoint includes prevalent or incident definitions
def_has_prevalent <- any(fields %in% c(
 valid_self_report, valid_icd9, valid_opcs3, valid_icd10[valid_icd10 %like% "prevalent"],
 valid_opcs4[valid_opcs4 %like% "prevalent"], "icd-10", "opcs-4"
))

def_has_incident <- any(fields %in% c(
  valid_icd10[valid_icd10 %like% "incident"], valid_opcs4[valid_opcs4 %like% "incident"], 
  "icd-10", "opcs-4"
))

if (!def_has_prevalent && !def_has_incident) {
  stop("No valid endpoint definition arguments found, nothing to curate")
}

# Function for converting codes and code ranges to individual codes in UK Biobank format.
# E.g.:
# 
#  I24.2 becomes I242
#  I20-I25 becomes I20, I21, I22, I23, I24, I25
#
parse_codelist <- function(x) {
  if (length(x) == 0) return(character(0))
  codes <- strsplit(split=", ", x)[[1]]
  chapters <- codes[grepl("-", codes)]
  codes <- codes[!grepl("-", codes)]
  if (length(chapters) > 0) {
    chapters <- strsplit(split="-", chapters)
    chapters <- data.table(start=sapply(chapters, `[`, 1), end=sapply(chapters, `[`, 2))
    chapters[, start := gsub("\\.", "", start)]
    chapters[, end := gsub("\\.", "", end)]
    chapters[, chapter_letter := substr(start, 0, 1)]
    chapters[!(tolower(chapter_letter) %in% letters), chapter_letter := ""]
    chapters[, start := gsub("[a-z]", "", tolower(start))]
    chapters[, end := gsub("[a-z]", "", tolower(end))]
    chapters[, group := .I]
    chapters[nchar(start) == 1, start := paste0(start, "0")]
    chapters[nchar(end) == 1, end := paste0(end, "0")]
    chapter_codes <- chapters[, .(codes=paste0(chapter_letter, sprintf("%02d", seq(start, end)))), by = group][, codes]
    return(sort(c(codes, chapter_codes)))
  } else {
    return(codes)
  }
}

# Get vector of codes for each field. Note fields may be repeated in the def file
opts <- lapply(unique(fields), function(ff) {
  this_opts <- opts[which(fields == ff)]
  if (!(ff %in% valid_options)) {
    this_opts <- lapply(this_opts, parse_codelist)
  }
  unique(unlist(this_opts))
})
names(opts) <- unique(fields)

if (args[["--verbose"]]) {
  message("Collated and simplified endpoint definition:")
  message(show(opts))
}

# Check fields with a single allowable entry:
n <- sapply(opts[valid_options], length)
names(n) <- valid_options
bad <- n[n > 1]
if (length(bad) > 0) {
  stop(sprintf("Multiple entries/options found for fields: %s", paste(paste0("'", names(bad), ":'"), collapse=", ")))
}

# Check max/min follow years are whole numbers
if (!is.null(opts[["max follow years"]]) && as.integer(opts[["max follow years"]]) != opts[["max follow years"]]) {
  stop("'max follow years:' must be given as an integer")
}

if (!is.null(opts[["min follow years"]]) && as.integer(opts[["min follow years"]]) != opts[["min follow years"]]) {
  stop("'min follow years:' must be given as an integer")
}

# Same for age
if (!is.null(opts[["max follow age"]]) && as.integer(opts[["max follow age"]]) != opts[["max follow age"]]) {
  stop("'max follow age:' must be given as an integer")
}

if (!is.null(opts[["min follow age"]]) && as.integer(opts[["min follow age"]]) != opts[["min follow age"]]) {
  stop("'min follow age:' must be given as an integer")
}

# Determine output directory
out_dir <- args[["--output"]]
if (is.null(out_dir)) {
  out_dir <- dirname(def_file)
} else {
  out_dir <- gsub("/$", "", out_dir) # remove trailing / for directories
}
if (args[["--verbose"]] || is.null(args[["--output"]])) {
  message(sprintf("Setting output directory as %s", paste0(out_dir, "/")))
}

########################################################################################
# Find events matching --def-file criteria
########################################################################################

# Function for matching ICD or OPCS codes
# Returns TRUE for each element of x if it starts
# with any pattern in code_list. E.g. returns TRUE for each
# hospital episode (x) if the icd10 code starts with any of the
# codes in the vec list.
# 
# E.g. I20 matches I20, I20.1, I20.2 and so on.
#
`%starts_with%` <- function(x, code_list) {
  if (length(code_list) == 0) {
    return(rep(FALSE, length(x)))
  }
  code_list <- gsub("\\.", "", code_list)
  code_list <- paste0("^", code_list)
  foreach(i = seq_along(code_list), .combine=`|`) %do% {
    tolower(x) %like% tolower(code_list[i])
  }
}

if (args[["--verbose"]]) {
  message("Loading follow-up information...")
}

# Load curated follow-up information
follow <- fread(sprintf("%s/data/curated/ukbiobank/followup/followup.txt", srcDir))

# Add in inferred date of birth for handling min/max follow-up age
anthro <- fread(sprintf("%s/data/curated/ukbiobank/anthropometrics/anthropometrics.txt", srcDir))
follow[anthro, on = .(eid), approx_birth_date := i.approx_birth_date]
rm(anthro)

if (args[["--verbose"]]) {
  message("Loading death records...")
}

# Load curated death records
deaths <- fread(sprintf("%s/data/curated/ukbiobank/deaths/deaths.txt", srcDir))

if (args[["--verbose"]]) {
  message("Filtering death records...")
}

# Filter to fatal events of interest (if any)
deaths <- deaths[
  cause_icd10 %starts_with% opts[["icd-10"]] |
  cause_icd10 %starts_with% opts[["incident icd-10"]] |
  cause_icd10 %starts_with% opts[["incident icd-10 primary cause only"]] & (primary_cause) |
  cause_icd10 %starts_with% opts[["fatal incident icd-10"]] |
  cause_icd10 %starts_with% opts[["fatal incident icd-10 primary cause only"]] & (primary_cause)
]

# Drop any exclusions
deaths <- deaths[!(
  cause_icd10 %starts_with% opts[["excluding icd-10"]] |
  cause_icd10 %starts_with% opts[["excluding incident icd-10"]] |
  cause_icd10 %starts_with% opts[["excluding incident icd-10 primary cause only"]] & (primary_cause) |
  cause_icd10 %starts_with% opts[["excluding fatal incident icd-10"]] |
  cause_icd10 %starts_with% opts[["excluding fatal incident icd-10 primary cause only"]] & (primary_cause)
)]

# Duplicate across UK Biobank visits so we can combine data later
deaths <- deaths[follow[,.(eid, visit_index)], on = .(eid), nomatch=0]

if (args[["--verbose"]]) {
  message("Loading hospital records...")
}

# Load curated hospital records
hes <- fread(sprintf("%s/data/curated/ukbiobank/hospital_records/hospital_records.txt", srcDir))

# Split into prevalent and incident events
prev_hes <- hes[follow[,.(eid, visit_index, assessment_date)], 
  on = .(eid, event_date <= assessment_date), nomatch=0, 
  .(eid, visit_index, assessment_date, event_type, event_date=x.event_date, episode_end, 
    hospital_nation, code, code_type, cause_type, fatal_event, fatality_in_episode)]

inci_hes <- hes[follow[,.(eid, visit_index, assessment_date)], 
  on = .(eid, event_date > assessment_date), nomatch=0, 
  .(eid, visit_index, assessment_date, event_type, event_date=x.event_date, episode_end, 
    hospital_nation, code, code_type, cause_type, fatal_event, fatality_in_episode)]

if (args[["--verbose"]]) {
  message("Filtering hospital records...")
}

# Filter to hospital records of interest (if any)
prev_hes <- prev_hes[
  code_type == "ICD-9" & code %starts_with% opts[["prevalent icd-9"]] |
  code_type == "ICD-9" & code %starts_with% opts[["prevalent icd-9 primary cause only"]] & cause_type == "primary" |
  code_type == "ICD-10" & code %starts_with% opts[["icd-10"]] |
  code_type == "ICD-10" & code %starts_with% opts[["prevalent icd-10"]] |
  code_type == "ICD-10" & code %starts_with% opts[["prevalent icd-10 primary cause only"]] & cause_type == "primary" |
  code_type == "OPCS-3" & code %starts_with% opts[["prevalent opcs-3"]] |
  code_type == "OPCS-3" & code %starts_with% opts[["prevalent opcs-3 primary cause only"]] & cause_type == "primary" |
  code_type == "OPCS-4" & code %starts_with% opts[["opcs-4"]] |
  code_type == "OPCS-4" & code %starts_with% opts[["prevalent opcs-4"]] |
  code_type == "OPCS-4" & code %starts_with% opts[["prevalent opcs-4 primary cause only"]] & cause_type == "primary" 
]

inci_hes <- inci_hes[
  code_type == "ICD-10" & code %starts_with% opts[["icd-10"]] |
  code_type == "ICD-10" & code %starts_with% opts[["incident icd-10"]] |
  code_type == "ICD-10" & code %starts_with% opts[["incident icd-10 primary cause only"]] & cause_type == "primary" |
  code_type == "ICD-10" & code %starts_with% opts[["non-fatal incident icd-10"]] & !(fatality_in_episode) |
  code_type == "ICD-10" & code %starts_with% opts[["non-fatal incident icd-10 primary cause only"]] & cause_type == "primary" & !(fatality_in_episode) |
  code_type == "OPCS-4" & code %starts_with% opts[["opcs-4"]] |
  code_type == "OPCS-4" & code %starts_with% opts[["incident opcs-4"]] |
  code_type == "OPCS-4" & code %starts_with% opts[["incident opcs-4 primary cause only"]] & cause_type == "primary" |
  code_type == "OPCS-4" & code %starts_with% opts[["non-fatal incident opcs-4"]] & !(fatality_in_episode) |
  code_type == "OPCS-4" & code %starts_with% opts[["non-fatal incident opcs-4 primary cause only"]] & cause_type == "primary" & !(fatality_in_episode)
]

# Remove any exclusions
prev_hes <- prev_hes[!(
  code_type == "ICD-9" & code %starts_with% opts[["excluding prevalent icd-9"]] |
  code_type == "ICD-9" & code %starts_with% opts[["excluding prevalent icd-9 primary cause only"]] & cause_type == "primary" |
  code_type == "ICD-10" & code %starts_with% opts[["excluding icd-10"]] |
  code_type == "ICD-10" & code %starts_with% opts[["excluding prevalent icd-10"]] |
  code_type == "ICD-10" & code %starts_with% opts[["excluding prevalent icd-10 primary cause only"]] & cause_type == "primary" |
  code_type == "OPCS-3" & code %starts_with% opts[["excluding prevalent opcs-3"]] |
  code_type == "OPCS-3" & code %starts_with% opts[["excluding prevalent opcs-3 primary cause only"]] & cause_type == "primary" |
  code_type == "OPCS-4" & code %starts_with% opts[["excluding opcs-4"]] |
  code_type == "OPCS-4" & code %starts_with% opts[["excluding prevalent opcs-4"]] |
  code_type == "OPCS-4" & code %starts_with% opts[["excluding prevalent opcs-4 primary cause only"]] & cause_type == "primary" 
)]

inci_hes <- inci_hes[!(
  code_type == "ICD-10" & code %starts_with% opts[["excluding icd-10"]] |
  code_type == "ICD-10" & code %starts_with% opts[["excluding incident icd-10"]] |
  code_type == "ICD-10" & code %starts_with% opts[["excluding incident icd-10 primary cause only"]] & cause_type == "primary" |
  code_type == "ICD-10" & code %starts_with% opts[["excluding non-fatal incident icd-10"]] & !(fatality_in_episode) |
  code_type == "ICD-10" & code %starts_with% opts[["excluding non-fatal incident icd-10 primary cause only"]] & cause_type == "primary" & !(fatality_in_episode) |
  code_type == "OPCS-4" & code %starts_with% opts[["excluding opcs-4"]] |
  code_type == "OPCS-4" & code %starts_with% opts[["excluding incident opcs-4"]] |
  code_type == "OPCS-4" & code %starts_with% opts[["excluding incident opcs-4 primary cause only"]] & cause_type == "primary" |
  code_type == "OPCS-4" & code %starts_with% opts[["excluding non-fatal incident opcs-4"]] & !(fatality_in_episode) |
  code_type == "OPCS-4" & code %starts_with% opts[["excluding non-fatal incident opcs-4 primary cause only"]] & cause_type == "primary" & !(fatality_in_episode)
)]

# Remove big hes table to free up memory
rm(hes)
invisible(gc())

if (args[["--verbose"]]) {
  message("Loading self-report medical history...")
}

# Load curated self-reported disease events
self_report <- fread(sprintf("%s/data/curated/ukbiobank/self_report/self_report_codes.txt", srcDir))

if (args[["--verbose"]]) {
  message("Filtering self-report medical history...")
}

# Get fields + codes of interest
numeric_self_report <- c(
 "self-report field 20127",  "self-report field 4609", "self-report field 4620",
 "self-report field 5375", "self-report field 5386"
)
                        
srmap <- foreach(ff = setdiff(valid_self_report, numeric_self_report), .combine=rbind) %do% {
  if (!is.null(opts[[ff]])) { 
		data.table(field_id = as.integer(gsub("self-report field ", "", ff)), code = as.integer(opts[[ff]]))
  } else {
    data.table(field_id = integer(0), code = integer(0))
  }
}

if (any(is.na(srmap$code))) {
  bad <- srmap[is.na(code), unique(field.id)]
  stop(sprintf("events for self report field(s) %s must be given as integer codes", paste(bad, collapse=", ")))
}

# For numeric fields, codes can be given as-is, as ranges (auto-expanded to single integers above),
# or using greater than (>) or less than (<) signs to indicate cut-offs
srnumranges <- fread(sprintf("%s/data/curated/ukbiobank/self_report/code_labels.txt", srcDir))
srnummap <- foreach(ff = numeric_self_report, .combine=rbind) %do% {
  if (!is.null(opts[[ff]])) {
    codes <- opts[[ff]]
    fid <- as.integer(gsub("self-report field ", "", ff))
    int_codes <- as.integer(codes[!grepl("(<)|(>)", codes)])
    gt_codes <- as.integer(gsub(">", "", codes[grepl("^>", codes)]))
    lt_codes <- as.integer(gsub("<", "", codes[grepl("^<", codes)]))
    fid_range <- srnumranges[field_id == fid & code >= 0]
    # Expand thresholds to cover full range of possibilities
    if (length(gt_codes) > 0) {
      gt_codes <- foreach(gt = gt_codes, .combine=c) %do% {
        fid_range[code > gt, code]
      }
    } 
    if (length(lt_codes) > 0) {
      lt_codes <- foreach(lt = lt_codes, .combine=c) %do% {
        fid_range[code < lt, code]
      }
    } 
    data.table(field_id = fid, code = c(gt_codes, lt_codes, int_codes))
  } else {
    data.table(field_id = integer(0), code = integer(0))
  }
}

if (any(is.na(srnummap$code))) {
  bad <- srnummap[is.na(code), unique(field.id)]
  stop(sprintf("events for self report field(s) %s must be given as integer codes, ranges, or < or > thresholds", paste(bad, collapse=", ")))
}

# Combine
srmap <- rbind(srmap, srnummap)

# Drop any duplicates that may arise
srmap <- unique(srmap)

# Get people with missing data answers (neither TRUE nor FALSE for a given condition)
# and thus cannot be used as controls.
#
# code == -1: "Do not know"
# code == -2: "Not applicable" (i.e. question not asked due to answer to another question)
# code == -3: "Prefer not to answer"
# code == -7 is "None of the above", which is basically code == 0 (FALSE) for multi-answer questions
#
missing <- self_report[(code == -1 | code == -3) & field_id %in% unique(srmap$field_id), .(field_id, eid, visit_index)]

# In these instances, in some cases the participant has data in another field, e.g. if a person
# Answers "Don't know" to the diabetes question (field 2443), this may actually be determined 
# in the verbal interview question (field 20002), in which case we don't want to introduce missing
# data for people we can otherwise determine have FALSE based on their answers to field 20002.
if (2443 %in% srmap$field_id && 20002 %in% srmap$field_id) {
  non_missing_20002 <- unique(self_report[field_id == 20002 & !(code == -1 | code == -3), .(eid, visit_index)])
  missing_2443_to_drop <- missing[field_id == 2443][non_missing_20002, on = .(eid, visit_index), nomatch=0]
  missing <- missing[!missing_2443_to_drop, on = .(field_id, eid, visit_index)]
}

if (6150 %in% srmap$field_id && 20002 %in% srmap$field_id) {
  non_missing_20002 <- unique(self_report[field_id == 20002 & !(code == -1 | code == -3), .(eid, visit_index)])
  missing_6150_to_drop <- missing[field_id == 6150][non_missing_20002, on = .(eid, visit_index), nomatch=0]
  missing <- missing[!missing_6150_to_drop, on = .(field_id, eid, visit_index)]
}

if (6152 %in% srmap$field_id && 20002 %in% srmap$field_id) {
  non_missing_20002 <- unique(self_report[field_id == 20002 & !(code == -1 | code == -3), .(eid, visit_index)])
  missing_6152_to_drop <- missing[field_id == 6152][non_missing_20002, on = .(eid, visit_index), nomatch=0]
  missing <- missing[!missing_6152_to_drop, on = .(field_id, eid, visit_index)]
}

if (2463 %in% srmap$field_id && 20004 %in% srmap$field_id) {
  non_missing_20004 <- unique(self_report[field_id == 20004 & !(code == -1 | code == -3), .(eid, visit_index)])
  missing_2463_to_drop <- missing[field_id == 2463][non_missing_20004, on = .(eid, visit_index), nomatch=0]
  missing <- missing[!missing_2463_to_drop, on = .(field_id, eid, visit_index)]
}

# Now we can filter to unique instances
missing <- unique(missing[,.(eid, visit_index)])

# Add in people not asked any of the self-report question - this doesn;t conflict with the
# above because no one who was asked 20002/20004 was not first asked one of the preceeding
# touchscreen questions listed above
if (length(srmap$field_id) > 0) {
	asked <- unique(self_report[field_id %in% unique(srmap$field_id), .(eid, visit_index)])
	not_asked <- follow[!asked, on = .(eid, visit_index), .(eid, visit_index)]
	missing <- rbind(missing, not_asked)
  missing <- unique(missing)
}

# Filter self report table
self_report <- self_report[srmap, on = .(field_id, code), nomatch=0]

########################################################################################
# Curate label tabel
########################################################################################
if (args[["--verbose"]]) {
  message("Loading and curating code labels...")
}

labels <- rbind(fill=TRUE,
  fread(sprintf("%s/data/curated/ukbiobank/self_report/code_labels.txt", srcDir))[,.(code_type=paste("self-report field", field_id), ukb_code=code, code, label)],
	rbind(idcol="code_type",
		"ICD-10"=fread(sprintf("%s/data/curated/ukbiobank/hospital_records/ICD10_codes.tsv", srcDir))[,.(ukb_code=coding, label=meaning)],
		"ICD-9"=fread(sprintf("%s/data/curated/ukbiobank/hospital_records/ICD9_codes.tsv", srcDir))[,.(ukb_code=coding, label=meaning)],
		"OPCS-4"=fread(sprintf("%s/data/curated/ukbiobank/hospital_records/OPCS4_codes.tsv", srcDir))[,.(ukb_code=coding, label=meaning)],
		"OPCS-3"=fread(sprintf("%s/data/curated/ukbiobank/hospital_records/OPCS3_codes.tsv", srcDir))[!(coding %like% "Chapter"),.(ukb_code=coding, label=meaning)]
	)
)
labels[, code := as.character(code)]

# UKB codes don't have periods in the appropriate places, fix.
labels[code_type == "ICD-10" & nchar(ukb_code) <= 3, code := ukb_code]
labels[code_type == "ICD-10" & nchar(ukb_code) > 3, 
  code := sprintf("%s.%s", substr(ukb_code, 1, 3), substr(ukb_code, 4, nchar(ukb_code)))]

labels[code_type == "ICD-9" & nchar(ukb_code) == 3, code := ukb_code]
labels[code_type == "ICD-9" & grepl("[A-Z]", ukb_code) & nchar(ukb_code) == 4, code := ukb_code]
labels[code_type == "ICD-9" & nchar(ukb_code) > 3 & !grepl("[A-Z]", ukb_code),
  code := sprintf("%s.%s", substr(ukb_code, 1, 3), substr(ukb_code, 4, nchar(ukb_code)))]
labels[code_type == "ICD-9" & nchar(ukb_code) > 4 & grepl("[A-Z]", ukb_code),
  code := sprintf("%s.%s", substr(ukb_code, 1, 4), substr(ukb_code, 5, nchar(ukb_code)))]

labels[code_type == "OPCS-4" & nchar(ukb_code) <= 3, code := ukb_code]
labels[code_type == "OPCS-4" & nchar(ukb_code) > 3, 
  code := sprintf("%s.%s", substr(ukb_code, 1, 3), substr(ukb_code, 4, nchar(ukb_code)))]

labels[code_type == "OPCS-3" & nchar(ukb_code) <= 3, code := ukb_code]
labels[code_type == "OPCS-3" & nchar(ukb_code) > 3, 
  code := sprintf("%s.%s", substr(ukb_code, 1, 3), substr(ukb_code, 4, nchar(ukb_code)))]

# Drop code from label
labels[code_type == "ICD-10", label := gsub("[A-Z][0-9][0-9]\\.?[0-9]* ", "", label)]
labels[code_type == "ICD-9", label := gsub("[A-Z]?[0-9][0-9][0-9]\\.?[0-9]* ", "", label)]
labels[code_type == "OPCS-4", label := gsub("[A-Z][0-9][0-9]\\.?[0-9]* ", "", label)]
labels[code_type == "OPCS-3", label := gsub("[0-9][0-9][0-9]\\.?[0-9]* ", "", label)]

########################################################################################
# Get first occurrences
########################################################################################
if (args[["--verbose"]]) {
  message("Getting first occurences of incident events...")
}

# Across all incident events, get the one that occurs first.
# Due to row ordering, the primary cause will always be selected over the secondary and so on
# if it is one of the codes of interest.
if (nrow(deaths) > 0 && nrow(inci_hes) > 0) {
	inci <- rbind(
		deaths[,.(eid, visit_index, event_type="death", record_source=death_source, code_type="ICD-10", code=cause_icd10,
							cause_type=ifelse(primary_cause, "primary", "secondary"), fatal_event=TRUE, event_date=date_of_death)],
		inci_hes[, .(eid, visit_index, event_type, record_source=paste("NHS records from hospitals in", hospital_nation),
								 code_type, code, cause_type, fatal_event, event_date)]
	)
} else if (nrow(deaths) > 0) {
  inci <- deaths[,.(eid, visit_index, event_type="death", record_source=death_source, code_type="ICD-10", code=cause_icd10,
                    cause_type=ifelse(primary_cause, "primary", "secondary"), fatal_event=TRUE, event_date=date_of_death)]
} else if (nrow(inci_hes) > 0) {
  inci <- inci_hes[, .(eid, visit_index, event_type, record_source=paste("NHS records from hospitals in", hospital_nation),
                       code_type, code, cause_type, fatal_event, event_date)]
} else {
  inci <- inci_hes[, .(eid, visit_index, event_type, record_source=character(0), code_type, code, cause_type, fatal_event, event_date)]
}
inci <- inci[, .SD[which.min(event_date)], by=.(eid, visit_index)]

# Update codes and add labels in incident table
inci[labels, on = .(code_type, code=ukb_code), c("code", "code_label") := .(i.code, i.label)]

# Rename columns to distinguish incident and prevalent events
setnames(inci, paste0("incident_", names(inci)))
setnames(inci, c("incident_eid", "incident_visit_index"), c("eid", "visit_index"))

if (args[["--verbose"]]) {
  message("Getting most recent occurences of prevalent events...")
}

# Across all prevalent events, get the one that occurs latest. This is a bit trickier than incident events
# because there's missing data in the self-report, where either age/date of event could not be inferred, or
# where a participant selected "don't know" or "prefer not to answer" for one or more fields.
if (nrow(prev_hes) > 0 && nrow(self_report) > 0) {
	prev <- rbind(
		prev_hes[, .(eid, visit_index, event_type, record_source=paste("NHS records from hospitals in", hospital_nation),
								 code_type, code, cause_type, event_date)],
		self_report[, .(eid, visit_index, event_type="self-reported", 
										record_source=ifelse(field_id %in% c(20001, 20002, 20004), "Verbal interview with nurse", "Touchscreen survey"),
										code_type=paste("self-report field", field_id), code, cause_type="self-reported", event_date=date)]
	)
} else if (nrow(prev_hes) > 0) {
  prev <- prev_hes[, .(eid, visit_index, event_type, record_source=paste("NHS records from hospitals in", hospital_nation),
                      code_type, code, cause_type, event_date)]
} else if (nrow(self_report) > 0) {
  prev <- self_report[, .(eid, visit_index, event_type="self-reported",
			 									  record_source=ifelse(field_id %in% c(20001, 20002, 20004), "Verbal interview with nurse", "Touchscreen survey"),
													code_type=paste("self-report field", field_id), code, cause_type="self-reported", event_date=date)]
} else {
  prev <- prev_hes[,.(eid, visit_index, event_type, record_source=character(0), code_type, code, cause_type, event_date)]
}

# First extract events with no event date attached, picking one representative event per participant/assessment visit
prev_nodate <- prev[is.na(event_date), .SD[1], by=.(eid, visit_index)] # ordered by priority in README.txt

# Next, get the most recent event for each participant/assessment among prevalent events with dates
prev_wdate <- prev[!is.na(event_date),.SD[which.max(event_date)], by=.(eid, visit_index)]

# latest event cannot be inferred where some events have missing dates, so self-report without event date will ultimately
# take priority in output
prev_wdate <- prev_wdate[!prev_nodate, on = .(eid, visit_index)]

# Recombine
prev <- rbind(prev_nodate, prev_wdate)

# Update codes and add labels in incident table
prev[labels, on = .(code_type, code=ukb_code), c("code", "code_label") := .(i.code, i.label)]

# Rename columns to distinguish incident and prevalent events
setnames(prev, paste0("prevalent_", names(prev)))
setnames(prev, c("prevalent_eid", "prevalent_visit_index"), c("eid", "visit_index"))

########################################################################################
# Merge in first occurence information to follow-up information table
########################################################################################

if (args[["--verbose"]]) {
  message("Merging with follow-up information...")
}

# Add incident events 
follow <- merge(follow, inci, by=c("eid", "visit_index"), all.x=TRUE)

# If endpoint includes hospital records and deaths, drop any events whose first event is a
# death occurring after maximum follow-up in hospital records
if (!is.null(opts[["icd-10"]]) || 
    !is.null(opts[["incident icd-10"]]) || 
    !is.null(opts[["incident icd-10 primary cause only"]]) || 
    !is.null(opts[["non-fatal incident icd-10"]]) || 
    !is.null(opts[["non-fatal incident icd-10 primary cause only"]]) || 
    !is.null(opts[["opcs-4"]]) || 
    !is.null(opts[["incident opcs-4"]]) || 
    !is.null(opts[["incident opcs-4 primary cause only"]]) || 
    !is.null(opts[["non-fatal incident opcs-4"]]) || 
    !is.null(opts[["non-fatal incident opcs-4 primary cause only"]])) {
  follow[incident_event_date > latest_hospital_date, 
		c("incident_event_type", "incident_record_source", "incident_code_type",
			"incident_code", "incident_cause_type", "incident_fatal_event",
			"incident_event_date", "incident_code_label") := 
    .(NA, NA, NA, NA, NA, NA, NA, NA)]
}

# Add prevalent events
follow <- merge(follow, prev, by=c("eid", "visit_index"), all.x=TRUE)

########################################################################################
# Truncate follow-up if requested
########################################################################################

if (args[["--verbose"]]) {
  message("Applying any truncations to follow-up time / date / age...")
}

# Function to compute years between two dates - note this preserves
# human notions of whole years, i.e:
#
#   2019-07-25 - 2009-07-25 = 10 years                                                                                                                                                                                                       #
# At the expense of (potentially) leading to small inaccuracies in                                                                                                                                                                           # rank ordering due to leap days, i.e. in pure terms of days,
#
#   2019-07-25 - 2009-07-25 = 9.9986 years                                                                                                                                                                                                   #
# As this makes it harder to accurately truncate follow-up (e.g.
# for testing model calibration)
years_between <- function(d1, d2) {
  year_diff <- as.period(interval(as.Date(d1), as.Date(d2))) / years(1)

  # When d1 or d2 is a feb 29th, we want to return a whole number for
  # number of years between d1 and d2 if the other date is (1) Feb 28th,
  # and (2) not also a leap year.
  leap_day_d1 <- which(month(d1) == 2 & day(d1) == 29)
  leap_day_d2 <- which(month(d2) == 2 & day(d1) == 29)
  feb_28_d1 <- which(month(d1) == 2 & day(d2) == 28 & !(leap_year(year(d1))))
  feb_28_d2 <- which(month(d2) == 2 & day(d2) == 28 & !(leap_year(year(d2))))

  to_correct <- union(intersect(leap_day_d1, feb_28_d2), intersect(feb_28_d1, leap_day_d2))
  year_diff[to_correct] <- year(d2[to_correct]) - year(d1[to_correct])
  return(year_diff)
}

# Likewise, add_years is consistent with the above, i.e.
# 
# 2009-07-25 + 10 years = 2019-07-25
#
# Note only works with whole years.
#
add_years <- function(d1, follow) {
  new_date <- as.IDate(as.Date(d1) + years(follow))
	# If the date input is Feb 29, and the corresponding year after
	# adding 'follow' is not a leap year, roll back to Feb 28.
  leap_day <- which(is.na(new_date) & month(d1) == 2 & day(d1) == 29)
  new_date[leap_day] <- as.IDate(as.Date(d1[leap_day]) - days(1) + years(follow)) 
  return(new_date)
}

# Apply any truncations to columns we will extract
if (!is.null(opts[["max follow years"]])) {
  max_years <- as.integer(opts[["max follow years"]])

  if (max_years <= 0) stop("'max follow years:' must be > 0 if provided")

  follow[latest_mortality_followup > max_years, 
    c("latest_mortality_followup", "latest_mortality_date", "all_cause_mortality") :=
    .(max_years, add_years(assessment_date, max_years), FALSE)]

  follow[latest_hospital_followup > max_years, 
    c("latest_hospital_followup", "latest_hospital_date") :=
    .(max_years, add_years(assessment_date, max_years))]

  follow[incident_event_date > add_years(assessment_date, max_years),
			c("incident_event_type", "incident_record_source", "incident_code_type",
				"incident_code", "incident_cause_type", "incident_fatal_event",
				"incident_event_date", "incident_code_label") :=
			.(NA, NA, NA, NA, NA, NA, NA, NA)]
}

if (!is.null(opts[["max follow date"]])) {
  max_date <- as.IDate(opts[["max follow date"]])
  if (as.character(max_date) != opts[["max follow date"]]) {
    warning(sprintf("Coerced 'max follow date:' from '%s' to '%s'. To avoid this warning use format YYYY-MM-DD", opts[["max follow date"]], max_date))
  }

  follow[latest_mortality_date > max_date, 
    c("latest_mortality_followup", "latest_mortality_date", "all_cause_mortality") :=
    .(years_between(assessment_date, max_date), max_date, FALSE)]

  follow[latest_hospital_date > max_date, 
    c("latest_hospital_followup", "latest_hospital_date") :=
    .(years_between(assessment_date, max_date), max_date)]
 
  follow[incident_event_date > max_date,
      c("incident_event_type", "incident_record_source", "incident_code_type",
        "incident_code", "incident_cause_type", "incident_fatal_event",
        "incident_event_date", "incident_code_label") :=
      .(NA, NA, NA, NA, NA, NA, NA, NA)]
}

if (!is.null(opts[["max follow age"]])) {
  max_age <- as.integer(opts[["max follow age"]])

  if (max_age <= 0) stop("'max follow age:' must be > 0 if provided")

  follow[latest_mortality_date > add_years(approx_birth_date, max_age),
    c("latest_mortality_followup", "latest_mortality_date", "all_cause_mortality") :=
    .(years_between(assessment_date, add_years(approx_birth_date, max_age)), add_years(approx_birth_date, max_age), FALSE)]

  follow[latest_hospital_date > add_years(approx_birth_date, max_age),
    c("latest_hospital_followup", "latest_hospital_date") :=
    .(years_between(assessment_date, add_years(approx_birth_date, max_age)), add_years(approx_birth_date, max_age))]

  follow[incident_event_date > add_years(approx_birth_date, max_age),
      c("incident_event_type", "incident_record_source", "incident_code_type",
        "incident_code", "incident_cause_type", "incident_fatal_event",
        "incident_event_date", "incident_code_label") :=
      .(NA, NA, NA, NA, NA, NA, NA, NA)]
}

if (!is.null(opts[["min follow years"]])) {
  min_years <- as.numeric(opts[["min follow years"]])

  if (min_years >= 0) stop("'min follow years:' must be < 0 if provided")

  follow[earliest_hospital_followup < min_years,
    c("earliest_hospital_followup", "earliest_hospital_date") :=
    .(min_years, add_years(assessment_date, min_years))]

  follow[prevalent_event_date > add_years(assessment_date, min_years),
      c("prevalent_event_type", "prevalent_record_source", "prevalent_code_type",
        "prevalent_code", "prevalent_cause_type", "prevalent_event_date", 
        "prevalent_code_label") := .(NA, NA, NA, NA, NA, NA, NA)]
}

if (!is.null(opts[["min follow date"]])) {
  min_date <- as.IDate(opts[["min follow date"]])
  if (as.character(min_date) != opts[["min follow date"]]) {
    warning(sprintf("Coerced 'min follow date:' from '%s' to '%s'. To avoid this warning use format YYYY-MM-DD", opts[["min follow date"]], min_date))
  }

  follow[earliest_hospital_date < min_date,
    c("earliest_hospital_followup", "earliest_hospital_date") :=
    .(years_between(assessment_date, min_date), min_date)]

  follow[prevalent_event_date > min_date,
      c("prevalent_event_type", "prevalent_record_source", "prevalent_code_type",
        "prevalent_code", "prevalent_cause_type", "prevalent_event_date", 
        "prevalent_code_label") := .(NA, NA, NA, NA, NA, NA, NA)]
}

if (!is.null(opts[["min follow age"]])) {
  min_age <- as.integer(opts[["min follow age"]])

  if (min_age <= 0) stop("'min follow age:' must be > 0 if provided")

  follow[earliest_hospital_date > add_years(approx_birth_date, min_age),
    c("earliest_hospital_followup", "earliest_hospital_date") :=
    .(years_between(assessment_date, add_years(approx_birth_date, min_age)), add_years(approx_birth_date, min_age))]

  follow[prevalent_event_date < add_years(approx_birth_date, min_age),
      c("prevalent_event_type", "prevalent_record_source", "prevalent_code_type",
        "prevalent_code", "prevalent_cause_type", "prevalent_event_date", 
        "prevalent_code_label") := .(NA, NA, NA, NA, NA, NA, NA)]
}

########################################################################################
# Curate prevalent disease information for survival analysis
########################################################################################

if (args[["--verbose"]]) {
  message("Curating retrospective follow-up time for prevalent events...")
}

follow[, prevalent_event := FALSE]

# Missing data in requested self-report answers become NA (if no other case definition met)
follow[missing, on = .(eid, visit_index), prevalent_event := NA]

# If prevalent disease definition includes hospital records set to 
# missing any person with withdrawn consent for record linkage 
# who does not have a prevalent event in the self-report data
if (!is.null(opts[["icd-10"]]) || 
    !is.null(opts[["prevalent icd-10"]]) || 
    !is.null(opts[["prevalent icd-10 primary cause only"]]) || 
    !is.null(opts[["prevalent icd-9"]]) || 
    !is.null(opts[["prevalent icd-9 primary cause only"]]) || 
    !is.null(opts[["opcs-4"]]) || 
    !is.null(opts[["prevalent opcs-4"]]) || 
    !is.null(opts[["prevalent opcs-4 primary cause only"]]) || 
    !is.null(opts[["prevalent opcs-3"]]) || 
    !is.null(opts[["prevalent opcs-3 primary cause only"]])) {
  follow[is.na(any_hospitalisations), prevalent_event := NA]
}

# Now that NAs have been handled, set any rows with prevalent events with known
# date to TRUE
follow[!is.na(prevalent_code), prevalent_event := TRUE]

# If minimum cut-offs on age, date, or year set, and the prevalent event is one
# with unknown date, set these to missing
if (!is.null(opts[["min follow years"]]) ||
    !is.null(opts[["min follow date"]]) ||
    !is.null(opts[["min follow age"]])) {
  follow[(prevalent_event) & is.na(prevalent_event_date), 
    c("prevalent_event_type", "prevalent_record_source", "prevalent_code_type", "prevalent_code", 
      "prevalent_cause_type", "prevalent_event_date", "prevalent_code_label", "prevalent_event") := NA]
} 

# Follow-up time is a bit tricky here, because hospital records don't go
# all the way back until birth, while self-report data does. Here we:
#
# 1. Set minimum follow-up time to (approximate) date of birth (exact date
#    is not made available by UKB) or to the date of the most recent event
# 2. Check if prevalent disease case definition includes hospital records,
#    if so, then minimum follow-up time is set to the minimum date in the
#    hospital records. 
#
follow[, prevalent_event_with_followup := prevalent_event]
follow[(prevalent_event) & is.na(prevalent_event_date), prevalent_event_with_followup := NA]
if (!is.null(opts[["icd-10"]]) || 
    !is.null(opts[["prevalent icd-10"]]) || 
    !is.null(opts[["prevalent icd-10 primary cause only"]]) || 
    !is.null(opts[["prevalent icd-9"]]) || 
    !is.null(opts[["prevalent icd-9 primary cause only"]]) || 
    !is.null(opts[["opcs-4"]]) || 
    !is.null(opts[["prevalent opcs-4"]]) || 
    !is.null(opts[["prevalent opcs-4 primary cause only"]]) || 
    !is.null(opts[["prevalent opcs-3"]]) || 
    !is.null(opts[["prevalent opcs-3 primary cause only"]])) {
  follow[prevalent_event_date < earliest_hospital_date, prevalent_event_with_followup := FALSE]
  follow[is.na(any_hospitalisations), prevalent_event_with_followup := FALSE]
  follow[!(prevalent_event_with_followup), prevalent_event_followup_date := earliest_hospital_date]
  follow[(prevalent_event_with_followup), prevalent_event_followup_date := prevalent_event_date]
} else {
  follow[!(prevalent_event), prevalent_event_followup_date := approx_birth_date]
  follow[(prevalent_event), prevalent_event_followup_date := prevalent_event_date]
}

# Compute retrospective follow-up time based on date of event or minimum follow-up
follow[, prevalent_event_followup := years_between(assessment_date, prevalent_event_followup_date)]

# Round follow-up for self-reported events. For these, each person responded with a
# calendar year, or age (whole number) at which the event occurred, from which a mid-point
# date was then determined by UK Biobank, so the precision of follow-up is at best dependent
# on the precision of the decimal age we have computed
follow[prevalent_event_type == "self-reported" & !is.na(prevalent_event_followup),
       prevalent_event_followup := round(prevalent_event_followup, digits=2)]

# Compute age at minimum retrospective follow-up (including age at event for events
# that fall within retrospective follow-up). Rounding is to reflect approximate nature
# of fractional age (accurate to within ~0.05 of a year, with maximum error of 16 days).
follow[, prevalent_event_followup_age := round(age_decimal + prevalent_event_followup, digits=2)]

# Compute age at event (includes events that occur prior to the minimum retrospective follow-up)
follow[, prevalent_event_age := age_decimal + years_between(assessment_date, prevalent_event_date)]
follow[, prevalent_event_age := round(prevalent_event_age, digits=2)]

########################################################################################
# Curate incident disease information for survival analysis
########################################################################################

if (args[["--verbose"]]) {
  message("Curating follow-up time for incident events...")
}

follow[, incident_event := FALSE]
follow[!is.na(incident_event_date), incident_event := TRUE]

# Set to missing any person with withdrawn consent for linkage to hospital/death records
follow[is.na(any_hospitalisations), incident_event := NA]

# Drop incident events where there is a prevalent event if requested
if (!is.null(opts[["incident only where not prevalent"]]) &&
     tolower(opts[["incident only where not prevalent"]]) == "true") {  
  follow[is.na(prevalent_event) | (prevalent_event),
      c("incident_event_type", "incident_record_source", "incident_code_type",
        "incident_code", "incident_cause_type", "incident_fatal_event",
        "incident_event_date", "incident_code_label", "incident_event") :=
      .(NA, NA, NA, NA, NA, NA, NA, NA, NA)]
}

# If endpoint includes a specific explicit request for fatal events, then
# we need to truncate maximum follow-up to maximum date in the death records,
# which is currently earlier than the maximum date available in the hospital 
# records for England and Scotland.
if (!is.null(opts[["fatal incident icd-10"]]) ||
    !is.null(opts[["fatal incident icd-10 primary cause only"]])) {
  follow[incident_event_date > latest_mortality_date, 
    c("incident_event_type", "incident_record_source", "incident_code_type",
      "incident_code", "incident_cause_type", "incident_fatal_event",
      "incident_event_date", "incident_code_label", "incident_event") :=
    .(NA, NA, NA, NA, NA, NA, NA, NA, FALSE)]
}

# If endpoint includes hospital records and deaths, drop any events whose first event is a
# death occurring after maximum follow-up in hospital records (currently only possible in Wales)
if (!is.null(opts[["icd-10"]]) || 
    !is.null(opts[["incident icd-10"]]) || 
    !is.null(opts[["incident icd-10 primary cause only"]]) || 
    !is.null(opts[["non-fatal incident icd-10"]]) || 
    !is.null(opts[["non-fatal incident icd-10 primary cause only"]]) || 
    !is.null(opts[["opcs-4"]]) || 
    !is.null(opts[["incident opcs-4"]]) || 
    !is.null(opts[["incident opcs-4 primary cause only"]]) || 
    !is.null(opts[["non-fatal incident opcs-4"]]) || 
    !is.null(opts[["non-fatal incident opcs-4 primary cause only"]])) {
  follow[incident_event_date > latest_hospital_date, 
		c("incident_event_type", "incident_record_source", "incident_code_type",
			"incident_code", "incident_cause_type", "incident_fatal_event",
			"incident_event_date", "incident_code_label", "incident_event") := 
    .(NA, NA, NA, NA, NA, NA, NA, NA, FALSE)]
}

# Follow-up date and time is either based on date of first event, or if
# no first event, the latest date in the hospital or death records (the
# latter if only fatal events are of interest).
follow[(incident_event), incident_event_followup_date := incident_event_date]
if (!is.null(opts[["icd-10"]]) || 
    !is.null(opts[["incident icd-10"]]) || 
    !is.null(opts[["incident icd-10 primary cause only"]]) || 
    !is.null(opts[["non-fatal incident icd-10"]]) || 
    !is.null(opts[["non-fatal incident icd-10 primary cause only"]]) || 
    !is.null(opts[["opcs-4"]]) || 
    !is.null(opts[["incident opcs-4"]]) || 
    !is.null(opts[["incident opcs-4 primary cause only"]]) || 
    !is.null(opts[["non-fatal incident opcs-4"]]) || 
    !is.null(opts[["non-fatal incident opcs-4 primary cause only"]])) {
  follow[!(incident_event), incident_event_followup_date := latest_hospital_date]
} else if (!is.null(opts[["fatal incident icd-10"]]) ||
           !is.null(opts[["fatal incident icd-10 primary cause only"]])) {
  follow[!(incident_event), incident_event_followup_date := latest_mortality_date]
}

follow[, incident_event_followup := years_between(assessment_date, incident_event_followup_date)]

# Set FALSE in the all_cause_mortality column any deaths that happened after the
# incident event
follow[(all_cause_mortality) & latest_mortality_date > incident_event_date, all_cause_mortality := FALSE]

# Compute age at follow-up, rounding to reflect approximate nature of age at assessment
# (accurate to within ~0.05 of a year, with maximum error of 16 days).
follow[, incident_event_followup_age := round(age_decimal + incident_event_followup, digits=2)]

########################################################################################
# Write out
########################################################################################

if (args[["--verbose"]]) {
  message("Writing out results...")
}

# Select and organise columns
if (def_has_prevalent && def_has_incident) {
	follow <- follow[, .(eid, sex, visit_index, assessment_date, assessment_centre, assessment_nation, 
		age=age_decimal, ehr_linkage_withdrawn = is.na(any_hospitalisations), earliest_hospital_date,
		earliest_hospital_nation, prevalent_event, prevalent_event_with_followup, prevalent_event_followup, 
		prevalent_event_followup_date, prevalent_event_followup_age, prevalent_event_date, prevalent_event_age,
    prevalent_event_type, prevalent_record_source, prevalent_cause_type, prevalent_code_type, 
    prevalent_code, prevalent_code_label, lost_to_followup, lost_to_followup_reason, lost_to_followup_date, 
    latest_hospital_date, latest_hospital_nation, latest_mortality_date, incident_event, 
    incident_event_followup, incident_event_followup_date, incident_event_followup_age,
    mortality_at_followup_date = all_cause_mortality, incident_event_type, incident_record_source, 
    incident_cause_type, incident_code_type, incident_code, incident_code_label)]
} else if (!def_has_prevalent) {
	follow <- follow[, .(eid, sex, visit_index, assessment_date, assessment_centre, assessment_nation, 
		age=age_decimal, ehr_linkage_withdrawn = is.na(any_hospitalisations), lost_to_followup, 
    lost_to_followup_reason, lost_to_followup_date, latest_hospital_date, latest_hospital_nation, 
    latest_mortality_date, incident_event, incident_event_followup, incident_event_followup_date, 
    incident_event_followup_age, mortality_at_followup_date = all_cause_mortality, incident_event_type, 
    incident_record_source, incident_cause_type, incident_code_type, incident_code, incident_code_label)]
} else if (!def_has_incident) {
	follow <- follow[, .(eid, sex, visit_index, assessment_date, assessment_centre, assessment_nation, 
		age=age_decimal, ehr_linkage_withdrawn = is.na(any_hospitalisations), earliest_hospital_date,
		earliest_hospital_nation, prevalent_event, prevalent_event_with_followup, prevalent_event_followup, 
		prevalent_event_followup_date, prevalent_event_followup_age, prevalent_event_date, prevalent_event_age,
    prevalent_event_type, prevalent_record_source, prevalent_cause_type, prevalent_code_type, 
    prevalent_code, prevalent_code_label)]
}

# Write out
fwrite(follow, sep="\t", quote=FALSE, file=sprintf("%s/events_and_followup.txt", out_dir))
if (dirname(def_file) != out_dir) {
	system(sprintf("cp %s %s/endpoint_definition.txt", def_file, out_dir), wait=TRUE)
}

if (args[["--verbose"]]) {
  message("Completed successfully!")
}

