########################################################################################
# Load R package dependencies
########################################################################################

# Set the R package location diroctory to one with the R packages pre-installed.
srcDir = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/GRS_resources/"
SoftwareDir = sprintf("%s/software_dependencies/", srcDir)
RpackageDir = sprintf("%s/Rpackages/", SoftwareDir)
RpackageDir = sprintf("%s/%s.%s/", RpackageDir, R.version$major, gsub("\\..*", "", R.version$minor))
.libPaths(RpackageDir)

suppressMessages(library("data.table"))
suppressMessages(library("foreach"))
suppressMessages(library("docopt"))
suppressMessages(library("bit64"))

########################################################################################
# Parse user input and sanity check
########################################################################################

# Parse input arguments
"Calculate the levels of a polygenic score in a group of samples

Note this script is intended to be run within a compute node on as part of an array job
submitted by calc_PS_lvls.sh. 

Usage:
  calc_PS_lvls.R --score-file <file> --work <directory> --mem <MB> [options]
  calc_PS_lvls.R -h | --help

Options:
  -h --help                   Show this screen.
  --score-file <file>         Path to polygenic score file, directory, or file containing list
                              of score files (see --type).
  --work <directory>          Working directory to store intermediate files and logs shared across all
                              scores for the duration of the run. 
  --mem <MB>                  Amount of memory to restrict plink to use.
  --type <type>               Type of filepath given to --score-file. 's': filepath points to
                              a single polygenic score, e.g. created by make_simple_PS.sh or
                              downloaded from the PGS catalog. 'd': filepath is a directory
                              containing multiple score files. 'l': filepath points to a file
                              containing a list of paths, one per line. [default: s]
  --score-rsid <col>          Name or number of the column in the polygenic score file corresponding to
                              the variant marker id to pass to plink. If column is not present, set to
                              'NULL'. [default: rsid]
  --score-chr <col>           Name or number of the column in the polygenic score file corresponding to
                              the variant's chromosome. If column is not present, set to 'NULL'.
                              If the score contains non-autosomal variants, then the chromosome field 
                              *must* contain one of X, Y, XY, or MT. [default: chr]
  --score-pos <col>           Name or number of the column in the polygenic score file corresponding to
                              the variant's position. If column is not present, set to 'NULL'. [default: pos]
  --score-EA <col>            Name or number of the column in the polygenic score file corresponding to
                              the variant's effect allele. [default: effect_allele]
  --score-EAF <col>           Name or number of the column in the polygenic score file corresponding to
                              the effect allele frequency. If column is not present, set to 'NULL'. [default: NULL]
  --score-OA <col>            Name or number of the column in the polygenic score file corresponding to
                              the variant's non-effect allele. If column is not present, set to 'NULL'. 
                              [default: other_allele]
  --score-weight <col>        Name or number of the column in the polygenic score file corresponding to
                              the effect allele's weight in the polygenic score. If the score file has
                              multiple weight columns for multiple scores, set this to 'm' to calculate
                              the levels of all these scores. In this case, the program will assume that
                              all columns not listed in the arguments above are weights columns. [default: weight]
  --score-dominant <col>      Column of TRUE/FALSE values indicating whether the effect for each variant
                              should be considered dominant (i.e. the weight is multiplied by the effect 
                              allele presence/absence rather than by the number of copies). [default: NULL]
  --score-recessive <col>     Column of TRUE/FALSE values indicating whether the effect for each variant
                              should be considered recessive (i.e. the weight is counted only when there 
                              are two copies of the effect allele). [default: NULL]
  --match-by-rsid             Flag to indicate that the score file should be matched to the
                              genotype data by the rsid column instead of by chromosome and
                              position when all three columns are provided.
  --cohort-name <name>        Name of the group of samples, used to name the output file name
                              and folder (if --out not provided). [default: UKBv3]
  --out <directory>           Directory to store the results in. By default, the results are stored
                              in a folder named <--cohort-name>_sample_levels/ in the same directory
                              as each input --score-file. If multiple score files are detected (files
                              ending in .txt.gz) in a score file's directory then a further sub-folder
                              is created <--cohort-name>_sample_levels/<score_name>/ for the score
                              being calculated by this program. [default: NULL]
  --single-out <name>         If provided, saves only a single file containing the levels of all scores
                              being calculated in a file with the given name (i.e. <--out>/<--single-out>.sscore.gz.
                              The default behaviour is to otherwise save the levels of each score into
                              separate files. [default: NULL]
  --genotype-prefix <prefix>  Path and prefix occurring before the chromosome number for the genotype
                              data for the samples you want to calculate the polygenic score levels in.
                              Defaults to UK Biobank, using the phase 3 release genotype data. For non-autosomal
                              chromosomes, the chromosomes are assumed to be 'X', 'Y', 'XY' and 'MT'.
                              [default: ~/rds/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/HRC_UK10K/plink_format/GRCh37/pgen/ukb_imp_v3_dedup_chr]
  --genotype-suffix <suffix>  Suffix for the filename occurring after the chromosome number but before the
                              .pgen/.pvar/.psam extension for the genotype data. [default: NULL]
  --genotype-format <format>  Format the genotype data is stored in, must correspond to one of the arguments to
                              plink, e.g. the default, 'pfile' is passed directly to plink as '--pfile'. To use
                              plink version 1 binary data (bed/bim/fam) set this as 'bfile'. [default: 'pfile']
  --single-geno               Flag to indicate that the genotype data is stored as a single file, not split across
                              multiple chromosomes. In this case, you can ignore the --genotype-suffix argument.
  --keep <file>               Optional, path to file to pass to plink2 --keep to subset to a given
                              set of samples when calculating the polygenic score levels. [default: NULL]
  --keep-ambiguous            Flag to force the program to keep variants with ambiguous alleles,
                              (e.g. A/T and G/C SNPs), which are normally excluded. In this case
                              the program proceeds assuming that the genotype data is on the
                              same strand as the GWAS whose summary statistics were used to
                              construct the score.
  --ambiguous-thresh <maf>    When provided, then variants with ambiguous alleles are kept only when their
                              minor allele frequency is below this threshold. Requires each score file to
                              have a column giving the effect allele's frequency. If multiple scores, then
                              any scores missing this column will exclude all variants with ambiguous alleles.
                              If using, we would suggest a threshold of 0.42. [default: NULL]
  --freqx-prefix <prefix>     Optional, path and prefix occurring before the chromosome number in the
                              filepaths for the plink1.9 --freqx reports containing the allele frequencies.
                              These may be used instead of those directly estimated from the genotype data,
                              e.g. when matching ambiguous alleles, or when imputing missing alleles when
                              using plink1 binary data. If --single-geno has been set, a single freqx file
                              is also assumed. [default: NULL]
  --freqx-suffix <suffix>     Optional, suffix for the above plink1.9 --freqx reports. [default: NULL]
  --remove-multiallelic       Flag that controls whether multi-allelic variants are kept or not. If given,
                              variants that have > 2 alleles in either the score file or genotype data are
                              discarded.
" -> doc

# Parse arguments and make additional sanity checks so errors are
# caught before long-running code fails
args <- docopt(doc)

# Set any "NULL" arguments to NULL
for (an in names(args)) {
  if (is.null(args[[an]]) || args[[an]] == "NULL") {
    args[[an]] <- NULL
  }
}

# Extract working directory
work_dir <- args[["--work"]]
stopifnot(dir.exists(work_dir))

# Save for easier debugging
saveRDS(args, file=sprintf("%s/args.rds", work_dir))

# Make sure we're running as part of an array job
if (Sys.getenv("SLURM_ARRAY_TASK_MAX") == "") {
  stop("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",
       "!!! Do not run this on the head node !!!\n",
       "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n",
       "This program must be run as part of an array job. Use calc_PS_lvls.sh to submit to slurm.\n\n")
}

# Check for --keep file if provided
if (!is.null(args[["--keep"]])) {
  if (!file.exists(args[["--keep"]])) {
    stop("File provided to --keep does not exist:\n", args[["--keep"]])
  }
}

# Check type
stopifnot(args[["type"]] %in% c('s', 'd', 'l'))

# Check thresholds
tryCatch({ args[["--mem"]] <- as.numeric(args[["--mem"]]) }, warning=function(w) { stop("--mem must be numeric") })
if (args[["--mem"]] <= 0) stop("--mem must be larger than 0")

if (!is.null(args[["--ambiguous-thresh"]])) {
	tryCatch({ args[["--ambiguous-thresh"]] <- as.numeric(args[["--ambiguous-thresh"]]) }, warning=function(w) { stop("--ambiguous-thresh must be numeric") })
	if ((args[["--ambiguous-thresh"]] < 0) || (args[["--ambiguous-thresh"]] > 0.5)) stop("--ambiguous-thresh must be between 0 and 0.5")
}

# Determine number of cores we can use:
ncores <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
if(is.na(ncores)) ncores <- 1
setDTthreads(ncores)

# Determine which chromosome we're running on:
taskIdx <- Sys.getenv("SLURM_ARRAY_TASK_ID")
chrIdx <- switch(taskIdx, "23" = "X", "24" = "Y", "25" = "XY", "26" = "MT", taskIdx)
taskIdx <- as.integer(taskIdx) # use to sleep processes to prevent race conditions on I/O
taskMax <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_MAX"))

# Determine chromosome specific file paths.
if (args[["--genotype-format"]] == "pfile") {
	varfile <- paste0(args[["--genotype-prefix"]], ifelse(args[["--single-geno"]], "", chrIdx), args[["--genotype-suffix"]], ".pvar")
	genofile <- paste0(args[["--genotype-prefix"]], ifelse(args[["--single-geno"]], "", chrIdx), args[["--genotype-suffix"]], ".pgen")
	samplefile <- paste0(args[["--genotype-prefix"]], ifelse(args[["--single-geno"]], "", chrIdx), args[["--genotype-suffix"]], ".psam")
} else if (args[["--genotype-format"]] == "bfile") {
	varfile <- paste0(args[["--genotype-prefix"]], ifelse(args[["--single-geno"]], "", chrIdx), args[["--genotype-suffix"]], ".bim")
	genofile <- paste0(args[["--genotype-prefix"]], ifelse(args[["--single-geno"]], "", chrIdx), args[["--genotype-suffix"]], ".bed")
	samplefile <- paste0(args[["--genotype-prefix"]], ifelse(args[["--single-geno"]], "", chrIdx), args[["--genotype-suffix"]], ".fam")
} else {
  stop("Unsuported genotype format: ", args[["--genotype-format"]])
}
freqxfile <- paste0(args[["--freqx-prefix"]], ifelse(args[["--single-geno"]], "", chrIdx), args[["--freqx-suffix"]])
if (freqxfile == "" || freqxfile == chrIdx) {
  freqxfile <- NULL
}

# Normalize paths
if (file.exists(varfile)) varfile <- normalizePath(varfile)
if (file.exists(genofile)) genofile <- normalizePath(genofile)
if (file.exists(samplefile)) samplefile <- normalizePath(samplefile)
if (!is.null(freqxfile)) freqxfile <- normalizePath(freqxfile)

# Check for genotype file at this point if a single file is provided
if (args[["--single-geno"]]) {
  stopifnot(file.exists(varfile))
  stopifnot(file.exists(genofile))
  stopifnot(file.exists(samplefile))
  stopifnot(is.null(freqxfile) || file.exists(freqxfile))
}

########################################################################################
# Define task safe I/O functions. These prevent multiple tasks from reading to or 
# writing to the same file. 
#
# If you run into issues with tasks inexplicably failing (and completing on subsequent
# runs) you might need to tweak the sleep times in these functions.
#
########################################################################################

# For a given file path, get the mount point and inode. We use this to create
# lock files for task safe I/O below. Inodes and mount points are used to prevent
# inappropriate locking, e.g. where many PGS files with the same name are read in
# (previous versions used the filename to create lockfiles). Mount point is included
# on the off chance files are being read/written from different filesystems in which
# case inodes might conflict
filestat <- function(fp) {
  fstat <- system(sprintf("stat -L -c '%%i %%m' %s", fp), intern=TRUE)
  mntpoint <- gsub(".* ", "", fstat)
  mntpoint <- gsub("/", "_", mntpoint)
  inode <- gsub(" .*", "", fstat)
  paste0(inode, mntpoint)
}

# Task safe readLines. Prevents readLines from running unless no other task is running
sreadLines <- function(con, ...) {
  lockfile <- sprintf("%s/00READLOCK-%s", work_dir, filestat(con))
  while(file.exists(lockfile)) {
    Sys.sleep(2) # file should be quick to read, so check frequently for read lock removal
  }
  cat(taskIdx, "\n", file=lockfile) # put task ID in the lockfile for debugging purposes
  on.exit({ system(sprintf("rm -f %s", lockfile), wait=TRUE) }) # on.exit means lock released even if readLines fails
  readLines(con, ...)
}

# fread can fall over when concurrently attempting to read the same score file from 26 processes.
# This function attempts to make this "task safe" by allowing only one call to fread to read from
# a given file at a time
sfread <- function(file, cmd, ...) {
  lockfile <- sprintf("%s/00FREADLOCK-%s", work_dir, filestat(file))
  while(file.exists(lockfile)) {
    Sys.sleep(5) 
  }
  cat(taskIdx, "\n", file=lockfile)
  on.exit({ system(sprintf("rm -f %s", lockfile), wait=TRUE) }) 
  if (!missing(cmd)) {
    fread(cmd=sprintf(cmd, file), tmpdir=work_dir, ...)
  } else {
    fread(file, tmpdir=work_dir, ...)
  }
}

# Task safe fwrite. Prevents fwrite from multiple tasks writing to the same file
sfwrite <- function(dt, file, ...) {
  fname <- basename(file)
  lockfile <- sprintf("%s/00FWRITELOCK-%s", work_dir, fname)
  while(file.exists(lockfile)) {
    Sys.sleep(5)
  }
  cat(taskIdx, "\n", file=lockfile) # put task ID in the lockfile for debugging purposes
  on.exit({ system(sprintf("rm -f %s", lockfile), wait=TRUE) }) # on.exit means lock released even if readLines fails
  fwrite(dt, file=file, ...)
  invisible(NULL)
}

########################################################################################
# Now we need to construct a single score file for all input scores by matching to the 
# designated genotype data. Each chromosome does this, this reduces memory load when
# loading the score variants (rather than loading them all for each chromosome).
########################################################################################

# construct table of score information - first step is to determine paths and
# column positions / names for each relevant column
is.integer <- function(value) {
  suppressWarnings(!is.na(as.integer(value)))
}
parsecolarg <- function(value) {
	if (is.null(value)) {
		return(NA)
	} else {
		return(value)
	}
}
if (args[["--type"]] == 's') {
	score_info <- data.table(path = args[["--score-file"]], 
													rsid = parsecolarg(args[["--score-rsid"]]),
													chr = parsecolarg(args[["--score-chr"]]),
													pos = parsecolarg(args[["--score-pos"]]),
													EA = parsecolarg(args[["--score-EA"]]),
													EAF = parsecolarg(args[["--score-EAF"]]),
													OA = parsecolarg(args[["--score-OA"]]),
													weight = parsecolarg(args[["--score-weight"]]),
													is_dom = parsecolarg(args[["--score-dominant"]]),
													is_rec = parsecolarg(args[["--score-recessive"]]),
													error = NA_character_)
} else if (args[["--type"]] == 'd') {
  score_info <- data.table(path = sprintf("%s/%s", args[["--score-file"]],
    setdiff(list.files(path=args[["--score-file"]]), list.dirs(path=args[["--score-file"]], full.names=FALSE, recursive=FALSE))),
													rsid = parsecolarg(args[["--score-rsid"]]),
													chr = parsecolarg(args[["--score-chr"]]),
													pos = parsecolarg(args[["--score-pos"]]),
													EA = parsecolarg(args[["--score-EA"]]),
													EAF = parsecolarg(args[["--score-EAF"]]),
													OA = parsecolarg(args[["--score-OA"]]),
													weight = parsecolarg(args[["--score-weight"]]),
													is_dom = parsecolarg(args[["--score-dominant"]]),
													is_rec = parsecolarg(args[["--score-recessive"]]),
													error = NA_character_)
	if (nrow(score_info) == 1 && is.na(score_info$path)) {
		stop("Directory is empty: ", args[["--score-file"]])
	}
} else if (args[["--type"]] == 'l') {
	tryCatch({ l <- sreadLines(args[["--score-file"]]) }, error=function(e) { 
		stop("Unable to read from file: ", args[["--score-file"]]) })
	score_info <- foreach(line = l, .combine=rbind) %do% {
		# Each line corresponds to a single score, optionally followed by 
		# 3-7 columns indicating the --score-<X> columns. Where this number
		# is < 7 reasonable defaults are inferred.
		fields <- strsplit(line, "\\s+")[[1]] # split on all whitespace.
		if (length(fields) == 1) { # Just a path to a score file, fill with defaults
			return(data.table(path = fields,
												rsid = parsecolarg(args[["--score-rsid"]]),
												chr = parsecolarg(args[["--score-chr"]]),
												pos = parsecolarg(args[["--score-pos"]]),
												EA = parsecolarg(args[["--score-EA"]]),
												EAF = parsecolarg(args[["--score-EAF"]]),
												OA = parsecolarg(args[["--score-OA"]]),
												weight = parsecolarg(args[["--score-weight"]]), 
												is_dom = parsecolarg(args[["--score-dominant"]]),
												is_rec = parsecolarg(args[["--score-recessive"]]),
												error = NA_character_))
		} else if (length(fields) == 10) { # All fields including dominant and recessive flags provided.
			return(data.table(path = fields[1], rsid = fields[2], chr = fields[3],
							pos = fields[4], EA = fields[5], EAF = fields[6],
							OA = fields[7], weight = fields[8], is_dom = fields[9], is_rec=fields[10],
							error = NA_character_))
		} else if (length(fields) == 9) { 
			return(data.table(path = fields[1], rsid = NA, chr = NA, pos = NA, EA = NA, EAF = NA, OA = NA, weight = NA, is_dom = NA, is_rec = NA,
							error = "Unable to infer use case when 9 fields provided."))
		} else if (length(fields) == 8) { # All fields provided, excluding dominant and recessive flags
			return(data.table(path = fields[1], rsid = fields[2], chr = fields[3], 
												pos = fields[4], EA = fields[5], EAF = fields[6],
												OA = fields[7], weight = fields[8], is_dom = NA, is_rec = NA, error = NA_character_))
		} else if (length(fields) == 7) {
			return(data.table(path = fields[1], rsid = fields[2], chr = fields[3],
							 pos = fields[4], EA = fields[5], EAF = NA, OA = fields[6], 
							 weight = fields[7], is_dom = NA, is_rec = NA, error = NA_character_))
		} else if (length(fields) == 6) {
			return(data.table(path = fields[1], rsid = NA, chr = fields[2],
							 pos = fields[3], EA = fields[4], EAF = NA, OA = fields[5], 
							 weight = fields[6], is_dom = NA, is_rec = NA, error = NA_character_))
		} else if (length(fields) == 5) {
			return(data.table(path = fields[1], rsid = NA, chr = fields[2],
							 pos = fields[3], EA = fields[4], EAF = NA, OA = NA, 
							 weight = fields[4], is_dom = NA, is_rec = NA, error = NA_character_))
		} else if (length(fields) == 4) {
			return(data.table(path = fields[1], rsid = fields[2], chr = NA,
							 pos = NA, EA = fields[3], EAF = NA, OA = NA, 
							 weight = fields[4], is_dom = NA, is_rec = NA, error = NA_character_))
		} else {
			return(data.table(path = fields[1], rsid = NA, chr = NA, pos = NA, EA = NA, EAF = NA, OA = NA, weight = NA, is_dom = NA, is_rec = NA,
												error = "At least three field to column mappings must be provided."))
		} 
	}
} else {
	stop("Internal error: allowed unknown --type")
}

# Get full path name
score_info[,path := normalizePath(path)]

# check all fields unique
field_check <- suppressWarnings(melt(score_info[!is.na(error)], id.vars="path"))
bad <- field_check[!is.na(value), .(fail=(length(unique(value)) == .N)), by=path][(fail)]
score_info[field_check, on = .(path), error := "The same column cannot be used for multiple fields."]
rm(bad, field_check)
invisible(gc())

# Check minimum set of fields present:
score_info[!is.na(error) & (is.na(EA) | is.na(weight) | (is.na(rsid) & (is.na(chr) | is.na(pos)))),
  error := "Must provide effect allele, weight, and either the rsid, or chromosome and position columns"]

# Give each score a unique number to uniquely identify it during computation
score_info[, compName := as.character(.I)]

# Reorganise
score_info = score_info[, .(path, compName, rsid, chr, pos, EA, EAF, OA, weight, is_dom, is_rec, error)]

# Check if each file exists
score_info[is.na(error), error := ifelse(file.exists(path), NA_character_, "score file does not exist")]

# At this point, if all scores have errored, we can stop all tasks.
if (all(!is.na(score_info$error))) {
  # Only one process needs to do this. Check for lock file or output file twice to prevent race conditions.
  if (!file.exists(sprintf("%s/00LOCK", work_dir)) && !file.exists("%s/score_summary.txt", work_dir)) {
    Sys.sleep(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))) 
  }
  if (!file.exists(sprintf("%s/00LOCK", work_dir)) && !file.exists("%s/score_summary.txt", work_dir)) {
    system(sprintf("touch %s/00LOCK", work_dir), wait=TRUE)
	  sfwrite(score_info, sep="\t", quote=FALSE, file=sprintf("%s/score_summary.txt", work_dir))
    system(sprintf("rm -f %s/00LOCK", work_dir), wait=TRUE)
  }
	stop("All score files had errors, see ", work_dir, "/score_summary.txt")
}

#################################################################################################
# Now, attempt to load each score file and map the columns. At this point we also detect whether
# they're PGS Catalog files and load appropriately
#################################################################################################

rbindf <- function(...) { rbind(..., fill=TRUE) }
scores <- foreach(idx = score_info[,.I], .combine=rbindf) %do% {
	# previously failed scores are skipped
	if (!is.na(score_info[idx, error])) {
		return(NULL)
	}

	# Read header line so we can determine if its a PGS Catalog file (and check we can read the file)
	tryCatch({ 
		line1 <- sreadLines(score_info[idx, path], 1)
	}, error = function(e) {
		score_info[idx, error := "Could not open score file"]
	})
	if (score_info[idx, !is.na(error)]) {
		return(NULL)
	}

	# Load the score with fread (and record an error if this fails)
	tryCatch({
    score <- sfread(file=score_info[idx, path], cmd="zgrep -v '^#' %s | zgrep -v '^[[:space:]]*$'")
	}, error = function(e) {
		score_info[idx, error := "Error when reading score with fread"]
	})
	if (score_info[idx, !is.na(error)]) {
		return(NULL)
	}

	pgs_catalog_file <- grepl("### ?PGS CATALOG SCORING FILE", line1)
	if (!pgs_catalog_file) {
		# Set the column names based on the information in the score_info field
		namecol <- function(score, old, new) {
			if (!is.na(old)) {
				if (is.integer(old)) {
          old <- as.integer(old)
					if (names(score)[old] != new && new %in% names(score)) {
						warning("Column named ", new, " found in score file (column ", 
						paste(which(names(score) == new), collapse=", "), " but using column ", 
						old, " [", names(score)[old], "] as --score-", new, " column instead")
						# rename extra columns so we don't pick the wrong one later
						which.new <- setdiff(which(names(score) == new), as.integer(old))
						names(score)[which.new] <- paste0(new, ".", seq_along(which.new))
					}
          score_info[idx, c(new) := names(score)[old]]
					setnames(score, names(score)[old], new)
				} else {
					if (sum(names(score) == old) > 1) {
						warning("Multiple columns named ", old, "in score file ", score_info[idx, path], " using first occurence as --score-", new, " column")
						# rename extra columns so we don't pick the wrong one later
						which.dup = which(names(score) == old)[-1]
						names(score)[which.dup] <- paste0(old, ".", seq_along(which.dup))
					} 
					if (new != old && new %in% names(score)) {
						warning("Column named ", new, " found in score file (column ", 
						paste(which(names(score) == new), collapse=", "), " but using column ", 
						old, " (column ", paste(which(names(score) == old), collapse=", "), 
						") as --score-", new, " column instead")
						# rename extra columns so we don't pick the wrong one later
						which.new <- which(names(score) == new)
						names(score)[which.new] <- paste0(new, ".", seq_along(which.new))
					}
 
					setnames(score, old, new, skip_absent=TRUE)
				}
			}
		}
		tryCatch({
			namecol(score, score_info[idx, rsid], "rsid")
			namecol(score, score_info[idx, chr], "chr")
			namecol(score, score_info[idx, pos], "pos")
			namecol(score, score_info[idx, EA], "EA")
			namecol(score, score_info[idx, EAF], "EAF")
			namecol(score, score_info[idx, OA], "OA")
			namecol(score, score_info[idx, is_dom], "is_dom")
			namecol(score, score_info[idx, is_rec], "is_rec")
			if (score_info[idx, weight != "m"]) { 
				namecol(score, score_info[idx, weight], "weight")
			}
		}, error=function(e) {
			score_info[idx, error := "Could not match provided --score-<X> columns to columns in score file"]
		})
    if (score_info[idx, !is.na(error)]) {
      return(NULL)
    }
 
    # Record how many variants are in the score file, and how many are on this chromosome
    score_info[idx, n_var := score[,.N]]
    if ("chr" %in% names(score)) {
			score = score[chr == chrIdx] # filter to this chromosome. 
    }

		# If multiple scores, transform to long format
		if (score_info[idx, weight == "m"]) {
			# first, make sure score name columns are unique, generate warning if necessary
			sn = data.table(name = names(score))
			dups = sn[,.N,by=name][N > 1, name]
			sn[name %chin% dups, new := paste0(name, ".", seq_len(.N)), by=name]
			sn[is.na(new), new := name]
			if (sn[, any(name != new)]) {
				warning("Non-unique score weight column names found in multi-score file", score_info[idx, path], 
								". Numbers have been appended to make each score's name unique")
			}
 
      # Give each sub-score a 'compName' 
      sn[!(name %chin% c("rsid", "chr", "pos", "EA", "EAF", "OA", "is_dom", "is_rec")),
           compName := gsub(" ", "0", paste0(score_info[idx, compName], ".", format(.I)))]
      sn[is.na(compName), compName := name]
      setnames(score, sn[,compName])
			
      # Melt to long format
			score <- melt(score, id.vars=intersect(names(score), c("rsid", "chr", "pos", "EA", "EAF", "OA", "is_dom", "is_rec")), 
										variable.name="compName", value.name="weight")

      # create score sub-table to append to score-info
      sub_info <- score_info[idx]
      sub_info[, compName := NULL]
      sub_info <- cbind(sn[!(name %chin% c("rsid", "chr", "pos", "EA", "EAF", "OA", "is_dom", "is_rec")), .(compName)], sub_info)
      sub_info[sn, on = .(compName), weight := new]
      score_info <- rbindf(score_info, sub_info)
      
      # Temporarily drop zero-weights to save on memory before we later cast back to wide format.
      if (score[, any(is.na(weight))]) {
        warning("Missing values found in one or more weight columns and converted to 0 for multi-score file ", score_info[idx, path])
        score[is.na(weight), weight := 0]
      }
      score <- score[weight != 0] 
		}

		# If we're working with a single score drop extra columns that might be in the file and add the compName
		if (score_info[idx, weight != "m"]) {
			score <- score[, which(names(score) %chin% c("rsid", "chr", "pos", "EA", "EAF", "OA", "is_dom", "is_rec", "weight")), with=FALSE]
      if (score[, any(is.na(weight))]) {
        warning("Missing values found in weight column and converted to 0 for score file ", score_info[idx, path])
        score[is.na(weight), weight := 0]
      }
      score[, compName := score_info[idx, compName]]
		}
    return(score)
	} else {
		# For PGS catalog files, we need to map score names to ones provided standard in the catalog.
    mapname <- function(old, new) {
      if (old %in% names(score)) {
        setnames(score, old, new)
        if (!is.character(score_info[[new]])) {
          score_info[, c(new) := as.character(score_info[[new]])]
        }
        score_info[idx, c(new) := old]
      } else {
        score_info[idx, c(new) := NA]
      }
    }
    tryCatch({
			mapname("rsID", "rsid")
			mapname("chr_name", "chr")
			mapname("chr_position", "pos")
			mapname("effect_allele", "EA")
			mapname("allelefrequency_effect", "EAF")
			mapname("effect_weight", "weight")
			mapname("is_dominant", "is_dom")
			mapname("is_recessive", "is_rec")

      # Old vs. new score file format has different names
      if ("reference_allele" %in% names(score)) {
				mapname("reference_allele", "OA")
      } else {
				mapname("other_allele", "OA")
      }
    }, error=function(e) {
      score_info[idx, error := "PGS Catalog file detected, but could not map column names"]
    })
    if (score_info[idx, !is.na(error)]) {
      return(NULL)
    }

    # Record how many variants are in the score file
    score_info[idx, n_var := score[,.N]]
    if ("chr" %in% names(score)) {
			score = score[chr == chrIdx] # filter to this chromosome. 
    }

		# Drop interaction terms
		if ("is_interaction" %in% names(score)) {
			int_var <- score[(is_interaction), .N,]
			score <- score[!(is_interaction)]
			score[, is_interaction := NULL]
			score_info[idx, n_interaction_skipped := int_var]
		}

    # Error if any variants have the same effect weight for both dominant and recessive effects
    if ("is_dom" %in% names(score) && "is_rec" %in% names(score)) {
      if (nrow(score[(is_dom) & (is_rec)]) > 0) {
        score_info[idx, error := "Entries with TRUE in both --score-dominant and --score-recessive columns"]
        return(NULL)
      }
    }

    # Add compName
    score[, compName := score_info[idx, compName]]

		# Drop any columns we won't use
	  score <- score[, which(names(score) %chin% c("rsid", "chr", "pos", "EA", "EAF", "OA", "is_dom", "is_rec", "weight", "compName")), with=FALSE]

    # Detect and convert missing values
		if (score[, any(is.na(weight))]) {
			warning("Missing values found in weight column and converted to 0 for score file ", score_info[idx, path])
			score[is.na(weight), weight := 0]
		}
    return(score)
	}
}
rm(score)
invisible(gc())

# At this point, if all scores have errored, we can stop all tasks.
if (all(!is.na(score_info$error))) {
  # Only one process needs to do this. Check for lock file or output file twice to prevent race conditions.
  if (!file.exists(sprintf("%s/00LOCK", work_dir)) && !file.exists("%s/score_summary.txt", work_dir)) {
    Sys.sleep(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))) 
  }
  if (!file.exists(sprintf("%s/00LOCK", work_dir)) && !file.exists("%s/score_summary.txt", work_dir)) {
    system(sprintf("touch %s/00LOCK", work_dir), wait=TRUE)
	  sfwrite(score_info, sep="\t", quote=FALSE, file=sprintf("%s/score_summary.txt", work_dir))
    system(sprintf("rm -f %s/00LOCK", work_dir), wait=TRUE)
  }
	stop("All score files had errors, see ", work_dir, "/score_summary.txt")
}

# make sure chromosome column is character:
if ("chr" %in% names(scores) && !is.character(scores$chr)) {
  scores[, chr := as.character(chr)]
}

# Fix any scores that use lower case alleles (WHY?)
scores[, EA := toupper(EA)]
if ("OA" %in% names(scores)) {
  scores[!is.na(OA), OA := toupper(OA)]
}

# Reorganise rows
score_info = score_info[order(as.numeric(compName))]
    
# There may be no variants on this chromosome 
if (nrow(scores) == 0) {
  system(sprintf("touch %s/00FIN-chr-%s", work_dir, chrIdx), wait=TRUE)
  quit(save="no")
}

################################################################################################
# Now each chromosome attempts to load the genotype variant information and orient variants 
# to the same allele in UK Biobank. Also at this point we want to filter out ambiguous alleles
# unless otherwise specified. When matching variants, default to matching by chromosome, 
# position, and alleles unless otherwise asked (i.e. --match-by-rsid, or score missing
# chr and pos columns).
################################################################################################

# Function for flipping the strand of an allele.
# Uses a series of gsub calls to replace A's with T's,
# G's with C's, and vice-versa. Also works for alleles
# with more than one nucleotide (e.g. indels).
flip_strand <- function(x) {
  # Swap each letter for a dummy, we need this intermediate
  # step so we can distinguish between alleles when swapping.
  # E.g if we did A -> T then T -> A we'd end up with all A's
  # and no T's. instead we do A -> V -> T and T -> X -> A.
  x <- gsub("A", "V", x)
  x <- gsub("T", "X", x)
  x <- gsub("C", "Y", x)
  x <- gsub("G", "Z", x)
  x <- gsub("V", "T", x)
  x <- gsub("X", "A", x)
  x <- gsub("Y", "G", x)
  x <- gsub("Z", "C", x)
  return(x)
}

# Attempt to load the variant information for this chromosome, exiting if not found
tryCatch({
  # Detect if pgen from vcf
  line1 <- sreadLines(varfile, 1)
  if (grepl("^##", line1)) {
    varinfo <- sfread(varfile, cmd="grep -v '^#' %s", sep="\t")
  } else {
		varinfo <- sfread(varfile)
  }
	if(args[["--genotype-format"]] == "pfile") {
    varinfo <- varinfo[,1:5] # drop extra columns from VCF pvar files
		setnames(varinfo, c("chr", "pos", "rsid", "ref", "alt"))
  } else if (args[["--genotype-format"]] == "bfile") {
		setnames(varinfo, c("chr", "rsid", "cm", "pos", "alt", "ref"))
  }
  if (!is.character(varinfo$chr)) {
    varinfo[,chr := as.character(chr)]
  }
}, error=function(e) {
  score_info[is.na(error) & compName %in% unique(scores[!is.na(chr), compName]), c("n_match", "chr_fail") := .(0, TRUE)]
  fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("%s/score_summary_%s.txt", work_dir, chrIdx))
  system(sprintf("touch %s/00FIN-chr-%s", work_dir, chrIdx), wait=TRUE)
  quit(save="no")
})

# At this point we can create the links to the genetic data for score computation.
# We need to give each variant a unique identifier in these files.
if (args[["--genotype-format"]] == "pfile") {
	system(sprintf("ln -s %s %s/chr%s.pgen", genofile, work_dir, chrIdx), wait=TRUE)
	system(sprintf("ln -s %s %s/chr%s.psam", samplefile, work_dir, chrIdx), wait=TRUE)
  fwrite(varinfo[, .(`#CHROM`=chr, POS=pos, ID=paste(chr, pos, alt, ref, sep=":"), REF=ref, ALT=alt)],
         sep="\t", quote=FALSE, file=sprintf("%s/chr%s.pvar", work_dir, chrIdx))
} else if (args[["--genotype-format"]] == "bfile") {
	system(sprintf("ln -s %s %s/chr%s.bed", genofile, work_dir, chrIdx), wait=TRUE)
	system(sprintf("ln -s %s %s/chr%s.fam", samplefile, work_dir, chrIdx), wait=TRUE)
  fwrite(varinfo[, .(chr, rsid=paste(chr, pos, alt, ref, sep=":"), cm, pos, alt, ref)],
         sep="\t", quote=FALSE, col.names=FALSE, file=sprintf("%s/chr%s.bim", work_dir, chrIdx))
}

# Note in the above, if there is a single genotype data file then this is linked N times,
# 1 per chromosome, and copies of varinfo are written out for each chromosome. This means
# all downstream computation can assume 1 file per chromosome.
# If this is the case, we also now need to check whether there are any variants on this
# chromosome (e.g. if the task is working on chromosomes X, Y, XY, or MT this may not be
# the case) and note error (only occurs if one or more score has a variant on these 
# chromosomes).
if (args[["--single-geno"]]) {
	varinfo <- varinfo[chr == chrIdx]
	if (nrow(varinfo[chr == chrIdx]) == 0) {
    score_info[compName %in% unique(scores$compName), c("n_match", "chr_fail") := .(0, TRUE)]
		fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("%s/score_summary_%s.txt", work_dir, chrIdx))
    system(sprintf("rm -f %s/chr%s*", work_dir, chrIdx), wait=TRUE) # remove temporary geno files.
		system(sprintf("touch %s/00FIN-chr-%s", work_dir, chrIdx), wait=TRUE)
		quit(save="no")
	}
}

# How many variants on this chromosome for each score?
if ("chr" %in% names(scores)) {
	score_info[is.na(error) & !is.na(chr), n_chr := 0]
	score_info[scores[!is.na(chr), .N, by=compName], on = .(compName), n_chr := i.N]
	if (!("chr" %in% names(scores)) && !("pos" %in% names(scores))) {
		score_info[scores[!is.na(chr), .N, by=.(compName=as.character(as.integer(as.vector(compName))), rsid, EA)][, .N, by=compName], on = .(compName), n_chr := i.N] # overall summary for multi-score file
	} else {
		score_info[scores[!is.na(chr), .N, by=.(compName=as.character(as.integer(as.vector(compName))), chr, pos, EA)][, .N, by=compName], on = .(compName), n_chr := i.N] # overall summary for multi-score file
	}
} else {
  score_info[, chr := NA_integer_]
}

# Obtain chromosome and position for variants missing this information,
# (or override if explictly asked to match by rsid)
# then subset again to just variants on this task's chromosome
if ("rsid" %in% names(scores)) {
  scores[varinfo[!is.na(rsid)], on = .(rsid), c("chr", "pos") :=
         .(ifelse((args[["--match-by-rsid"]] & !is.na(rsid)) | !("chr" %in% names(scores)) | is.na(chr), i.chr, chr),
           ifelse((args[["--match-by-rsid"]] & !is.na(rsid)) | !("pos" %in% names(scores)) | is.na(pos), i.pos, pos))]
}
scores <- scores[chr == chrIdx]

# can exit if no variants on this chromosome
if (nrow(scores) == 0) {
  if (args[["--genotype-format"]] == "pfile") {
    system(sprintf("rm -f %s/chr%s.pvar", work_dir, chrIdx), wait=TRUE)
    system(sprintf("rm -f %s/chr%s.psam", work_dir, chrIdx), wait=TRUE)
    system(sprintf("rm -f %s/chr%s.pgen", work_dir, chrIdx), wait=TRUE)
  } else if (args[["--genotype-format"]] == "bfile") {
    system(sprintf("rm -f %s/chr%s.bim", work_dir, chrIdx), wait=TRUE)
    system(sprintf("rm -f %s/chr%s.fam", work_dir, chrIdx), wait=TRUE)
    system(sprintf("rm -f %s/chr%s.bed", work_dir, chrIdx), wait=TRUE)
  }
  fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("%s/score_summary_%s.txt", work_dir, chrIdx))
  system(sprintf("touch %s/00FIN-chr-%s", work_dir, chrIdx), wait=TRUE)
  quit(save="no")
}

# Before matching variants by alleles, check for field separators in
# the other allele column. The matching process won't work with these,
# so we can just set those other allele fields to NA to match by effect
# allele only
if ("OA" %in% names(scores)) {
	scores[!grepl("^(A|G|C|T)+$", OA), OA := NA]
}

# Check if we can match by chromosome, position, and alleles
scores[, match := FALSE]
if ("OA" %in% names(scores)) {
  # Where both effect and other allele provided, must match entry by both
  scores[varinfo, on = .(chr, pos, EA=alt, OA=ref), match := !is.na(OA)]
  scores[varinfo, on = .(chr, pos, EA=ref, OA=alt), match := !is.na(OA)]
  # If no other allele provided, match just by effect allele
  scores[varinfo, on = .(chr, pos, EA=alt), match := match | is.na(OA)]
  scores[varinfo, on = .(chr, pos, EA=ref), match := match | is.na(OA)]
} else {
  scores[varinfo, on = .(chr, pos, EA=alt), match := TRUE]
  scores[varinfo, on = .(chr, pos, EA=ref), match := TRUE]
}

# For variants that did not match above, but can match by chromosome and position,
# flip the strand of the alleles and attempt match again
if ("OA" %in% names(scores)) {
  # Flip alleles
  scores[varinfo, on = .(chr, pos), c("EA", "OA") := 
         .(ifelse(match, EA, flip_strand(EA)), ifelse(match, OA, flip_strand(OA)))]
  # Where both effect and other allele provided, must match entry by both
  scores[varinfo, on = .(chr, pos, EA=alt, OA=ref), match := ifelse(match, match, !is.na(OA))]
  scores[varinfo, on = .(chr, pos, EA=ref, OA=alt), match := ifelse(match, match, !is.na(OA))]
  # If no other allele provided, match just by effect allele
  scores[varinfo, on = .(chr, pos, EA=alt), match := ifelse(match, match, match | is.na(OA))]
  scores[varinfo, on = .(chr, pos, EA=ref), match := ifelse(match, match, match | is.na(OA))]
} else {
  scores[varinfo, on = .(chr, pos), c("EA") := .(ifelse(match, EA, flip_strand(EA)))]
  scores[varinfo, on = .(chr, pos, EA=alt), match := ifelse(match, match, TRUE)]
  scores[varinfo, on = .(chr, pos, EA=ref), match := ifelse(match, match, TRUE)]
}

# Drop score variants that could not be matched to the genotype data
scores <- scores[(match)]
scores[, match := NULL]
invisible(gc())

# How many variants could we match to the genotype data?
score_info[is.na(error), n_match := 0]
score_info[scores[, .N, by=compName], on = .(compName), n_match := i.N]
score_info[scores[, .N, by=.(compName=as.character(as.integer(as.vector(compName))), chr, pos)][, .N, by=compName], on = .(compName), n_match := i.N] # overall summary for multi-score file

# can exit if no variants matched
if (nrow(scores) == 0) {
  if (args[["--genotype-format"]] == "pfile") {
    system(sprintf("rm -f %s/chr%s.pvar", work_dir, chrIdx), wait=TRUE)
    system(sprintf("rm -f %s/chr%s.psam", work_dir, chrIdx), wait=TRUE)
    system(sprintf("rm -f %s/chr%s.pgen", work_dir, chrIdx), wait=TRUE)
  } else if (args[["--genotype-format"]] == "bfile") {
    system(sprintf("rm -f %s/chr%s.bim", work_dir, chrIdx), wait=TRUE)
    system(sprintf("rm -f %s/chr%s.fam", work_dir, chrIdx), wait=TRUE)
    system(sprintf("rm -f %s/chr%s.bed", work_dir, chrIdx), wait=TRUE)
  }
  fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("%s/score_summary_%s.txt", work_dir, chrIdx))
  system(sprintf("touch %s/00FIN-chr-%s", work_dir, chrIdx), wait=TRUE)
  quit(save="no")
}

# Identify and flag variants that are multi-allelic in either the score file or in the genotype data
genomulti <- varinfo[,.N,by=pos][N > 1, .(pos)]
multi <- unique(rbind(
  scores[genomulti, on=.(pos), nomatch=0, .(pos, compName)],
  scores[,.N,by=.(pos,compName)][N > 1, .(pos, compName)]))

# Collate number of multi-allelic variants
score_info[is.na(error), n_multiallele := 0]
score_info[multi[,.N, by=compName], on = .(compName), n_multiallele := i.N]
score_info[multi[, .N, by=.(compName=as.character(as.integer(as.vector(compName))), pos)][, .N, by=compName], on = .(compName), n_multiallele := i.N] # overall summary for multi-score file

# Remove them if requested
if (args[["--remove-multiallelic"]]) {
  scores <- scores[!multi, on = .(pos)]
  score_info[, n_multiallele_removed := n_multiallele]
} else {
  score_info[, n_multiallele_removed := ifelse(n_multiallele > 0, 0, NA)]
}

# Check at this point whether each score has only 1 effect weight per variant.
bad <- scores[,.N,by=.(pos, compName)][N > 1, .(compName)]
score_info[bad, on = .(compName), error := "Score has multiple effect weights for the same variant/position"]
scores <- scores[!bad, on = .(compName)]

# can exit if all errors
if (nrow(scores) == 0) {
  if (args[["--genotype-format"]] == "pfile") {
    system(sprintf("rm -f %s/chr%s.pvar", work_dir, chrIdx), wait=TRUE)
    system(sprintf("rm -f %s/chr%s.psam", work_dir, chrIdx), wait=TRUE)
    system(sprintf("rm -f %s/chr%s.pgen", work_dir, chrIdx), wait=TRUE)
  } else if (args[["--genotype-format"]] == "bfile") {
    system(sprintf("rm -f %s/chr%s.bim", work_dir, chrIdx), wait=TRUE)
    system(sprintf("rm -f %s/chr%s.fam", work_dir, chrIdx), wait=TRUE)
    system(sprintf("rm -f %s/chr%s.bed", work_dir, chrIdx), wait=TRUE)
  }
  fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("%s/score_summary_%s.txt", work_dir, chrIdx))
  system(sprintf("touch %s/00FIN-chr-%s", work_dir, chrIdx), wait=TRUE)
  quit(save="no")
}

# Now that we've done the matching, give the variants in the score file
# the same unique identifier they have in the variant information file.
# For multi-allelic variants we want to avoid double counting alleles,
# e.g. where the effect allele is the one in common across the multiple
# rows.
scores[varinfo, on = .(pos, EA=alt), rsid := paste(i.chr, i.pos, i.alt, i.ref, sep=":")]
scores[varinfo, on = .(pos, EA=ref), rsid := paste(i.chr, i.pos, i.alt, i.ref, sep=":")]

# Now identify variants whose effect allele is ambiguous (e.g. A/T or G/C SNPs). The 
# following works for biallelic sites and some multi-allelic sites
scores[varinfo, on = .(pos, EA=alt), ambig := EA == flip_strand(ref)]
scores[varinfo, on = .(pos, EA=ref), ambig := EA == flip_strand(alt)]

# For multi-allelic sites that otherwise weren't identified as ambiguous:
# 1. if the OA column is provided, check if the effect allele is strand ambiguous
#    with any other allele at that position.
# 2. Remove the OA column.
# 3. Add in all OA from the varinfo table, then check again
mult_ambig <- scores[multi, on = .(pos, compName)][!(ambig)]
if ("OA" %in% names(mult_ambig)) {
  mult_ambig[, ambig := any(EA == flip_strand(OA)), by=pos]
  scores[mult_ambig[(ambig)], on = .(pos, EA), ambig := TRUE]
  mult_ambig <- mult_ambig[!(ambig)]
  mult_ambig[, OA := NULL]
}
mult_ambig <- rbind(
  varinfo[, .(pos, EA=alt, OA=ref)][mult_ambig, on = .(pos, EA), nomatch=0],
  varinfo[, .(pos, EA=ref, OA=alt)][mult_ambig, on = .(pos, EA), nomatch=0])
mult_ambig[, ambig := any(EA == flip_strand(OA)), by=pos]
mult_ambig <- unique(mult_ambig[,.(pos, EA, ambig)])
scores[mult_ambig[(ambig)], on = .(pos, EA), ambig := TRUE]
rm(mult_ambig)

# Extract ambiguous variants so we can do further checks by allele frequency if requested
ambig <- scores[(ambig)]
scores <- scores[!(ambig)]
scores[, ambig := NULL]
ambig[, ambig := NULL]
invisible(gc())

# How many ambiguous alleles in each score?
score_info[is.na(error), n_ambig := 0]
if (nrow(ambig) > 0) {
  score_info[ambig[,.N,by=compName], on = .(compName), n_ambig := i.N]
  score_info[ambig[,.N, by=.(compName=as.character(as.integer(as.vector(compName))), chr, pos)][, .N, by=compName], on = .(compName), n_ambig := i.N] # overall summary for multi-score file
}

# If the freqx file has been provided (either for allele frequency checking of ambiguous SNPs,
# or for providing directly to plink2 to skip allele frequency calculations) we need to add 
# the missing position information and give it the same unique variant identifiers as the
# variant info file.
tryCatch({
  if (!is.null(freqxfile)) {
	  freqx <- sfread(freqxfile)
  }
}, error=function(e) {
  warning("Error when attempting to read freqx file: ", freqxfile, " computing from genotype data instead")
})
# Parse the freqx file. Depends on the version of plink and format requested:
# Parse the freqx file. Depends on the version of plink and format requested:
if (exists("freqx")) {
  if ("C(HOM A1)" %in% names(freqx)) { # Plink 1.9 --freqx
    if (inherits(freqx$CHR, "integer")) freqx[, CHR := as.character(CHR)]
    freqx[varinfo, on = .(SNP=rsid, CHR=chr, A2=ref, A1=alt), POS := i.pos]
    freqx[, SNP := paste(CHR, POS, A1, A2, sep=":")]
    freqx[, POS := NULL]
    fwrite(freqx, sep="\t", quote=FALSE, file=sprintf("%s/chr%s.frqx", work_dir, chrIdx))
    freqx_ext <- "frqx"
  } else if ("MAF" %in% names(freqx)) { # Plink 1.9 --freq
    if (inherits(freqx$CHR, "integer")) freqx[, CHR := as.character(CHR)]
    freqx[varinfo, on = .(SNP=rsid, CHR=chr, A2=ref, A1=alt), POS := i.pos]
    freqx[, SNP := paste(CHR, POS, A1, A2, sep=":")]
    freqx[, POS := NULL]
    fwrite(freqx, sep="\t", quote=FALSE, file=sprintf("%s/chr%s.frq", work_dir, chrIdx))
    freqx_ext <- "frq"
  } else if ("ALT" %in% names(freqx)) { # Plink 2 --freq
    setnames(freqx, "#CHROM", "CHROM")
    if (inherits(freqx$CHROM, "integer")) freqx[, CHROM := as.character(CHROM)]
    freqx[varinfo, on = .(ID=rsid, CHROM=chr, ALT=alt, REF=ref), POS := i.pos]
    freqx[, ID := paste(CHROM, POS, ALT, REF, sep=":")]
    freqx[, POS := NULL]
    fwrite(freqx, sep="\t", quote=FALSE, file=sprintf("%s/chr%s.afreq", work_dir, chrIdx))
    freqx_ext <- "afreq"
  } else {
    rm(freqx)
    warning("Unrecognised freqx format, computing from genotype data instead")
  }
}

# Now we can give the loaded varinfo table the appropriate unique identifiers
# and subset to just the score variants
varinfo[, rsid := paste(chr, pos, alt, ref, sep=":")]
varinfo <- varinfo[pos %in% unique(c(scores$pos, ambig$pos))] # we want to preserve row order, so filter instead of join.

# Handle ambiguous variants depending of various input options and collate
if (nrow(ambig) > 0) {
	if (!args[["--keep-ambiguous"]]) {
		ambig[, keep := FALSE]
	} else if (args[["--keep-ambiguous"]] && is.null(args[["--ambiguous-thresh"]])) {
		ambig[, keep := TRUE]
	} else if (!("EAF" %in% names(scores))) {
		# Requested to match by effect allele frequency, but no scores provided this.
		ambig[, keep := FALSE] 
	} else { # keep if can match by effect allele frequency
		if (!exists("freqx")) {
			# Need to make sure we give variants unique identifiers when writing out
			fwrite(ambig[,.(unique(rsid))], col.names=FALSE, quote=FALSE,
						 file=sprintf("%s/ambig_freqx_extract_chr%s", work_dir, chrIdx))

      # If there's only a single genotype file, we want to prevent multiple tasks 
      # reading from the same file.
      if (args[["--single-geno"]]) {
        lockfile = sprintf("%s/00READLOCK-GENO", work_dir)
        while(file.exists(lockfile)) {
          Sys.sleep(5)
        }
        system(sprintf("touch %s", lockfile), wait=TRUE)
      }

      # Do frequency calculation for ambiguous variants
			cmd <- normalizePath(sprintf("%s/plink2", SoftwareDir))
			cmd[2] <- sprintf("--%s %s/chr%s", args[["--genotype-format"]], work_dir, chrIdx)
			cmd[3] <- sprintf("--extract %s/ambig_freqx_extract_chr%s", work_dir, chrIdx)
			cmd[4] <- sprintf("--threads %s --memory 1024 --silent", ncores)
			cmd[5] <- sprintf("--freq --out %s/ambig_freqx_extract_chr%s", work_dir, chrIdx)
			if (!is.null(args[["--keep"]])) cmd[6] <- paste("--keep", args[["--keep"]])
			errorcode <- system(paste(cmd, collapse=" "), wait=TRUE)

			if (errcode != 0) {
				logfile <- sprintf("%s/ambig_freqx_extract_chr%s", work_dir, chrIdx)
				score_info[, error := sprintf("%sError from plink2 when computing frequencies of ambiguous alleles, see log file: %s", 
					ifelse(is.na(error) | error == "", "", paste0(error, ".")), logfile)]
				fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("%s/score_summary_%s.txt", work_dir, chrIdx))
				system(sprintf("touch %s/00FIN-chr-%s", work_dir, chrIdx), wait=TRUE)
				quit(save="no")
			}

      # Release readlock
      if (args[["--single-geno"]]) {
        system(sprintf("rm -f %s", lockfile), wait=TRUE)
      }

      # load allele frequencies
			freqx <- sfread(sprintf("%s/ambig_freqx_extract_chr%s.afreq", work_dir, chrIdx))
		  freqx <- freqx[, .(rsid=ID, EA=ALT, EAF=ALT_FREQS)]
		} else {
      # Need to extract EAF from loaded freqx file
      if (freqx_ext == "frqx") {
        freqx <- freqx[, .(rsid=SNP, EA=A1, EAF=(`C(HOM A1)`*2+`C(HET)`)/(`C(HOM A1)`*2+`C(HET)`+`C(HOM A2)`*2))]
      } else if (freqx_ext == "frq") {
        freqx <- freqx[, .(rsid=SNP, EA=A1, EAF=MAF)]
      } else if (freqx_ext == "afreq") {
        freqx <- freqx[, .(rsid=ID, EA=ALT, EAF=ALT_FREQS)]
      }
    }

    # Flag multi-allelic SNPs in freqx file - for these we can only match the the allele whose frequency has been calculated
    freqx[, pos := as.numeric(gsub(":.*", "", gsub("^.*?:", "", rsid)))]
    freqx[varinfo[,.N,by=pos], on=.(pos), multi := N > 1]

		# Where EAF column not provided, set keep to FALSE
		ambig[, keep := TRUE]
		ambig[is.na(EAF), keep := FALSE]

		# Attempt to match ambiguous variants by EAF below threshold:
    ambig[freqx, on = .(rsid, EA), # effect allele is the one we have frequency for in freqx.
          keep := ifelse(keep &    # ignore score variants with no EAF information
            ((EAF < args[["--ambiguous-thresh"]] & i.EAF < args[["--ambiguous-thresh"]]) |    # effect allele is minor allele, and frequency below threshold
             (EAF > (1 - args[["--ambiguous-thresh"]]) & i.EAF > (1 - args[["--ambiguous-thresh"]]))), # effect allele is major allele, and frequency below threshold
            TRUE, FALSE)]
    # The same as above, but instead the effect allele is the allele whose frequency was not calculated.
    # In this case we can take the effect allele's EAF in the genotype data as 1 - the calculated EAF,
    # provided we're dealing only with bi-allelic SNPs.
    ambig[freqx[!(multi)], on = .(rsid, OA=EA),
          keep := ifelse(keep &   
            ((EAF < args[["--ambiguous-thresh"]] & (1 - i.EAF) < args[["--ambiguous-thresh"]]) |    # effect allele is minor allele, and frequency below threshold
             (EAF > (1 - args[["--ambiguous-thresh"]]) & (1 - i.EAF) > (1 - args[["--ambiguous-thresh"]]))), # effect allele is major allele, and frequency below threshold
            TRUE, FALSE)]
    
    # Remove no longer needed freqx object
    rm(freqx)
	}
  
  # Filter to just ambiguous variants we want to keep
  ambig <- ambig[(keep)]
  ambig[, keep := NULL]
  invisible(gc())
}

# How many ambiguous alleles are kept?
score_info[is.na(error), n_ambig_kept := 0]
if (nrow(ambig) > 0) {
	score_info[ambig[,.N,by=compName], on = .(compName), n_ambig_kept := i.N]
	score_info[ambig[,.N, by=.(compName=as.character(as.integer(as.vector(compName))), chr, pos)][, .N, by=compName], on = .(compName), n_ambig_kept := i.N] # overall summary for multi-score file
}

# Add kept variants back to score and collate
scores <- rbind(scores, ambig)
rm(ambig)
invisible(gc())
score_info[, n_used := n_match - n_ambig + n_ambig_kept]

# can exit if no variants left (i.e. all matched were ambiguous)
if (nrow(scores) == 0) {
  if (args[["--genotype-format"]] == "pfile") {
    system(sprintf("rm -f %s/chr%s.pvar", work_dir, chrIdx), wait=TRUE)
    system(sprintf("rm -f %s/chr%s.psam", work_dir, chrIdx), wait=TRUE)
    system(sprintf("rm -f %s/chr%s.pgen", work_dir, chrIdx), wait=TRUE)
  } else if (args[["--genotype-format"]] == "bfile") {
    system(sprintf("rm -f %s/chr%s.bim", work_dir, chrIdx), wait=TRUE)
    system(sprintf("rm -f %s/chr%s.fam", work_dir, chrIdx), wait=TRUE)
    system(sprintf("rm -f %s/chr%s.bed", work_dir, chrIdx), wait=TRUE)
  }
  fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("%s/score_summary_%s.txt", work_dir, chrIdx))
  system(sprintf("touch %s/00FIN-chr-%s", work_dir, chrIdx), wait=TRUE)
  quit(save="no")
}

# Split scores into those for linear, dominant, and recessive effects,
# transform to wide format, write out, and remove.
sdcast <- function(...) {
  tryCatch({
    dcast(...)
  }, warning=function(w) {
    score_info[!is.na(error), error := paste0("Critical error collating score weights on chromosome ", 
                                              chrIdx, ": encountered multiple weights for the same effect allele.")]
    fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("%s/score_summary_%s.txt", work_dir, chrIdx))
    system(sprintf("touch %s/00FIN-chr-%s", work_dir, chrIdx), wait=TRUE)
    quit(save="no")
  })
}

# If different scores use different effect alleles, we will need to create more than one file since
# plink can handle only one entry per variant
make_score_files <- function(dt, model_name, file_prefix) {
  var_ea <- unique(dt[,.(rsid, EA)])
  var_ea[, group := seq_len(.N), by=rsid]
  sfinfo <- data.table(model = model_name, group = unique(var_ea$group), ncol=0) 
  for (group_id in unique(var_ea$group)) {
    subdt <- dt[var_ea[group == group_id], on = .(rsid, EA)]
    wide <- sdcast(subdt, rsid + EA ~ compName, value.var="weight", fill=0)
    wide <- wide[varinfo[, .(rsid)], on = .(rsid), nomatch=0] # row order to match that of genotype data.
    fwrite(wide, sep="\t", quote=FALSE, file=sprintf("%s_group%s.txt", file_prefix, group_id))
    sfinfo[model == model_name & group == group_id, ncol := ncol(wide)]
  }
  return(sfinfo)
}

plink_input_info <- data.table()
if ("is_dom" %in% names(scores)) {
  dominant <- score[(is_dom)]
  scores <- scores[!(is_dom)]

  if (nrow(dominant) > 0) {
    plink_input_info <- rbind(fill=TRUE, plink_input_info, 
        make_score_files(dominant, "dominant", sprintf("%s/collated_scores_dominant_chr%s", work_dir, chrIdx)))
  }

  rm(dominant)
}

if ("is_rec" %in% names(scores)) {
  recessive <- score[(is_rec)]
  scores <- scores[!(is_rec)]

  if (nrow(recessive) > 0) {
    plink_input_info <- rbind(fill=TRUE, plink_input_info, 
        make_score_files(recessive, "recessive", sprintf("%s/collated_scores_recessive_chr%s", work_dir, chrIdx)))
  }

  rm(recessive)
}


if (nrow(scores) > 0) {
  plink_input_info <- rbind(fill=TRUE, plink_input_info, 
      make_score_files(scores, "linear", sprintf("%s/collated_scores_linear_chr%s", work_dir, chrIdx)))
  rm(scores)
}

# save score_info and remove uncessary objects in memory while plink2 runs
fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("%s/score_summary_%s.txt", work_dir, chrIdx))
rm(varinfo)
invisible(gc())

##############################################################################################
# Now that we have collated the score files, we can compute the score
##############################################################################################

# For hard call genotype data, mean-imputation of missing alleles only works when > 50 samples
no_mean_imp <- FALSE
if (args[["--genotype-format"]] == "bfile") {
  n_samples <- as.numeric(system(sprintf("wc -l %s | cut -f 1 -d ' '", samplefile), intern=TRUE))
  if (n_samples < 50) {
		warning("Less than 50 samples, no-mean-imputation flag to --score has been set.") 
		no_mean_imp <- TRUE
  }
}

# Calculate the score levels for each model. Also extract subcohort if asked.
for (idx in plink_input_info[,.I]) {
  # If there's only a single genotype file, we want to prevent multiple tasks 
  # reading from the same file.
  if (args[["--single-geno"]]) {
    lockfile = sprintf("%s/00READLOCK-GENO", work_dir)
    while(file.exists(lockfile)) {
      Sys.sleep(5)
    }
    system(sprintf("touch %s", lockfile), wait=TRUE)
  }

  # Get set of variants to extract
  system(sprintf("tail -n +2 %s/collated_scores_%s_chr%s_group%s.txt | cut -f 1 > %s/collated_scores_%s_variants_chr%s_group%s.txt", 
         work_dir, plink_input_info[idx, model], chrIdx, plink_input_info[idx, group],
         work_dir, plink_input_info[idx, model], chrIdx, plink_input_info[idx, group]), wait=TRUE)

  # Calculate the score sums
  cmd <- normalizePath(sprintf("%s/plink2", SoftwareDir))
  cmd[2] <- sprintf("--%s %s/chr%s", args[["--genotype-format"]], work_dir, chrIdx)
  cmd[3] <- sprintf("--out %s/collated_scores_%s_chr%s_group%s", work_dir, plink_input_info[idx, model], chrIdx, plink_input_info[idx, group])
  cmd[4] <- sprintf("--threads %s --memory %s --silent", ncores, args[["--mem"]])
  cmd[5] <- sprintf("--extract %s/collated_scores_%s_variants_chr%s_group%s.txt", work_dir, plink_input_info[idx, model], chrIdx, plink_input_info[idx, group])
  cmd[6] <- sprintf("--score %s/collated_scores_%s_chr%s_group%s.txt", work_dir, plink_input_info[idx, model], chrIdx, plink_input_info[idx, group])
  cmd[7] <- sprintf("'header-read' 'ignore-dup-ids' 'cols=scoresums'")
  if (no_mean_imp) cmd[8] <- paste("no-mean-imputation")
  if (plink_input_info[idx, model] != "linear") cmd[9] <- sprintf("'%s' %s", plink_input_info[idx, model])
  if (plink_input_info[idx, ncol] > 3) cmd[10] <- sprintf("--score-col-nums 3-%s", plink_input_info[idx, ncol])
  if (exists("freqx_ext") && args[["--genotype-format"]] == "bfile") 
    cmd[11] <- sprintf("--read-freq %s/chr%s.%s", work_dir, chrIdx, freqx_ext)
  if (!is.null(args[["--keep"]])) cmd[12] <- paste("--keep", args[["--keep"]])
  errcode <- system(paste(na.omit(cmd), collapse=" "), wait=TRUE)

  if (errcode != 0) {
    logfile <- sprintf("%s/collated_scores_%s_chr%s_group%s", work_dir, plink_input_info[idx, model], chrIdx, plink_input_info[idx, group])
    score_info[, error := sprintf("%sError from plink2 when computing scores, see log file: %s", 
      ifelse(is.na(error) | error == "", "", paste0(error, ".")), logfile)]
    fwrite(score_info, sep="\t", quote=FALSE, file=sprintf("%s/score_summary_%s.txt", work_dir, chrIdx))
    system(sprintf("touch %s/00FIN-chr-%s", work_dir, chrIdx), wait=TRUE)
    quit(save="no")
  }
 
  # Release readlock
  if (args[["--single-geno"]]) {
    system(sprintf("rm -f %s", lockfile), wait=TRUE)
  }
}

# Create file to let others know this chromosome has finished processing
system(sprintf("touch %s/00FIN-chr-%s", work_dir, chrIdx), wait=TRUE)

########################################################################################
# Collate the results and move to target output directory(s)
########################################################################################
finfiles <- list.files(path=work_dir, pattern="00FIN-chr-*")
if (length(finfiles) <  taskMax) {
  quit(save="no") # other tasks still running, finish.
}

lockfile <- sprintf("%s/00LOCK-collate", work_dir)
if (!file.exists(lockfile)) {
  # If two or more tasks reach this point at the same time, each one sleeps for a different
  # amount of time then checks for lock file again to prevent multiple tasks doing the collation
  # step and messing up the output.
  Sys.sleep(taskIdx) 
}

if (file.exists(lockfile)) {
  quit(save="no") # Another task that finished at the same time and has taken over the collation process
}

# Create lock file and enter the task ID for debugging purposes.
cat(taskIdx, "\n", file=lockfile) 

# remove no longer needed 00FIN files
for (chr in c(1:22, "X", "Y", "XY", "MT")) {
  system(sprintf("rm -f %s/00FIN-chr-%s", work_dir, chr), wait=TRUE)
}

# Collate the plink logs
logfile <- sprintf("%s/collated_plink_logs.txt", work_dir)
appendifexists <- function(path) {
  if (!file.exists(logfile)) {
    tryCatch({
			system(sprintf("touch %s", logfile))
    }, warning=function(w) {
      stop("Could not write to working directory ", work_dir)
    })
  }
  if (file.exists(path)) {
    tryCatch({
			system(sprintf("cat %s >> %s", path, logfile), wait=TRUE)
    }, warning=function(w) {
      stop("Write failure when collating plink2 log files in working directory ", work_dir)
    })
    system(sprintf("rm -f %s", path), wait=TRUE)
  }
}
system(sprintf("touch %s/collated_plink_logs.txt", work_dir), wait=TRUE)
for (chr in c(1:22, "X", "Y", "XY", "MT")) {
  appendifexists(sprintf("%s/ambig_freqx_extract_chr%s.log", work_dir, chr))
}
scoringlogs <- list.files(path=work_dir, pattern="collated_scores_.*.log", full.names=TRUE)
for (ff in scoringlogs) {
  system(sprintf("cat %s >> %s/collated_plink_logs.txt", ff, work_dir), wait=TRUE)
  system(sprintf("rm -f %s", ff), wait=TRUE)
}

# Collate a list of matched variants
rbindu <- function(...) { unique(rbind(...)) }
for(chr in c(1:22, "X", "Y", "XY", "MT")) {
  # Load subset of variants that were matched to any score
  varmatchfiles <- list.files(path=work_dir, pattern=sprintf("collated_scores_.*_variants_chr%s_.*.txt", chr), full.names=TRUE)
  if (length(varmatchfiles) == 0) next # nothing on this chromosome
  varmatch <- foreach(ff = varmatchfiles, .combine=rbindu) %do% {
	  sfread(ff, header=FALSE)
  }
  setnames(varmatch, "match_id")

  # Load full variant information for this chromosome
  if (args[["--genotype-format"]] == "pfile") {
    varinfo <- sfread(sprintf("%s/chr%s.pvar", work_dir, chr))
    setnames(varinfo, c("chromosome", "position", "match_id", "ref_allele", "alt_allele"))
  } else if (args[["--genotype-format"]] == "bfile") {
    varinfo <- sfread(sprintf("%s/chr%s.bim", work_dir, chr))
    setnames(varinfo, c("chromosome", "match_id", "centimorgan", "position", "minor_allele", "major_allele"))
  }

  
  # Obtain the rsid in the original data
  origfile <- paste0(args[["--genotype-prefix"]], ifelse(args[["--single-geno"]], "", chr), args[["--genotype-suffix"]])  
  if (args[["--genotype-format"]] == "pfile") {
    origfile <- sprintf("%s.pvar", origfile)
    line1 <- sreadLines(origfile, 1)
		if (grepl("^##", line1)) {
			varinfo[, rsid := sfread(origfile, cmd="grep -v '^#' %s", sep="\t")[["V3"]]]
    } else {
			varinfo[, rsid := sfread(origfile, cmd="tail -n +2 %s | cut -f 3", header=FALSE)[[1]]]
    }
  } else if (args[["--genotype-format"]] == "bfile") {
    varinfo[, rsid := sfread(origfile, cmd="cut -f 2 %s.bim", header=FALSE)[[1]]]
  }

  # Obtain the effect allele frequency, if applicable
  denovofreqxfile <- sprintf("%s/ambig_freqx_extract_chr%s.afreq", work_dir, chr)
  if (file.exists(denovofreqxfile)) {
    freqx <- sfread(denovofreqxfile)
    freqx <- freqx[, .(match_id=ID, alt_frequency=ALT_FREQS)]
    varinfo[freqx, on = .(match_id), alt_frequency := i.alt_frequency]
  } else if (exists("freqx_ext")) {
    freqx <- sfread(sprintf("%s/chr%s.%s", work_dir, chr, freqx_ext))
		if (freqx_ext == "frqx") {
			freqx <- freqx[, .(match_id=SNP, MAF=(`C(HOM A1)`*2+`C(HET)`)/(`C(HOM A1)`*2+`C(HET)`+`C(HOM A2)`*2))]
      varinfo[freqx, on = .(match_id), MAF := i.MAF]
		} else if (freqx_ext == "frq") {
			freqx <- freqx[, .(match_id=SNP, MAF)]
      varinfo[freqx, on = .(match_id), MAF := i.MAF]
		} else if (freqx_ext == "afreq") {
			freqx <- freqx[, .(match_id=ID, alt_frequency=ALT_FREQS)]
      varinfo[freqx, on = .(match_id), alt_frequency := i.alt_frequency]
		}
  }
  
  # Extract the subset of matched variants
  varinfo <- varinfo[varmatch, on = .(match_id)]
  
  # reorganise columns (not all will be present, hence intersect())
  varinfo <- varinfo[, intersect(c("match_id", "chromosome", "position", 
     "alt_allele", "ref_allele", "minor_allele", "major_allele", "rsid", 
     "alt_frequency", "MAF"), 
     names(varinfo)), with=FALSE]

  # write out
  matchvarfile <- sprintf("%s/matched_variants.txt", work_dir)
  if (!file.exists(matchvarfile)) {
    cat("# This file contains variants that matched the genotype data across all input score files\n", file=matchvarfile)
    cat("# which may include scores whose levels are not saved in this folder (see run_score_file_list.txt)\n", file=matchvarfile, append=TRUE)
    cat("# Note different input scores may use either pair of alleles as the effect allele,\n", file=matchvarfile, append=TRUE)
    cat("# and may be oriented to the opposite strand.\n", file=matchvarfile, append=TRUE)
    fwrite(varinfo[0], sep="\t", quote=FALSE, file=sprintf("%s/header.txt", work_dir))
    system(sprintf("cat %s/header.txt >> %s", work_dir, matchvarfile), wait=TRUE)
    system(sprintf("rm -f %s/header.txt", work_dir), wait=TRUE)
  }
  fwrite(varinfo, sep="\t", quote=FALSE, append=TRUE, file=matchvarfile)
 
  # clean up big objects before next loop
  rm(varinfo, varmatch)
  if(exists("freqx")) rm(freqx)
  invisible(gc())

  # Remove files no longer needed
  system(sprintf("rm -f %s", denovofreqxfile), wait=TRUE) 
  if (exists("freqx_ext")) {
    system(sprintf("rm -f %s/chr%s.%s", work_dir, chr, freqx_ext), wait=TRUE)
    system(sprintf("rm -f %s/ambig_freqx_extract_chr%s", work_dir, chr), wait=TRUE)
  }
  if (args[["--genotype-format"]] == "pfile") {
    system(sprintf("rm -f %s/chr%s.pvar", work_dir, chr), wait=TRUE)
    system(sprintf("rm -f %s/chr%s.psam", work_dir, chr), wait=TRUE)
    system(sprintf("rm -f %s/chr%s.pgen", work_dir, chr), wait=TRUE)
  } else if (args[["--genotype-format"]] == "bfile") {
    system(sprintf("rm -f %s/chr%s.bim", work_dir, chr), wait=TRUE)
    system(sprintf("rm -f %s/chr%s.fam", work_dir, chr), wait=TRUE)
    system(sprintf("rm -f %s/chr%s.bed", work_dir, chr), wait=TRUE)
  }
  for (model in c("linear", "dominant", "recessive")) {
    system(sprintf("rm -f %s/collated_scores_%s_chr%s_*.txt", work_dir, model, chr))
    system(sprintf("rm -f %s/collated_scores_%s_variants_chr%s_*.txt", work_dir, model, chr))
  }
}
system(sprintf("gzip %s/matched_variants.txt", work_dir), wait=TRUE)


# Load all the score sum files and total across chromosomes and models. We need to 
# progressively load and sum to avoid having to load all scores across all chromosomes
# into memory at once. 
sscorefiles <- list.files(path=work_dir, pattern="*.sscore", full.names=TRUE)
sscores <- sfread(sscorefiles[1])
setnames(sscores, names(sscores)[1], "IID")
setnames(sscores, names(sscores)[-1], gsub("_SUM", "", names(sscores)[-1]))
cn <- names(sscores)[-1] # keep track of scores loaded so far (may not be all if score has no variants on this chromosome)
integer_scores <- setdiff(which(sapply(sscores, class) == "integer"), 1L)
for (ic in integer_scores) sscores[[ic]] <- as.numeric(sscores[[ic]])
sscores <- melt(sscores, id.vars=c("IID"), variable.name="compName", value.name="score_sum")
if (length(sscorefiles) > 1) {
  for (sidx in 2:length(sscorefiles)) {
    sscores2 <- sfread(sscorefiles[sidx])
    setnames(sscores2, names(sscores2)[1], "IID")
    setnames(sscores2, names(sscores2)[-1], gsub("_SUM", "", names(sscores2)[-1]))
    cn_new <- setdiff(names(sscores2)[-1], cn)
		integer_scores <- setdiff(which(sapply(sscores2, class) == "integer"), 1L)
		for (ic in integer_scores) sscores2[[ic]] <- as.numeric(sscores2[[ic]])
    sscores2 <- melt(sscores2, id.vars=c("IID"), variable.name="compName", value.name="score_sum")
    sscores[sscores2, on = .(compName, IID), score_sum := .(score_sum + i.score_sum)]
    sscores <- rbind(sscores, sscores2[compName %in% cn_new]) # need to add in scores that did not have variants on any previous chromosome
    cn <- c(cn, cn_new)
    rm(sscores2)
    invisible(gc())
  }
}

# Load all score_info files and collate
sinfofiles <- list.files(path=work_dir, pattern="score_summary_*")
score_info <- lapply(sprintf("%s/%s", work_dir, sinfofiles), sfread, colClasses=c("compName"="character", "error"="character"))
names(score_info) <- gsub("score_summary_", "", gsub(".txt", "", sinfofiles))
score_info <- rbindlist(score_info, idcol="chromosome", fill=TRUE, use.names=TRUE)

# Collate numbers that are always present
nasum <- function(...) { base::sum(..., na.rm=TRUE) }
collated_info <- score_info[, .(n_chr=nasum(n_chr), n_match=nasum(n_match), n_multiallele=nasum(n_multiallele),
                                n_multiallele_removed=nasum(n_multiallele_removed),
                                n_ambig=nasum(n_ambig), n_ambig_kept=nasum(n_ambig_kept), n_used=nasum(n_used)),
                            by=.(path, compName, rsid, chr, pos, EA, EAF, OA, weight, is_dom, is_rec, n_var)]

# Collate numbers that are sometimes present:
if ("n_interaction_skipped" %in% names(score_info)) {
  col2 <- score_info[, .(n_interaction_skipped=nasum(n_interaction_skipped)), by=compName]
  collated_info <- collated_info[col2, on = .(compName)]
  collated_info <- collated_info[, .(path, compName, rsid, chr, pos, EA, EAF, OA, weight, is_dom, is_rec,
                                     n_var, n_chr, n_interaction_skipped, n_match, n_ambig, n_ambig_kept,
                                     n_used)]
} 

# Collate errors: first those that are the same across all chromosomes (should always be the case,
# unless theres' been I/O blocking shenanigans (like the disk quota running out partway through a run).
all_errors <- score_info[,.(error=unique(error)), by=compName]
n_errors <- all_errors[,.(N=length(unique(error))), by=compName]
errors <- all_errors[n_errors[N == 1, .(compName)], on=.(compName)]

# Collate errors that are different on each chromosome into a single error message per score
diff_errors <- score_info[n_errors[N > 1, .(compName)], on = .(compName), .(compName, chromosome, error)][!is.na(error)]
diff_errors <- diff_errors[, .(msg = sprintf("For chromosomes %s:", paste(chromosome, collapse=", "))), by=.(compName, error)]
diff_errors <- diff_errors[, .(error = paste(sprintf("%s %s.", msg, error), collapse=" ")), by=.(compName)]
errors[diff_errors, on=.(compName), error := i.error] # Add to global error table

# Collate warnings generated when the score has variants on a chromosome for which there
# was no genotype data (or where that file could not be opened).
if ("chr_fail" %in% names(score_info)) {
  chr_fail <- score_info[(chr_fail), 
    .(error = sprintf("Variants on chromosome(s) %s not counted due to failure to locate or open genotype data for those chromosome(s).", 
                      paste(chromosome, collapse=", "))),
     by=compName]

  errors[chr_fail, on = .(compName), error := paste(error, i.error, sep=". ")]
  errors[, error := gsub("^NA. ", "", error)]
}

# Add errors to collated_info
collated_info[errors, on = .(compName), error := error]

# clean up a bit
rm(score_info, all_errors, n_errors, errors, diff_errors)
if (exists("col2")) rm(col2)
if (exists("chr_fail")) rm(chr_fail)
invisible(gc())

# Give a name to each score that we may use in the output
collated_info[, name := gsub("\\.txt(\\..*)?", "", basename(path))]
collated_info[as.integer(compName) != as.numeric(compName), name := weight]

# We also want the integer version of the compName to handle multi-score files
collated_info[, compNameInt := as.integer(compName)]

# Determine path for score file output(s)
if (is.null(args[["--out"]]) && is.null(args[["--single-out"]])) {
  # The default. The levels of each score are stored in the same location in which the score-file are
  # stored.
  if (is.null(args[["--cohort-name"]])) {
    collated_info[, outpath := sprintf("%s/sample_levels", dirname(path))]
    fname <- collated_info[compName == compNameInt][,.(compNameInt, outpath = sprintf("%s/%s.sscore.gz", outpath, name))]
    collated_info[fname, on = .(compNameInt), outpath := i.outpath]
    rm(fname)
  } else {
    collated_info[, outpath := sprintf("%s/%s_sample_levels", dirname(path), args[["--cohort-name"]])]
    fname <- collated_info[compName == compNameInt][,.(compNameInt, outpath = sprintf("%s/%s_%s.sscore.gz", outpath, name, args[["--cohort-name"]]))]
    collated_info[fname, on = .(compNameInt), outpath := i.outpath]
    rm(fname)
  }

  # if there are multiple score files in that directory create new sub folders named for the score.
  n_files <- collated_info[compName == compNameInt][, .(N=length(list.files(path=dirname(path), pattern="*\\..*"))), by=.(compNameInt)]
  newpath <- collated_info[compName == compNameInt][n_files[N > 1], on=.(compNameInt), 
     .(compNameInt, outpath = sprintf("%s/%s/%s", dirname(outpath), name, basename(outpath)))]
  collated_info[newpath, on=.(compNameInt), outpath := i.outpath]
} else if (!is.null(args[["--out"]]) && is.null(args[["--single-out"]])) {
  # Single root output directory given, but still want to store each input score-file separately 
  if (is.null(args[["--cohort-name"]])) {
    collated_info[, outpath := sprintf("%s/sample_levels", args[["--out"]])]
    fname <- collated_info[compName == compNameInt][,.(compNameInt, outpath = sprintf("%s/%s/%s.sscore.gz", outpath, name, name))]
    collated_info[fname, on = .(compNameInt), outpath := i.outpath]
    rm(fname)
  } else {
    collated_info[, outpath := sprintf("%s/%s_sample_levels", args[["--out"]], args[["--cohort-name"]])]
    fname <- collated_info[compName == compNameInt][,.(compNameInt, outpath = sprintf("%s/%s/%s_%s.sscore.gz", outpath, name, name, args[["--cohort-name"]]))]
    collated_info[fname, on = .(compNameInt), outpath := i.outpath]
    rm(fname)
  }
} else if (is.null(args[["--out"]]) && !is.null(args[["--single-out"]])) {
  # Single output file requested, but output directory not given (default to directory --score-file is in)
  if (args[["--type"]] == "d") {
    root_dir <- args[["--score-file"]]
  } else {
    root_dir <- dirname(args[["--score-file"]])
  }
  if (is.null(args[["--cohort-name"]])) {
    collated_info[, outpath := sprintf("%s/sample_levels/%s.sscore.gz", root_dir, args[["--single-out"]])]
  } else {
    collated_info[, outpath := sprintf("%s/%s_sample_levels/%s.sscore.gz", root_dir, args[["--cohort-name"]], args[["--single-out"]])]
  }
} else {
  # both --out and --single-out given
  collated_info[, outpath := sprintf("%s/%s.sscore.gz", args[["--out"]], args[["--single-out"]])]
}

# if there are files with only one score, we will name these as "score_sum" in the score-file output.
n_scores <- collated_info[,.N,by=outpath]
collated_info[n_scores[N == 1], on = .(outpath), name := "score_sum"]

# conversely, if there are any duplicated names (for example when requesting a single output file),
# make them unique:
dup_names <- collated_info[weight != "m",.N,by=.(name, outpath)][N > 1]
collated_info[dup_names, on = .(name, outpath), name := paste0(name, ".", seq_len(.N)), by=.(name, outpath)]

# rename some of the columns and drop any unused ones:
setnames(collated_info, c("rsid", "chr", "pos", "EA", "EAF", "OA", "weight", "is_dom", "is_rec"),
  c("rsid_column", "chr_column", "pos_column", "effect_allele_column", "EAF_column", "other_allele_column",
    "weight_column", "is_dominant_column", "is_recessive_column"))
collated_info[, compNameInt := NULL]
setnames(collated_info, "name", "score_name")
setnames(collated_info, "path", "score_path")
collated_info <- collated_info[, which(!sapply(collated_info, function(x) { all(is.na(x)) })), with=FALSE] # drop all columns that are filled with all NA
if (!("error" %in% names(collated_info))) collated_info[, error := NA_character_] # need to keep this column for later errors
collated_info[weight_column == "m", c("weight_column", "score_name") := .(NA, "(Totals for input multi-score):")]
collated_info <- collated_info[, c("score_path", "score_name", setdiff(names(collated_info), c("score_path", "score_name"))), with=FALSE]

# At this point, we can save the collated files (temporarily) and clean up the working directory a bit
sfwrite(collated_info, sep="\t", quote=FALSE, file=sprintf("%s/score_summary.txt", work_dir))
for (chr in c(1:22, "X", "Y", "XY", "MT")) {
  system(sprintf("rm -f %s/score_summary_%s.txt", work_dir, chr), wait=TRUE)
}
sfwrite(sscores, sep="\t", quote=FALSE, compress="gzip", file=sprintf("%s/collated_scores.sscore.gz", work_dir))
system(sprintf("rm -f %s/collated_scores_*_chr*.sscore", work_dir), wait=TRUE)

# We can now remove the args.rds and lock file
system(sprintf("rm -f %s/args.rds", work_dir), wait=TRUE)
system(sprintf("rm -f %s", lockfile))

# if there is only a single output requested, and its the same
# as the working directory then all we need to do is dcast the
# score file and clean up.
outfiles <- unique(collated_info$outpath)
if (length(outfiles) == 1 && normalizePath(dirname(outfiles), mustWork=FALSE) == normalizePath(work_dir)) {
  outfile <- unique(outfiles)
  out_dir <- dirname(outfile)
  # give the scores their names
  sscores <- sscores[collated_info[, .(compName, score_name)], on = .(compName), nomatch=0, .(IID, score_name, score_sum)]
  sscores <- dcast(sscores, IID ~ score_name, value.var="score_sum")
  sscores <- sscores[, intersect(c("IID", collated_info$score_name), names(sscores)), with=FALSE] # preserve score order
  tryCatch({
    fwrite(sscores, sep="\t", quote=FALSE, compress="gzip", file=outfile)
  }, error = function(e) {
    collated_info[, error := paste0(error, ". Error when trying to write wide-format collated scores to", outfile)]
  })
  if (basename(outfile) != "collated_scores.sscore.gz") {
		system(sprintf("rm -f %s/collated_scores.sscore.gz", out_dir), wait=TRUE)
  }

  # Clean up collated info
  dropnames <- c("compName", "outpath", "score_fail")
  if (collated_info[,all(is.na(error))]) dropnames <- c(dropnames, "error")
  if (collated_info[,all(score_name == "score_sum")]) dropnames <- c(dropnames, "score_name")
  fwrite(collated_info[, setdiff(names(collated_info), dropnames), with=FALSE],
         sep="\t", quote=FALSE, file=sprintf("%s/score_summary.txt", out_dir))

  # Copy across information about scores run:
  if (args[["--type"]] == 'l') {
    system(sprintf("cp %s %s/run_score_file_list.txt", args[["--score-file"]], out_dir))
  } else {
    cat(args[["--score-file"]], "\n", file=sprintf("%s/run_score_file_list.txt", out_dir))
  }

  # Copy this script across to freeze version
	r_cmd <- commandArgs()
	sfile <- grep("--file", r_cmd, value=TRUE)
	sfile <- gsub("--file=", "", sfile)
	system(sprintf("cp %s %s", sfile, out_dir), wait=TRUE)
  
  # Freeze files:
  system(sprintf("chmod -x %s/calc_PS_lvls.*", out_dir), wait=TRUE)
  system(sprintf("chmod -w %s", out_dir), wait=TRUE)

  quit(save="no")
}

# For each output directory, attempt to:
#
# (1) save the sscore file
# (2) save the relevant subset of the information file
# (3) copy the collated plink logs
#
# We only do this if the file does not exist, otherwise we preserve
# the working directory and error out.
for (outIdx in unique(collated_info$outpath)) {
  out_dir <- dirname(outIdx)
  if (dir.exists(out_dir)) {
    collated_info[outpath == outIdx, error := paste0(error, ". Output directory ", out_dir, " already exists!")] 
    next
  }
  tryCatch({
		dir.create(out_dir, recursive=TRUE)
  }, error = function(e) {
    collated_info[outpath == outIdx, error := paste0(error, ". Could not create output directory ", out_dir, "!")]
    next
  })
 
  # Extract subset of scores to write out and cast to wide format
  this_info <- collated_info[outpath == outIdx]
  this_info[, score_fail := FALSE]
  this_info[!(compName %in% unique(sscores$compName)) || is.na(weight_column), score_fail := TRUE]
  this_sscore <- sscores[this_info[, .(compName, score_name)], on = .(compName), nomatch=0, .(IID, score_name, score_sum)]
  this_sscore <- dcast(this_sscore, IID ~ score_name, value.var="score_sum")
  this_sscore <- this_sscore[, intersect(c("IID", collated_info$score_name), names(this_sscore)), with=FALSE] # preserve score order

  tryCatch({
    fwrite(this_sscore, sep="\t", quote=FALSE, compress="gzip", file=outIdx)
  }, error = function(e) {
    collated_info[outpath == outIdx, error := paste0(error, ". Error when trying to write collated scores to", outIdx)]
    next
  })
  
  # Remove successfull transferred scores from the collated scores
  sscores <- sscores[!this_info[, .(compName)], on = .(compName)]

  # Write out subset of the score info
  this_info <- this_info[!(score_fail)]
  tryCatch({
    dropnames <- c("compName", "outpath", "score_fail")
    if (this_info[,all(is.na(error))]) dropnames <- c(dropnames, "error")
    if (this_info[,all(score_name == "score_sum")]) dropnames <- c(dropnames, "score_name")
    fwrite(this_info[, setdiff(names(this_info), dropnames), with=FALSE],
           sep="\t", quote=FALSE, file=sprintf("%s/score_summary.txt", out_dir))
  }, error = function(e) {
    collated_info[outpath == outIdx, error := paste0(error, ". Error when trying to write out score information to ", out_dir, "/score_summary.txt")]
    next
  })
  
  # Copy across shared files and logs
  tryCatch({
    copy_files <- c("collated_plink_logs.txt", "slurm_logs/", "calc_PS_lvls.sh", "command_log.txt", "matched_variants.txt.gz")
    copy_cmds <- sprintf("cp -r %s/%s %s", work_dir, copy_files, out_dir)
    invisible(sapply(copy_cmds, system, wait=TRUE))
  }, error=function(e) {
    collated_info[outpath == outIdx, error := paste0(error, ". Error when trying to copy ", paste(copyfiles, collapse=","), " from working directory ", work_dir)]
  })

  # Copy across information about scores run:
  tryCatch({
    if (args[["--type"]] == 'l') {
      system(sprintf("cp %s %s/run_score_file_list.txt", args[["--score-file"]], out_dir))
    } else {
      cat(args[["--score-file"]], "\n", file=sprintf("%s/run_score_file_list.txt", out_dir))
    }
  }, error=function(e) {
    collated_info[outpath == outIdx, error := paste0(error, ". Error when trying to log input score files for this run ")]
  })

  # Copy this script across to freeze version
  tryCatch({
		r_cmd <- commandArgs()
		sfile <- grep("--file", r_cmd, value=TRUE)
		sfile <- gsub("--file=", "", sfile)
		system(sprintf("cp %s %s", sfile, out_dir), wait=TRUE)
  }, error=function(e) {
    collated_info[outpath == outIdx, error := paste0(error, ". Error when trying to copy this program from ", sfile)]
  })
  
  # Freeze files:
  system(sprintf("chmod -x %s/calc_PS_lvls.*", out_dir), wait=TRUE)
  system(sprintf("chmod -w %s", out_dir), wait=TRUE)

  # Remove successfully transferred scores from the collated information
  collated_info <- collated_info[!this_info, on=.(compName)]
}
collated_info[, error := gsub("^NA\\. ", "", error)]
collated_info[, error := gsub("^\\. ", "", error)]
collated_info[, error := gsub("^ \\. ", "", error)]

if (nrow(collated_info) > 0) {
  # Update collated scores and score summary files in working dir to just the failed scores
  tryCatch({
    fwrite(collated_info, sep="\t", quote=FALSE, file=sprintf("%s/score_summary.txt", work_dir))
    fwrite(sscores, sep="\t", quote=FALSE, compress="gzip", file=sprintf("%s/collated_scores.sscore.gz", work_dir))
  }, error = function(e) {
    warning("Error updating collated score and information files to include only scores with write errors")
  })
  stop("Some scores had errors: see ", work_dir, "/score_summary.txt")
} 

# Final cleanup
system(sprintf("rm -f %s/collated_scores.sscore.gz", work_dir), wait=TRUE)
system(sprintf("rm -f %s/score_summary.txt", work_dir), wait=TRUE)
system(sprintf("rm -f %s/collated_plink_logs.txt", work_dir), wait=TRUE)
system(sprintf("rm -f %s/calc_PS_lvls.sh", work_dir), wait=TRUE)
system(sprintf("rm -f %s/command_log.txt", work_dir), wait=TRUE)
system(sprintf("rm -f %s/matched_variants.txt.gz", work_dir), wait=TRUE)
system(sprintf("rm -rf %s/slurm_logs", work_dir), wait=TRUE) # will this cause the job to fail?
system(sprintf("rmdir %s", work_dir), wait=TRUE)

