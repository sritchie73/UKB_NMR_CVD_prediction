#!/usr/bin/env bash

src_dir=$HOME/rds/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/GRS_resources/

eval "$($src_dir/software_dependencies/docopts -h - : "$@" <<EOF
Calculate the levels of a polygenic score in a group of samples

By default, this program calculates score levels in UK Biobank
participants using the phase 3 release of UK Biobank's genotype 
data (stored as probablistic dosages in plink 2 binary format)
by matching variants by chromosome, position, and alleles between
the genotype data (build GRCh37/hg19) and the score file (must
be the same genome build). Variants with ambiguous alleles (A/T or
G/C SNPs are excluded with warning). Please read through the list
of options to modify these default behaviours.

This script submits a job to slurm, sensible defaults are provided.
If you are not a member of the inouyelab, you will need to manually
specify the --account to use when submitting to slurm.

NOTE: when calculating multiple scores, please submit one job using
the multiple score functionality. Launching many parallel jobs to
calculate single score files with this program will cause all these
jobs to be very slow as they compete to read from the same genotype
data.

Usage:
  calc_PS_lvls.sh --score-file <file> [options]
  calc_PS_lvls.sh -h | --help

Options:
  -h --help                   Show this screen.
  --score-file <file>         Path to polygenic score file, directory, or file containing list
                              of score files (see --type).
  --type <type>               Type of filepath given to --score-file. 's': filepath points to
                              a single polygenic score, e.g. created by make_simple_PS.sh or
                              downloaded from the PGS catalog. 'd': filepath is a directory 
                              containing multiple score files. 'l': filepath points to a file
                              containing a list of paths, one per line. Each line may optionally
                              be followed by three to ten column position or column name arguments to 
                              specify for each file the --score-X options (described below) to 
                              handle different locations, names, or presence/absence of the respective 
                              fields in the score files. If not provided, the columns for that particular
                              score will be loaded as provided to the respective program arguments.
                              Reasonable combinations can be inferred when listing only a subset of columns 
                              (assuming they are in the same order as documented below) unless you are using 
                              the ambiguous MAF threshold option, or there are columns indicating that some 
                              or all effect weights should be taken as either dominant or recessive. 
                              The score_summary.txt output will always report which columns have been used 
                              for which field. Note that files from the PGS Catalog are automatically
                              detected and handled, so you do not need to specify these manually. [default: s]
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
                              the effect allele frequency (i.e. if using --keep-ambiguous and --ambiguous-thresh).
                              If column is not present, set to 'NULL'. [default: NULL]
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
  --work <directory>          Working directory to store intermediate files and logs shared across all 
                              scores for the duration of the run. The default is to store this in a 
                              directory called 'work/' in the folder in which the input --score-file
                              is scored, then all logs will be copied to each score output folder under
                              'logs/' at the end of the run. [default: NULL]
  --genotype-prefix <prefix>  Path and prefix occurring before the chromosome number for the genotype
                              data for the samples you want to calculate the polygenic score levels in.
                              Defaults to UK Biobank, using the phase 3 release genotype data.
                              [default: $HOME/rds/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/HRC_UK10K/plink_format/GRCh37/pgen/ukb_imp_v3_dedup_chr]
  --genotype-suffix <suffix>  Suffix for the filename occurring after the chromosome number but before the
                              .pgen/.pvar/.psam extension for the genotype data. [default: NULL]
  --single-geno               Flag to indicate that the genotype data is stored as a single file, not split across
                              multiple chromosomes. In this case, you can ignore the --genotype-suffix argument.
  --genotype-format <format>  Format the genotype data is stored in, must correspond to one of the arguments to
                              plink, e.g. the default, 'pfile' is passed directly to plink as '--pfile'. To use
                              plink version 1 binary data (bed/bim/fam) set this as 'bfile'. [default: 'pfile']
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
  --time <runtime>            Passed to sbatch. [default: 6:0:0]
  --account <names>           Passed to sbatch. [default: INOUYE-SL3-CPU]
  --partition <names>         Passed to sbatch. [default: skylake,skylake-himem,cclake,cclake-himem]
  --mem-per-chr <MB>          Memory to allocate to each task (26 tasks, 1 per chromosome).
                              This also dictates the number of cores allocated to each task,
                              as on CSD3 the memory requested dictats the number of CPUs 
                              allocated (12GB = 1 CPU on cardio and skylake-himem partions,
                              6GB = 1CPU on skylake partition). [default: 24000]
  --parsable                  Prints the job id and output directory name to stdout enabling embedding in other scripts.
  --dependency <depends>      Passed to sbatch. Allows you to set dependencies for the slurm job, see the man page for
                              sbatch for usage. [default: afterany:1]
EOF
)"

# Check type switch
if ! [[ $type = 's' || $type = 'd' || $type = 'l' ]]; then
  echo "--type must be one of 's', 'd', or 'l'. See --help." >&2
  exit 1
fi

# Check for logging/working directory
if [[ $work = "NULL" ]]; then
  indir=$(dirname $score_file)
  work=$indir/logs
fi
if [[ -d $work ]]; then
  echo "Working directory $work already exists. Overwrite? (y/n)" 1>&2
  read ans
  while true; do
    if [[ $ans = "y" || $ans = "Y" || $ans = "Yes" || $ans = "YES" || $ans = "yes" ]]; then
      rm -rf $work
      if [[ $? -ne 0 ]]; then
        echo "Could do not overwrite $work" >&2
        exit 1
      fi
      break
    elif [[ $ans = "n" || $ans = "N" || $ans = "No" || $ans = "NO" || $ans = "no" ]]; then
      exit 0
    else
     echo "Unrecognised user input. Please answer 'y' or 'n'." 1>&2
     read ans
    fi
  done
fi
mkdir -p $work/slurm_logs
if [[ $? -ne 0 ]]; then
  echo "Could not create $work" >&2
  exit 1
fi

echo "Working and temporary logging directory is: $work" 1>&2

# Log submitted command
arg_string=$@
echo "Batch command:" > $work/command_log.txt
echo "------------------------------------------------------------------------" >> $work/command_log.txt
echo "$src_dir/calc_PS_lvls.sh $arg_string" >> $work/command_log.txt
echo "" >> $work/command_log.txt

# Copy across this script file
cp $src_dir/calc_PS_lvls.sh $work
chmod -x $work/calc_PS_lvls.sh

# Determine memory to allow plink to use. This is less
# than the specified memory as R, used to dispatch and collate,
# will have some overhead.
plink_mem=$(( $mem_per_chr - 512 )) # reserve 512MB for R
if [[ $plink_mem -lt 0 ]]; then
  echo "Insufficient memory allocated" >&2
  exit 1
fi

# build command string
cmd[0]="Rscript --vanilla $src_dir/calc_PS_lvls.R"
cmd[1]="--score-file $score_file"
cmd[2]="--mem $plink_mem"
cmd[3]="--work $work"
cmd[4]="--type $type"
cmd[5]="--score-rsid $score_rsid"
cmd[6]="--score-chr $score_chr"
cmd[7]="--score-pos $score_pos"
cmd[8]="--score-EA $score_EA"
if [[ $score_EAF != "NULL" ]]; then  cmd[9]="--score-EAF $score_EAF"; fi
cmd[10]="--score-OA $score_OA"
cmd[11]="--score-weight $score_weight"
if [[ $score_dominant != "NULL" ]]; then  cmd[12]="--score-dominant $score_dominant"; fi
if [[ $score_recessive != "NULL" ]]; then  cmd[13]="--score-recessive $score_recessive"; fi
if $match_by_rsid; then  cmd[14]="--match-by-rsid"; fi
cmd[15]="--cohort-name $cohort_name"
if [[ $out != "NULL" ]]; then  cmd[16]="--out $out"; fi
if [[ $single_out != "NULL" ]]; then  cmd[17]="--single-out $single_out"; fi
cmd[18]="--genotype-prefix $genotype_prefix"
if [[ $genotype_suffix != "NULL" ]]; then  cmd[19]="--genotype-suffix $genotype_suffix"; fi
cmd[20]="--genotype-format $genotype_format"
if $single_geno; then  cmd[21]="--single-geno"; fi
if [[ $keep != "NULL" ]]; then  cmd[22]="--keep $keep"; fi
if $keep_ambiguous; then  cmd[23]="--keep-ambiguous"; fi
if [[ $ambiguous_thresh != "NULL" ]]; then  cmd[24]="--ambiguous-thresh $ambiguous_thresh"; fi
if [[ $freqx_prefix != "NULL" ]]; then  cmd[25]="--freqx-prefix $freqx_prefix"; fi
if [[ $freqx_suffix != "NULL" ]]; then  cmd[26]="--freqx-suffix $freqx_suffix"; fi
if $remove_multiallelic; then  cmd[27]="--remove-multiallelic"; fi

cmd_string=${cmd[@]}

# submit to slurm
job=$(sbatch --job-name "calculate polygenic score levels" \
             --parsable --dependency $dependency \
             --kill-on-invalid-dep yes \
             --array 1-26 --mem $mem_per_chr \
						 --time $time --nodes 1 \
						 --partition $partition \
						 --output $work/slurm_logs/slurm-%A_%a.out \
						 --error $work/slurm_logs/slurm-%A_%a.err \
             --account $account \
						 --wrap "$cmd_string")


if $parsable; then
  echo "$job"
else
  echo "Submitted batch job $job" # mimic sbatch behaviour when not given --parsable
fi
