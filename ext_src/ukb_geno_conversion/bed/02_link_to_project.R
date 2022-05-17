# Maps to project specific sample IDs and outputs traits
# in relevant folders
library(data.table)

rds_dir="/rds/project/asb38/rds-asb38-ceu-ukbiobank"
common_dir=sprintf("%s/genetics/COMMON/post_qc_data/imputed/HRC_UK10K/plink_format/GRCh37/bed", rds_dir)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 0 && length(args) != 2) {
  stop("Must provide the project number (starting with 'P') and path to autosomal .sample file as arguments")
}
if (length(args) == 0) {
  projects = c("P7439", "P20480")
  sample_files = c(
    "P7439"=sprintf("%s/genetics/P7439/post_qc_data/imputed/HRC_UK10K/ukb_adiposity_imp_v3.sample", rds_dir),
    "P20480"=sprintf("%s/genetics/P20480/post_qc_data/imputed/HRC_UK10K/ukb_BP_imp_v3.sample", rds_dir)
  )
} else {
  projects = args[1]
  sample_files = args[2]
  names(sample_files) = args[1]
}

for (proj_id in projects) {
  target_dir = sprintf("%s/genetics/%s/post_qc_data/imputed/HRC_UK10K/plink_format/GRCh37/bed", rds_dir, proj_id)
  dir.create(target_dir, recursive=TRUE, showWarnings=FALSE)
  sample = fread(sample_files[proj_id])[-1]
  sample[, row := .I]
  for (chr_id in c(1:22, "X", "XY")) {
    system(sprintf("ln -s %s/ukb_imp_v3_dedup_chr%s.bed %s/", common_dir, chr_id, target_dir), wait=TRUE)
    system(sprintf("ln -s %s/ukb_imp_v3_dedup_chr%s.bim %s/", common_dir, chr_id, target_dir), wait=TRUE)
    system(sprintf("ln -s %s/ukb_imp_v3_dedup_chr%s.gcount %s/", common_dir, chr_id, target_dir), wait=TRUE)
    fam = fread(sprintf("%s/ukb_imp_v3_dedup_chr%s.fam", common_dir, chr_id), header=FALSE)
    setnames(fam, c("FID", "IID", "PID", "MID", "SEX", "PHEN"))
    fam[sample, on = .(IID=row), c("FID", "IID") := .(ID_1, ID_2)]
    fwrite(fam, sep="\t", quote=FALSE, col.names=FALSE, file=sprintf("%s/ukb_imp_v3_dedup_chr%s.fam", target_dir, chr_id))
  }
  system(sprintf("ln -s %s/README.txt %s", common_dir, target_dir), wait=TRUE)
}

