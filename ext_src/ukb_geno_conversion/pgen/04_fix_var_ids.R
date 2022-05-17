library(data.table)

out_dir = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/COMMON/post_qc_data/imputed/HRC_UK10K/plink_format/GRCh37/pgen"
ref_dir = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/COMMON/post_qc_data/imputed/reference_files"

chr_id = Sys.getenv("SLURM_ARRAY_TASK_ID")
if (chr_id == "23") chr_id = "X"
if (chr_id == "24") chr_id = "XY"

# Remove extra row identifier tacked onto the variant IDs
pvar = fread(sprintf("%s/ukb_imp_v3_dedup_chr%s.pvar", out_dir, chr_id))
pvar[, ID := gsub("\\.[0-9]*$", "", ID)]
fwrite(pvar, sep="\t", quote=FALSE, file=sprintf("%s/ukb_imp_v3_dedup_chr%s.pvar", out_dir, chr_id))

