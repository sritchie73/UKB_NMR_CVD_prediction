library(data.table)

out_dir = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/COMMON/post_qc_data/imputed/HRC_UK10K/plink_format/GRCh37/pgen"
ref_dir = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/COMMON/post_qc_data/imputed/reference_files"

chr_id = Sys.getenv("SLURM_ARRAY_TASK_ID")
if (chr_id == "23") chr_id = "X"
if (chr_id == "24") chr_id = "XY"

# Load variant information
pvar = fread(sprintf("%s/ukb_imp_v3_chr%s.pvar", out_dir, chr_id))

# Load in variant statistics so we can use INFO scores to flag 
# which of each pair of duplicates to remove
snpstats = fread(sprintf("%s/ukb_impv3_chr%s_snpstats.txt", ref_dir, chr_id), skip=8)
pvar[, INFO := snpstats$impute_info]

# Get sorted combination of alleles to aid duplicate identification
pvar[,row := .I]
pvar[, sorted_alleles := paste(sort(c(REF, ALT)), collapse=":"), by=.(row)]

# Flag variants that are duplicates by rsID, keeping the one with the highest INFO
dup_by_rsid = pvar[, .N, by=.(ID, sorted_alleles)][N > 1]
dup_by_rsid = pvar[dup_by_rsid[,.(ID, sorted_alleles)], on = .(ID, sorted_alleles)]
keep_by_rsid = dup_by_rsid[, .SD[which.max(INFO)], by = .(ID, sorted_alleles)]
rm_by_rsid = dup_by_rsid[!keep_by_rsid, on = .(row)]

# Separately, do the same for variants by position
dup_by_pos = pvar[, .N, by=.(POS, sorted_alleles)][N > 1]
dup_by_pos = pvar[dup_by_pos[,.(POS, sorted_alleles)], on = .(POS, sorted_alleles)]
dup_by_pos = dup_by_pos[!dup_by_rsid, on = .(row)] # in case a duplicated rsID has the wrong position, conflicting with a different variant
keep_by_pos = dup_by_pos[, .SD[which.max(INFO)], by = .(POS, sorted_alleles)]
rm_by_pos = dup_by_pos[!keep_by_pos, on = .(row)]

# Make sure every row in the pvar file has a unique ID so we can accurately flag
# rows to remove
pvar[, ID := paste0(ID, ".", row)]

# Write out pvar file with new IDs
fwrite(pvar[,.(`#CHROM`, POS, ID, REF, ALT)], sep="\t", quote=FALSE,
       file=sprintf("%s/ukb_imp_v3_chr%s.pvar", out_dir, chr_id))

# Write out list of variants to remove:
if (nrow(dup_by_rsid) > 0 && nrow(dup_by_pos) > 0) {
  dups = rbind(dup_by_rsid[, .(ID = paste0(ID, ".", row))],
               dup_by_pos[, .(ID = paste0(ID, ".", row))])
  fwrite(dups, quote=FALSE, col.names=FALSE, file=sprintf("%s/chr%s_duplicates.txt", out_dir, chr_id))
} else if (nrow(dup_by_rsid) > 0) {
  fwrite(dup_by_rsid[, .(ID = paste0(ID, ".", row))], quote=FALSE, col.names=FALSE, file=sprintf("%s/chr%s_duplicates.txt", out_dir, chr_id))
} else if (nrow(dup_by_pos) > 0) {
  fwrite(dup_by_pos[, .(ID = paste0(ID, ".", row))], quote=FALSE, col.names=FALSE, file=sprintf("%s/chr%s_duplicates.txt", out_dir, chr_id))
}

