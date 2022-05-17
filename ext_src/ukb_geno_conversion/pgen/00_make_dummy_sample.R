library(data.table)

ref_dir = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/HRC_UK10K/"
out_dir = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/COMMON/post_qc_data/imputed/HRC_UK10K/plink_format/GRCh37/pgen"

# the X and XY chromosome files have different .sample files, but we want to make just one row to 
# sample ID mapping to the common file
dummy = fread(sprintf("%s/ukb_adiposity_imp_v3.sample", ref_dir))
dummyX = fread(sprintf("%s/ukb7439_imp_chrX_v3_s486645.sample", ref_dir))
dummyXY = fread(sprintf("%s/ukb7439_imp_chrXY_v3_s486445.sample", ref_dir))

dummy[, row := .I-1]
dummyX[dummy, on = .(ID_1, ID_2), row := i.row]
dummyXY[dummy, on = .(ID_1, ID_2), row := i.row]

dummy[ID_1 > 0, c("ID_1", "ID_2") := row]
dummy[, row := NULL]

dummyX[ID_1 > 0, c("ID_1", "ID_2") := row]
dummyX[, row := NULL]

dummyXY[ID_1 > 0, c("ID_1", "ID_2") := row]
dummyXY[, row := NULL]

# write out dummy sample files
fwrite(dummy, file=sprintf("%s/dummy.sample", out_dir), sep="\t", quote=FALSE)
fwrite(dummyX, file=sprintf("%s/dummyX.sample", out_dir), sep="\t", quote=FALSE)
fwrite(dummyXY, file=sprintf("%s/dummyXY.sample", out_dir), sep="\t", quote=FALSE)

# Write out sample exclusion lists
fwrite(dummy[ID_1 < 0, .(`#FID`=ID_1, IID=ID_2)], sep="\t", quote=FALSE, file=sprintf("%s/dummy.remove", out_dir))
fwrite(dummyX[ID_1 < 0, .(`#FID`=ID_1, IID=ID_2)], sep="\t", quote=FALSE, file=sprintf("%s/dummyX.remove", out_dir))
fwrite(dummyXY[ID_1 < 0, .(`#FID`=ID_1, IID=ID_2)], sep="\t", quote=FALSE, file=sprintf("%s/dummyXY.remove", out_dir))

