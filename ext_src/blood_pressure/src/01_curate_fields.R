library(data.table)

system("mkdir -p data/raw/ukbiobank/extracted/", wait=TRUE)

info <- rbind(use.names=FALSE,
  data.table(field.id=93, name="manual SBP"),
  data.table(field.id=4080, name="automatic SBP"),
  data.table(field.id=94, name="manual DBP"),
  data.table(field.id=4079, name="automatic DBP"),
  data.table(field.id=95, name="pulse rate (during manual blood pressure readings)"),
  data.table(field.id=102, name="pulse rate (during autmatic blood pressure readings)")
)

fwrite(info, sep="\t", quote=FALSE, file="data/raw/ukbiobank/extracted/field_info.txt")
fwrite(info[,.(field.id)], quote=FALSE, col.names=FALSE, file="data/raw/ukbiobank/extracted/fields.txt")

