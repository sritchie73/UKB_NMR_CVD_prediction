library(data.table)

system("mkdir -p data/raw/ukbiobank/extracted/", wait=TRUE)

info <- rbind(use.names=FALSE,
  data.table(field.id=20107, name="Illnesses of father", var="father"),
  data.table(field.id=20110, name="Illnesses of mother", var="mother"),
  data.table(field.id=20111, name="Illnesses of siblings", var="siblings")
)

fwrite(info, sep="\t", quote=FALSE, file="data/raw/ukbiobank/extracted/field_info.txt")
fwrite(info[,.(field.id)], quote=FALSE, col.names=FALSE, file="data/raw/ukbiobank/extracted/fields.txt")
