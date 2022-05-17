library(data.table)

system("mkdir -p data/raw/ukbiobank/extracted/followup", wait=TRUE)

info <- rbind(use.names=FALSE,
  data.table(field.id=190, name="Reason lost to follow-up", var="lost_to_followup_reason"),
  data.table(field.id=191, name="Date lost to follow-up", var="lost_to_followup_date")
)

fwrite(info, sep="\t", quote=FALSE, file="data/raw/ukbiobank/extracted/followup/field_info.txt")
fwrite(info[,.(field.id)], quote=FALSE, col.names=FALSE, file="data/raw/ukbiobank/extracted/followup/fields.txt")
