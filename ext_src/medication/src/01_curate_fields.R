library(data.table)

system("mkdir -p data/raw/ukbiobank/extracted/", wait=TRUE)

info <- rbind(use.names=FALSE,
  data.table(field.id=6177, name="Medication for cholesterol, blood pressure or diabetes"),
  data.table(field.id=6153, name="Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones"),
  data.table(field.id=2492, name="Taking other prescription medications"),
  data.table(field.id=137, name="Number of treatments/medications taken"),
  data.table(field.id=20003, name="Treatment/medication code"),
  data.table(field.id=6671, name="Number of antibiotics taken in last 3 months"),
  data.table(field.id=20199, name="Antibiotic codes for last 3 months")
)

fwrite(info, sep="\t", quote=FALSE, file="data/raw/ukbiobank/extracted/field_info.txt")
fwrite(info[,.(field.id)], quote=FALSE, col.names=FALSE, file="data/raw/ukbiobank/extracted/fields.txt")
