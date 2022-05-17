library(data.table)

system("mkdir -p data/raw/ukbiobank/extracted/", wait=TRUE)

info <- rbind(use.names=FALSE,
  data.table(field.id=54, name="Assessment centre", var="assessment_centre"),
  data.table(field.id=53, name="Date of attending assessment centre", var="assessment_date"),
  data.table(field.id=189, name="Townsend deprivation index at recruitment", var="townsend"),
  data.table(field.id=31, name="Sex", var="sex"),
  data.table(field.id=21003, name="Age at assessment", var="age"),
  data.table(field.id=34, name="Year of birth", var="birth_year"),
  data.table(field.id=52, name="Month of birth", var="birth_month"),
  data.table(field.id=21002, name="Weight", var="weight"),
  data.table(field.id=50, name="Standing height", var="height"),
  data.table(field.id=21001, name="Body Mass Index", var="bmi"),
  data.table(field.id=48, name="Waist circumference", var="waist"),
  data.table(field.id=49, name="Hip circumference", var="hip"),
  data.table(field.id=21000, name="Ethnicity", var="ethnicity"),
  data.table(field.id=3140, name="Pregnant", var="pregnant")
)

fwrite(info, sep="\t", quote=FALSE, file="data/raw/ukbiobank/extracted/field_info.txt")
fwrite(info[,.(field.id)], quote=FALSE, col.names=FALSE, file="data/raw/ukbiobank/extracted/fields.txt")
