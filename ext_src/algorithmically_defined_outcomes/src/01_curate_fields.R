library(data.table)

system("mkdir -p data/raw/ukbiobank/extracted/", wait=TRUE)

info <- rbind(use.names=FALSE,
  data.table(field.id=42014, name="Date of asthma report", var="asthma_date"),
  data.table(field.id=42015, name="Source of asthma report", var="asthma_source"),
  data.table(field.id=42016, name="Date of COPD report", var="copd_date"),
  data.table(field.id=42017, name="Source of COPD report", var="copd_source"),
  data.table(field.id=42018, name="Date of all-cause dementia report", var="dementia_date"),
  data.table(field.id=42019, name="Source of all-cause dementia report", var="dementia_source"),
  data.table(field.id=42020, name="Date of alzheimers report", var="alzheimers_date"),
  data.table(field.id=42021, name="Source of alzheimers report", var="alzheimers_source"),
  data.table(field.id=42022, name="Date of vascular dementia report", var="dementia_vasc_date"),
  data.table(field.id=42023, name="Source of vascular dementia report", var="dementia_vasc_source"),
  data.table(field.id=42024, name="Date of frontotemporal dementia report", var="dementia_ft_date"),
  data.table(field.id=42025, name="Source of frontotemporal dementia report", var="dementia_ft_source"),
  data.table(field.id=42026, name="Date of ESRD report", var="esrd_date"),
  data.table(field.id=42027, name="Source of ESRD report", var="esrd_source"),
  data.table(field.id=42028, name="Date of motor neurone disease report", var="mnd_date"),
  data.table(field.id=42029, name="Source of motor neurone disease report", var="mnd_source"),
  data.table(field.id=42000, name="Date of myocardial infarction report", var="mi_date"),
  data.table(field.id=42001, name="Source of myocardial infarction report", var="mi_source"),
  data.table(field.id=42002, name="Date of STEMI report", var="stemi_date"),
  data.table(field.id=42003, name="Source of STEMI report", var="stemi_source"),
  data.table(field.id=42004, name="Date of NSTEMI report", var="nstemi_date"),
  data.table(field.id=42005, name="Source of NSTEMI report", var="nstemi_source"),
  data.table(field.id=42030, name="Date of all-cause parkinsonism report", var="parkinsonism_date"),
  data.table(field.id=42031, name="Source of all-cause parkinsonism report", var="parkinsonism_source"),
  data.table(field.id=42032, name="Date of parkinsons report", var="parkinsons_date"),
  data.table(field.id=42033, name="Source of parkinsons report", var="parkinsons_source"),
  data.table(field.id=42034, name="Date of progressive supranuclear palsy report", var="palsy_ps_date"),
  data.table(field.id=42035, name="Source of progressive supranuclear palsy report", var="palsy_ps_source"),
  data.table(field.id=42036, name="Date of multiple system atrophy report", var="atrophy_ms_date"),
  data.table(field.id=42037, name="Source of multiple system atrophy report", var="atrophy_ms_source"),
  data.table(field.id=42006, name="Date of stroke report", var="stroke_date"),
  data.table(field.id=42007, name="Source of stroke report", var="stroke_source"),
  data.table(field.id=42008, name="Date of ischaemic stroke report", var="stroke_is_date"),
  data.table(field.id=42009, name="Source of ischaemic stroke report", var="stroke_is_source"),
  data.table(field.id=42010, name="Date of intracerebral haemorrhage report", var="stroke_ih_date"),
  data.table(field.id=42011, name="Source of intracerebral haemorrhage report", var="stroke_ih_source"),
  data.table(field.id=42012, name="Date of subarachnoid haemorrhage report", var="stroke_sh_date"),
  data.table(field.id=42013, name="Source of subarachnoid haemorrhage report", var="stroke_sh_source")
)

fwrite(info, sep="\t", quote=FALSE, file="data/raw/ukbiobank/extracted/field_info.txt")
fwrite(info[,.(field.id)], quote=FALSE, col.names=FALSE, file="data/raw/ukbiobank/extracted/fields.txt")
