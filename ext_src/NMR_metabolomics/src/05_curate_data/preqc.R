library(data.table)
library(ukbnmr)

# make output directory
system("mkdir -p data/curated/preqc", wait=TRUE)

# Load and curate biomarker measurements
nmr <- fread("output/basic_participant_QC/biomarker_measurements.txt")
nmr[, visit := ifelse(visit == "Main Phase", 0L, 1L)]
setnames(nmr, "visit", "visit_index")
setnames(nmr, "eid_30418", "eid")
fwrite(nmr, sep="\t", quote=FALSE, file="data/curated/preqc/biomarker_measurements.txt")

# Curate biomarker information
nmr_info <- ukbnmr::nmr_info[(Nightingale), .(Biomarker, Description, Units, Group, Sub.Group, UKB.Field.ID, QC.Flag.Field.ID)]
fwrite(nmr_info, sep="\t", quote=FALSE, file="data/curated/preqc/biomarker_information.txt")

# Load and curate biomarker QC flags
flags <- fread("output/basic_participant_QC/biomarker_QC_tags.txt")
flags[, visit := ifelse(visit == "Main Phase", 0L, 1L)]
setnames(flags, "visit", "visit_index")
setnames(flags, "eid_30418", "eid")
fwrite(flags, sep="\t", quote=FALSE, file="data/curated/preqc/measurement_qc_flags.txt")

# load and curate relevant sample meta-data
sinfo <- fread("output/basic_participant_QC/sample_information.txt")
sinfo <- sinfo[!(removed), .(eid=eid_30418, visit_index=ifelse(visit == "Main Phase", 0L, 1L), 
  Shipment.Plate=as.character(plate_id), 
  Spectrometer=sprintf("Spectrometer %s", as.integer(factor(spectrometer))), 
  High.Lactate=ifelse(tag %like% "High_lactate", "Yes", ""),
  High.Pyruvate=ifelse(tag %like% "High_pyruvate", "Yes", ""),
  Low.Glucose=ifelse(tag %like% "Low_glucose", "Yes", ""),
  Low.Protein=ifelse(tag %like% "Low_protein", "Yes", ""),
  Sample.Measured.Date.and.Time=sprintf("%sT%sZ", date_MEASURED_AT, time_MEASURED_AT),
  Sample.Prepared.Date.and.Time=sprintf("%sT%sZ", date_PREPARED_AT, time_PREPARED_AT),
  Well.Position.Within.Plate=plate_position, Shipment.Batch=batch)]

fwrite(sinfo, sep="\t", quote=FALSE, file="data/curated/preqc/sample_processing_information.txt")

