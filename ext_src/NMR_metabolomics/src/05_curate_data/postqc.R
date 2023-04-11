library(data.table)
library(ukbnmr)

# make output directory
system("mkdir -p data/curated/postqc", wait=TRUE)

# Load and curate biomarker measurements
nmr <- fread("output/techqc_drift_batch_fix/biomarker_measurements.txt")
nmr[, visit := ifelse(visit == "Main Phase", 0L, 1L)]
setnames(nmr, "visit", "visit_index")
setnames(nmr, "eid_30418", "eid")
fwrite(nmr, sep="\t", quote=FALSE, file="data/curated/postqc/biomarker_measurements.txt")

# Record biomarker information
fwrite(ukbnmr::nmr_info, sep="\t", quote=FALSE, file="data/curated/postqc/biomarker_information.txt")

# Record log offset diagnostic information
system("cp output/techqc_drift_batch_fix/log_offset.txt data/curated/postqc/log_offset_diagnostic_information.txt")

# Load and curate biomarker QC flags
flags <- fread("output/techqc_drift_batch_fix/biomarker_QC_tags.txt")
flags[, visit := ifelse(visit == "Main Phase", 0L, 1L)]
setnames(flags, "visit", "visit_index")
setnames(flags, "eid_30418", "eid")
fwrite(flags, sep="\t", quote=FALSE, file="data/curated/postqc/measurement_qc_flags.txt")

# Record outlier plate limit detection
system("cp output/techqc_drift_batch_fix/outlier_plate_detection_limits.txt data/curated/postqc/outlier_plate_detection_information.txt")

# load and curate relevant sample meta-data
sinfo <- fread("output/techqc_drift_batch_fix/sample_information.txt")
sinfo <- sinfo[!(removed), .(eid=eid_30418, visit_index=ifelse(visit == "Main Phase", 0L, 1L), 
  Shipment.Plate=as.character(plate_id), 
  Spectrometer=sprintf("Spectrometer %s", as.integer(factor(spectrometer))), 
  High.Lactate=ifelse(tag %like% "High_lactate", "Yes", ""),
  High.Pyruvate=ifelse(tag %like% "High_pyruvate", "Yes", ""),
  Low.Glucose=ifelse(tag %like% "Low_glucose", "Yes", ""),
  Low.Protein=ifelse(tag %like% "Low_protein", "Yes", ""),
  Sample.Measured.Date.and.Time=sprintf("%sT%sZ", date_MEASURED_AT, time_MEASURED_AT),
  Sample.Prepared.Date.and.Time=sprintf("%sT%sZ", date_PREPARED_AT, time_PREPARED_AT),
  Well.Position.Within.Plate=plate_position, Well.Row=well_row, Well.Column=well_column,
  Sample.Measured.Date=date_MEASURED_AT, Sample.Measured.Time=time_MEASURED_AT,
  Sample.Prepared.Date=date_PREPARED_AT, Sample.Prepared.Time=time_PREPARED_AT,
  Prep.to.Measure.Duration=prep_to_measured, 
  Plate.Measured.Date=plate_MEASURED_DATE,
  Spectrometer.Date.Bin=spectrometer_plate_bin,
  Shipment.Batch=batch)]

fwrite(sinfo, sep="\t", quote=FALSE, file="data/curated/postqc/sample_processing_information.txt")

