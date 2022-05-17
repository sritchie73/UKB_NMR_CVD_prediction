library(data.table)
library(ukbnmr)

# Load raw data extracted from decoded UKB
raw <- fread("data/raw/ukbiobank/extracted/nmr.csv")

# Remove technical variation
clean <- remove_technical_variation(raw) # takes about 15 minutes

# Extract fasting time (not handled by ukbnmr)
fasting <- raw[,.SD,.SDcols=c("eid", names(raw)[names(raw) %like% "^74-"])]
fasting <- melt(fasting, id.vars="eid", na.rm=TRUE)
fasting[, visit.repeat := gsub(".*-", "", variable)]
fasting[, visit_index := as.numeric(gsub("\\..*", "", visit.repeat))]
fasting <- fasting[,.(eid, visit_index, fasting_time=value)]
fasting <- fasting[clean$biomarkers[,.(eid, visit_index)], on = .(eid, visit_index)] # filter to samples present with NMR data


# Write out 
system("mkdir -p output", wait=TRUE)

fwrite(clean$biomarkers, sep="\t", quote=FALSE, file="output/nmr_techadj.txt")
fwrite(clean$biomarker_qc_flags, sep="\t", quote=FALSE, file="output/measurement_qc_flags.txt")
fwrite(clean$sample_processing, sep="\t", quote=FALSE, file="output/sample_processing_information.txt")
fwrite(clean$outlier_plate_detection, sep="\t", quote=FALSE, file="output/outlier_plate_detection_information.txt")
fwrite(clean$log_offset, sep="\t", quote=FALSE, file="output/log_offset_diagnostic_information.txt")
fwrite(ukbnmr::nmr_info, sep="\t", quote=FALSE, file="output/biomarker_information.txt")
fwrite(fasting, sep="\t", quote=FALSE, file="output/fasting_time_hours.txt")

