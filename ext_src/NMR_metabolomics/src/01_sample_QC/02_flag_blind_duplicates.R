library(data.table)

# Load data passing basic sample QC
sinfo <- fread("output/basic_sample_QC/sample_information.txt")
tags <- fread("output/basic_sample_QC/biomarker_QC_tags.txt")
nmr <- fread("output/basic_sample_QC/biomarker_measurements.txt")

# identify and flag blind duplicates
bd <- sinfo[!(removed), .N, by=.(eid_30418, visit)][N > 1]
sinfo[!(removed), has_blind_duplicate := FALSE]
sinfo[bd, on = .(eid_30418, visit), has_blind_duplicate := TRUE]

# From our previous QC of the Phase 1 data we can identify the blind duplicates
# that UKB dropped from their public release to also drop here
phase1_sinfo <- fread("data/Phase1_NMR_QC/data/raw/processed/harmonized_sample_information.txt")
phase1_bd_drop <- phase1_sinfo[in_ukb_raw == "No biomarkers" & (blind_duplicate) & !(removed)]
phase1_bd_drop <- phase1_bd_drop[!sinfo[(removed)], on = .(sample_id, visit)]
sinfo[phase1_bd_drop, on = .(sample_id, visit), c("removed", "removal_reason") := .(TRUE, "Blind duplicate sample (phase 1 sample excluded by UKB)")]

# For phase 2, we can try to take a slightly more principled approach by keeping
# the sample with the smallest time between sample prep and measurement, although
# the exact row kept will likely change with the UKB public release (I have no
# insight into their process for choosing which sample to keep from among blind
# duplicates, it appears to be random).
phase2_bd <- sinfo[!(removed), .N, by=.(eid_30418, visit)][N > 1]
phase2_bd <- phase2_bd[sinfo[!(removed)], on = .(eid_30418, visit), nomatch=0, .(uid, eid_30418, visit, prep_to_measured)]
phase2_bd[, keep := FALSE]
phase2_bd_keep <- phase2_bd[,.SD[which.min(prep_to_measured)], by=.(eid_30418, visit)]
phase2_bd[phase2_bd_keep, on = .(uid), keep := TRUE]
sinfo[phase2_bd[!(keep)], on = .(uid), c("removed", "removal_reason") := .(TRUE, "Blind duplicate sample (phase 2 sample, second or more by sample degradation time)")]

# Drop removed samples
nmr <- nmr[!sinfo[(removed)], on = .(sample_id, visit)]
tags <- tags[!sinfo[(removed)], on = .(sample_id, visit)]

# Replace sample ID with eid
nmr[sinfo[!(removed)], on = .(sample_id, visit), sample_id := eid_30418]
setnames(nmr, "sample_id", "eid_30418")

tags[sinfo[!(removed)], on = .(sample_id, visit), sample_id := eid_30418]
setnames(tags, "sample_id", "eid_30418")

# Write out 
if (!dir.exists("output/basic_participant_QC")) dir.create("output/basic_participant_QC", recursive=TRUE)

fwrite(nmr, sep="\t", quote=FALSE, file="output/basic_participant_QC/biomarker_measurements.txt")
fwrite(sinfo[visit != "Internal Control"], sep="\t", quote=FALSE, file="output/basic_participant_QC/sample_information.txt")
fwrite(tags, sep="\t", quote=FALSE, file="output/basic_participant_QC/biomarker_QC_tags.txt")

