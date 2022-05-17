library(data.table)
library(lubridate)

system("mkdir -p data/curated/ukbiobank/self_report/", wait=TRUE)

raw <- fread("data/raw/ukbiobank/extracted/self_report/self_report.csv")
info <- fread("data/raw/ukbiobank/extracted/self_report/field_info.txt")

# Convert to long format per field - we split this out into a list, one per 
# field, as fields are a mix of codes (integer), age of event (numeric), and
# date of event.
raw <- lapply(info$field.id, function(fid) {
  dt <- raw[,.SD,.SDcols=c("eid", names(raw)[names(raw) %like% sprintf("^%s-", fid)])]
  suppressWarnings(dt <- melt(dt, id.vars="eid", na.rm=TRUE)) # don't warn on coercion
  dt[, visit_index := gsub("(.*-)|(\\..*)", "", variable)]
  dt[, repeat_index := gsub(".*\\.", "", variable)]
  dt <- dt[, .(eid, visit_index, repeat_index, value)]
  return(dt)
})
names(raw) <- as.character(info$field.id)

# Curate code lists
codes <- rbind(idcol="field.id", "6150" = raw[["6150"]], "6152" = raw[["6152"]], 
  "2443" = raw[["2443"]], "4041" = raw[["4041"]], "2986" = raw[["2986"]], 
  "2463" = raw[["2463"]], "6151" = raw[["6151"]], "3005" = raw[["3005"]], 
  "4728" = raw[["4728"]], "5452" = raw[["5452"]], "5463" = raw[["5463"]], 
  "5474" = raw[["5474"]], "5485" = raw[["5485"]], "5496" = raw[["5496"]], 
  "5507" = raw[["5507"]], "5518" = raw[["5518"]], "5529" = raw[["5529"]], 
  "5540" = raw[["5540"]], "20002" = raw[["20002"]], "20001" = raw[["20001"]], 
  "20004" = raw[["20004"]], "20126" = raw[["20126"]], "20122" = raw[["20122"]], 
  "20127" = raw[["20127"]], "20124" = raw[["20124"]], "20125" = raw[["20125"]], 
  "20123" = raw[["20123"]], "1920" = raw[["1920"]], "1930" = raw[["1930"]], 
  "1940" = raw[["1940"]], "1950" = raw[["1950"]], "1960" = raw[["1960"]], 
  "1970" = raw[["1970"]], "1980" = raw[["1980"]], "1990" = raw[["1990"]], 
  "2000" = raw[["2000"]], "2010" = raw[["2010"]], "2020" = raw[["2020"]], 
  "2030" = raw[["2030"]], "2040" = raw[["2040"]], "4526" = raw[["4526"]], 
  "4537" = raw[["4537"]], "4548" = raw[["4548"]], "4559" = raw[["4559"]], 
  "4570" = raw[["4570"]], "4581" = raw[["4581"]], "2050" = raw[["2050"]], 
  "2060" = raw[["2060"]], "2070" = raw[["2070"]], "2080" = raw[["2080"]], 
  "2090" = raw[["2090"]], "2100" = raw[["2100"]], "4598" = raw[["4598"]], 
  "4609" = raw[["4609"]], "4620" = raw[["4620"]], "4631" = raw[["4631"]], 
  "5375" = raw[["5375"]], "5386" = raw[["5386"]], "4642" = raw[["4642"]], 
  "4653" = raw[["4653"]], "6156" = raw[["6156"]], "5663" = raw[["5663"]], 
  "5674" = raw[["5674"]], "6145" = raw[["6145"]])

# Add in fields that encode missing data for verbal interview questions
verbal_missing <- rbind(idcol="field.id", 
  "20001" = raw[["2453"]][value < 0], 
  "20002" = raw[["2473"]][value < 0],
  "20004" = raw[["2415"]][value < 0],
  "20004" = raw[["2844"]][value < 0]
)

# Some people who answered "Don't know" or "Prefer not to answer" had entries
# in the verbal interview corresponding to diseases/operations. We don't want
# to drop these answers as samples with missing data downstream, so we remove
# these from the 'verbal_missing' table before adding it to the 'codes' table.
verbal_missing <- verbal_missing[!codes, on = .(field.id, eid, visit_index)]

codes <- rbind(codes, verbal_missing)

# Create function for backfilling missing data, e.g. setting "None of the above"
# for verbal interview questions where the participant has reported no illness in
# the preceeding touchscreen survey questions, or "Not applicable" for specific 
# questions asked as follow-ups in touchscreen survey, or "No" for mental health
# question fields which only contain entries for "Yes" answers.
backfill <- function(backfill_field, backfill_code, fields_with_all_participants) {
  stopifnot(length(backfill_field) == 1)
  stopifnot(length(backfill_code) == 1)
  answered <- unique(codes[field.id == backfill_field, .(eid, visit_index)])
  all_participants <- unique(rbindlist(raw[fields_with_all_participants], idcol="field.id")[, .(eid, visit_index)])
  to_backfill <- fsetdiff(all_participants, answered)
  if (nrow(to_backfill) == 0) { 
    return(NULL)
  }
  to_backfill[, .(field.id=backfill_field, eid, visit_index, repeat_index=0, value=backfill_code)]
}

codes <- rbind(codes, 
  # None of the above for verbal interview questions where participant has reported 
  # No cancers, No illnesses, or No operations
  backfill("20001", -7, "2453"), 
  backfill("20002", -7, c("6150", "6152", "2443", "2473")),
  backfill("20004", -7, c("2415", "2844")),
  # Not applicable for touchscreen survey questions asked as follow-ups to Yes to 
  # fracturing a bone within last 5 years
  backfill("6151", -2, "2463"),
  backfill("3005", -2, "2463"),
  # Not applicable for touchscreen survey questions asked as follow-ups to Yes to 
  # diabetes diagnosis
  backfill("2986", -2, "2443"),
  backfill("4041", -2, "2443"),
  # Not applicable for touchscreen survey questions asked as follow-ups to Yes to
  # leg pain on walking
  backfill("5452", -2, "4728"),
  backfill("5463", -2, "4728"), 
  backfill("5474", -2, "4728"),
  backfill("5485", -2, "4728"),
  backfill("5496", -2, "4728"),
  backfill("5507", -2, "4728"),
  backfill("5518", -2, "4728"),
  backfill("5529", -2, "4728"),
  backfill("5540", -2, "4728"),
  # Mental health survey questions only added to touchscreen survey question after April 2009.
  # Many of the curated fields here only contain "Yes" answers, so we backfill "No" for those
  # asked these questions (captured by field 20126), 
  backfill("20124", 0, "20126"),
  backfill("20125", 0, "20126"),
  backfill("20123", 0, "20126"),
  backfill("20122", 0, "20126")
)

# Rename value column to code
setnames(codes, "value", "code")

# Curate label sets
labels <- rbind(idcol="field.id",
  "6150" = rbind(
    data.table(code=1, label="Heart attack"),
    data.table(code=2, label="Angina"),
    data.table(code=3, label="Stroke"),
    data.table(code=4, label="High blood pressure"),
    data.table(code=-7, label="None of the above"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "6152" = rbind(
    data.table(code=5, label="Blood clot in the leg (DVT)"),
    data.table(code=7, label="Blood clot in the lung"),
    data.table(code=6, label="Emphysema/chronic bronchitis"),
    data.table(code=8, label="Asthma"),
    data.table(code=9, label="Hayfever, allergic rhinitis or eczema"),
    data.table(code=-7, label="None of the above"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "2443" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "4041" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-2, label="Not applicable"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "2986" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-2, label="Not applicable"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "2463" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "6151" = rbind(
    data.table(code=1, label="Ankle"),
    data.table(code=2, label="Leg"),
    data.table(code=3, label="Hip"),
    data.table(code=4, label="Spine"),
    data.table(code=5, label="Wrist"),
    data.table(code=6, label="Arm"),
    data.table(code=7, label="Other bones"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-2, label="Not applicable"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "3005" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-2, label="Not applicable"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "4728" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "5452" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-2, label="Not applicable"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "5463" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-2, label="Not applicable"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "5474" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-2, label="Not applicable"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "5485" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-2, label="Not applicable"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "5496" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-2, label="Not applicable"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "5507" = rbind(
    data.table(code=1, label="Stop"),
    data.table(code=2, label="Slow down"),
    data.table(code=3, label="Continue at same pace"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-2, label="Not applicable"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "5518" = rbind(
    data.table(code=1, label="Pain usually continues for more than 10 minutes"),
    data.table(code=2, label="Pain usually disappears in less than 10 minutes"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-2, label="Not applicable"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "5529" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-2, label="Not applicable"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "5540" = rbind(
    data.table(code=0, label="No"),
    data.table(code=1, label="Yes, toes"),
    data.table(code=2, label="Yes, leg below the knee"),
    data.table(code=3, label="Yes, leg above the knee"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-2, label="Not applicable"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "20002" = rbind(
    data.table(code=1065, label="hypertension"),
    data.table(code=1066, label="heart/cardiac problem"),
    data.table(code=1067, label="peripheral vascular disease"),
    data.table(code=1068, label="venous thromboembolic disease"),
    data.table(code=1072, label="essential hypertension"),
    data.table(code=1073, label="gestational hypertension/pre-eclampsia"),
    data.table(code=1074, label="angina"),
    data.table(code=1075, label="heart attack/myocardial infarction"),
    data.table(code=1076, label="heart failure/pulmonary odema"),
    data.table(code=1077, label="heart arrhythmia"),
    data.table(code=1078, label="heart valve problem/heart murmur"),
    data.table(code=1079, label="cardiomyopathy"),
    data.table(code=1080, label="pericardial problem"),
    data.table(code=1081, label="stroke"),
    data.table(code=1082, label="transient ischaemic attack (tia)"),
    data.table(code=1083, label="subdural haemorrhage/haematoma"),
    data.table(code=1086, label="subarachnoid haemorrhage"),
    data.table(code=1087, label="leg claudication/ intermittent claudication"),
    data.table(code=1088, label="arterial embolism"),
    data.table(code=1093, label="pulmonary embolism +/- dvt"),
    data.table(code=1094, label="deep venous thrombosis (dvt)"),
    data.table(code=1111, label="asthma"),
    data.table(code=1112, label="chronic obstructive airways disease/copd"),
    data.table(code=1113, label="emphysema/chronic bronchitis"),
    data.table(code=1114, label="bronchiectasis"),
    data.table(code=1115, label="interstitial lung disease"),
    data.table(code=1117, label="other respiratory problems"),
    data.table(code=1120, label="asbestosis"),
    data.table(code=1121, label="pulmonary fibrosis"),
    data.table(code=1122, label="fibrosing alveolitis/unspecified alveolitis"),
    data.table(code=1123, label="sleep apnoea"),
    data.table(code=1124, label="respiratory failure"),
    data.table(code=1125, label="pleurisy"),
    data.table(code=1126, label="spontaneous pneumothorax/recurrent pneumothorax"),
    data.table(code=1134, label="oesophageal disorder"),
    data.table(code=1135, label="stomach disorder"),
    data.table(code=1136, label="liver/biliary/pancreas problem"),
    data.table(code=1137, label="other abdominal problem"),
    data.table(code=1138, label="gastro-oesophageal reflux (gord) / gastric reflux"),
    data.table(code=1139, label="oesophagitis/barretts oesophagus"),
    data.table(code=1140, label="oesophageal stricture"),
    data.table(code=1141, label="oesophageal varicies"),
    data.table(code=1142, label="gastric/stomach ulcers"),
    data.table(code=1143, label="gastritis/gastric erosions"),
    data.table(code=1154, label="irritable bowel syndrome"),
    data.table(code=1155, label="hepatitis"),
    data.table(code=1156, label="infective/viral hepatitis"),
    data.table(code=1157, label="non-infective hepatitis"),
    data.table(code=1158, label="liver failure/cirrhosis"),
    data.table(code=1159, label="bile duct disease"),
    data.table(code=1160, label="bile duct obstruction/ascending cholangitis"),
    data.table(code=1161, label="gall bladder disease"),
    data.table(code=1162, label="cholelithiasis/gall stones"),
    data.table(code=1163, label="cholecystitis"),
    data.table(code=1164, label="pancreatic disease"),
    data.table(code=1165, label="pancreatitis"),
    data.table(code=1190, label="peritonitis"),
    data.table(code=1191, label="gastrointestinal bleeding"),
    data.table(code=1192, label="renal/kidney failure"),
    data.table(code=1193, label="renal failure requiring dialysis"),
    data.table(code=1194, label="renal failure not requiring dialysis"),
    data.table(code=1196, label="urinary tract infection/kidney infection"),
    data.table(code=1197, label="kidney stone/ureter stone/bladder stone"),
    data.table(code=1200, label="ureteric obstruction/hydronephrosis"),
    data.table(code=1201, label="bladder problem (not cancer)"),
    data.table(code=1202, label="urinary frequency / incontinence"),
    data.table(code=1207, label="prostate problem (not cancer)"),
    data.table(code=1210, label="scrotal problem (not cancer)"),
    data.table(code=1214, label="testicular problems (not cancer)"),
    data.table(code=1220, label="diabetes"),
    data.table(code=1221, label="gestational diabetes"),
    data.table(code=1222, label="type 1 diabetes"),
    data.table(code=1223, label="type 2 diabetes"),
    data.table(code=1224, label="thyroid problem (not cancer)"),
    data.table(code=1225, label="hyperthyroidism/thyrotoxicosis"),
    data.table(code=1226, label="hypothyroidism/myxoedema"),
    data.table(code=1228, label="thyroid radioablation therapy"),
    data.table(code=1229, label="parathyroid gland problem (not cancer)"),
    data.table(code=1230, label="parathyroid hyperplasia/adenoma"),
    data.table(code=1232, label="disorder of adrenal gland"),
    data.table(code=1233, label="adrenal tumour"),
    data.table(code=1234, label="adrenocortical insufficiency/addison's disease"),
    data.table(code=1235, label="hyperaldosteronism/conn's syndrome"),
    data.table(code=1236, label="phaeochromocytoma"),
    data.table(code=1237, label="disorder or pituitary gland"),
    data.table(code=1238, label="pituitary adenoma/tumour"),
    data.table(code=1239, label="cushings syndrome"),
    data.table(code=1240, label="neurological injury/trauma"),
    data.table(code=1242, label="eye/eyelid problem"),
    data.table(code=1243, label="psychological/psychiatric problem"),
    data.table(code=1244, label="infection of nervous system"),
    data.table(code=1245, label="brain abscess/intracranial abscess"),
    data.table(code=1246, label="encephalitis"),
    data.table(code=1247, label="meningitis"),
    data.table(code=1248, label="spinal abscess"),
    data.table(code=1249, label="cranial nerve problem/palsy"),
    data.table(code=1250, label="bell's palsy/facial nerve palsy"),
    data.table(code=1251, label="spinal cord disorder"),
    data.table(code=1252, label="paraplegia"),
    data.table(code=1254, label="peripheral nerve disorder"),
    data.table(code=1255, label="peripheral neuropathy"),
    data.table(code=1256, label="acute infective polyneuritis/guillain-barre syndrome"),
    data.table(code=1257, label="trapped nerve/compressed nerve"),
    data.table(code=1258, label="chronic/degenerative neurological problem"),
    data.table(code=1259, label="motor neurone disease"),
    data.table(code=1260, label="myasthenia gravis"),
    data.table(code=1261, label="multiple sclerosis"),
    data.table(code=1262, label="parkinsons disease"),
    data.table(code=1263, label="dementia/alzheimers/cognitive impairment"),
    data.table(code=1264, label="epilepsy"),
    data.table(code=1265, label="migraine"),
    data.table(code=1266, label="head injury"),
    data.table(code=1267, label="spinal injury"),
    data.table(code=1274, label="eye infection"),
    data.table(code=1275, label="retinal problem"),
    data.table(code=1276, label="diabetic eye disease"),
    data.table(code=1277, label="glaucoma"),
    data.table(code=1278, label="cataract"),
    data.table(code=1279, label="eye trauma"),
    data.table(code=1281, label="retinal detachment"),
    data.table(code=1282, label="retinal artery/vein occlusion"),
    data.table(code=1286, label="depression"),
    data.table(code=1287, label="anxiety/panic attacks"),
    data.table(code=1288, label="nervous breakdown"),
    data.table(code=1289, label="schizophrenia"),
    data.table(code=1290, label="deliberate self-harm/suicide attempt"),
    data.table(code=1291, label="mania/bipolar disorder/manic depression"),
    data.table(code=1293, label="bone disorder"),
    data.table(code=1294, label="back problem"),
    data.table(code=1295, label="joint disorder"),
    data.table(code=1297, label="muscle/soft tissue problem"),
    data.table(code=1308, label="osteomyelitis"),
    data.table(code=1309, label="osteoporosis"),
    data.table(code=1310, label="paget's disease"),
    data.table(code=1311, label="spine arthritis/spondylitis"),
    data.table(code=1312, label="prolapsed disc/slipped disc"),
    data.table(code=1313, label="ankylosing spondylitis"),
    data.table(code=1322, label="myositis/myopathy"),
    data.table(code=1327, label="low platelets/platelet disorder"),
    data.table(code=1328, label="haemophilia"),
    data.table(code=1330, label="iron deficiency anaemia"),
    data.table(code=1331, label="pernicious anaemia"),
    data.table(code=1332, label="aplastic anaemia"),
    data.table(code=1339, label="sickle cell disease"),
    data.table(code=1340, label="thalassaemia"),
    data.table(code=1344, label="stevens johnson syndrome"),
    data.table(code=1345, label="pemphigoid/pemphigus"),
    data.table(code=1348, label="gynaecological disorder (not cancer)"),
    data.table(code=1349, label="ovarian cyst or cysts"),
    data.table(code=1350, label="polycystic ovaries/polycystic ovarian syndrome"),
    data.table(code=1351, label="uterine fibroids"),
    data.table(code=1352, label="uterine polyps"),
    data.table(code=1353, label="vaginal prolapse/uterine prolapse"),
    data.table(code=1364, label="breast disease (not cancer)"),
    data.table(code=1366, label="fibrocystic disease"),
    data.table(code=1367, label="breast cysts"),
    data.table(code=1371, label="sarcoidosis"),
    data.table(code=1372, label="vasculitis"),
    data.table(code=1373, label="connective tissue disorder"),
    data.table(code=1374, label="allergy/hypersensitivity/anaphylaxis"),
    data.table(code=1376, label="giant cell/temporal arteritis"),
    data.table(code=1377, label="polymyalgia rheumatica"),
    data.table(code=1378, label="wegners granulmatosis"),
    data.table(code=1379, label="microscopic polyarteritis"),
    data.table(code=1380, label="polyartertis nodosa"),
    data.table(code=1381, label="systemic lupus erythematosis/sle"),
    data.table(code=1382, label="sjogren's syndrome/sicca syndrome"),
    data.table(code=1383, label="dermatopolymyositis"),
    data.table(code=1384, label="scleroderma/systemic sclerosis"),
    data.table(code=1385, label="allergy or anaphylactic reaction to food"),
    data.table(code=1386, label="allergy or anaphylactic reaction to drug"),
    data.table(code=1387, label="hayfever/allergic rhinitis"),
    data.table(code=1394, label="peripheral nerve injury"),
    data.table(code=1396, label="enlarged prostate"),
    data.table(code=1397, label="other demyelinating disease (not multiple sclerosis)"),
    data.table(code=1398, label="pneumonia"),
    data.table(code=1400, label="peptic ulcer"),
    data.table(code=1402, label="endometriosis"),
    data.table(code=1403, label="female infertility"),
    data.table(code=1404, label="male infertility"),
    data.table(code=1405, label="other renal/kidney problem"),
    data.table(code=1406, label="muscle or soft tissue injuries"),
    data.table(code=1407, label="burns"),
    data.table(code=1408, label="alcohol dependency"),
    data.table(code=1409, label="opioid dependency"),
    data.table(code=1410, label="other substance abuse/dependency"),
    data.table(code=1411, label="lung abscess"),
    data.table(code=1412, label="bronchitis"),
    data.table(code=1413, label="nasal/sinus disorder"),
    data.table(code=1414, label="throat or larynx disorder"),
    data.table(code=1415, label="ear/vestibular disorder"),
    data.table(code=1416, label="chronic sinusitis"),
    data.table(code=1417, label="nasal polyps"),
    data.table(code=1418, label="chronic laryngitis"),
    data.table(code=1419, label="vocal cord polyp"),
    data.table(code=1420, label="otosclerosis"),
    data.table(code=1421, label="meniere's disease"),
    data.table(code=1425, label="cerebral aneurysm"),
    data.table(code=1426, label="myocarditis"),
    data.table(code=1427, label="polycystic kidney"),
    data.table(code=1428, label="thyroiditis"),
    data.table(code=1429, label="acromegaly"),
    data.table(code=1430, label="hypopituitarism"),
    data.table(code=1431, label="hyperprolactinaemia"),
    data.table(code=1432, label="carcinoid syndrome/tumour"),
    data.table(code=1433, label="cerebral palsy"),
    data.table(code=1434, label="other neurological problem"),
    data.table(code=1435, label="optic neuritis"),
    data.table(code=1436, label="headaches (not migraine)"),
    data.table(code=1437, label="myasthenia gravis"),
    data.table(code=1438, label="polycythaemia vera"),
    data.table(code=1445, label="clotting disorder/excessive bleeding"),
    data.table(code=1446, label="anaemia"),
    data.table(code=1447, label="pancytopenia"),
    data.table(code=1448, label="neutropenia/lymphopenia"),
    data.table(code=1449, label="myeloproliferative disorder"),
    data.table(code=1450, label="monoclonal gammopathy/not myeloma"),
    data.table(code=1451, label="hereditary/genetic haematological disorder"),
    data.table(code=1452, label="eczema/dermatitis"),
    data.table(code=1453, label="psoriasis"),
    data.table(code=1454, label="blistering/desquamating skin disorder"),
    data.table(code=1455, label="chronic skin ulcers"),
    data.table(code=1456, label="malabsorption/coeliac disease"),
    data.table(code=1457, label="duodenal ulcer"),
    data.table(code=1458, label="diverticular disease/diverticulitis"),
    data.table(code=1459, label="colitis/not crohns or ulcerative colitis"),
    data.table(code=1460, label="rectal or colon adenoma/polyps"),
    data.table(code=1461, label="inflammatory bowel disease"),
    data.table(code=1462, label="crohns disease"),
    data.table(code=1463, label="ulcerative colitis"),
    data.table(code=1464, label="rheumatoid arthritis"),
    data.table(code=1465, label="osteoarthritis"),
    data.table(code=1466, label="gout"),
    data.table(code=1467, label="other joint disorder"),
    data.table(code=1468, label="diabetic neuropathy/ulcers"),
    data.table(code=1469, label="post-traumatic stress disorder"),
    data.table(code=1470, label="anorexia/bulimia/other eating disorder"),
    data.table(code=1471, label="atrial fibrillation"),
    data.table(code=1472, label="emphysema"),
    data.table(code=1473, label="high cholesterol"),
    data.table(code=1474, label="hiatus hernia"),
    data.table(code=1475, label="sclerosing cholangitis"),
    data.table(code=1476, label="sciatica"),
    data.table(code=1477, label="psoriatic arthropathy"),
    data.table(code=1478, label="cervical spondylosis"),
    data.table(code=1479, label="rheumatic fever"),
    data.table(code=1480, label="dermatomyositis"),
    data.table(code=1481, label="polymyositis"),
    data.table(code=1482, label="chronic fatigue syndrome"),
    data.table(code=1483, label="atrial flutter"),
    data.table(code=1484, label="wolff parkinson white / wpw syndrome"),
    data.table(code=1485, label="irregular heart beat"),
    data.table(code=1486, label="sick sinus syndrome"),
    data.table(code=1487, label="svt / supraventricular tachycardia"),
    data.table(code=1488, label="mitral valve prolapse"),
    data.table(code=1489, label="mitral stenosis"),
    data.table(code=1490, label="aortic stenosis"),
    data.table(code=1491, label="brain haemorrhage"),
    data.table(code=1492, label="aortic aneurysm"),
    data.table(code=1493, label="other venous/lymphatic disease"),
    data.table(code=1494, label="varicose veins"),
    data.table(code=1495, label="lymphoedema"),
    data.table(code=1496, label="alpha-1 antitrypsin deficiency"),
    data.table(code=1497, label="pneumothorax"),
    data.table(code=1498, label="empyema"),
    data.table(code=1499, label="labyrinthitis"),
    data.table(code=1500, label="vertigo"),
    data.table(code=1501, label="pyloric stenosis"),
    data.table(code=1502, label="appendicitis"),
    data.table(code=1503, label="anal problem"),
    data.table(code=1504, label="anal fissure"),
    data.table(code=1505, label="haemorrhoids / piles"),
    data.table(code=1506, label="primary biliary cirrhosis"),
    data.table(code=1507, label="haemochromatosis"),
    data.table(code=1508, label="jaundice (unknown cause)"),
    data.table(code=1509, label="gastroenteritis/dysentry"),
    data.table(code=1510, label="dyspepsia / indigestion"),
    data.table(code=1511, label="abdominal hernia"),
    data.table(code=1512, label="umbilical hernia"),
    data.table(code=1513, label="inguinal hernia"),
    data.table(code=1514, label="cystitis"),
    data.table(code=1515, label="pyelonephritis"),
    data.table(code=1516, label="bph / benign prostatic hypertrophy"),
    data.table(code=1517, label="prostatitis"),
    data.table(code=1518, label="erectile dysfunction / impotence"),
    data.table(code=1519, label="kidney nephropathy"),
    data.table(code=1520, label="iga nephropathy"),
    data.table(code=1521, label="diabetes insipidus"),
    data.table(code=1522, label="grave's disease"),
    data.table(code=1523, label="trigemminal neuralgia"),
    data.table(code=1524, label="spina bifida"),
    data.table(code=1525, label="benign / essential tremor"),
    data.table(code=1526, label="polio / poliomyelitis"),
    data.table(code=1527, label="retinitis pigmentosa"),
    data.table(code=1528, label="macular degeneration"),
    data.table(code=1529, label="dry eyes"),
    data.table(code=1530, label="iritis"),
    data.table(code=1531, label="post-natal depression"),
    data.table(code=1532, label="disc problem"),
    data.table(code=1533, label="disc degeneration"),
    data.table(code=1534, label="back pain"),
    data.table(code=1535, label="scoliosis"),
    data.table(code=1536, label="spinal stenosis"),
    data.table(code=1537, label="joint pain"),
    data.table(code=1538, label="arthritis (nos)"),
    data.table(code=1540, label="plantar fascitis"),
    data.table(code=1541, label="carpal tunnel syndrome"),
    data.table(code=1542, label="fibromyalgia"),
    data.table(code=1544, label="dupuytren's contracture"),
    data.table(code=1545, label="neck problem/injury"),
    data.table(code=1546, label="essential thrombocytosis"),
    data.table(code=1548, label="acne/acne vulgaris"),
    data.table(code=1549, label="lichen planus"),
    data.table(code=1550, label="lichen sclerosis"),
    data.table(code=1551, label="ovarian problem"),
    data.table(code=1552, label="uterine problem"),
    data.table(code=1553, label="cervical problem"),
    data.table(code=1554, label="cervical intra-epithelial neoplasia (cin) / precancerous cells cervix"),
    data.table(code=1555, label="cervical polyps"),
    data.table(code=1556, label="menorrhagia (unknown cause)"),
    data.table(code=1557, label="pelvic inflammatory disease/ pid"),
    data.table(code=1558, label="ectopic pregnancy"),
    data.table(code=1559, label="miscarriage"),
    data.table(code=1560, label="breast fibroadenoma"),
    data.table(code=1561, label="raynaud's phenomenon/disease"),
    data.table(code=1562, label="food intolerance"),
    data.table(code=1563, label="urticaria"),
    data.table(code=1564, label="antiphospholipid syndrome"),
    data.table(code=1578, label="hepatitis a"),
    data.table(code=1579, label="hepatitis b"),
    data.table(code=1580, label="hepatitis c"),
    data.table(code=1581, label="hepatitis d"),
    data.table(code=1582, label="hepatitis e"),
    data.table(code=1583, label="ischaemic stroke"),
    data.table(code=1584, label="mitral valve disease"),
    data.table(code=1585, label="mitral regurgitation / incompetence"),
    data.table(code=1586, label="aortic valve disease"),
    data.table(code=1587, label="aortic regurgitation / incompetence"),
    data.table(code=1588, label="hypertrophic cardiomyopathy (hcm / hocm)"),
    data.table(code=1589, label="pericarditis"),
    data.table(code=1590, label="pericardial effusion"),
    data.table(code=1591, label="aortic aneurysm rupture"),
    data.table(code=1592, label="aortic dissection"),
    data.table(code=1593, label="varicose ulcer"),
    data.table(code=1594, label="respiratory infection"),
    data.table(code=1595, label="pleural plaques (not known asbestosis)"),
    data.table(code=1596, label="pleural effusion"),
    data.table(code=1597, label="tinnitus / tiniitis"),
    data.table(code=1598, label="tonsiltis"),
    data.table(code=1599, label="constipation"),
    data.table(code=1600, label="bowel / intestinal perforation"),
    data.table(code=1601, label="bowel / intestinal infarction"),
    data.table(code=1602, label="bowel / intestinal obstruction"),
    data.table(code=1603, label="rectal prolapse"),
    data.table(code=1604, label="alcoholic liver disease / alcoholic cirrhosis"),
    data.table(code=1605, label="femoral hernia"),
    data.table(code=1606, label="incisional hernia"),
    data.table(code=1607, label="diabetic nephropathy"),
    data.table(code=1608, label="nephritis"),
    data.table(code=1609, label="glomerulnephritis"),
    data.table(code=1610, label="thyroid goitre"),
    data.table(code=1611, label="hyperparathyroidism"),
    data.table(code=1613, label="blepharitis / eyelid infection"),
    data.table(code=1614, label="stress"),
    data.table(code=1615, label="obsessive compulsive disorder (ocd)"),
    data.table(code=1616, label="insomnia"),
    data.table(code=1617, label="osteopenia"),
    data.table(code=1618, label="soft tissue inflammation"),
    data.table(code=1619, label="tendonitis / tendinitis / tenosynovitis"),
    data.table(code=1620, label="bursitis"),
    data.table(code=1621, label="synovitis"),
    data.table(code=1622, label="epicondylitis"),
    data.table(code=1623, label="tennis elbow / lateral epicondylitis"),
    data.table(code=1624, label="housemaid's knee (prepatellar bursitis)"),
    data.table(code=1625, label="cellulitis"),
    data.table(code=1626, label="fracture skull / head"),
    data.table(code=1627, label="fracture jaw"),
    data.table(code=1628, label="fracture nose"),
    data.table(code=1629, label="fracture face / orbit / eye socket"),
    data.table(code=1630, label="fracture neck / cervical fracture"),
    data.table(code=1631, label="fracture clavicle / collar bone"),
    data.table(code=1632, label="fracture shoulder / scapula"),
    data.table(code=1633, label="fracture upper arm / humerus / elbow"),
    data.table(code=1634, label="fracture forearm / wrist"),
    data.table(code=1635, label="fracture radius"),
    data.table(code=1636, label="fracture ulna"),
    data.table(code=1637, label="fracture wrist / colles fracture"),
    data.table(code=1638, label="fracture hand"),
    data.table(code=1639, label="fracture finger"),
    data.table(code=1640, label="fracture thumb"),
    data.table(code=1644, label="fracture rib"),
    data.table(code=1645, label="fracture sternum"),
    data.table(code=1646, label="fracture vertebra / crush fracture / vertebral collapse"),
    data.table(code=1647, label="fracture pelvis"),
    data.table(code=1648, label="fracture neck of femur / hip"),
    data.table(code=1649, label="fracture shaft of femur"),
    data.table(code=1650, label="fracture patella / knee"),
    data.table(code=1651, label="fracture lower leg / ankle"),
    data.table(code=1652, label="fracture tibia"),
    data.table(code=1653, label="fracture fibula"),
    data.table(code=1654, label="fracture foot"),
    data.table(code=1655, label="fracture metatarsal"),
    data.table(code=1656, label="fracture toe"),
    data.table(code=1657, label="septicaemia / sepsis"),
    data.table(code=1658, label="myelofibrosis"),
    data.table(code=1659, label="meningioma / benign meningeal tumour"),
    data.table(code=1660, label="rosacea"),
    data.table(code=1661, label="vitiligo"),
    data.table(code=1662, label="cervical erosion"),
    data.table(code=1663, label="abnormal smear (cervix)"),
    data.table(code=1664, label="dysmenorrhoea / dysmenorrhea"),
    data.table(code=1665, label="menopausal symptoms / menopause"),
    data.table(code=1666, label="benign breast lump"),
    data.table(code=1667, label="alopecia / hair loss"),
    data.table(code=1668, label="allergy to house dust mite"),
    data.table(code=1669, label="contact dermatitis"),
    data.table(code=1670, label="allergy to elastoplast"),
    data.table(code=1671, label="allergy to nickel"),
    data.table(code=1439, label="hiv/aids"),
    data.table(code=1568, label="measles / morbillivirus"),
    data.table(code=1570, label="rubella / german measles"),
    data.table(code=1674, label="varicella zoster virus"),
    data.table(code=1571, label="chickenpox"),
    data.table(code=1573, label="shingles"),
    data.table(code=1575, label="herpes simplex"),
    data.table(code=1567, label="infectious mononucleosis / glandular fever / epstein barr virus (ebv)"),
    data.table(code=1569, label="mumps / epidemic parotitis"),
    data.table(code=1566, label="mrsa / methicillin resistant staphylococcus aureus"),
    data.table(code=1442, label="helicobacter pylori"),
    data.table(code=1440, label="tuberculosis (tb)"),
    data.table(code=1441, label="malaria"),
    data.table(code=1443, label="schistosomiasis/bilharzia"),
    data.table(code=1576, label="dengue fever"),
    data.table(code=1577, label="typhoid fever"),
    data.table(code=1675, label="giardia / giardiasis"),
    data.table(code=1574, label="diphtheria"),
    data.table(code=1572, label="whooping cough / pertussis"),
    data.table(code=1676, label="yellow fever"),
    data.table(code=1677, label="scarlet fever / scarlatina"),
    data.table(code=1678, label="chlamydia"),
    data.table(code=1679, label="undescended testicle"),
    data.table(code=1680, label="bowen's disease"),
    data.table(code=1681, label="hydatiform mole"),
    data.table(code=1682, label="benign insulinoma"),
    data.table(code=1683, label="benign neuroma"),
    data.table(code=99999, label="unclassifiable"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer"),
    data.table(code=-7, label="None of the above")
  ),
  "20001" = rbind(
    data.table(code=1001, label="lung cancer"),
    data.table(code=1002, label="breast cancer"),
    data.table(code=1003, label="skin cancer"),
    data.table(code=1004, label="cancer of lip/mouth/pharynx/oral cavity"),
    data.table(code=1005, label="salivary gland cancer"),
    data.table(code=1006, label="larynx/throat cancer"),
    data.table(code=1007, label="nasal cavity cancer"),
    data.table(code=1008, label="ear cancer"),
    data.table(code=1009, label="sinus cancer"),
    data.table(code=1010, label="lip cancer"),
    data.table(code=1011, label="tongue cancer"),
    data.table(code=1012, label="gum cancer"),
    data.table(code=1015, label="parotid gland cancer"),
    data.table(code=1016, label="other salivary gland cancer"),
    data.table(code=1017, label="oesophageal cancer"),
    data.table(code=1018, label="stomach cancer"),
    data.table(code=1019, label="small intestine/small bowel cancer"),
    data.table(code=1020, label="large bowel cancer/colorectal cancer"),
    data.table(code=1021, label="anal cancer"),
    data.table(code=1022, label="colon cancer/sigmoid cancer"),
    data.table(code=1023, label="rectal cancer"),
    data.table(code=1024, label="liver/hepatocellular cancer"),
    data.table(code=1025, label="gallbladder/bile duct cancer"),
    data.table(code=1026, label="pancreas cancer"),
    data.table(code=1027, label="small cell lung cancer"),
    data.table(code=1028, label="non-small cell lung cancer"),
    data.table(code=1029, label="peripheral nerve/autonomic nerve cancer"),
    data.table(code=1030, label="eye and/or adnexal cancer"),
    data.table(code=1031, label="meningeal cancer / malignant meningioma"),
    data.table(code=1032, label="brain cancer / primary malignant brain tumour"),
    data.table(code=1033, label="spinal cord or cranial nerve cancer"),
    data.table(code=1034, label="kidney/renal cell cancer"),
    data.table(code=1035, label="bladder cancer"),
    data.table(code=1036, label="other cancer of urinary tract"),
    data.table(code=1037, label="female genital tract cancer"),
    data.table(code=1038, label="male genital tract cancer"),
    data.table(code=1039, label="ovarian cancer"),
    data.table(code=1040, label="uterine/endometrial cancer"),
    data.table(code=1041, label="cervical cancer"),
    data.table(code=1042, label="vaginal cancer"),
    data.table(code=1043, label="vulval cancer"),
    data.table(code=1044, label="prostate cancer"),
    data.table(code=1045, label="testicular cancer"),
    data.table(code=1046, label="penis cancer"),
    data.table(code=1047, label="lymphoma"),
    data.table(code=1048, label="leukaemia"),
    data.table(code=1050, label="multiple myeloma"),
    data.table(code=1051, label="myelofibrosis or myelodysplasia"),
    data.table(code=1052, label="hodgkins lymphoma / hodgkins disease"),
    data.table(code=1053, label="non-hodgkins lymphoma"),
    data.table(code=1055, label="chronic lymphocytic"),
    data.table(code=1056, label="chronic myeloid"),
    data.table(code=1058, label="other haematological malignancy"),
    data.table(code=1059, label="malignant melanoma"),
    data.table(code=1060, label="non-melanoma skin cancer"),
    data.table(code=1061, label="basal cell carcinoma"),
    data.table(code=1062, label="squamous cell carcinoma"),
    data.table(code=1063, label="primary bone cancer"),
    data.table(code=1064, label="mesothelioma"),
    data.table(code=1065, label="thyroid cancer"),
    data.table(code=1066, label="parathyroid cancer"),
    data.table(code=1067, label="adrenal cancer"),
    data.table(code=1068, label="sarcoma/fibrosarcoma"),
    data.table(code=1070, label="malignant lymph node, unspecified"),
    data.table(code=1071, label="metastatic cancer (unknown primary)"),
    data.table(code=1072, label="cin/pre-cancer cells cervix"),
    data.table(code=1073, label="rodent ulcer"),
    data.table(code=1074, label="acute myeloid leukaemia"),
    data.table(code=1075, label="retinoblastoma"),
    data.table(code=1076, label="kaposis sarcoma"),
    data.table(code=1077, label="mouth cancer"),
    data.table(code=1078, label="tonsil cancer"),
    data.table(code=1079, label="oropharynx / oropharyngeal cancer"),
    data.table(code=1080, label="trachea cancer"),
    data.table(code=1081, label="thymus cancer / malignant thymoma"),
    data.table(code=1082, label="heart / mediastinum cancer"),
    data.table(code=1084, label="respiratory / intrathoracic cancer"),
    data.table(code=1085, label="bone metastases / bony secondaries"),
    data.table(code=1086, label="appendix cancer"),
    data.table(code=1087, label="fallopian tube cancer"),
    data.table(code=1088, label="malignant insulinoma"),
    data.table(code=99999, label="unclassifiable"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer"),
    data.table(code=-7, label="None of the above")
  ),
  "20004" = rbind(
    data.table(code=1069, label="heart surgery"),
    data.table(code=1070, label="coronary angioplasty (ptca) +/- stent"),
    data.table(code=1071, label="other arterial surgery/revascularisation procedures"),
    data.table(code=1095, label="coronary artery bypass grafts (cabg)"),
    data.table(code=1096, label="pacemaker/defibrillator insertion"),
    data.table(code=1097, label="heart valve surgery"),
    data.table(code=1098, label="heart transplant"),
    data.table(code=1099, label="aortic valve repair/replacement"),
    data.table(code=1100, label="mitral valve repair/replacement"),
    data.table(code=1101, label="other valve repair/replacement"),
    data.table(code=1102, label="fem-pop bypass/leg artery bypass"),
    data.table(code=1103, label="leg artery aneurysm repair"),
    data.table(code=1104, label="aortic aneurysm/repair or stent"),
    data.table(code=1105, label="carotid artery surgery/endarterectomy"),
    data.table(code=1106, label="cerebral artery aneurysm surgery or clipping"),
    data.table(code=1107, label="non-coronary artery angioplasty +/- stent"),
    data.table(code=1108, label="leg artery angioplasty +/- stent"),
    data.table(code=1109, label="carotid artery angioplasty +/- stent"),
    data.table(code=1110, label="renal artery angioplasty +/- stent"),
    data.table(code=1118, label="ent surgery"),
    data.table(code=1127, label="ear surgery"),
    data.table(code=1128, label="nasal/sinus/nose surgery"),
    data.table(code=1129, label="throat/larynx surgery (incl tracheostomy)"),
    data.table(code=1130, label="mouth/salivary gland surgery"),
    data.table(code=1195, label="renal/kidney transplant"),
    data.table(code=1197, label="percutaneous/open kidney stone surgery/lithotripsy"),
    data.table(code=1199, label="ureteric surgery"),
    data.table(code=1207, label="prostate operation"),
    data.table(code=1208, label="radical prostatectomy"),
    data.table(code=1209, label="transurethral resection of prostate (turp)"),
    data.table(code=1216, label="testicular/scrotal operation"),
    data.table(code=1217, label="removal of testicle/orchidectomy"),
    data.table(code=1218, label="vasectomy"),
    data.table(code=1219, label="other urological surgery"),
    data.table(code=1224, label="thyroid surgery"),
    data.table(code=1228, label="thyroid radioablation therapy"),
    data.table(code=1229, label="parathyroidectomy"),
    data.table(code=1232, label="adrenal surgery"),
    data.table(code=1296, label="bone surgery/joint surgery"),
    data.table(code=1297, label="muscle/soft tissue surgery"),
    data.table(code=1318, label="hip replacement/revision"),
    data.table(code=1319, label="knee replacement/revision"),
    data.table(code=1320, label="spine or back surgery"),
    data.table(code=1354, label="gynaecological surgery"),
    data.table(code=1355, label="bilateral oophorectomy"),
    data.table(code=1356, label="unilateral oophorectomy"),
    data.table(code=1357, label="hysterectomy"),
    data.table(code=1358, label="hysterectomy with oophorectomy"),
    data.table(code=1359, label="hysterectomy with cervical sparing"),
    data.table(code=1360, label="endometrial ablation"),
    data.table(code=1361, label="vaginal prolapse/colposuspension"),
    data.table(code=1362, label="sterilisation"),
    data.table(code=1363, label="intrauterine contraceptive device insertion/removal"),
    data.table(code=1365, label="breast surgery"),
    data.table(code=1368, label="mastectomy"),
    data.table(code=1369, label="lumpectomy"),
    data.table(code=1370, label="mammoplasty/cosmetic operation on breast"),
    data.table(code=1393, label="venous surgery/procedures"),
    data.table(code=1402, label="hernia surgery"),
    data.table(code=1403, label="inguinal/femoral hernia repair"),
    data.table(code=1404, label="umbilical hernia repair"),
    data.table(code=1405, label="incisional hernia repair"),
    data.table(code=1427, label="bladder surgery"),
    data.table(code=1428, label="transurethral resection bladder tumour (turbt)"),
    data.table(code=1429, label="cystectomy"),
    data.table(code=1431, label="urethral surgery"),
    data.table(code=1432, label="thyroidectomy/partial thyroidectomy"),
    data.table(code=1434, label="eye surgery"),
    data.table(code=1435, label="cataract extraction/lens implant"),
    data.table(code=1436, label="glaucoma surgery/trabeculectomy"),
    data.table(code=1437, label="retinal operation/vitrectomy"),
    data.table(code=1438, label="bone marrow transplant"),
    data.table(code=1439, label="reduction or fixation of bone fracture"),
    data.table(code=1440, label="amputation of leg"),
    data.table(code=1441, label="amputation of foot"),
    data.table(code=1442, label="amputation of toe"),
    data.table(code=1443, label="other amputation"),
    data.table(code=1444, label="lung removal/pneumonectomy/lobectomy"),
    data.table(code=1445, label="bullaectomy/bullae removal"),
    data.table(code=1446, label="lung transplant"),
    data.table(code=1447, label="other thoracic surgery"),
    data.table(code=1448, label="stomach surgery"),
    data.table(code=1449, label="oesophageal surgery"),
    data.table(code=1450, label="liver surgery"),
    data.table(code=1451, label="biliary surgery"),
    data.table(code=1452, label="pancreas surgery"),
    data.table(code=1453, label="splenectomy"),
    data.table(code=1454, label="liver transplant"),
    data.table(code=1455, label="cholecystectomy/gall bladder removal"),
    data.table(code=1456, label="bile duct surgery"),
    data.table(code=1457, label="pancreatic transplant"),
    data.table(code=1458, label="appendicectomy"),
    data.table(code=1459, label="bowel resection"),
    data.table(code=1460, label="anal surgery"),
    data.table(code=1461, label="large bowel resection +/- colostomy"),
    data.table(code=1462, label="small bowel resection"),
    data.table(code=1463, label="rectal or colon polypectomy"),
    data.table(code=1464, label="colectomy/hemicolectomy"),
    data.table(code=1465, label="rectal/sigmoid resection"),
    data.table(code=1466, label="other bowel surgery"),
    data.table(code=1467, label="brain surgery"),
    data.table(code=1468, label="intracranial haematoma drainage"),
    data.table(code=1469, label="cerebroventricular shunt/shunt for hydrocephalus"),
    data.table(code=1470, label="spinal cord surgery"),
    data.table(code=1471, label="peripheral nerve surgery"),
    data.table(code=1472, label="pituitary surgery"),
    data.table(code=1473, label="pituitary transplant"),
    data.table(code=1474, label="skin operation or plastic surgery"),
    data.table(code=1475, label="fertility treatment procedures"),
    data.table(code=1476, label="fistula for dialysis"),
    data.table(code=1477, label="caesarean section / caesarian section"),
    data.table(code=1478, label="tonsillectomy +/- adenoids"),
    data.table(code=1479, label="varicose vein surgery"),
    data.table(code=1480, label="wisdom teeth surgery"),
    data.table(code=1481, label="maxillo-facial surgery"),
    data.table(code=1482, label="oesphageal fundoplication/hiatus hernia surgery"),
    data.table(code=1483, label="haemorroidectomy / piles surgery/ banding of piles"),
    data.table(code=1484, label="pilonidal sinus surgery (anal)"),
    data.table(code=1485, label="abdominal/pelvic adhesion surgery"),
    data.table(code=1486, label="laparotomy nos"),
    data.table(code=1487, label="nephrectomy/kidney removed"),
    data.table(code=1488, label="anterior/posterior repair bladder/uterus"),
    data.table(code=1489, label="testicular hydrocoele surgery / drainage"),
    data.table(code=1490, label="male circumcision"),
    data.table(code=1491, label="squint correction"),
    data.table(code=1492, label="spinal laminectomy"),
    data.table(code=1493, label="arthroscopy nos"),
    data.table(code=1494, label="foot surgery"),
    data.table(code=1495, label="lower limb surgery"),
    data.table(code=1496, label="knee surgery (not replacement)"),
    data.table(code=1497, label="hip surgery (not replacement)"),
    data.table(code=1498, label="upper limb surgery"),
    data.table(code=1499, label="shoulder surgery"),
    data.table(code=1500, label="elbow surgery"),
    data.table(code=1501, label="carpal tunnel surgery"),
    data.table(code=1502, label="cruciate ligament surgery"),
    data.table(code=1503, label="bunion/hallus valgus surgery"),
    data.table(code=1504, label="removal of mole/skin lesion"),
    data.table(code=1505, label="ovarian/tubal surgery"),
    data.table(code=1506, label="ovarian cyst removal/surgery"),
    data.table(code=1507, label="ectopic pregnancy surgery"),
    data.table(code=1508, label="uterine surgery"),
    data.table(code=1509, label="myomectomy/fibroids removed"),
    data.table(code=1510, label="dilatation and curettage/d+c"),
    data.table(code=1511, label="cervix surgery"),
    data.table(code=1512, label="cone biopsy"),
    data.table(code=1513, label="breast cyst/abscess removal"),
    data.table(code=1514, label="coronary angiogram"),
    data.table(code=1515, label="angiogram not coronary"),
    data.table(code=1516, label="echocardiogram/transoesophageal echo"),
    data.table(code=1517, label="endoscopy / gastroscopy"),
    data.table(code=1518, label="ercp/endoscopic retrograde cholangiopancreatogram"),
    data.table(code=1519, label="colonoscopy/sigmoidoscopy"),
    data.table(code=1520, label="laparoscopy"),
    data.table(code=1521, label="hysteroscopy +/- biopsy"),
    data.table(code=1522, label="colposcopy +/- biopsy"),
    data.table(code=1523, label="triple heart bypass"),
    data.table(code=1524, label="cardioversion"),
    data.table(code=1525, label="mastoidectomy/mastoid surgery"),
    data.table(code=1526, label="stapedectomy"),
    data.table(code=1527, label="dental/tooth surgery"),
    data.table(code=1528, label="gallstones removed"),
    data.table(code=1529, label="duodenal surgery"),
    data.table(code=1530, label="ileostomy/colostomy surgery"),
    data.table(code=1531, label="discectomy"),
    data.table(code=1532, label="ankle surgery"),
    data.table(code=1533, label="achilles tendon repair"),
    data.table(code=1534, label="hand/finger surgery"),
    data.table(code=1535, label="dupuytren's contracture surgery"),
    data.table(code=1536, label="lipoma removed / excision of lipoma"),
    data.table(code=1537, label="lymph node surgery"),
    data.table(code=1538, label="salpingectomy"),
    data.table(code=1539, label="uterine polypectomy/uterine polyps removed"),
    data.table(code=1540, label="termination of pregnancy/top"),
    data.table(code=1541, label="cervical polyps removed"),
    data.table(code=1542, label="laser treatment cervix"),
    data.table(code=1543, label="pelvic floor surgery"),
    data.table(code=1544, label="breast biopsy"),
    data.table(code=1545, label="ecg/electrocardiogram"),
    data.table(code=1546, label="bronchoscopy"),
    data.table(code=1547, label="cystoscopy"),
    data.table(code=1548, label="pacemaker insertion"),
    data.table(code=1549, label="pacemaker battery change"),
    data.table(code=1550, label="defibrillator/icd insertion"),
    data.table(code=1551, label="defibrillator/icd battery change"),
    data.table(code=1552, label="electrophysiological studies (eps)"),
    data.table(code=1553, label="cardiac ablation"),
    data.table(code=1554, label="pericardial surgery"),
    data.table(code=1555, label="femoral/popliteal/iliac aneurysm repair"),
    data.table(code=1556, label="tympanic membrane surgery / ear drum repair"),
    data.table(code=1557, label="rhinoplasty / nose surgery"),
    data.table(code=1558, label="sinus surgery"),
    data.table(code=1559, label="nasal polyp surgery / nasal polypectomy"),
    data.table(code=1560, label="parotid surgery / parotidectomy"),
    data.table(code=1561, label="tonsillectomy / tonsil surgery"),
    data.table(code=1562, label="adenoid surgery / adenoidectomy"),
    data.table(code=1563, label="inguinal hernia repair"),
    data.table(code=1564, label="femoral hernia repair"),
    data.table(code=1565, label="gastrectomy"),
    data.table(code=1566, label="peptic ulcer surgery"),
    data.table(code=1567, label="gastric ulcer surgery"),
    data.table(code=1568, label="duodenal ulcer surgery"),
    data.table(code=1570, label="ileostomy surgery"),
    data.table(code=1571, label="ileostomy formation"),
    data.table(code=1572, label="ileostomy reversal"),
    data.table(code=1573, label="surgery for pyloric stenosis"),
    data.table(code=1574, label="colostomy surgery"),
    data.table(code=1575, label="colostomy formation"),
    data.table(code=1576, label="colostomy reversal"),
    data.table(code=1577, label="anterior resection"),
    data.table(code=1578, label="abdomino-perineal resection / ap resection"),
    data.table(code=1579, label="anal fissure surgery"),
    data.table(code=1580, label="dialysis access surgery"),
    data.table(code=1581, label="haemodialysis access / fistula surgery"),
    data.table(code=1582, label="peritoneal dialysis (capd) access surgery"),
    data.table(code=1583, label="urethral stricture surgery / dilatation"),
    data.table(code=1584, label="testicular lump removal"),
    data.table(code=1586, label="varicocoele / varicocele surgery"),
    data.table(code=1587, label="removal of epididymal cyst"),
    data.table(code=1588, label="orchidopexy"),
    data.table(code=1589, label="reversal of vasectomy / vasectomy reversal"),
    data.table(code=1590, label="spinal fusion"),
    data.table(code=1591, label="eyelid surgery"),
    data.table(code=1592, label="nasolacrimal / tear duct surgery"),
    data.table(code=1593, label="removal of malignant melanoma"),
    data.table(code=1594, label="removal of benign skin lesion"),
    data.table(code=1595, label="removal of squamous cell carcinoma (scc)"),
    data.table(code=1596, label="removal of rodent ulcer / basal cell carcinoma (bcc)"),
    data.table(code=1597, label="skin graft"),
    data.table(code=1598, label="nail surgery / nail removal"),
    data.table(code=1599, label="exercise ecg"),
    data.table(code=1600, label="mri / magnetic resonance imaging"),
    data.table(code=1601, label="laryngoscopy"),
    data.table(code=1602, label="nasal endoscopy"),
    data.table(code=1603, label="scan / x-ray"),
    data.table(code=1604, label="cardiac investigation"),
    data.table(code=1605, label="gastrointestinal investigation"),
    data.table(code=1606, label="genito-urinary investigation"),
    data.table(code=1607, label="ct scan"),
    data.table(code=1608, label="bone scan"),
    data.table(code=1609, label="ultrasound scan"),
    data.table(code=1610, label="pet scan"),
    data.table(code=1611, label="barium meal / barium swallow"),
    data.table(code=1612, label="barium enema"),
    data.table(code=1613, label="ct colonoscopy"),
    data.table(code=1614, label="biopsy"),
    data.table(code=1615, label="bone marrow biopsy"),
    data.table(code=1616, label="liver biopsy"),
    data.table(code=1617, label="prostate biopsy"),
    data.table(code=1618, label="renal biopsy / kidney biopsy"),
    data.table(code=1619, label="lung biopsy"),
    data.table(code=1620, label="lumbar puncture"),
    data.table(code=1621, label="cystoscopy"),
    data.table(code=99999, label="unclassifiable"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer"),
    data.table(code=-7, label="None of the above")
  ),
  "20126" = rbind(
		data.table(code=0, label="No Bipolar or Depression"),
		data.table(code=1, label="Bipolar I Disorder"),
		data.table(code=2, label="Bipolar II Disorder"),
		data.table(code=3, label="Probable Recurrent major depression (severe)"),
		data.table(code=4, label="Probable Recurrent major depression (moderate)"),
		data.table(code=5, label="Single Probable major depression episode")
  ), 
  "20122" = rbind(
    data.table(code=1, label="Bipolar Type I (Mania)"),
    data.table(code=2, label="Bipolar Type II (Hypomania)"),
    data.table(code=0, label="No Bipolar")
  ), 
  "20124" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No") 
  ),
  "20125" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No") 
  ),
  "20123" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No") 
  ),
  "1920" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "1930" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "1940" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "1950" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "1960" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "1970" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "1980" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "1990" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "2000" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "2010" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "2020" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "2030" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "2040" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "4526" = rbind(
    data.table(code=1, label="Extremely happy"),
    data.table(code=2, label="Very happy"),
    data.table(code=3, label="Moderately happy"),
    data.table(code=4, label="Moderately unhappy"),
    data.table(code=5, label="Very unhappy"),
    data.table(code=6, label="Extremely unhappy"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "4537" = rbind(
    data.table(code=1, label="Extremely happy"),
    data.table(code=2, label="Very happy"),
    data.table(code=3, label="Moderately happy"),
    data.table(code=4, label="Moderately unhappy"),
    data.table(code=5, label="Very unhappy"),
    data.table(code=6, label="Extremely unhappy"),
    data.table(code=7, label="I am not employed"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "4548" = rbind(
    data.table(code=1, label="Extremely happy"),
    data.table(code=2, label="Very happy"),
    data.table(code=3, label="Moderately happy"),
    data.table(code=4, label="Moderately unhappy"),
    data.table(code=5, label="Very unhappy"),
    data.table(code=6, label="Extremely unhappy"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "4559" = rbind(
    data.table(code=1, label="Extremely happy"),
    data.table(code=2, label="Very happy"),
    data.table(code=3, label="Moderately happy"),
    data.table(code=4, label="Moderately unhappy"),
    data.table(code=5, label="Very unhappy"),
    data.table(code=6, label="Extremely unhappy"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "4570" = rbind(
    data.table(code=1, label="Extremely happy"),
    data.table(code=2, label="Very happy"),
    data.table(code=3, label="Moderately happy"),
    data.table(code=4, label="Moderately unhappy"),
    data.table(code=5, label="Very unhappy"),
    data.table(code=6, label="Extremely unhappy"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "4581" = rbind(
    data.table(code=1, label="Extremely happy"),
    data.table(code=2, label="Very happy"),
    data.table(code=3, label="Moderately happy"),
    data.table(code=4, label="Moderately unhappy"),
    data.table(code=5, label="Very unhappy"),
    data.table(code=6, label="Extremely unhappy"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "2050" = rbind(
    data.table(code=1, label="Not at all"),
    data.table(code=2, label="Several days"),
    data.table(code=3, label="More than half the days"),
    data.table(code=4, label="Nearly every day"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ), 
  "2060" = rbind(
    data.table(code=1, label="Not at all"),
    data.table(code=2, label="Several days"),
    data.table(code=3, label="More than half the days"),
    data.table(code=4, label="Nearly every day"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ), 
  "2070" = rbind(
    data.table(code=1, label="Not at all"),
    data.table(code=2, label="Several days"),
    data.table(code=3, label="More than half the days"),
    data.table(code=4, label="Nearly every day"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ), 
  "2080" = rbind(
    data.table(code=1, label="Not at all"),
    data.table(code=2, label="Several days"),
    data.table(code=3, label="More than half the days"),
    data.table(code=4, label="Nearly every day"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ), 
  "2090" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "2100" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "4598" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "4631" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  
  "4642" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "4653" = rbind(
    data.table(code=1, label="Yes"),
    data.table(code=0, label="No"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "6156" = rbind(
    data.table(code=11, label="I was more active than usual"),
    data.table(code=12, label="I was more talkative than usual"),
    data.table(code=13, label="I needed less sleep than usual"),
    data.table(code=14, label="I was more creative or had more ideas than usual"),
    data.table(code=15, label="All of the above"),
    data.table(code=-7, label="None of the above")
  ), 
  "5663" = rbind(
    data.table(code=11, label="At least two days, but less than a week"),
    data.table(code=12, label="Less than a week"),
    data.table(code=13, label="A week or more"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "5674" = rbind(
    data.table(code=11, label="No problems"),
    data.table(code=12, label="Needed treatment or caused problems with work, relationships, finances, the law or other aspects of life"),
    data.table(code=-1, label="Do not know"),
    data.table(code=-3, label="Prefer not to answer")
  ),
  "6145" = rbind(
    data.table(code=1, label="Serious illness, injury or assault to yourself"),
    data.table(code=2, label="Serious illness, injury or assault of a close relative"),
    data.table(code=3, label="Death of a close relative"),
    data.table(code=4, label="Death of a spouse or partner"),
    data.table(code=5, label="Marital separation/divorce"),
    data.table(code=6, label="Financial difficulties"),
    data.table(code=-7, label="None of the above"),
    data.table(code=-3, label="Prefer not to answer")
  )
)

# Some are numeric/integer with special codes, populate
populate_numeric <- function(fid, specials) {
  possible <- unique(codes[field.id == fid, .(field.id=fid, code=code, label=as.character(code))])
  possible <- possible[order(code)]
	if (any(possible$code < 0)) {
		possible <- rbind(possible[code >= 0], possible[code < 0][.N:1])
  }
  if (!missing(specials)) {
		specials <- data.table(code=as.integer(names(specials)), label=as.vector(specials))
		possible[specials, on = .(code), label := i.label]
  }
  possible
}

labels <- rbind(labels,
  populate_numeric(20127),
  populate_numeric(4609, c("-1"="Do not know", "-3"="Prefer not to answer")),
  populate_numeric(4620, c("-1"="Do not know", "-3"="Prefer not to answer")),
  populate_numeric(5375, c("-1"="Do not know", "-3"="Prefer not to answer")),
  populate_numeric(5386, c("-1"="Do not know", "-3"="Prefer not to answer"))
)

# Add labels
codes[labels, on = .(field.id, code), label := i.label]

# Curate data.table of age at diagnosis
touchscreen_age_at_diagnosis <- rbind(
  cbind(field.id = "6150", code = 1, label = "Heart attack", raw[["3894"]]),
  cbind(field.id = "6150", code = 2, label = "Angina", raw[["3627"]]),
  cbind(field.id = "6150", code = 3, label = "Stroke", raw[["4056"]]),
  cbind(field.id = "6150", code = 4, label = "High blood pressure", raw[["2966"]]),
  cbind(field.id = "6152", code = 5, label = "Blood clot in the leg (DVT)", raw[["4012"]]),
  cbind(field.id = "6152", code = 7, label = "Blood clot in the lung", raw[["4022"]]),
  cbind(field.id = "6152", code = 6, label = "Emphysema/chronic bronchitis", raw[["3992"]]),
  cbind(field.id = "6152", code = 8, label = "Asthma", raw[["3786"]]),
  cbind(field.id = "6152", code = 9, label = "Hayfever, allergic rhinitis or eczema", raw[["3761"]]),
  cbind(field.id = "2443", code = 1, label = "Yes", raw[["2976"]]),
  cbind(field.id = "4041", code = 1, label = "Yes", raw[["2976"]]) # note no matches despite what UKB showcase implies
)
codes[touchscreen_age_at_diagnosis, on = .(eid, visit_index, field.id, code, label), age := as.numeric(i.value)]

interview_age_at_diagnosis <- rbind(
  cbind(field.id = "20002", raw[["20009"]]),
  cbind(field.id = "20001", raw[["20007"]]),
  cbind(field.id = "20004", raw[["20011"]])
)
codes[interview_age_at_diagnosis, on = .(eid, visit_index, repeat_index, field.id), age := i.value]

codes[age < 0, age := NA] # Dont know / Prefer not to answer set to NA

# Curate date of diagnosis
interview_year_at_diagnosis <- rbind(
  cbind(field.id = "20002", raw[["20008"]]),
  cbind(field.id = "20001", raw[["20006"]]),
  cbind(field.id = "20004", raw[["20010"]])
)
interview_year_at_diagnosis[value < 0, value := NA]
interview_year_at_diagnosis[, value := as.IDate(date_decimal(value))]
codes[interview_year_at_diagnosis, on = .(eid, visit_index, repeat_index, field.id), date := i.value]

# Add in age and date at assessment so we can impute dates and other missing age data
assessment <- fread(
  file = "data/curated/ukbiobank/anthropometrics/anthropometrics.txt", 
  select = c("eid", "visit_index", "age_decimal", "assessment_date"),
  colClasses = c("visit_index"="character")
) 
codes[assessment, on = .(eid, visit_index), age_assessment := i.age_decimal] # age decimal incorporates birth month, with specific day imputed as 15th of month
codes[assessment, on = .(eid, visit_index), date_assessment := i.assessment_date]

# Follow verbal interview method for interpolating date of assessment:
#
#  - If the participant gave their age then the value presented is the fractional 
#    year corresponding to the mid-point of that age. For example, if the participant 
#    said they were 30 years old then the value is the date at which they were 
#    30 years + 6 months.
#  - Interpolated values before the date of birth were truncated forwards to that time.
#  - Interpolated values after the time of data acquisition were truncated back to that time.
#
# E.g. see Notes field of https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20008
#
# Note precise date of birth is not known (censored for privacy reasons), however, there are no
# age of event within the first year of life for touchscreen survey questions (for which we're 
# imputing dates here)
codes[!is.na(age) & is.na(date), age := age + 0.5]
codes[!is.na(age) & is.na(date), date := as.IDate(date_decimal(decimal_date(date_assessment) - (age_assessment - age)))]
codes[date > date_assessment, date := date_assessment]
codes[age > age_assessment, age := age_assessment]

# Impute age/date insulin started (within 1 year of diabetes diagnosis)
diabetes <- codes[field.id == "2443" & label == "Yes"]
diabetes[, field.id := "2986"] # now matches to start insulin treatment field
codes[diabetes, on = .(field.id, eid, visit_index, label), age := i.age + 0.5] # set at mid-point of 1 year followup
codes[field.id == "2986" & label == "Yes" & !is.na(age), date := as.IDate(date_decimal(decimal_date(date_assessment) - (age_assessment - age)))]
codes[field.id == "2986" & label == "Yes" & date > date_assessment, date := date_assessment]
codes[field.id == "2986" & label == "Yes" & age > age_assessment, age := age_assessment]

# Fields 2463, 6151, and 3005 relate to the question 2463 "Fractured bones in the last 5 years".
# In the absence of a precise date, impute as the mid-point (i.e. 2.5 years prior to assessment).
codes[field.id %in% c("2463", "6151", "3005") & code > 0, age := age_assessment - 2.5]
codes[field.id %in% c("2463", "6151", "3005") & code > 0, date := as.IDate(date_decimal(decimal_date(date_assessment) - (age_assessment - 2.5)))]

# Fields 2050, 2060, 2070, 2080 relate to mood in previous two weeks, set to midpoint of 7 days prior to assessment
codes[field.id %in% c("2050", "2060", "2070", "2080") & code > 1, age := age_assessment - days(7)/years(1)]
codes[field.id %in% c("2050", "2060", "2070", "2080") & code > 1, date := as.IDate(date_decimal(decimal_date(date_assessment) - (age_assessment - days(7)/years(1))))]

# Field 6145 relates to stress in previous two years, set to midpoint of 1 year prior to assessment
codes[field.id == "6145" & code > 0, age := age_assessment - 1]
codes[field.id == "6145" & code > 0, date := as.IDate(date_decimal(decimal_date(date_assessment) - (age_assessment - 1)))]

# Ages and dates for other fields cannot be reasonably imputed.

# Filter codes table to fields we want to write out
codes <- codes[, .(eid, visit_index, field_id=field.id, code, age, date)]

# Set field order so that when multiple events are matched, the most informative
# one is reported by 'curate_endpoints.R' (see top-level README.txt)
codes[, field_id := factor(field_id, levels=c(
  # Verbal interview with nurse,
	20001,  # Cancer-related illnesses
	20002,  # Non-cancer related illnesses
	20004,  # Medical operations

	# Touchscreen survey questions relating to diabetes
	2443,   # "Has a doctor ever told you that you have diabetes?"
	2986,   # "Did you start insulin within one year of your diagnosis of diabetes?"
	4041,   # "Did you only have diabetes during pregnancy?"

	# Touchscreen survey questions relating to vascular/heart problems diagnosed by doctor
	6150,   # "Has a doctor ever told you that you have had any of the following conditions? (You can select more than one answer)"

	# Touchscreen survey questions relating to blood clot, DVT, bronchitis, emphysema, asthma, rhinitis, eczema, allergy
	# diagnosed by doctor
	6152,   # "Has a doctor ever told you that you have had any of the following conditions? (You can select more than one answer)"

	# Touchscreen survey questions relating to claudication and peripheral artery disease
	4728,   # "Do you get a pain in either leg on walking?"
	5452,   # "Does this pain ever begin when you are standing still or sitting?"
	5463,   # "Do you get this pain in your calf (calves)?"
	5474,   # "Do you get pain when you walk uphill or hurry?"
	5485,   # "Do you get pain when you walk at an ordinary pace on the level?"
	5496,   # "Does the pain you get while walking ever disappear when you continue walking?"
	5507,   # "What do you do if you get pain when you are walking?"
	5518,   # "What happens to the pain you get while walking if you stand still?"
	5529,   # "Have you ever had surgery on the arteries of your legs (other than for varicose veins)?"
	5540,   # "Have you ever had surgery to remove any of the following?"

	# Touchscreen survey questions relating to bone fractures,
	2463,   # "Have you fractured/broken any bones in the last 5 years?"
	3005,   # "Did the fracture result from a simple fall (i.e. from standing height)?"
	6151,   # "Which bones did you fracture/break? (You can select more than one answer)"

  # Touchscreen survey questions relating to mental health:
  20126,  # "Bipolar and major depression status"
  20122,  # "Bipolar disorder status"
  20127,  # "Neuroticism score"
  20124,  # "Probable recurrent major depression (moderate)"
  20125,  # "Probable recurrent major depression (severe)"
  20123,  # "Single episode of probable major depression"
  1920,   # "Mood swings"
  1930,   # "Miserableness"
  1940,   # "Irritability"
  1950,   # "Sensitivity / hurt feelings"
  1960,   # "Fed-up feelings"
  1970,   # "Nervous feelings"
  1980,   # "Worrier / anxious feelings"
  1990,   # "Tense / 'highly strung'"
  2000,   # "Worry too long after embarrassment"
  2010,   # "Suffer from 'nerves'"
  2020,   # "Loneliness, isolation"
  2030,   # "Guilty feelings"
  2040,   # "Risk taking"
  4526,   # "Happiness"
  4537,   # "Work/job satisfaction"
  4548,   # "Health satisfaction"
  4559,   # "Family relationship satisfaction"
  4570,   # "Friendships satisfaction"
  4581,   # "Financial situation satisfaction"
  2050,   # "Frequency of depressed mood in last 2 weeks"
  2060,   # "Frequency of unenthusiasm / disinterest in last 2 weeks"
  2070,   # "Frequency of tenseness / restlessness in last 2 weeks"
  2080,   # "Frequency of tiredness / lethargy in last 2 weeks"
  2090,   # "Seen doctor (GP) for nerves, anxiety, tension or depression"
  2100,   # "Seen a psychiatrist for nerves, anxiety, tension or depression"
  4598,   # "Ever depressed for a whole week"
  4609,   # "Longest period of depression"
  4620,   # "Number of depression episodes"
  4631,   # "Ever unenthusiastic/disinterested for a whole week"
  5375,   # "Longest period of unenthusiasm / disinterest"
  5386,   # "Number of unenthusiastic/disinterested episodes"
  4642,   # "Ever manic/hyper for 2 days"
  4653,   # "Ever highly irritable/argumentative for 2 days"
  6156,   # "Manic/hyper symptoms"
  5663,   # "Length of longest manic/irritable episode"
  5674,   # "Severity of manic/irritable episodes"
  6145    # "Illness, injury, bereavement, stress in last 2 years"
))]

codes <- codes[order(code)][order(field_id)][order(visit_index)][order(eid)]

# Write out
fwrite(codes, sep="\t", quote=FALSE, file="data/curated/ukbiobank/self_report/self_report_codes.txt")

# Write out code labels
setnames(labels, "field.id", "field_id")
labels[, field_id := factor(field_id, levels=levels(codes$field_id))]
labels <- labels[order(field_id)]
fwrite(labels, sep="\t", quote=FALSE, file="data/curated/ukbiobank/self_report/code_labels.txt")

# Curate detailed field information
info <- rbind(
	data.table(field_id=20001, codes="Cancer code (verbal interview), missing data curated from field #2453",
						 age="From field #20007", date="From field #20006"),
	data.table(field_id=20002, codes="Non-cancer illness code (verbal interview), missing data curated from field #2473",
						 age="From field #20009", date="From field #20008"),
	data.table(field_id=20004, codes="Operation code (verbal interview), missing data curated from fields #2415 and #2844",
						 age="From field #20011", date="From field #20010"),
	data.table(field_id=2443, codes="Diabetes diagnosed by doctor (touchscreen)",
						 age="From field #2976 plus half a year", date="Imputed from age at event"),
	data.table(field_id=2986, codes="Started insulin within one year diagnosis of diabetes (touchscreen)",
						 age="Age of diabetes diagnosis (see entry for field_id 4041) plus half a year", date="Imputed from age at event"),
	data.table(field_id=4041, codes="Gestational diabetes only (touchscreen)",
						 age="No age at event recorded in field #2976", date="None available"),
	data.table(field_id=6150, codes="Vascular/heart problems diagnosed by doctor (touchscreen)",
						 age="From fields #3894, #3627, #4056, and #2966 plus half a year", date="Imputed from age at event"),
	data.table(field_id=6152, codes="Blood clot, DVT, bronchitis, emphysema, asthma, rhinitis, eczema, allergy diagnosed by doctor (touchscreen)",
						 age="From fields #2966, #4012, #4022, #3992, #3786, and #3761 plus half a year", date="Imputed from age at event"),
	data.table(field_id=4728, codes="Leg pain on walking (touchscreen)",
						 age="Could not be imputed", date="None available"),
	data.table(field_id=5452, codes="Leg pain when standing still or sitting (touchscreen)",
						 age="Could not be imputed", date="None available"),
	data.table(field_id=5463, codes="Leg pain in calf/calves (touchscreen)",
						 age="Could not be imputed", date="None available"),
	data.table(field_id=5474, codes="Leg pain when walking uphill or hurrying (touchscreen)",
						 age="Could not be imputed", date="None available"),
	data.table(field_id=5485, codes="Leg pain when walking normally (touchscreen)",
						 age="Could not be imputed", date="None available"),
	data.table(field_id=5496, codes="Leg pain when walking ever disappears while walking (touchscreen)",
						 age="Could not be imputed", date="None available"),
	data.table(field_id=5507, codes="Leg pain on walking : action taken (touchscreen)",
						 age="Could not be imputed", date="None available"),
	data.table(field_id=5518, codes="Leg pain on walking : effect of standing still (touchscreen)",
						 age="Could not be imputed", date="None available"),
	data.table(field_id=5529, codes="Surgery on leg arteries (other than for varicose veins) (touchscreen)",
						 age="Could not be imputed", date="None available"),
	data.table(field_id=5540, codes="Surgery/amputation of toe or leg (touchscreen)",
						 age="Could not be imputed", date="None available"),
	data.table(field_id=2463, codes="Fractured/broken bones in last 5 years (touchscreen)",
						 age="Age of assessment minus 2.5 years", date="Imputed from age at event"),       
	data.table(field_id=3005, codes="Fracture resulting from simple fall (touchscreen)",
						 age="Age of assessment minus 2.5 years", date="Imputed from age at event"),       
	data.table(field_id=6151, codes="Fractured bone site(s) (touchscreen)",
						 age="Age of assessment minus 2.5 years", date="Imputed from age at event"),
  data.table(field_id=20126, codes="Bipolar and major depression status (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=20122, codes="Bipolar disorder status (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=20127, codes="Neuroticism score (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=20124, codes="Probable recurrent major depression (moderate) (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=20125, codes="Probable recurrent major depression (severe) (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=20123, codes="Single episode of probable major depression (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=1920, codes="Mood swings (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=1930, codes="Miserableness (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=1940, codes="Irritability (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=1950, codes="Sensitivity / hurt feelings (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=1960, codes="Fed-up feelings (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=1970, codes="Nervous feelings (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=1980, codes="Worrier / anxious feelings (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=1990, codes="Tense / 'highly strung' (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=2000, codes="Worry too long after embarrassment (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=2010, codes="Suffer from 'nerves' (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=2020, codes="Loneliness, isolation (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=2030, codes="Guilty feelings (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=2040, codes="Risk taking (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=4526, codes="Happiness (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=4537, codes="Work/job satisfaction (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=4548, codes="Health satisfaction (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=4559, codes="Family relationship satisfaction (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=4570, codes="Friendships satisfaction (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=4581, codes="Financial situation satisfaction (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=2050, codes="Frequency of depressed mood in last 2 weeks (touchscreen)",
             age="Age of assessment minus 7 days", date="Imputed from age at event"),
  data.table(field_id=2060, codes="Frequency of unenthusiasm / disinterest in last 2 weeks (touchscreen)",
             age="Age of assessment minus 7 days", date="Imputed from age at event"),
  data.table(field_id=2070, codes="Frequency of tenseness / restlessness in last 2 weeks (touchscreen)",
             age="Age of assessment minus 7 days", date="Imputed from age at event"),
  data.table(field_id=2080, codes="Frequency of tiredness / lethargy in last 2 weeks (touchscreen)",
             age="Age of assessment minus 7 days", date="Imputed from age at event"),
  data.table(field_id=2090, codes="Seen doctor (GP) for nerves, anxiety, tension or depression (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=2100, codes="Seen a psychiatrist for nerves, anxiety, tension or depression (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=4598, codes="Ever depressed for a whole week (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=4609, codes="Longest period of depression (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=4620, codes="Number of depression episodes (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=4631, codes="Ever unenthusiastic/disinterested for a whole week (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=5375, codes="Longest period of unenthusiasm / disinterest (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=5386, codes="Number of unenthusiastic/disinterested episodes (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=4642, codes="Ever manic/hyper for 2 days (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=4653, codes="Ever highly irritable/argumentative for 2 days (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=6156, codes="Manic/hyper symptoms (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=5663, codes="Length of longest manic/irritable episode (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=5674, codes="Severity of manic/irritable episodes (touchscreen)",
             age="Could not be imputed", date="None available"),
  data.table(field_id=6145, codes="Illness, injury, bereavement, stress in last 2 years (touchscreen)",
             age="Age of assessment minus 1 year", date="Imputed from age at event")
)

fwrite(info, sep="\t", quote=FALSE, file="data/curated/ukbiobank/self_report/field_information.txt")



