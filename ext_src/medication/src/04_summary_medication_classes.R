library(data.table)

meds <- fread("output/detailed_medications_field_20003.txt")

med_classes <- meds[, by=.(eid, visit_index), .(
	# Curate data on lipid lowering medications
	#
	# List of drugs from BNF chapter 2.12: Lipid-regulating drugs
	# https://openprescribing.net/bnf/0212/
	#
	# Note not all drugs listed in BNF are found in the list of UKB
	# medication codes. All drugs have also been cross-checked against
	# drugbank to ensure they are capturing the right condition (i.e.
	# lipid lowering to treat hypercholesterolemia and prevent CVD)
	#
	# E.g., while Ispaghula husk is listed, is primary indication is for constipation,
	# not lipid lowering, so its not included.
	# https://bnf.nice.org.uk/drug/ispaghula-husk.html
	#
	# Conversely, colestyramine is listed, even though it has multiple indications,
	# because one of its main indications is for primary prevention of coronary heart disease
	# https://bnf.nice.org.uk/drug/colestyramine.html
	#
	lipid_lowering_medication = ifelse(
		any(medication_code == 1140861892) |  # acipimox
		any(medication_code == 1141146234) |  # atorvastatin
		any(medication_code == 1140861924) |  # bezafibrate
		any(medication_code == 1141157260) |  # bezafibrate product
		any(medication_code == 1140862026) |  # ciprofibrate
		any(medication_code == 1140888590) |  # colestipol
		any(medication_code == 1140909780) |  # colestyramine
		any(medication_code == 1141180734) |  # colestyramine product
		any(medication_code == 1141180722) |  # colestyramine+aspartame 4g/sachet powder
		any(medication_code == 1141192736) |  # ezetimibe
		any(medication_code == 1140861954) |  # fenofibrate
		any(medication_code == 1140888594) |  # fluvastatin
		any(medication_code == 1140861856) |  # gemfibrozil
		any(medication_code == 1141157262) |  # gemfibrozil product
		any(medication_code == 1140861868) |  # nicotinic acid product
		any(medication_code == 1140888648) |  # pravastatin
		any(medication_code == 1141192410) |  # rosuvastatin
		any(medication_code == 1140861958),   # simvastatin
		TRUE, FALSE),

	# Curate data on treated hypertension
	#
	# List of drugs from BNF chapter 2.5: Hypertension and heart failure
	# https://openprescribing.net/bnf/0205/
	#
	# Cross-checked with drug indication at https://bnf.nice.org.uk/drug to ensure
	# those specifically for hypertension (e.g. instead of heart failure) are curated
	#
	# Note not all drugs listed in BNF are found in the list of UKB medication codes.
	#
	hypertension_medication = ifelse(
		any(medication_code == 1141186674) |  # bosentan
		any(medication_code == 1140888686) |  # hydralazine
		any(medication_code == 1140860532) |  # minoxidil
		any(medication_code == 1141168936) |  # sildenafil
		any(medication_code == 1141187810) |  # tadalafil
		any(medication_code == 1140883468) |  # clonidine
		any(medication_code == 1140871986) |  # clonidine hydrochloride 25micrograms tablet
		any(medication_code == 1140860470) |  # methyldopa
		any(medication_code == 1140860562) |  # methyldopa+hydrochlorothiazide 250mg/15mg tablet
		any(medication_code == 1140910606) |  # alpha methyldopa
		any(medication_code == 1140928284) |  # moxonidine
		any(medication_code == 1140888536) |  # guanethidine
		any(medication_code == 1140879778) |  # doxazosin
		any(medication_code == 1140879782) |  # indoramin
		any(medication_code == 1141157490) |  # indoramin product
		any(medication_code == 1140879794) |  # prazosin
		any(medication_code == 1140879798) |  # terazosin
		any(medication_code == 1141156836) |  # candesartan cilexetil
		any(medication_code == 1140860750) |  # captopril
		any(medication_code == 1140860764) |  # captopril+hydrochlorothiazide 25mg/12.5mg tablet
		any(medication_code == 1140860882) |  # cilazapril
		any(medication_code == 1141181186) |  # co-zidocapt 25mg/12.5mg tablet
		any(medication_code == 1140888552) |  # enalapril
		any(medication_code == 1140860790) |  # enalapril maleate+hydrochlorothiazide 20mg/12.5mg tablet
		any(medication_code == 1141171336) |  # eprosartan
		any(medication_code == 1140888556) |  # fosinopril
		any(medication_code == 1141164148) |  # imidapril hydrochloride
		any(medication_code == 1141152998) |  # irbesartan
		any(medication_code == 1141172682) |  # irbesartan+hydrochlorothiazide 150mg/12.5mg tablet
		any(medication_code == 1140860696) |  # lisinopril
		any(medication_code == 1140864952) |  # lisinopril+hydrochlorothiazide 10mg/12.5mg tablet
		any(medication_code == 1140916356) |  # losartan
		any(medication_code == 1141151016) |  # losartan potassium+hydrochlorothiazide 50mg/12.5mg tablet
		any(medication_code == 1140923712) |  # moexipril
		any(medication_code == 1141193282) |  # olmesartan
		any(medication_code == 1140879802) |  # amlodipine
		any(medication_code == 1140888560) |  # perindopril
		any(medication_code == 1141180592) |  # perindopril+indapamide
		any(medication_code == 1140860728) |  # quinapril
		any(medication_code == 1140860806) |  # ramipril
		any(medication_code == 1141165470) |  # felodipine+ramipril
		any(medication_code == 1141166006) |  # telmisartan
		any(medication_code == 1141187788) |  # telmisartan+hydrochlorothiazide 40mg/12.5mg tablet
		any(medication_code == 1140860904) |  # trandolapril
		any(medication_code == 1141153328) |  # trandolapril+verapamil hydrochloride
		any(medication_code == 1140888510) |  # verapamil
		any(medication_code == 1141145660) |  # valsartan
		any(medication_code == 1141201038),   # valsartan+hydrochlorothiazide 80mg/12.5mg tablet
		TRUE, FALSE),

	# Curate data on insulin treatment
  # 
  # From BNF chapter 6.1.1
  # https://openprescribing.net/bnf/060101/
	insulin_medication = ifelse(
		any(medication_code == 1140883066),   # insulin product
		TRUE, FALSE)
)]

fwrite(med_classes, sep="\t", quote=FALSE, file="output/detailed_medications_summarised.txt")

