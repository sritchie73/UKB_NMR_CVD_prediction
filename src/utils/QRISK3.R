library(data.table)

# Adaptation of QRISK3 algorithm (https://qrisk.org/src.php), which splits
# out the linear predictor from the absolute risk % calculation
qrisk3 <- function(
  sex, age, ethnicity, townsend, smoking, sbp, sbp_sd, weight, height, chol_hdl_ratio, famhist_young_mi,
  atrial_fibrillation, CKD, erectile_dysfunction, mental_illness, migraine, rheumatoid_arthritis, SLE, type_1_diabetes, type_2_diabetes,
  atypical_antipsychotics, bp_treatment, corticosteroids, type="absolute risk"
) {
  # Basic error checking
  type <- match.arg(type, c("linear predictor", "absolute risk"))

  # combine inputs into a work in progress table
  wip <- data.table(
		sex, age, ethnicity, townsend, smoking, sbp, sbp_sd, weight, height, chol_hdl_ratio, famhist_young_mi,
		atrial_fibrillation, CKD, erectile_dysfunction, mental_illness, migraine, rheumatoid_arthritis, SLE, type_1_diabetes, type_2_diabetes,
		atypical_antipsychotics, bp_treatment, corticosteroids
  )

  # Compute common fractional polynomial transformations and derived variables
  wip[, age := age/10]
  wip[, height := height/100] 
  wip[, bmi := weight / height^2]
  wip[, bmi := bmi/10]
  wip[, bmi_1 := bmi^-2]
  wip[, bmi_2 := bmi^-2 * log(bmi)]

  # Apply sex-specific fractional polynomial transforms (includes scaling) and centering of continuous variabels
  wip[, age_1 := ifelse(sex == "Female", age^-2, age^-1)]
  wip[, age_2 := ifelse(sex == "Female", age, age^3)]
  
  wip[, age_1 := age_1 - ifelse(sex == "Female", 0.053274843841791, 0.234766781330109)]
  wip[, age_2 := age_2 - ifelse(sex == "Female", 4.332503318786621, 77.284080505371094)]
  wip[, bmi_1 := bmi_1 - ifelse(sex == "Female", 0.154946178197861, 0.149176135659218)]
  wip[, bmi_2 := bmi_2 - ifelse(sex == "Female", 0.144462317228317, 0.141913309693336)]
  wip[, chol_hdl_ratio := chol_hdl_ratio - ifelse(sex == "Female", 3.476326465606690, 4.300998687744141)]
  wip[, sbp := sbp - ifelse(sex == "Female", 123.130012512207030, 128.571578979492190)]
  wip[, sbp_sd := sbp_sd - ifelse(sex == "Female", 9.002537727355957, 8.756621360778809)]
  wip[, townsend := townsend - ifelse(sex == "Female", 0.392308831214905, 0.526304900646210)]

  # Compute sex-specific linear predictors
  wip[sex == "Female", QRISK3 :=
	  fcase(
      ethnicity == "White or not stated", 0, 
      ethnicity == "Indian", 0.2804031433299542500000000,
      ethnicity == "Pakistani", 0.5629899414207539800000000,
      ethnicity == "Bangladeshi", 0.2959000085111651600000000,
      ethnicity == "Other Asian", 0.0727853798779825450000000,
      ethnicity == "Black Caribbean", -0.1707213550885731700000000,
      ethnicity == "Black African", -0.3937104331487497100000000,
      ethnicity == "Chinese", -0.3263249528353027200000000,
      ethnicity == "Other ethnic group", -0.1712705688324178400000000
    ) +
    fcase(
      smoking == "non-smoker", 0,
      smoking == "ex-smoker", 0.1338683378654626200000000,
      smoking == "light smoker", 0.5620085801243853700000000,
      smoking == "moderate smoker", 0.6674959337750254700000000, 
      smoking == "heavy smoker", 0.8494817764483084700000000
    ) +

    age_1 * -8.1388109247726188000000000 +
    age_2 * 0.7973337668969909800000000 +
    bmi_1 * 0.2923609227546005200000000 +
    bmi_2 * -4.1513300213837665000000000 +
    chol_hdl_ratio * 0.1533803582080255400000000 +
    sbp * 0.0131314884071034240000000 +
    sbp_sd * 0.0078894541014586095000000 +
    townsend * 0.0772237905885901080000000 +

    atrial_fibrillation * 1.5923354969269663000000000 +
    atypical_antipsychotics * 0.2523764207011555700000000 +
    corticosteroids * 0.5952072530460185100000000 +
    migraine * 0.3012672608703450000000000 +
    rheumatoid_arthritis * 0.2136480343518194200000000 +
    CKD * 0.6519456949384583300000000 +
    mental_illness * 0.1255530805882017800000000 +
    SLE * 0.7588093865426769300000000 +
    bp_treatment * 0.5093159368342300400000000 +
    type_1_diabetes * 1.7267977510537347000000000 +
    type_2_diabetes * 1.0688773244615468000000000 +
    famhist_young_mi * 0.4544531902089621300000000 +

    age_1 * fcase(
      smoking == "non-smoker", 0,
      smoking == "ex-smoker", -4.7057161785851891000000000,
      smoking == "light smoker", -2.7430383403573337000000000,
      smoking == "moderate smoker", -0.8660808882939218200000000,
      smoking == "heavy smoker", 0.9024156236971064800000000
    ) +
    age_1 * atrial_fibrillation * 19.9380348895465610000000000 +
    age_1 * corticosteroids * -0.9840804523593628100000000 +
    age_1 * migraine * 1.7634979587872999000000000 +
    age_1 * CKD * -3.5874047731694114000000000 +
    age_1 * SLE * 19.6903037386382920000000000 +
    age_1 * bp_treatment * 11.8728097339218120000000000 +
    age_1 * type_1_diabetes * -1.2444332714320747000000000 +
    age_1 * type_2_diabetes * 6.8652342000009599000000000 +
    age_1 * bmi_1 * 23.8026234121417420000000000 +
    age_1 * bmi_2 * -71.1849476920870070000000000 +
    age_1 * famhist_young_mi * 0.9946780794043512700000000 +
    age_1 * sbp * 0.0341318423386154850000000 +
    age_1 * townsend * -1.0301180802035639000000000 +

    age_2 * fcase(
      smoking == "non-smoker", 0,
      smoking == "ex-smoker", -0.0755892446431930260000000,
      smoking == "light smoker", -0.1195119287486707400000000,
      smoking == "moderate smoker", -0.1036630639757192300000000,
      smoking == "heavy smoker", -0.1399185359171838900000000
    ) + 
    age_2 * atrial_fibrillation * -0.0761826510111625050000000 +
    age_2 * corticosteroids * -0.1200536494674247200000000 +
    age_2 * migraine * -0.0655869178986998590000000 +
    age_2 * CKD * -0.2268887308644250700000000 +
    age_2 * SLE * 0.0773479496790162730000000 +
    age_2 * bp_treatment * 0.0009685782358817443600000 +
    age_2 * type_1_diabetes * -0.2872406462448894900000000 +
    age_2 * type_2_diabetes * -0.0971122525906954890000000 +
    age_2 * bmi_1 * 0.5236995893366442900000000 +
    age_2 * bmi_2 * 0.0457441901223237590000000 +
    age_2 * famhist_young_mi * -0.0768850516984230380000000 +
    age_2 * sbp * -0.0015082501423272358000000 +
    age_2 * townsend * -0.0315934146749623290000000

  ]

  wip[sex == "Male", QRISK3 :=  
    fcase(
      ethnicity == "White or not stated", 0,
      ethnicity == "Indian", 0.2771924876030827900000000,
      ethnicity == "Pakistani", 0.4744636071493126800000000,
      ethnicity == "Bangladeshi", 0.5296172991968937100000000, 
      ethnicity == "Other Asian", 0.0351001591862990170000000, 
      ethnicity == "Black Caribbean", -0.3580789966932791900000000,
      ethnicity == "Black African", -0.4005648523216514000000000,
      ethnicity == "Chinese", -0.4152279288983017300000000,
      ethnicity == "Other ethnic group", -0.2632134813474996700000000
    ) + 
    fcase(
      smoking == "non-smoker", 0, 
      smoking == "ex-smoker", 0.1912822286338898300000000,
      smoking == "light smoker", 0.5524158819264555200000000, 
      smoking == "moderate smoker", 0.6383505302750607200000000,
      smoking == "heavy smoker", 0.7898381988185801900000000
    ) + 

    age_1 * -17.8397816660055750000000000 +
    age_2 * 0.0022964880605765492000000 +
    bmi_1 * 2.4562776660536358000000000 +
    bmi_2 * -8.3011122314711354000000000 +
    chol_hdl_ratio * 0.1734019685632711100000000 +
    sbp * 0.0129101265425533050000000 +
    sbp_sd * 0.0102519142912904560000000 +
    townsend * 0.0332682012772872950000000 +

    atrial_fibrillation * 0.8820923692805465700000000 +
    atypical_antipsychotics * 0.1304687985517351300000000 +
    corticosteroids * 0.4548539975044554300000000 +
    erectile_dysfunction * 0.2225185908670538300000000 +
    migraine * 0.2558417807415991300000000 +
    rheumatoid_arthritis * 0.2097065801395656700000000 +
    CKD * 0.7185326128827438400000000 +
    mental_illness * 0.1213303988204716400000000 +
    SLE * 0.4401572174457522000000000 +
    bp_treatment * 0.5165987108269547400000000 +
    type_1_diabetes * 1.2343425521675175000000000 +
    type_2_diabetes * 0.8594207143093222100000000 +
    famhist_young_mi * 0.5405546900939015600000000 +
  
    age_1 * fcase(
      smoking == "non-smoker", 0,
      smoking == "ex-smoker", -0.2101113393351634600000000,
      smoking == "light smoker", 0.7526867644750319100000000,
      smoking == "moderate smoker", 0.9931588755640579100000000,
      smoking == "heavy smoker", 2.1331163414389076000000000
    ) + 
    age_1 * atrial_fibrillation * 3.4896675530623207000000000 +
    age_1 * corticosteroids * 1.1708133653489108000000000 +
    age_1 * erectile_dysfunction * -1.5064009857454310000000000 +
    age_1 * migraine * 2.3491159871402441000000000 +
    age_1 * CKD * -0.5065671632722369400000000 +
    age_1 * bp_treatment * 6.5114581098532671000000000 +
    age_1 * type_1_diabetes * 5.3379864878006531000000000 +
    age_1 * type_2_diabetes * 3.6461817406221311000000000 +
    age_1 * bmi_1 * 31.0049529560338860000000000 +
    age_1 * bmi_2 * -111.2915718439164300000000000 +
    age_1 * famhist_young_mi * 2.7808628508531887000000000 +
    age_1 * sbp * 0.0188585244698658530000000 +
    age_1 * townsend * -0.1007554870063731000000000 +

    age_2 * fcase(
      smoking == "non-smoker", 0,
      smoking == "ex-smoker", -0.0004985487027532612100000,
      smoking == "light smoker", -0.0007987563331738541400000,
      smoking == "moderate smoker", -0.0008370618426625129600000,
      smoking == "heavy smoker", -0.0007840031915563728900000
    ) + 
    age_2 * atrial_fibrillation * -0.0003499560834063604900000 +
    age_2 * corticosteroids * -0.0002496045095297166000000 +
    age_2 * erectile_dysfunction * -0.0011058218441227373000000 +
    age_2 * migraine * 0.0001989644604147863100000 +
    age_2 * CKD * -0.0018325930166498813000000 +
    age_2 * bp_treatment * 0.0006383805310416501300000 +
    age_2 * type_1_diabetes * 0.0006409780808752897000000 +
    age_2 * type_2_diabetes * -0.0002469569558886831500000 +
    age_2 * bmi_1 * 0.0050380102356322029000000 +
    age_2 * bmi_2 * -0.0130744830025243190000000 +
    age_2 * famhist_young_mi * -0.0002479180990739603700000 +
    age_2 * sbp * -0.0000127187419158845700000 +
    age_2 * townsend * -0.0000932996423232728880000
  ]

  if (type == "absolute risk") {
    wip[, QRISK3 := 100 * 1 - ifelse(sex == "Female", 0.988876402378082, 0.977268040180206)^exp(QRISK3)]
  }
  return(wip$QRISK3)
}

QRISK3_absrisk <- function(sex, linear_predictor) {
  1 - ifelse(sex == "Female", 0.988876402378082, 0.977268040180206)^exp(linear_predictor)
}

