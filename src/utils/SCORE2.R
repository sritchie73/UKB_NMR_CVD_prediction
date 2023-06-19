library(data.table)

# Implementation of SCORE2 prediciton algorithm as outlined in the Supplementary Methods of
# Heart Journal, E., and 2021 (2021). SCORE2 risk prediction algorithms: new models to estimate 
# 10-year risk of cardiovascular disease in Europe. Eur. Heart J. 42, 2439–2454.
#
# Supplementary Methods Table 2 gives the coefficients and normalizations to apply to each variable
# and interaction term
#
# Supplementary Methods Table 4 part 1 illustrates how to use these to compute a linear predictor
# based on the variables
#
# Supplementary Methods Table 4 part 2 shows how to compute absolute risk from these linear predictors,
# using baseline hazards from Supplementary Methods Table 2
#
score2 <- function(sex, age, smoking, sbp, tchol, hdl, type="absolute risk", calibration=TRUE, risk_region="low") {
  # Basic error checking
  type <- match.arg(type, c("linear predictor", "absolute risk"))
  stopifnot(length(calibration) == 1 || !is.na(calibration))
  risk_region <- match.arg(risk_region, c("low", "medium", "high", "very high"))

  stopifnot(age >= 40 & age <= 70) 
  stopifnot(is.logical(smoking))
  stopifnot(all(sex %in% c("Male", "Female")))
  stopifnot(all(!is.na(tchol)))
  stopifnot(all(!is.na(hdl)))
  stopifnot(all(!is.na(sbp)))
  stopifnot(all(!is.na(age)))
  stopifnot(all(!is.na(sex)))

  if (length(unique(length(sex), length(age), length(tchol), length(hdl), length(sbp), length(smoking))) != 1) {
    stop("Arguments 'sex', 'age', 'tchol', 'hdl', 'sbp', and 'smoking' must all be the same length")
  }

  # Apply standardisations as per Supplementary Methods Table 2
  cage <- (age - 60)/5
  smoking <- as.integer(smoking)
  smoking[is.na(smoking)] <- 0
  csbp <- (sbp - 120)/20
  ctchol <- (tchol - 6)/1
  chdl <- (hdl - 1.3)/0.5

  # Combine into work-in-progress table
  wip <- data.table(sex, cage, smoking, csbp, ctchol, chdl)

  # Calculate linear predictor as per Supplmentary Methods Table 4
  # using coefficients in Supplementary Methods Table 2
  wip[sex == "Male", linear_predictor :=
    cage * 0.3742 + 
    smoking * 0.6012 + 
    csbp * 0.2777 +
    ctchol * 0.1458 +
    chdl * -0.2698 +
    cage * smoking * -0.0755 +
    cage * csbp * -0.0255 +
    cage * ctchol * -0.0281 +
    cage * chdl * 0.0426]
  wip[sex == "Female", linear_predictor :=
    cage * 0.4648 + 
    smoking * 0.7744 + 
    csbp * 0.3131 +
    ctchol * 0.1002 +
    chdl * -0.2606 +
    cage * smoking * -0.1088 +
    cage * csbp * -0.0277 +
    cage * ctchol * -0.0226 +
    cage * chdl * 0.0613]

  # If all we want is the linear predictor, we return that now
  if (type == "linear predictor") {
    return(wip$linear_predictor)
  }

  # Compute un-calibrated absolute risk as per Supplementary Methods Table 4 
  # using baseline hazards in Supplementary Methods Table 2
  wip[, uncalibrated_absrisk := score2_absrisk(sex, linear_predictor)] # moved into separate function below

  # If all we want is uncalibrated absolute 10 year risk, return that now
  if (!calibration) {
    return(wip$uncalibrated_absrisk)
  }

  # Calibrate absolute risk based on chosen risk region as per Supplementary Methods Table 4
  # using scaling factors using scaling factors in Supplementary Methods Table 3
  wip[, absrisk := score2_recalibration(sex, uncalibrated_absrisk, risk_region)]
  return(wip$absrisk)
}

# convert a linear predictor from a Cox model into absolute 10-year risk using baseline hazards obtained
# from Score 2 paper in Supplementary Methods Table 2
score2_absrisk <- function(sex, linear_predictor) {
  stopifnot(all(sex %in% c("Male", "Female")))
  stopifnot(length(sex) == length(linear_predictor))

  wip <- data.table(sex, linear_predictor)
  wip[sex == "Male", uncalibrated_absrisk := 1-0.9605^exp(linear_predictor)]
  wip[sex == "Female", uncalibrated_absrisk := 1-0.9776^exp(linear_predictor)]

  return(wip$uncalibrated_absrisk) 
}

# recalibrate absolute risk to different European risk regions using SCORE2 paper scaling factors as
# detailed in Supplementary Methods Table 3 (scaling factors) and Supplementary Methods Table 4 
# (their application)
score2_recalibration <- function(sex, absrisk, risk_region="low") {
  stopifnot(all(sex %in% c("Male", "Female")))
  stopifnot(length(sex) == length(absrisk))
  risk_region <- match.arg(risk_region, c("low", "medium", "high", "very high"))
   
  wip <- data.table(sex, uncalibrated_absrisk=absrisk)

  if (risk_region == "low") {
    wip[sex == "Male", absrisk := 1-exp(-exp(-0.5699 + 0.7476 * log(-log(1-uncalibrated_absrisk))))]
    wip[sex == "Female", absrisk := 1-exp(-exp(-0.7380 + 0.7019 * log(-log(1-uncalibrated_absrisk))))]
  } else if (risk_region == "medium") {
    wip[sex == "Male", absrisk := 1-exp(-exp(-0.1565 + 0.8009 * log(-log(1-uncalibrated_absrisk))))]
    wip[sex == "Female", absrisk := 1-exp(-exp(-0.3143 + 0.7701 * log(-log(1-uncalibrated_absrisk))))]
  } else if (risk_region == "high") {
    wip[sex == "Male", absrisk := 1-exp(-exp(0.3207 + 0.9360 * log(-log(1-uncalibrated_absrisk))))]
    wip[sex == "Female", absrisk := 1-exp(-exp(0.5710 + 0.9369 * log(-log(1-uncalibrated_absrisk))))]
  } else if (risk_region == "very high") {
    wip[sex == "Male", absrisk := 1-exp(-exp(0.5836 + 0.8294 * log(-log(1-uncalibrated_absrisk))))]
    wip[sex == "Female", absrisk := 1-exp(-exp(0.9412 + 0.8329 * log(-log(1-uncalibrated_absrisk))))]
  }

  return(wip$absrisk)
}

# Validation as per Supplementary Methods Table 4 - note some digits off-by-one as table assumes rounding
# to 4 digits at each step, whereas in the function no rounding is performed
stopifnot(all.equal(score2("Male", 50, TRUE, 140, 6.3, 1.4, type="linear predictor"), 0.3221))
stopifnot(all.equal(score2("Female", 50, TRUE, 140, 6.3, 1.4, type="linear predictor"), 0.39788))
stopifnot(all.equal(round(score2("Male", 50, TRUE, 140, 6.3, 1.4, calibration=FALSE), digits=4), 0.0541))
stopifnot(all.equal(round(score2("Female", 50, TRUE, 140, 6.3, 1.4, calibration=FALSE), digits=4), 0.0332))
stopifnot(all.equal(round(score2("Male", 50, TRUE, 140, 6.3, 1.4, risk_region="low"), digits=4), 0.0631))
stopifnot(all.equal(round(score2("Female", 50, TRUE, 140, 6.3, 1.4, risk_region="low"), digits=4), 0.0433)) # off-by-one
stopifnot(all.equal(round(score2("Male", 50, TRUE, 140, 6.3, 1.4, risk_region="medium"), digits=4), 0.0811))
stopifnot(all.equal(round(score2("Female", 50, TRUE, 140, 6.3, 1.4, risk_region="medium"), digits=4), 0.0523))
stopifnot(all.equal(round(score2("Male", 50, TRUE, 140, 6.3, 1.4, risk_region="high"), digits=4), 0.0881))
stopifnot(all.equal(round(score2("Female", 50, TRUE, 140, 6.3, 1.4, risk_region="high"), digits=4), 0.0713))
stopifnot(all.equal(round(score2("Male", 50, TRUE, 140, 6.3, 1.4, risk_region="very high"), digits=4), 0.1506))
stopifnot(all.equal(round(score2("Female", 50, TRUE, 140, 6.3, 1.4, risk_region="very high"), digits=4), 0.1413)) # off-by-one

