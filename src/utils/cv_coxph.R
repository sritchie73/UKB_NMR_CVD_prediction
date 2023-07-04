require(survival)
require(data.table)
require(foreach)
source('src/utils/absrisk.R')

# Workhorse function for others in this script: fits a Cox proportional hazards model
# in cross-validation, returning a list of coxph objects, one per test fold (with the
# coxph being fitted on the remaining folds). The partitioning of the data into test
# folds is expected to have been done in advance.
cv.coxph <- function(formula, data, partition_col) {
  if (typeof(formula) == "character") {
    formula <- as.formula(formula)
  }

  # Determine if the model is the null model (intercept only, or just strata terms)
  formula_rhs <- attr(terms(formula), "term.labels")
  null_model <- length(formula_rhs) == 0 || all(grepl("^strata\\(", formula_rhs))
  if (null_model) {
    warning("Nothing to return for Null model (model with intercept only, or strata terms only)")
    return(NULL)
  }

  # Temporarily rename columns for programmatic access
  setnames(data, partition_col, "partition_col")
  on.exit(setnames(data, "partition_col", partition_col))

  test_folds <- unique(data[["partition_col"]])
  message(sprintf("Fitting Cox proportional hazards model in %s-fold cross-validation:", length(test_folds)))
  message("  Test fold:", appendLF=FALSE)
  cx_list <- foreach(tidx = 1:length(test_folds), .inorder=TRUE) %dopar% {
     message(sprintf(" %s...", tidx), appendLF=FALSE)
     coxph(formula, data=data[partition_col != test_folds[tidx]], x=TRUE)
  }
  names(cx_list) <- test_folds
  message("")

  return(cx_list)
}

# Fit Cox proportional hazards model in cross-validation and extract linear predictors
# in each test set to be used downstream for C-index and absolute risk predictions.
# Uses a predefinied assignment of samples to cross-validation folds in 'partition_col'.
cv.linear.predictor <- function(formula, data, partition_col, cx_list=NULL) {
  if (typeof(formula) == "character") {
    formula <- as.formula(formula)
  }

  # Fit the Cox proportional hazards models in each training set if not provided by the user
  if (is.null(cx_list)) {
    cx_list <- cv.coxph(formula, data, partition_col)
  }

  # Temporarily rename columns for programmatic access
  setnames(data, partition_col, "partition_col")
  on.exit(setnames(data, "partition_col", partition_col))

  # Build data.table to house results, making sure row order matches input data
  pred_dat <- data[,.SD,.SDcols="partition_col"]
  pred_dat[, linear.predictor := NA_real_]
  pred_dat[, rn := .I]

  # Predict linear predictor in each test set using Cox proportional hazards model fit
  # in the remaining data, and insert in the correct places in 'pred_dat'
  message("Extracting linear predictors...")
  for (test_fold in names(cx_list)) {
     train_cx <- cx_list[[test_fold]]
     test_data <- data[partition_col == test_fold] 
  
     # If the model formula has offset terms, we need to do some wrangling as predict.coxph
     # doesn't natively handle these. In summary, we need to remove them from the formula 
     # embedded in the train_cx object, then after predicting the linear predictor, add the
     # offset term to the result.
     offset_index <- attr(terms(formula), "offset")
     has_offset <- !is.null(offset_index)
     if (has_offset) {
       offset_terms <- rownames(attr(terms(formula), "factors"))[offset_index]
       offset_vars <- gsub("\\)$", "", gsub("offset\\(", "", offset_terms))

       surv_term <- as.character(formula)[2]
       depend_vars <- attr(terms(formula), "term.labels")
       new_formula <- as.formula(sprintf("%s ~ %s", surv_term, paste(depend_vars, collapse=" + ")))  

       new_terms <- coxph(new_formula, test_data)$terms # discard model fit, just want correctly formatted terms object
       train_cx$terms <- new_terms       
     }
    
     test_lp <- predict(train_cx, test_data, type = "lp", reference="zero")

     if (has_offset) {
       test_offset <- rowSums(test_data[as.integer(names(test_lp)), .SD, .SDcols=offset_vars])
       test_lp <- test_offset + test_lp
     }
  
     this_pred_dat <- pred_dat[partition_col == test_fold]
     this_pred_dat[as.integer(names(test_lp)), linear.predictor := test_lp]
     pred_dat[this_pred_dat, on = .(rn), linear.predictor := i.linear.predictor]
  }

  # The linear.predictors computed inside coxph are equivalent to predict(formula, data, reference="sample").
  # Above, we've used reference = "zero" to avoid doing per-test fold coefficient centering of the linear
  # predictors. Now, we apply that centering dataset-wide so that the returned linear.predictors are
  # functionally equivalent for downstream purposes.
  message("Centering linear predictors...")
  full_cx <- coxph(formula, data=data, x=TRUE)
  centering <- full_cx$linear.predictors - predict(full_cx, type="lp", reference="zero")
  pred_dat[as.integer(rownames(full_cx$x)), linear.predictor := linear.predictor + centering]

  # return as vector with names indicate rows (i.e. vector is shorter where missing data omitted)
  pred_dat[!is.na(linear.predictor), structure(linear.predictor, names=rn)]
}

# Fit C-index in cross-validation
cv.cindex <- function(formula, data, partition_col, cx_list=NULL, linear.predictor=NULL) {
  # process formula to extract relevant components needed for concordancefit
  if (typeof(formula) == "character") {
    formula <- as.formula(formula)
  }

  # Determine if the model is the null model (intercept only, or just strata terms)
  formula_rhs <- attr(terms(formula), "term.labels")
  null_model <- length(formula_rhs) == 0 || all(grepl("^strata\\(", formula_rhs))
  if (null_model) {
    warning("Nothing to return for Null model (model with intercept only, or strata terms only)")
    return(NULL)
  }

  # Extract survival term
  surv_term <- as.character(formula)[2]
  follow_col <- gsub("Surv\\(", "", gsub(",.*", "", surv_term))
  event_col <- gsub(".*, ?", "", gsub("\\)", "", surv_term))

  # Extract strata columns
  strata_terms <- formula_rhs[formula_rhs %like% "^strata\\("]
  strata_vars <- gsub("strata\\(", "", gsub("\\)", "", strata_terms))

  # Extract data columns later required for concordance fit
  pred_dat <- data[,.SD,.SDcols=c(event_col, follow_col, strata_vars)]

  # Build integer column that splits data along all strata variables
  if (length(strata_vars) > 0) {
    strata_num <- pred_dat[,.GRP, by=strata_vars]
    pred_dat[strata_num, on = c(strata_vars), strata_num := i.GRP]
  }

  # Get linear predictor, computing in cross-validation if not provided by the user
  if (is.null(linear.predictor) && is.null(cx_list)) {
    linear.predictor <- cv.linear.predictor(formula, data, partition_col)
  } else if (is.null(linear.predictor)) {
    linear.predictor <- cv.linear.predictor(formula, data, partition_col, cx_list=cx_list)
  }
  pred_dat[as.integer(names(linear.predictor)), c("linear.predictor") := linear.predictor]

  # Drop samples with incomplete data
  pred_dat <- pred_dat[!is.na(linear.predictor)]

  # Compute C-index
  message("Computing C-index...")
  if (length(strata_vars) > 0) {
    cf <- survival::concordancefit(
      y = Surv(pred_dat[[follow_col]], pred_dat[[event_col]]),
      x = pred_dat[["linear.predictor"]],
      strata = pred_dat[["strata_num"]],
      reverse = TRUE,
      timefix = FALSE
    )
  } else {
    cf <- survival::concordancefit(
      y = Surv(pred_dat[[follow_col]], pred_dat[[event_col]]),
      x = pred_dat[["linear.predictor"]],
      reverse = TRUE,
      timefix = FALSE
    )
  }

  # Extract point estimate, SE, and 95% CI
  cindex = cf$concordance
  cindex.se = sqrt(cf$var)
  cindex.l95 = cindex - qnorm(0.975)*cindex.se
  cindex.u95 = cindex + qnorm(0.975)*cindex.se

  # collate results
  data.table(C.index=cindex, SE=cindex.se, L95=cindex.l95, U95=cindex.u95,
             Samples=pred_dat[,.N], Cases=sum(pred_dat[[event_col]]),
             Missing=data[,.N]-pred_dat[,.N],
             Missing.Cases=sum(data[[event_col]]) - sum(pred_dat[[event_col]]))
}

# Predict absolute risk in cross-validation.
cv.absrisk <- function(data, partition_col, years=10, formula=NULL, cx_list=NULL, linear.predictor=NULL) {
  # Fit Cox proportional hazards models cross-validation if not provided by the user
  if (is.null(cx_list)) {
    if (is.null(formula)) {
      stop("'formula' must be provided if 'cx_list' is not")
    }
    cx_list <- cv.coxph(formula, data, partition_col)
  }

  # Extract linear predictors if not provided by the user
  if (is.null(linear.predictor)) {
    if (is.null(formula)) {
      stop("'formula' must be provided if 'linear.predictor' is not")
    }
    linear.predictor <- cv.linear.predictor(formula, data, partition_col, cx_list=cx_list)
  }  

  # Temporarily rename columns for programmatic access
  setnames(data, partition_col, "partition_col")
  on.exit(setnames(data, "partition_col", partition_col))

  # Build data.table to house results, making sure row order matches input data
  pred_dat <- data[,.SD,.SDcols="partition_col"]
  pred_dat[as.integer(names(linear.predictor)), c("linear.predictor") := linear.predictor]
  pred_dat[, absolute.risk := NA_real_]
  pred_dat[, rn := .I]

  # Predict absolute risk in each test fold using the model fit in the remaining folds.
  # Code in the for loop adapted from Aki Havulinna's Coxar function in the absrisk.R file
  message(sprintf("Predicting absolute %s-year risk...", years))
  for (test_fold in names(cx_list)) {
     train_cx <- cx_list[[test_fold]]
     test_data <- data[partition_col == test_fold] 
     this_pred_dat <- pred_dat[partition_col == test_fold]

     # Extract baseline hazards fit in the training folds
     train_bh <- basehaz(train_cx, centered=TRUE) # basehaz() is part of survival package

     # Extract a vector coding of the response term in the test fold
     surv_term <- train_cx$formula[[2]]
     test_surv <- eval(parse(text=deparse(surv_term)), envir=test_data)

     # Extract a data.table containing the strata in the test data
     test_strata <- getStrata(train_cx, newdata=test_data, warn=FALSE) # getStrata is a function in absrisk.R
   
     # Extract linear predictors predicted in the test data
     test_lp <- this_pred_dat[, linear.predictor]

     # Drop NAs
     non_missing <- this_pred_dat[, which(!is.na(linear.predictor))]
     test_surv <- test_surv[non_missing]
     test_strata <- test_strata[non_missing]
     test_lp <- test_lp[non_missing]
 
     # predict risk using absrisk function in absrisk.R
     test_absrisk <- absrisk(surv=test_surv, bh=train_bh, lp=test_lp, stratum=test_strata, years=years)

     # Add back into overall result matching row-order
     names(test_absrisk) <- non_missing
     this_pred_dat[as.integer(names(test_absrisk)), absolute.risk := test_absrisk]
     pred_dat[this_pred_dat, on = .(rn), absolute.risk := i.absolute.risk]
  }

  # return as vector with names indicate rows (i.e. vector is shorter where missing data omitted)
  pred_dat[!is.na(absolute.risk), structure(absolute.risk, names=rn)]
}

 
