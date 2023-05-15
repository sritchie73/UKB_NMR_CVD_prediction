library(survival)

# Function to extract all relevant details from cox models
cox.test = function(formula, event_col, data, partition_col) {
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

  cx = coxph(formula, data=data, x=TRUE) # fit the cox model
  ci = confint(cx) # get 95% confidence intervals
  cz = cox.zph(cx, terms=FALSE) # test proportional hazards assumption

  # Get C-index and its 95% CI (in cross-validation if partition_col is given)
  if (!missing(partition_col)) {
    cf <- cv.cindex(formula, partition_col, data)
    cindex <- cf$C.index
    cindex.se <- cf$SE
    cindex.l95 <- cf$L95
    cindex.u95 <- cf$U95
  } else {
    cindex = summary(cx)$concordance[1]
    cindex.se = summary(cx)$concordance[2]
    cindex.l95 = cindex - qnorm(0.975)*cindex.se
    cindex.u95 = cindex + qnorm(0.975)*cindex.se
  }

  # collate information about all model coefficients
  cxdt = as.data.table(coef(summary(cx)), keep.rownames="coefficient")
  cidt = as.data.table(ci, keep.rownames="coefficient")
  czdt = as.data.table(cz$table, keep.rownames="coefficient")
  dt = cxdt[cidt, on = .(coefficient)][czdt, on = .(coefficient), nomatch=0]


  dt = dt[, .(coefficient, logHR = coef, SE=`se(coef)`, HR = `exp(coef)`,
              L95=exp(`2.5 %`), U95=exp(`97.5 %`), Pvalue=`Pr(>|z|)`,
              Proportionality.chisq=chisq, Proportionality.df=df, Proportionality.Pvalue=p)]

  # collate whole model stats
  wm = data.table(C.index=cindex, SE=cindex.se, L95=cindex.l95, U95=cindex.u95,
                  Proportionality.chisq=czdt[coefficient == "GLOBAL", chisq],
                  Proportionality.df=czdt[coefficient == "GLOBAL", df],
                  Proportionality.Pvalue=czdt[coefficient == "GLOBAL", p],
                  Samples=cx$n, Cases=cx$nevent, Missing=data[,.N]-cx$n, 
                  Missing.Cases=sum(data[[event_col]]) - cx$nevent)
  return(list(coefficients=dt, model_fit=wm))
}

# Fit Cox regression in cross-validation (predefined assignment of samples
# to cross-validation folds in 'partition_col') to estimate C-index
cv.cindex = function(formula, partition_col, data) {
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

  # Temporarily rename columns for programmatic access
  setnames(data, partition_col, "partition_col") 
  on.exit(setnames(data, "partition_col", partition_col))

  # Extract data columns later required for concordance fit
  pred_dat <- data[,.SD,.SDcols=c(event_col, follow_col, strata_vars, "partition_col")]
  pred_dat[, linear.predictor := NA_real_]
  pred_dat[, rn := .I]
  
  # Build integer column that splits data along all strata variables
  strata_num <- pred_dat[,.GRP, by=strata_vars]
  pred_dat[strata_num, on = c(strata_vars), strata_num := i.GRP]

  # For each test fold, fit cox model in the remaining data and obtain the linear predictors
  # in the test fold
  test_folds <- unique(data[["partition_col"]])
  message(sprintf("Calculating C-index in %s-fold cross-validation...", length(test_folds)))

  for (tidx in 1:length(test_folds)) {
     message(sprintf(" Test fold %s...", tidx))
     this_cx <- coxph(formula, data=data[partition_col != test_folds[tidx]], x=TRUE)
   
     this_lp <- predict(this_cx, newdata=data[partition_col == test_folds[tidx]], type = "lp")

     this_pred_dat <- pred_dat[partition_col == test_folds[tidx]]
     this_pred_dat[as.integer(names(this_lp)), linear.predictor := this_lp]
     pred_dat[this_pred_dat, on = .(rn), linear.predictor := i.linear.predictor]
  }

  # Drop samples with incomplete data
  pred_dat <- pred_dat[!is.na(linear.predictor)]

  # Compute C-index
  cf <- survival::concordancefit(
    y = Surv(pred_dat[[follow_col]], pred_dat[[event_col]]),
    x = pred_dat[["linear.predictor"]],
    strata = pred_dat[["strata_num"]],
    reverse = TRUE,
    timefix = FALSE
  )

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
