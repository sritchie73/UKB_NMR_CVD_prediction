library(survival)

# Formula must mirror coxph formula
score_cindex <- function(formula, data) {
  if (typeof(formula) == "character") {
    formula <- as.formula(formula)
  }

  # Extract survival term
  surv_term <- as.character(formula)[2]
  follow_col <- gsub("Surv\\(", "", gsub(",.*", "", surv_term))
  event_col <- gsub(".*, ?", "", gsub("\\)", "", surv_term))

  # Extract strata columns
  formula_rhs <- attr(terms(formula), "term.labels")
  strata_terms <- formula_rhs[formula_rhs %like% "^strata\\("]
  strata_vars <- gsub("strata\\(", "", gsub("\\)", "", strata_terms))

  # Extract dependent variables
  depend_vars <- formula_rhs[!(formula_rhs %like% "^strata\\(")]

  # Combine scores as linear predictor
  lp <- rowSums(data[,.SD,.SDcols=depend_vars])

  # Build integer column that splits data along all strata variables
  if (length(strata_vars) > 0) {
    strata_num <- data[,.GRP, by=strata_vars]
    strata_num <- data[strata_num, on = c(strata_vars), i.GRP]
  }

  # Compute C-index
  if (length(strata_vars) > 0) {
    cf <- survival::concordancefit(y = Surv(data[[follow_col]], data[[event_col]]), x = lp, strata = strata_num, reverse = TRUE, timefix = FALSE)
  } else {
    cf <- survival::concordancefit(y = Surv(data[[follow_col]], data[[event_col]]), x = lp, reverse = TRUE, timefix = FALSE)
  }

  # Extract point estimate, SE, and 95% CI
  cindex = cf$concordance
  cindex.se = sqrt(cf$var)
  cindex.l95 = cindex - qnorm(0.975)*cindex.se
  cindex.u95 = cindex + qnorm(0.975)*cindex.se

  # collate results
  data.table(C.index=cindex, SE=cindex.se, L95=cindex.l95, U95=cindex.u95)
}
