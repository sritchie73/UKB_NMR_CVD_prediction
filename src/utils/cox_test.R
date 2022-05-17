library(survival)

# Function to extract all relevant details from cox models
cox.test = function(formula, event_col, data) {
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

  # Get C-index and its 95% CI
  cindex = summary(cx)$concordance[1]
  cindex.se = summary(cx)$concordance[2]
  cindex.l95 = cindex - qnorm(0.975)*cindex.se
  cindex.u95 = cindex + qnorm(0.975)*cindex.se

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
