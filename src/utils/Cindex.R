library(survival)

# Estimate the C-index for a model or score
#
# 'model' can be any of the following:
# 
#  1. A Cox proportional hazards model fit using coxph()
#  2. A model formula to be fit using coxph()
#  3. A linear predictor (e.g. a risk score)
#
# 'data' must be supplied if 'model' is a formula to fit with coxph()
# 'y' must be supplied if 'model' is linear predictors and must have been constructed with Surv()
# 'strata.term' can be optionally supplied when 'model' is a linear predictor to stratify by some factor (e.g. sex)
# 'strata.term' is otherwise ignored, i.e. you don't need to supply it if 'model' is a coxph object that is already stratified,
#    or when 'model' is a formula that includes strata terms. 
# 
# Returns a data.frame containing columns:
#  
#  'n': the number of bservations with complete data 
#  'nevent': the number of events with complete data
#  'C.index': the C-index for the model
#  'SE': the standard error for the C-index, computed using the usual infinitesiminal jackknife method
#  'L95': the lower bound of the 95% confidence interval for the C-index
#  'U95': the upper bound of the 95% confidence interval for the C-index
#
cindex <- function(model, y, strata.term, data) {
  if (length(model) == 1 && "character" %in% class(model)) {
    if (missing(data)) stop("'data' must be provided if 'model' is a formula")
    if (length(model) == 1 && "character" %in% class(model)) model <- as.formula(model)
    return(cindex(model, data=data))
  } else if ("formula" %in% class(model)) { 
    if (missing(data)) stop("'data' must be provided if 'model' is a formula")
    if ("formula" %in% class(model)) model <- coxph(model, data=data)
    return(cindex(model))
  } else if ("coxph" %in% class(model)) {
    # Extract information directly from cox model
    C.index <- model$concordance["concordance"]
    SE <- model$concordance["std"]
    n <- model$n
    nevent <- model$nevent
  } else if ("numeric" %in% class(model)) {
    if(missing(y)) stop("'y' must be provided if 'model' is a linear predictor")
    if(!("Surv" %in% class(y))) stop("'y' must be an survival object constructed with Surv()")

    # Filter to complete cases
    if (any(is.na(model))  || any(is.na(y)) || (!missing(strata.term) && any(is.na(strata.term)))) {
      non_missing_m <- which(!is.na(model))
      non_missing_y <- which(!is.na(y))
      if (!missing(strata.term)) non_missing_s <- which(!is.na(strata.term))

      non_missing <- intersect(non_missing_y, non_missing_m)
      if (!missing(strata.term)) non_missing <- intersect(non_missing, non_missing_s)

      model <- model[non_missing]
      y <- y[non_missing]
      if (!missing(strata.term)) strata.term <- strata.term[non_missing]
    }

    # Get information of complete cases
    n <- length(y[,1])
    nevent <- sum(y[,2])

    # Compute concordance based on score
    if (!missing(strata.term)) {
      conc_fit <- concordance(y ~ model + strata(strata.term), reverse=TRUE)
    } else {
      conc_fit <- concordance(y ~ model, reverse=TRUE)
    }
 
    # Extract C-index and SE
    C.index <- conc_fit$concordance
    SE <- sqrt(conc_fit$var)
  } else {
    stop("'model' must either be a model fit by 'coxph()', a formula to fit that models, or a linear predictor (e.g. risk score)")
  }
 
  # Compute 95% CI
  L95 <- C.index - qnorm(0.975) * SE
  U95 <- C.index + qnorm(0.975) * SE

  # collate and return the results
  res <- data.table("n"=n, "nevent"=nevent, "C.index"=C.index, "SE"=SE, "L95"=L95, "U95"=U95)
  return(res)
}

# Estimate the change in C-index between two models or scores
#
# 'model1' and 'model2' can be any combination of:
# 
#  1. A Cox proportional hazards models fit using coxph()
#  2. A model formula to be fit using coxph()
#  3. A linear predictor (e.g. a risk score)
#
# 'data' must be supplied if either 'model1' or 'model2' are model formula to fit with coxph()
# 'y' must be supplied if either 'model1' or 'model2' are linear predictors and must be constructed with Surv()
# 'strata.term' can be optionally supplied when either 'model1' or 'model2' are linear predictors to stratify by some factor (e.g. sex)
# 'strata.term' is otherwise ignored, i.e. you don't need to supply it if 'model1' or 'model2' are coxph objects that are already stratified,
#    or when 'model1' or 'model2' are formula that include strata terms. 
# 
# Returns a data.frame containing columns:
#  
#  'n': the number of shared observations with complete data for both models
#  'nevent': the number of events with complete data for both models
#  'm1.C.index': the C-index for model1, computed on the complete data shared by both models
#  'm1.SE': the standard error for the C-index for model1, computed using the usual infinitesiminal jackknife method
#  'm1.L95': the lower bound of the 95% confidence interval for the C-index for model1
#  'm1.U95': the upper bound of the 95% confidence interval for the C-index for model1
#  'm2.C.index': the C-index for model2, computed on the complete data shared by both models
#  'm2.SE': the standard error for the C-index for model2, computed using the usual infinitesiminal jackknife method
#  'm2.L95': the lower bound of the 95% confidence interval for the C-index for model2
#  'm2.U95': the upper bound of the 95% confidence interval for the C-index for model2
#  'deltaC': the change in C-index from model1 to model2
#  'deltaC.SE': the standard error of the change in C-index, accounting for the positive correlation between the two models
#  'deltaC.L95': the lower bound of the 95% confidence interval for the change in C-index
#  'deltaC.U95': the upper bound of the 95% confidence interval for the change in C-index
#  'deltaC.pval': the two-sided p-value that the change in C-index is significantly different to 0
#
# Special thanks to Joel T. Gibson for working out the covariance modelling steps required for fast computation of correlated standard errors
# 
delta_cindex <- function(model1, model2, y, strata.term, data) {
  if (length(model1) == 1 && "character" %in% class(model1) || length(model2) == 1 && "character" %in% class(model2)) {
    if (missing(data)) stop("'data' must be provided if either 'model1' or 'model2' are formulas")

    if (length(model1) == 1 && "character" %in% class(model1)) model1 <- as.formula(model1)
    if (length(model2) == 1 && "character" %in% class(model2)) model2 <- as.formula(model2)
     
    return(delta_cindex(model1, model2, data=data))
  } else if ("formula" %in% class(model1) || "formula" %in% class(model2)) { # If provided formula(s) fit models
    if (missing(data)) stop("'data' must be provided if either 'model1' or 'model2' are formulas")

    if ("formula" %in% class(model1)) model1 <- coxph(model1, data=data, x=TRUE) 
    if ("formula" %in% class(model2)) model2 <- coxph(model2, data=data, x=TRUE) 
    
    return(delta_cindex(model1, model2))
  } else if("coxph" %in% class(model1) && "coxph" %in% class(model2)) { # Compare two fitted cox models
    # Check models are fit on the same data, and if not (i.e. due to differential missingness) attempt to fix
    if (model1$n != model2$n || model1$nevent != model2$nevent || any(rownames(model2$y) != rownames(model1$y))) {
      if(!("x" %in% names(model1)) || !("x" %in% names(model2))) {
        stop("Unable to handle differential missing data: 'x=TRUE' must be set when running coxph() on the provided models")
      }

      # get intersection of complete cases
      non_missing_1 <- as.integer(rownames(model1$y))
      non_missing_2 <- as.integer(rownames(model2$y))
      non_missing <- intersect(non_missing_1, non_missing_2)

      # Define function for refitting a cox model on a given set of samples 
      refit_coxph <- function(cx, samples) {
        new <- as.data.frame(cx$x)
        if ("offset" %in% names(cx)) new <- cbind(new, cx$offset)
        new_xvar <- paste0("x", 1:ncol(new))
        colnames(new) <- new_xvar
        new <- cbind("time"=cx$y[,1], "status"=cx$y[,2], new)
        if ("strata" %in% names(cx)) new <- cbind(new, s=cx$strata)
        new <- new[as.character(samples),]

        mf <- paste("Surv(time, status) ~", paste(new_xvar, collapse=" + "))
        if ("strata" %in% names(cx)) mf <- paste(mf, "+ strata(s)")

        coxph(as.formula(mf), data=new)
      }
       
      model1 <- refit_coxph(model1, non_missing)
      model2 <- refit_coxph(model2, non_missing)
    }

    # Also sanity check the response vectors are identical (i.e. we're not models fit on different datasets), concordance
    # only throws a warning, we want to error
    if (any(model1$y[,1] != model2$y[,1]) || any(model1$y[,2] != model2$y[,2])) {
      stop("Models have different response vectors, cannot compare models fit on different datasets")
    }

    # Get information of shared observations and events
    shared_n <- model1$n
    shared_events <- model2$nevent

    # Compute C-index information for both models now on the shared complete data
    c1 <- cindex(model1)
    c2 <- cindex(model2)

    # Now compute the concordance between the two models from which to calculate the change in C-index
    c_compare <- concordance(model1, model2)
  } else if ("numeric" %in% class(model1) && "numeric" %in% class(model2)) { # Compare two scores (linear predictors)
    if(missing(y)) stop("'y' must be provided if both model1 and model2 are linear predictors")
    if(!("Surv" %in% class(y))) stop("'y' must be an survival object constructed with Surv()")

    # Filter to complete cases
    if (any(is.na(model1)) || any(is.na(model2)) || any(is.na(y)) || (!missing(strata.term) && any(is.na(strata.term)))) {
      non_missing_1 <- which(!is.na(model1))
      non_missing_2 <- which(!is.na(model2))
      non_missing_y <- which(!is.na(y))
      if (!missing(strata.term)) non_missing_s <- which(!is.na(strata.term))

      non_missing <- intersect(non_missing_y, intersect(non_missing_1, non_missing_2))
      if (!missing(strata.term)) non_missing <- intersect(non_missing, non_missing_s)
      
      model1 <- model1[non_missing]
      model2 <- model2[non_missing]
      y <- y[non_missing]
      if (!missing(strata.term)) strata.term <- strata.term[non_missing]
    }

    # Get information of shared observations and events
    shared_n <- length(y[,1])
    shared_events <- sum(y[,2])

    # Compute C-index information for both models now on the shared complete data
    if (!missing(strata.term)) {
			c1 <- cindex(model1, y, strata.term)
			c2 <- cindex(model2, y, strata.term)
    } else {
      c1 <- cindex(model1, y)
      c2 <- cindex(model2, y)
    }

    # Now compute the concordance between the two models from which to calculate the change in C-index
    if (!missing(strata.term)) {
			# Strata terms require fitting new cox models to make the comparison as concordance() does not natively work with 
			# multiple predictors and a strata term (see: https://github.com/therneau/survival/issues/249)
      cx1 <- coxph(y ~ strata(strata.term) + model1) 
      cx2 <- coxph(y ~ strata(strata.term) + model2)
      c_compare <- concordance(cx1, cx2)
    } else {
			c_compare <- concordance(y ~ model1 + model2, reverse=TRUE)
    }
  } else if ( # Compare a cox model to a linear predictor
    ("coxph" %in% class(model1) || "coxph" %in% class(model2)) &&
    ("numeric" %in% class(model1) || "numeric" %in% class(model2))
  ) {
    if(missing(y)) stop("'y' must be provided if either model1 or model2 are linear predictors")
    if(!("Surv" %in% class(y))) stop("'y' must be an survival object constructed with Surv()")

    # Conver the linear predictor to a fitted cox model, then re-call deltaC 
    if ("numeric" %in% class(model1)) {
      m1_lp <- model1 # prevent namespace collision when calling delta_cindex on result
      if (!missing(strata.term)) {
        model1 <- coxph(y ~ strata(strata.term) + m1_lp, x=TRUE)
      } else {
        model1 <- coxph(y ~ m1_lp, x=TRUE)
      }
    } else if ("numeric" %in% class(model2)) {
      m2_lp <- model2 # prevent namespace collision when calling delta_cindex on result
      if (!missing(strata.term)) {
        model2 <- coxph(y ~ strata(strata.term) + m2_lp, x=TRUE)
      } else {
        model2 <- coxph(y ~ m2_lp, x=TRUE)
      }
    }
    return(delta_cindex(model1, model2))
  } else {
    stop("'model1' and 'model2' must either be models fit by 'coxph()', formulae to fit those models, or linear predictors (e.g. risk scores) to compare")
  }

  # Compute change in C-index and its standard error accounting for positive correlation between the two models
  deltaC <- as.vector(c(-1, 1) %*% coef(c_compare))                   # Courtesy of Joel T. Gibson
  SE <- as.vector(sqrt(c(-1, 1) %*% vcov(c_compare) %*% c(-1, 1)))    # Courtest of Joel T. Gibson
  L95 <- deltaC - qnorm(0.975) * SE
  U95 <- deltaC + qnorm(0.975) * SE
  pval <- pmin(1, pnorm(abs(deltaC / SE), lower.tail = FALSE) * 2)

  # collated and return the results
  res <- data.table(
    "shared.n"=shared_n, "shared.nevent"=shared_events, 
    "m1.C.index"=c1$C.index, "m1.SE"=c1$SE, "m1.L95"=c1$L95, "m1.U95"=c1$U95,
    "m2.C.index"=c2$C.index, "m2.SE"=c2$SE, "m2.L95"=c2$L95, "m2.U95"=c2$U95,
    "deltaC"=deltaC, "deltaC.SE"=SE, "deltaC.L95"=L95, "deltaC.U95"=U95, "deltaC.pval"=pval
  )
  return(res)
}

