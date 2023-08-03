# Aggregate multiple P-values measured from separate tests in the same samples
# (e.g. multiple aptamers targetting the same protein) into a single P-value
# by taking the mean Z-score, and converting back into a two-sided p-value.
mean_pvalue <- function(pval, beta, avg.func=mean) {
  if(length(pval) == 1) {
    # When there is only test, simply return the p-value.
    return(pval)
  } else {
    # First, we need to convert two-sided p-values to one-sided p-values,
    # which we'll orient to the direction of effect of the smallest PRS to
    # aptamer p-value for this protein.
    one.sided <- pval/2
    to.invert <- sign(beta) != sign(beta[which.min(pval)])
    one.sided[to.invert] <- 1 - one.sided[to.invert]

    # Convert to Z-scores, average, then convert back to p-value
    zs <- qnorm(one.sided, lower.tail=FALSE)
    meta.one.sided <- pnorm(avg.func(zs), lower.tail=FALSE)

    # Convert to two-sided p-value
    meta.two.sided <- ifelse(meta.one.sided > 0.5, (1-meta.one.sided)*2, meta.one.sided*2)

    return(meta.two.sided)
  }
}

