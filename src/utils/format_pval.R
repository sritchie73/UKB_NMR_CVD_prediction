# Custom p-value formatting
format_pval <- function(p) {
  sapply(p, function(x) {
		if (is.na(x)) {
			return("")
		} else if (x >= 0.06) {
			sprintf("%.2f", round(x, digits=2))
		} else if (x >= 0.001) {
			sprintf("%.3f", round(x, digits=3))
		} else {
			format(x, digits=1, scientific=TRUE)
		}
  })
}
