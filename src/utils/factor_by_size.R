library(data.table)

# Function for converting a set of numbers to a factor, ordering levels by group size
factor_by_size <- function(x) {
  group_sizes <- as.data.table(table(x))
  group_sizes <- group_sizes[order(-N)]
  factor(x, levels=group_sizes$x)
}
