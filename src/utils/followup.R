library(lubridate)

# Get number of years between two sets of dates
years_between <- function(d1, d2) {
  time_length(as.Date(d2) - as.Date(d1), unit="years")
}

# Get a date after adding some number of years to it.
# Note non-integer numbers given to 'follow' sometimes
# leads to off-by-one errors in the day of the returned
# event where the numeric precision leads to < 1 days worth
# of seconds.
add_years_to_date <- function(d1, follow) {
  as.IDate(date_decimal(decimal_date(d1) + follow))
}
