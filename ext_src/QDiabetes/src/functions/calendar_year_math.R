library(data.table)
library(lubridate)

# Function to compute years between two dates - note this preserves
# human notions of whole years, i.e:
#
#   2019-07-25 - 2009-07-25 = 10 years
#
# At the expense of (potentially) leading to small inaccuracies in
# rank ordering due to leap days, i.e. in pure terms of days,
#
#   2019-07-25 - 2009-07-25 = 9.9986 years
#
# As this makes it harder to accurately truncate follow-up (e.g.
# for testing model calibration)
years_between <- function(d1, d2) {
  year_diff <- as.period(interval(as.Date(d1), as.Date(d2))) / years(1)

  # When d1 or d2 is a feb 29th, we want to return a whole number for
  # number of years between d1 and d2 if the other date is (1) Feb 28th,
  # and (2) not also a leap year.
  leap_day_d1 <- which(month(d1) == 2 & day(d1) == 29)
  leap_day_d2 <- which(month(d2) == 2 & day(d1) == 29)
  feb_28_d1 <- which(month(d1) == 2 & day(d2) == 28 & !(leap_year(year(d1))))
  feb_28_d2 <- which(month(d2) == 2 & day(d2) == 28 & !(leap_year(year(d2))))

  to_correct <- union(intersect(leap_day_d1, feb_28_d2), intersect(feb_28_d1, leap_day_d2))
  year_diff[to_correct] <- year(d2[to_correct]) - year(d1[to_correct])
  return(year_diff)
}

# Likewise, add_years is consistent with the above, i.e.
#
# 2009-07-25 + 10 years = 2019-07-25
#
# Note only works with whole years.
#
add_years <- function(d1, follow) {
  new_date <- as.IDate(as.Date(d1) + years(follow))
  # If the date input is Feb 29, and the corresponding year after
  # adding 'follow' is not a leap year, roll back to Feb 28.
  leap_day <- which(is.na(new_date) & month(d1) == 2 & day(d1) == 29)
  new_date[leap_day] <- as.IDate(as.Date(d1[leap_day]) - days(1) + years(follow))
  return(new_date)
}
