library(data.table)

# Load results
blanket <- fread("analyses/public_health_modelling/blanket_screening/population_screening.txt")
targeted <- fread("analyses/public_health_modelling/targeted_screening/population_screening.txt")

# Combine
dat <- rbind(idcol="strategy", "blanket"=blanket, "targeted"=targeted)

# Format columns of interest
shownum <- function(num) { format(num, big.mark=",", trim=TRUE) }
dat <- dat[,.(
  strategy, model, sex, 
  high_risk = sprintf("%s (%s–%s)", shownum(high_risk), shownum(high_risk_L95), shownum(high_risk_U95)),
  cases = sprintf("%s (%s–%s)", shownum(high_risk_cases), shownum(high_risk_cases_L95), shownum(high_risk_cases_U95)),
  non_cases = sprintf("%s (%s–%s)", shownum(high_risk_non_cases), shownum(high_risk_non_cases_L95), shownum(high_risk_non_cases_U95)),
  prevented = sprintf("%s (%s–%s)", shownum(events_prevented), shownum(events_prevented_L95), shownum(events_prevented_U95)),
  NNS = sprintf("%s (%s–%s)", shownum(NNS), shownum(NNS_L95), shownum(NNS_U95)),
  NNT = sprintf("%s (%s–%s)", shownum(NNT), shownum(NNT_L95), shownum(NNT_U95))
)]

# Impose ordering
dat[, model := factor(model, levels=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
dat[, sex := factor(sex, levels=c("Males", "Females"))]
dat[, strategy := factor(strategy, levels=c("blanket", "targeted"))]
dat <- dat[order(sex)][order(model)][order(strategy)]

dat <- dcast(dat, strategy + model ~ sex, value.var=c("high_risk", "cases", "non_cases", "prevented", "NNS", "NNT"))
dat <- dat[, .(strategy, model, high_risk_Males, cases_Males, non_cases_Males, prevented_Males, NNS_Males, NNT_Males,
  high_risk_Females, cases_Females, non_cases_Females, prevented_Females, NNS_Females, NNT_Females)]

# Write out
fwrite(dat, sep="\t", quote=FALSE, file="analyses/public_health_modelling/screening_summary.txt")

