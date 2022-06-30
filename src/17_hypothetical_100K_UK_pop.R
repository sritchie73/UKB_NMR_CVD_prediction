library(openxlsx)
library(data.table)
library(foreach)

# Create output directory
out_dir <- "analyses/public_health_modelling/UK_population_generalised"
system(sprintf("mkdir -p %s", out_dir))

# Load mid-2020 population estimates for the UK downloaded from ONS:
# https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/bulletins/annualmidyearpopulationestimates/mid2020
ons_pop <- read.xlsx("data/ONS/datadownload.xlsx", sheet="2020")
setDT(ons_pop)
ons_pop <- ons_pop[geogname == "UNITED KINGDOM"] # filter to UK level summary
ons_pop[, c("variable", "geogcode") := NULL]
ons_pop <- melt(ons_pop, id.vars="geogname", value.name="N")
ons_pop[, geogname := NULL]
ons_pop <- ons_pop[!(variable %like% "_al$")] # Drop total counts (across all ages) for each sex
ons_pop[, age := as.integer(gsub(".*_", "", variable))]
ons_pop[, sex := ifelse(variable %like% "^m_", "Male", "Female")]

# Get total numbers in each five year age group by sex for age groups 40-70
ons_pop <- ons_pop[age >= 40 & age < 70]
ons_pop[, age_group := age %/% 5 * 5]
ons_pop <- ons_pop[, .(N=sum(N)), by=.(sex, age_group)]

# Standardise population to 100,000 individuals
total <- ons_pop[,sum(N)]
ons_pop[, N := N/total * 100000] 

# Note, hypothetical number of samples, cases, controls, and those allocated
# to risk strata are all fractional throughout - if rounded to whole numbers
# numbers don't add up to totals due to imbalance of rounding up/rounding down 
# across groups. Final numbers should only be rounded when preparing flowchart
# or inline text numbers - note this rounding sometimes will result in 
# samples/cases/controls being off by one when comparing different steps of the
# flowchart - you will have to manually add or subtract cases/controls at the
# relevant steps to harmonise numbers across all steps of the flowchart (i.e. 
# doing this at the minimum number of steps to harmonize numbers)

# Age- and sex-specific incidence rates (per 1000 person-years) of CVD
# among CPRD participants (n=2.1 million) from the Appendix table of
# Sun L, *et al*. Polygenic risk scores in cardiovascular risk
# prediction: A cohort study and modelling analyses. *PLOS Med* (2021).
CPRD = data.table(
	age_at_risk_start = seq(40, 75, by=5),
	age_at_risk_end = seq(44, 79, by=5),
	male_rate_per_1000py = c(1.578, 2.964, 5.141, 7.561, 10.848, 13.731, 19.271, 25.667),
	female_rate_per_1000py = c(0.726, 1.309, 2.165, 2.969, 4.560, 7.007, 10.840, 17.007)
)

# Calculate annual CVD incidence in each age-at-risk group from the
# 1000 person-year CVD rates:
CPRD[, male_annual_incidence := male_rate_per_1000py / 1000]
CPRD[, female_annual_incidence := female_rate_per_1000py / 1000]

# The expected 10-year CVD risk is calculated for each age group
# based on the mid-point of the next interval ahead, e.g. for the
# 40-44 year age-group the expected 10-year risk is calculated based
# on the annual incidence rates in the 45-49 year olds in CPRD:
CPRD[, age_group_start := age_at_risk_start - 5]
CPRD[, age_group_end := age_at_risk_end - 5]

# The expected risk is calculated assuming exponential survival (i.e.
# constant hazard) in the 10-years ahead.
CPRD[, male_expected_risk := 1 - exp(-male_annual_incidence * 10)]
CPRD[, female_expected_risk := 1 - exp(-female_annual_incidence * 10)]

# Use expected risks to estimate number of cases in each age group and sex
# In our standardised population
CPRD <- melt(CPRD, id.vars="age_group_start", measure.vars=c("male_expected_risk", "female_expected_risk"),
             value.name="expected_risk")
CPRD[, sex := ifelse(variable %like% "^male", "Male", "Female")]

# Estimate number of cases in ONS population
ons_pop[CPRD, on = .(sex, age_group=age_group_start), cases := N * expected_risk]
ons_pop[, controls := N - cases]

# Summarise all to population totals
ons_pop_summary <- ons_pop[, .(N=round(sum(N)), cases=round(sum(cases)), controls=round(sum(controls)))]

# Write out hypothetical population
fwrite(ons_pop, sep="\t", quote=FALSE, file=sprintf("%s/ONS_hypothetical_100k_pop_by_age_sex.txt", out_dir))
fwrite(ons_pop_summary, sep="\t", quote=FALSE, file=sprintf("%s/ONS_hypothetical_100k_pop.txt", out_dir))

