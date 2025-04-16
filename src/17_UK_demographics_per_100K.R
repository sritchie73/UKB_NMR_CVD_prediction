library(openxlsx)
library(data.table)
library(foreach)

# Create output directory
out_dir <- "analyses/public_health_modelling"
system(sprintf("mkdir -p %s", out_dir))

# Load mid-2023 population estimates for the UK downloaded from ONS:
# https://www.ons.gov.uk/file?uri=/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland/mid2023/mye23tablesuk.xlsx
ons_pop <- rbind(idcol="sex",
  "Male"=as.data.table(read.xlsx("data/ONS/mye23tablesuk.xlsx", sheet="MYE2 - Males", startRow=8)),
  "Female"=as.data.table(read.xlsx("data/ONS/mye23tablesuk.xlsx", sheet="MYE2 - Females", startRow=8))
)
ons_pop <- ons_pop[Name == "UNITED KINGDOM"] # filter to UK level summary
 
# Get total numbers in each five year age group by sex for age groups 40-70
ons_pop <- melt(ons_pop, id.vars="sex", measure.vars=as.character(40:69), variable.name="age", value.name="N")
ons_pop[, age := as.integer(as.character(age))] # convert to character first so we dont get integer representation of factor levels
ons_pop[, age_group := sprintf("%s-%s", age %/% 5 * 5, age %/% 5 * 5 + 4)]
ons_pop <- ons_pop[, .(N=sum(N)), by=.(sex, age_group)]

# Standardise population to 100,000 individuals
ons_pop[, N := N/sum(N) * 100000]

# Age- and sex-specific incidence rates (per 1000 person-years) of CVD 
# among CPRD participants (n=3,117,544) from the Appendix table of 
# Sun L, *et al*. Use of polygenic risk scores and other molecular 
# markers to enhance cardiovascular risk prediction: prospective 
# cohort study and modelling analysis. *bioRxiv* (2019). 
# doi: 10.1101/744565.
CPRD <- data.table(
  age_at_risk_start = seq(40, 75, by=5),
  age_at_risk_end = seq(44, 79, by=5),
  male_rate_per_1000py = c(1.578, 2.964, 5.141, 7.561, 10.848, 13.731, 19.271, 25.667),
  female_rate_per_1000py = c(0.726, 1.309, 2.165, 2.969, 4.560, 7.007, 10.840, 17.007)
)

# The expected 10-year CVD risk is calculated for each age group
# based on the mid-point of the next interval ahead, e.g. for the 
# 40-44 year age-group the expected 10-year risk is calculated based
# on the annual incidence rates in the 45-49 year olds in CPRD:
CPRD[, age_group := sprintf("%s-%s", age_at_risk_start - 5, age_at_risk_end - 5)]

# Melt so that sex is a column
CPRD <- melt(CPRD, id.vars=c("age_group"), measure.vars=c("male_rate_per_1000py", "female_rate_per_1000py"), value.name="rate_per_1000py", variable.name="sex")
CPRD[, sex := gsub("_rate_per_1000py", "", sex)]
CPRD[, sex := ifelse(sex == "male", "Male", "Female")]

# Calculate annual CVD incidence in each age-at-risk group from the 
# 1000 person-year CVD rates:
CPRD[, annual_incidence := rate_per_1000py / 1000]

# The expected risk is calculated assuming exponential survival (i.e.
# constant hazard) in the 10-years ahead.
CPRD[, expected_risk := 1 -exp(-annual_incidence * 10)]

# Estimate number of cases in ONS population
ons_pop[CPRD, on = .(sex, age_group), cases := N * expected_risk]
ons_pop[, controls := N - cases]

# Summarise all to population totals
ons_pop_summary <- ons_pop[, .(N=sum(N), cases=sum(cases), controls=sum(controls)), by=sex]

# Write out hypothetical population
fwrite(ons_pop, sep="\t", quote=FALSE, file=sprintf("%s/ONS_hypothetical_100k_pop_by_age_sex.txt", out_dir))
fwrite(ons_pop_summary, sep="\t", quote=FALSE, file=sprintf("%s/ONS_hypothetical_100k_pop.txt", out_dir))

