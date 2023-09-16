library(openxlsx)
library(data.table)
library(foreach)

# Create output directory
out_dir <- "analyses/public_health_modelling"
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
ons_pop[, age_group := sprintf("%s-%s", age %/% 5 * 5, age %/% 5 * 5 + 4)]
ons_pop <- ons_pop[, .(N=sum(N)), by=.(sex, age_group)]

# Standardise population to 100,000 individuals
# Note number will be slightly less to avoid fractional people
total <- ons_pop[,sum(N)]
ons_pop[, N := floor(N/total * 100000)]

# Load in CRPD incidences previously computed
CPRD <- fread("analyses/risk_recalibration/CPRD_incidence_and_risk.txt")

# Estimate number of cases in ONS population
ons_pop[CPRD, on = .(sex, age_group), cases := floor(N * expected_risk)]
ons_pop[, controls := N - cases]

# Summarise all to population totals
ons_pop_summary <- ons_pop[, .(N=sum(N), cases=sum(cases), controls=sum(controls)), by=sex]

# Write out hypothetical population
fwrite(ons_pop, sep="\t", quote=FALSE, file=sprintf("%s/ONS_hypothetical_100k_pop_by_age_sex.txt", out_dir))
fwrite(ons_pop_summary, sep="\t", quote=FALSE, file=sprintf("%s/ONS_hypothetical_100k_pop.txt", out_dir))

