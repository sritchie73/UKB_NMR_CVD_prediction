library(data.table)
library(ggplot2)
library(ggforce)
library(scales)

# Load information about case identification
caseid <- rbind(idcol="strategy",
  blanket = fread("analyses/public_health_modelling/UK_population_generalised/blanket_screening/public_health_statistics_by_treatment_reason.txt"),
  targeted = fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening/public_health_statistics_by_treatment_reason.txt"),
  targeted_above_PRS = fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening_above_PRS/public_health_statistics_by_treatment_reason.txt")
)

# Harmonize treatment reason
caseid[, treatment_reason := fcase(
  treatment_reason %like% "High risk", "High risk in blanket screening",
  treatment_reason %like% "reclassified", "High risk in targeted screening",
  default = "Medium risk but diabetic or LDL-C >= 5 mmol/L"
)]

# Enforce row ordering of models
caseid[, name := factor(name, levels=c("Conventional RF", "NMR", "CRP", "NMR + Assays", "Assays"))]
caseid[, PGS := factor(PGS, levels=c("FALSE", "TRUE"))]
caseid[, treatment_reason := factor(treatment_reason, levels=rev(c("High risk in blanket screening", "Medium risk but diabetic or LDL-C >= 5 mmol/L", "High risk in targeted screening")))]
caseid <- caseid[order(name)][order(PGS)]
caseid[, long_name := factor(long_name, levels=rev(unique(long_name)))]
caseid <- caseid[order(strategy)][order(guidelines)]

# Load ONS pop so we know total cases
ons_pop <- fread("analyses/public_health_modelling/UK_population_generalised/ONS_hypothetical_100k_pop.txt")

# Load information about number needed to screen/treat
phs <- rbind(idcol="strategy",
  blanket = fread("analyses/public_health_modelling/UK_population_generalised/blanket_screening/public_health_statistics.txt"),
  targeted = fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening/public_health_statistics.txt"),
  targeted_above_PRS = fread("analyses/public_health_modelling/UK_population_generalised/targeted_screening_above_PRS/public_health_statistics.txt")
)

# Enforce row ordering of models
phs[, name := factor(name, levels=c("Conventional RF", "NMR", "CRP", "NMR + Assays", "Assays"))]
phs[, PGS := factor(PGS, levels=c("FALSE", "TRUE"))]
phs <- phs[order(name)][order(PGS)]
phs[, long_name := factor(long_name, levels=rev(unique(long_name)))]
phs <- phs[order(strategy)][order(guidelines)]


# Generate separate plots depending on risk thresholds
for (ii in c("NICE.2014", "ACC.AHA.2019")) {
  # Plot number of cases identified
  ggdt <- caseid[guidelines == ii & lambda != "lambda.1se"]
  xlim <- c( 
    ggdt[treatment_reason == "High risk in blanket screening", min(cases_treated)],
    ggdt[, .(total_cases=sum(cases_treated)), by=.(strategy, guidelines, long_name)][, max(total_cases)]
  )
  xlim <- c(xlim[1] - diff(xlim)*0.1, xlim[2] + diff(xlim)*0.1)
  xlim <- c(floor(xlim[1]/100)*100 - 30, ceiling(xlim[2]/100)*100 + 10)
  g <- ggplot(ggdt) +
    aes(x=cases_treated, y=long_name, fill=treatment_reason) +
    geom_bar(position="stack", stat="identity", color="black") +
    scale_fill_manual(name="", values=c(
      "High risk in blanket screening"="#bd0026",
      "Medium risk but diabetic or LDL-C >= 5 mmol/L"="#fc4e2a",
      "High risk in targeted screening"="#feb24c"
    )) + 
    facet_col(facets = vars(strategy), space="free", scales="free_y") +
    scale_x_continuous(name="CVD cases caught", limits=xlim, oob=oob_keep, expand=c(0,0),
      sec.axis=sec_axis(~ . / ons_pop$cases * 100, name="% CVD cases caught")) + 
    ylab("") + 
    theme_bw() +
    theme(axis.text=element_text(size=7), axis.title=element_text(size=8), 
          panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
          legend.position="bottom", legend.text=element_text(size=8))
  ggsave(g, width=10, height=6, file=sprintf("analyses/public_health_modelling/UK_population_generalised/%s_cases_identified.pdf", ii))

  # Plot NNS
  ggdt <- phs[guidelines == ii & lambda != "lambda.1se"]
  xlim <- range(ggdt$NNS)
  xlim <- c(xlim[1] - diff(xlim)*0.1, xlim[2] + diff(xlim)*0.1)
  xlim <- c(floor(xlim[1]/10)*10 - 3, ceiling(xlim[2]/10)*10 + 1)
  g <- ggplot(ggdt) +
    aes(x=NNS, y=long_name) +
    geom_col() + 
    geom_vline(linetype=2, color="black", xintercept=ggdt[strategy == "blanket" & name == "Conventional RF" & PGS == "FALSE", NNS]) +
    geom_vline(linetype=2, color="red", xintercept=ggdt[strategy == "blanket" & name == "NMR" & PGS == "FALSE", NNS]) +
    facet_col(facets = vars(strategy), space="free", scales="free_y") +
    scale_x_continuous(name="Number needed to screen", limits=xlim, oob=squish, expand=c(0,0)) +
    ylab("") + 
    theme_bw() +
    theme(axis.text=element_text(size=7), axis.title=element_text(size=8), 
          panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank())
  ggsave(g, width=10, height=5, file=sprintf("analyses/public_health_modelling/UK_population_generalised/%s_NNS.pdf", ii))

  # Plot NNT
  ggdt <- phs[guidelines == ii & lambda != "lambda.1se"]
  xlim <- range(ggdt$NNT)
  xlim <- c(xlim[1] - diff(xlim)*0.05, xlim[2] + diff(xlim)*0.05)
  xlim <- c(floor(xlim[1]/5)*5 - 3, ceiling(xlim[2]/5)*5 + 1)
  g <- ggplot(ggdt) +
    aes(x=NNT, y=long_name) +
    geom_col() + 
    geom_vline(linetype=2, color="black", xintercept=ggdt[strategy == "blanket" & name == "Conventional RF" & PGS == "FALSE", NNT]) +
    geom_vline(linetype=2, color="red", xintercept=ggdt[strategy == "blanket" & name == "NMR" & PGS == "FALSE", NNT]) +
    facet_col(facets = vars(strategy), space="free", scales="free_y") +
    scale_x_continuous(name="Number needed to treat", limits=xlim, oob=squish, expand=c(0,0)) +
    ylab("") + 
    theme_bw() +
    theme(axis.text=element_text(size=7), axis.title=element_text(size=8), 
          panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank())
  ggsave(g, width=10, height=5, file=sprintf("analyses/public_health_modelling/UK_population_generalised/%s_NNT.pdf", ii))
}

 

