#open GWAS https://opengwas.io/datasets/

#Get authentication for ieugwasr
#https://mrcieu.github.io/ieugwasr/articles/guide.html#authentication
Sys.getenv("R_ENVIRON_USER")
Sys.getenv(c("HOME", "R_ENVIRON_USER"))
file.edit('~/.Renviron')
readRenviron("~/.Renviron")

#Load packages
library("ieugwasr")
ieugwasr::get_opengwas_jwt()
ieugwasr::user()
api_status()
library("openxlsx")
library("TwoSampleMR")
library("dplyr")

#Load exposure data
exposures_format<-read.csv("data/exposures_format_LST.csv")

phenotypes<- read.xlsx(xlsxFile = "data/replication/phenotype_opengwas.xlsx")
outcomes<-c("ukb-b-20208","ukb-b-13019")

#Get outcome data
#Set proxies=TRUE
out<-outcomes[i]
outcomes_format<-extract_outcome_data(
  snps = exposures_format$SNP, 
    outcomes=out,
  proxies= TRUE
)

write.csv(outcomes_format,paste0("data/replication/opengwas/",out,"_outcome.csv"), row.names = FALSE)

#Load number of cases and controls
sel<-filter(phenotypes, phenotype==out)
ncase<-sel$n_cases
ncase
ncontrol<-sel$n_controls
ncontrol
u<-ncase/(ncase+ncontrol)
uu<-u*(1-u)

#Transform beta se
outcomes_format$beta.outcome<-outcomes_format$beta.outcome/uu
outcomes_format$beta.outcome[1]
outcomes_format$se.outcome<-outcomes_format$se.outcome/uu
outcomes_format$se.outcome[1]

#Run MR IVW
dat <- harmonise_data(exposures_format, outcomes_format, action = 2)
res <- mr(dat, method_list = c("mr_wald_ratio", "mr_ivw"))

#Transform to or
res_or<-generate_odds_ratios(res)
res_or$or_ivwci<-sprintf("%.2f (%.2f, %.2f)", res_or$or, res_or$or_lci95, res_or$or_uci95)
write.csv(res_or, paste0("results/replication/opengwas/",out,"_LST.csv"), row.names = FALSE)
