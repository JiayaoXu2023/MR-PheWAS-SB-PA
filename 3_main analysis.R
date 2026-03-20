# Install necessary packages
library("dplyr")
library("openxlsx")
library("TwoSampleMR")
data_endpoint<-read.xlsx("data/endpoints.xlsx")

#Load exposure and outcome data
exposures_format<-read.csv("data/exposures_format.csv")
outcomes_format<-read.csv("data/outcome/outcomes_proxy_num.csv")

#Harmonize data sets     
dat <- harmonise_data(exposures_format, outcomes_format, action = 2)
write.csv(dat, "data/expout_data.csv", row.names = FALSE)

#Run main analyses 
res <- mr(dat, method_list = c("mr_wald_ratio", "mr_ivw"))

#Transform to or
res_or<-generate_odds_ratios(res)
res_or$or_ci<-sprintf("%.2f (%.2f, %.2f)", res_or$or, res_or$or_lci95, res_or$or_uci95)

#Calculate p value
names(res)
res_or$pval_bonferroni<- ifelse(res_or$pval <= 0.05/144, 1, 0)

#Calculate FDR 5%
res_or$pval_fdr <- p.adjust(res_or$pval, method = "BH", n = length(res_or$pval))
