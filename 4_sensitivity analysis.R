# Load packages
library("dplyr")
library("TwoSampleMR")
library("openxlsx")

#Run sensitivity analyses for LST
#Load harmonized data sets 
res<-read.csv("results/main_LST_pval.csv")
dat<-read.csv("data/expout_data_LST.csv")

#Select and save results passing Bonferroni
res_bon <- filter(res, pval_bonferroni==1)
outcomes<-res_bon$outcome
dat_bon<-dat[dat$outcome %in% outcomes, ]
  
#Add p-q for ivw hetero 
het_b <- mr_heterogeneity(dat_bon)
head(het_b)
het_b<-subset(het_b, method=="Inverse variance weighted")
p_het_b<-select(het_b,exposure, outcome, Q, Q_pval)
head(p_het_b)
write.csv(p_het_b,"results/cochran_q_LST.csv", row.names = FALSE)

#Run MR egger & wm
res <- mr(dat_bon, method_list = c("mr_egger_regression", "mr_weighted_median"))
#Transform to or
res_or<-generate_odds_ratios(res)
res_or$or_ci<-sprintf("%.2f (%.2f, %.2f)", res_or$or, res_or$or_lci95, res_or$or_uci95)
write.csv(res_or, "results/res_sens_LST_or.csv", row.names = FALSE)

#Get p-egger intercept to egger res
plt_b <- mr_pleiotropy_test(dat_bon)
head(plt_b)
write.csv(plt_b,"results/egger_inter_LST.csv", row.names = FALSE)

#Run Steiger filtering
#Load number of cases and controls
outcome_info<-read.csv("data/endpoint_cate_total.csv")
names(outcome_info)
outcome_info<-select(outcome_info, phenocode, num_cases, num_cases_prev,num_controls)
colnames(outcome_info)[1]<-"outcome"
dat_bon<-left_join(dat_bon, outcome_info, by= "outcome")
summary(dat_bon$num_cases)

#Get r from bse
dat_bon$r_exp <- get_r_from_bsen(
  b =  dat_bon$beta.exposure, 
  se = dat_bon$se.exposure, 
  n =  dat_bon$samplesize.exposure)

#Get r from lor
dat_bon$r_out <- get_r_from_lor(
  lor =  dat_bon$beta.outcome,
  af =  dat_bon$eaf.outcome,
  ncase =  dat_bon$num_controls,
  ncontrol =  dat_bon$num_cases,
  prevalence =   dat_bon$num_cases_prev/(dat_bon$num_cases_prev+dat_bon$num_controls),
  model = "logit",
  correction = FALSE
)

dat_bon$correct_causal_direction<-NA
for (i in 1:9309)  { 
  steiger_res <- mr_steiger(
    p_exp = dat_bon$pval.exposure[i], 
    p_out = dat_bon$pval.outcome[i], 
    n_exp = dat_bon$samplesize.exposure[i], 
    n_out = dat_bon$samplesize.outcome[i], 
    r_exp = dat_bon$r_exp[i], 
    r_out = dat_bon$r_out[i])
  
  dat_bon$correct_causal_direction[i]<-steiger_res$correct_causal_direction}

write.xlsx(dat_bon, file = "results/expout_data_LST_steiger.xlsx", overwrite=F)

#Select SNPs passing Steiger filtering
dat_filtered <-dat_bon %>% filter(correct_causal_direction == TRUE)

#Run MR IVW
res <- mr(dat_filtered , method_list = c("mr_wald_ratio","mr_ivw"))

#Transform to or
res_or<-generate_odds_ratios(res)
res_or$or_ci<-sprintf("%.2f (%.2f, %.2f)", res_or$or, res_or$or_lci95, res_or$or_uci95)
write.csv(res_or, "results/steiger_LST_or.csv", row.names = F)

#Run sensitivity analyses for MPVA
#Load harmonized data sets 
res<-read.csv("results/main_MVPA_pval.csv")
dat<-read.csv("data/expout_data_MVPA.csv")

#Select and save results passing Bonferroni
res_bon <- filter(res, pval_bonferroni==1)
names(res_bon)
outcomes<-res_bon$outcome
head(dat)
dat_bon<-dat[dat$outcome %in% outcomes, ]

#Add p-q for ivw hetero 
het_b <- mr_heterogeneity(dat_bon)
head(het_b)
het_b<-subset(het_b, method=="Inverse variance weighted")
p_het_b<-select(het_b,exposure, outcome, Q, Q_pval)
head(p_het_b)
write.csv(p_het_b,"results/cochran_q_MVPA.csv", row.names = FALSE)

#Run MR egger & wm
res <- mr(dat_bon, method_list = c("mr_egger_regression", "mr_weighted_median"))
#Transform to or
res_or<-generate_odds_ratios(res)
res_or$or_ci<-sprintf("%.2f (%.2f, %.2f)", res_or$or, res_or$or_lci95, res_or$or_uci95)
write.csv(res_or, "results/res_sens_MVPA_or.csv", row.names = FALSE)

#Get p-egger intercept to egger res
plt_b <- mr_pleiotropy_test(dat_bon)
head(plt_b)
write.csv(plt_b,"results/egger_inter_MVPA.csv", row.names = FALSE)

#Run Steiger filtering
#Get r from lor
head(dat_bon)
dat_bon$ncase_exp <- 299659
dat_bon$ncontrol_exp <- 325299
dat_bon$prev_exp <- 0.48

dat_bon<-left_join(dat_bon, outcome_info, by= "outcome")

dat_bon$r_exp <- get_r_from_lor(
  lor =  dat_bon$beta.exposure,
  af =  dat_bon$eaf.exposure,
  ncase =  dat_bon$ncase_exp,
  ncontrol =  dat_bon$ncontrol_exp,
  prevalence =   dat_bon$prev_exp,
  model = "logit",
  correction = TRUE
)

#get r from lor
dat_bon$r_out <- get_r_from_lor(
  lor =  dat_bon$beta.outcome,
  af =  dat_bon$eaf.outcome,
  ncase =  dat_bon$num_controls,
  ncontrol =  dat_bon$num_cases,
  prevalence =   dat_bon$num_cases_prev/(dat_bon$num_cases_prev+dat_bon$num_controls),
  model = "logit",
  correction = FALSE
)

head(dat_bon)

dat_bon$correct_causal_direction<-NA
head(dat_bon)
for (i in 1:96)  { 
  steiger_res <- mr_steiger(
    p_exp = dat_bon$pval.exposure[i], 
    p_out = dat_bon$pval.outcome[i], 
    n_exp = dat_bon$samplesize.exposure[i], 
    n_out = dat_bon$samplesize.outcome[i], 
    r_exp = dat_bon$r_exp[i], 
    r_out = dat_bon$r_out[i])
  
  dat_bon$correct_causal_direction[i]<-steiger_res$correct_causal_direction}

write.xlsx(dat_bon, file = "results/expout_data_MVPA_steiger.xlsx", overwrite=F)

#Select SNPs passing Steiger filtering
dat_filtered <-dat_bon %>% filter(correct_causal_direction == TRUE)

#Run MR IVW
res <- mr(dat_filtered , method_list = c("mr_wald_ratio","mr_ivw"))

#Transform to or
res_or<-generate_odds_ratios(res)
res_or$or_ci<-sprintf("%.2f (%.2f, %.2f)", res_or$or, res_or$or_lci95, res_or$or_uci95)
write.csv(res_or, "results/steiger_MVPA_or.csv", row.names = FALSE)
