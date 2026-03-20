#Neale lab https://www.nealelab.is/uk-biobank
#Load packages
library("openxlsx")
library("data.table")
library("dplyr")
library("readr")
library("TwoSampleMR")

#Load information
outcome<-read.xlsx("data/replication/nealeoutcomes.xlsx")

#LST
outcome_LST<-outcome[outcome$Exposure=="LST",]
outlist<-outcome_LST$Outcome
risd_LST<-readRDS("data/replication/neale/rsidinf_correct_LST.rds")

#Prepare outcome data
for (i in 1:51){
  print(i)
  out<-outlist[i]
  print(out)
  phedat<-fread(paste0("data/replication/neale/",out,".gwas.imputed_v3.both_sexes.tsv.bgz"))  
  #separate for female-only outcomes
  #phedat<-fread(paste0("data/replication/neale/",out,".gwas.imputed_v3.female.tsv.bgz")) 
  outcomes<-left_join(risd_LST, phedat, by = c("variant"))
  outcomes$af<-ifelse (outcomes$alt==outcomes$minor_allele,outcomes$minor_AF, 
                       ifelse (outcomes$ref==outcomes$minor_allele, 1-outcomes$minor_AF, NA))
  outcomes_gwas <- outcomes[,c("rsid","alt","ref","af","beta","se","pval","n_complete_samples", "chr", "pos")]
  write.csv(outcomes_gwas, paste0("data/replication/neale/", out, "_outcomeall.csv"), row.names = FALSE)
}

#Run analysis
exposures_format<-read.csv("data/exposures_format_LST.csv")

phenotypes<-read_delim(file = "data/replication/neale/phenotypes.both_sexes.tsv.bgz",delim="\t")
f_phenotypes<-read_delim(file = "data/replication/neale/phenotypes.female.tsv.bgz",delim="\t")

res_all<-data.frame()
for (i in 1:51){
  print(i)
  out<-outlist[i]
  print(out)
  outcomes_gwas<-read.csv(paste0("data/replication/neale/", out, "_outcomeall.csv"))
  outcomes_format<-format_data(outcomes_gwas, 
                               type="outcome", 
                               header=TRUE,  
                               snp= exposures_format$SNP,
                               snp_col = "rsid",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "af",
                               effect_allele_col = "alt",
                               other_allele_col = "ref",
                               pval_col = "pval",
                               chr_col = "chr",
                               pos_col = "pos")
  
  sel<-filter(phenotypes, phenotype==out)
  #sel<-filter(f_phenotypes, phenotype==out)
  
  #Load number of cases and controls
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
  #All SNPs have been harmonised so no need for proxy SNPs
  res <- mr(dat, method_list = c("mr_wald_ratio", "mr_ivw"))
  
  #Transform to or
  res_or<-generate_odds_ratios(res)
  res_or$or_ivwci<-sprintf("%.2f (%.2f, %.2f)", res_or$or, res_or$or_lci95, res_or$or_uci95)
  res_or$outcome<-out
  res_or[1,16]<-ncase
  res_or[1,17]<-ncontrol
  
  res_all<-rbind(res_all,res_or)
}
write.csv(res_all, "results/replication/neale/neale_LST_all.csv", row.names = FALSE)

#MVPA
risd_MVPA<-readRDS("data/replication/neale/rsidinf_correct_MVPA.rds")
outcome_MVPA<-outcome[outcome$Exposure=="Lower MVPA",]
outlist<-outcome_MVPA$Outcome

#Prepare outcome data
for (i in 1:2){
  print(i)
  out<-outlist[i]
  print(out)
  phedat<-fread(paste0("data/replication/neale/",out,".gwas.imputed_v3.both_sexes.tsv.bgz"))  
  #separate for female-only outcomes
  #phedat<-fread(paste0("data/replication/neale/",out,".gwas.imputed_v3.female.tsv.bgz")) 
  outcomes<-left_join(risd_MVPA, phedat, by = c("variant"))
  
  outcomes$af<-ifelse (outcomes$alt==outcomes$minor_allele,outcomes$minor_AF, 
                       ifelse (outcomes$ref==outcomes$minor_allele, 1-outcomes$minor_AF, NA))
  outcomes_gwas <- outcomes[,c("rsid","alt","ref","af","beta","se","pval","n_complete_samples", "chr", "pos")]
  write.csv(outcomes_gwas, paste0("data/replication/neale/", out, "_outcomeall_MVPA.csv"), row.names = FALSE)
}

#Run analysis
exposures_format<-read.csv("data/exposures_format_MVPA.csv")

phenotypes<-read_delim(file = "data/replication/neale/phenotypes.both_sexes.tsv.bgz",delim="\t")
f_phenotypes<-read_delim(file = "data/replication/neale/phenotypes.female.tsv.bgz",delim="\t")

res_all<-data.frame()
for (i in 1:2){
  print(i)
  out<-outlist[i]
  print(out)
  outcomes_gwas<-read.csv(paste0("data/replication/neale/", out, "_outcomeall_MVPA.csv"))
  outcomes_format<-format_data(outcomes_gwas, 
                               type="outcome", 
                               header=TRUE,  
                               snp= exposures_format$SNP,
                               snp_col = "rsid",
                               beta_col = "beta",
                               se_col = "se",
                               eaf_col = "af",
                               effect_allele_col = "alt",
                               other_allele_col = "ref",
                               pval_col = "pval",
                               chr_col = "chr",
                               pos_col = "pos")
  
  sel<-filter(phenotypes, phenotype==out)
  #sel<-filter(f_phenotypes, phenotype==out)
  
  #Load number of cases and controls
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
  #All SNPs have been harmonised so no need for proxy SNPs
  res <- mr(dat, method_list = c("mr_wald_ratio", "mr_ivw"))
  
  #Transform to or
  res_or<-generate_odds_ratios(res)
  res_or$or_ivwci<-sprintf("%.2f (%.2f, %.2f)", res_or$or, res_or$or_lci95, res_or$or_uci95)
  res_or$outcome<-out
  res_or[1,16]<-ncase
  res_or[1,17]<-ncontrol

  res_all<-rbind(res_all,res_or)
}
write.csv(res_all, "results/replication/neale/neale_MVPA_all.csv", row.names = FALSE)