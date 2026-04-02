#Load packages
library("vroom")
library("dplyr")
library("TwoSampleMR")
library("openxlsx")

#Phenotype one by one
data_endpoint<-vroom("data/endpoint_filter.csv")

#Load exposure data
exposures_format<-read.csv("data/exposures_format.csv")

#Extract SNPs
for (i in 1:1719) {
  out<-data_endpoint$phenocode[i]
  outcomes<-vroom(paste0("~/finngen_R11_",out,".gz"))
  #Transform outcomes
  outcomes_format<-format_data(outcomes, 
                               type="outcome", 
                               header=TRUE,  
                               snp= exposures_format$SNP,
                               snp_col = "rsids",
                               beta_col = "beta",
                               se_col = "sebeta",
                               eaf_col = "af_alt",
                               effect_allele_col = "alt",
                               other_allele_col = "ref",
                               pval_col = "pval",
                               min_pval = 1e-1000,
                               chr_col = "#chrom",
                               pos_col = "pos",
                               log_pval = FALSE)
  
  outcomes_format$outcome <- out
  write.csv(outcomes_format, paste0("data/outcome/outcome_",out,".csv"), row.names = FALSE)
  gc()}

#Combine extracted SNPs
file_list<-list.files(path="~", pattern=".csv$", full.names=T)
df_list<-lapply(file_list, read_csv)
outcomes_all<-bind_rows(df_list)
write.csv(outcomes_all,"data/outcome/outcomes.csv")

#Check proxy SNPs
table<-LDproxy(snp, #missing SNPs
               pop = "EUR", 
               r2d = "r2", 
               token = "~", 
               file = FALSE,
               genome_build = "grch38",
               win_size = "500000")
#proxy rs4863692 for rs5862438
#proxy rs1188869 for rs200307517
#proxy rs35584753 for rs58541850

#Extract proxy SNPs
for (i in 1:1719)
{ 
print(i)
out<-data_endpoint$Phenocode[i]
print(out)
#access data
url <- paste0("https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_",out,".gz")
data <- fread(url, select = c(5, 9, 10))
data_snps <- data[c("rs4863692", "rs1188869", "rs35584753"), on = "rsids"]
write.csv(data_snps, paste0("data/outcome/outcome_proxy_",out,".csv"), row.names = FALES)
rm(data)
}

#Load combined extracted SNPs
outcomes_origin<-read.csv("data/outcome/outcomes.csv")
outcomes_format<-outcomes_origin

#Add proxy SNPs
for (i in 1:1719)
{
  out<-data_endpoint$Phenocode[i]
  print(out)
  info = read.csv(paste0("data/outcome/outcome_proxy_",out,".csv"))
  
  #rsids=rs5862438 
  outcomes_format <- outcomes_format %>%
    add_row(
      chr.outcome=4, 
      pos.outcome=165752075, 
      other_allele.outcome="ATT", 
      effect_allele.outcome="A", 
      SNP="rs5862438", 
      beta.outcome= info[info$rsids=="rs4863692",2], 
      se.outcome= info[info$rsids=="rs4863692",3], 
      eaf.outcome= 0.3594, 
      outcome=out, 
      mr_keep.outcome=TRUE, 
      id.outcome=outcomes_format$id.outcome[outcomes_format$outcome==out][1])
  
  #rsids=rs200307517
  outcomes_format <- outcomes_format %>%
    add_row(
      chr.outcome=6, 
      pos.outcome=165752075, 
      other_allele.outcome="C", 
      effect_allele.outcome="CT", 
      SNP="rs200307517", 
      beta.outcome= info[info$rsids=="rs1188869",2], 
      se.outcome= info[info$rsids=="rs1188869",3], 
      eaf.outcome= 0.0585, 
      outcome=out, 
      mr_keep.outcome=TRUE, 
      id.outcome=outcomes_format$id.outcome[outcomes_format$outcome==out][1])
  
  #rsids=rs58541850 
  outcomes_format <- outcomes_format %>%
    add_row(
      chr.outcome=6, 
      pos.outcome=165752075, 
      other_allele.outcome="G", 
      effect_allele.outcome="A", 
      SNP="rs58541850", 
      beta.outcome= info[info$rsids=="rs35584753",9], 
      se.outcome= info[info$rsids=="rs35584753",10], 
      eaf.outcome= 0.0585, 
      outcome=out, 
      mr_keep.outcome=TRUE, 
      id.outcome=outcomes_format$id.outcome[outcomes_format$outcome==out][1])
}

write.csv(outcomes_format, "data/outcome/outcomes_proxy.csv", row.names = FALES)

#Add case number
outcome<-read.xlsx("data/outcome/outcome_info.xlsx")
outcome$samplesize.outcome<- outcome$Num_cases+outcome$Num_controls
samplesize<-select(outcome, Phenocode,samplesize.outcome)
colnames(samplesize)[1]<-"outcome"
outcomes_origin<-read.csv("data/outcome/outcomes_proxy.csv")
outcomes_format<-left_join(outcomes_format, samplesize, by= "outcome")
write.csv(outcomes_format, "data/outcome/outcomes_proxy_num.csv", row.names = FALSE)

