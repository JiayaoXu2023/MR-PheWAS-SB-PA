#geneatlas http://geneatlas.roslin.ed.ac.uk/downloads/

# Load packages
library("vroom")
library("dplyr")
library(SNPlocs.Hsapiens.dbSNP155.GRCh37)
library(GenomicRanges)
library(readr)
library("TwoSampleMR")
libaray("openxlsx")

#Clean SNPs rsid for each chrom
for (chrom in c(1:20, 22)){
data<-vroom(paste0("data/replication/geneatlas/snps.imputed.chr", chrom, ".csv.gz"))
data_select<-data %>% filter(!grepl("^rs", SNP))
data_select$CHR<-chrom
gr <- GRanges(seqnames = as.character(data_select$CHR),
              ranges = IRanges(start = data_select$Position, end = data_select$Position))
hits <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh37, gr)
hits_df <- as.data.frame(hits)
rsid<-hits_df[,c(2,4)]
colnames(rsid)[2]<-"rsid"
data_rsid<-left_join(data, rsid, by= c("Position"="pos") )
idx<-grepl("^rs", data_rsid$SNP)
data_rsid$rsid[idx]<-data_rsid$SNP[idx]
write.csv(data_rsid,paste0("data/replication/geneatlas/snps.imputed.chr",chrom,"_rsid.csv"), row.names = F)}

#Load instrumental variables for LST
exposures_format<-read.csv("data/exposures_format_LST.csv")
#Load outcome names
outcome_list<-c(
  "clinical_c_Block_J40-J47", 
  "clinical_c_Block_J30-J39",
  "clinical_c_Block_K20-K31",
  "clinical_c_Block_K00-K14",
  "clinical_c_Block_K55-K64",
  "clinical_c_Block_M50-M54",
  "clinical_c_Block_M05-M14",
  "clinical_c_Block_N10-N16"
  )

#Load female-only outcome names
outcome_list<-c(
"clinical_c_Block_N80-N98",
"clinical_c_Block_N70-N77")

#Select SNPs in each chrom and format data
for (i in 1:8){
  print(i)
  out<-outcome_list[i]
  print(out)
  outcome_combine<-data.frame()
for (chrom in c(1:20, 22)){
  print(chrom)
  url <- paste0("http://static.geneatlas.roslin.ed.ac.uk/gwas/allWhites/imputed/data.copy/imputed.allWhites.",out,".chr",chrom,".csv.gz")
  #url <- paste0("http://static.geneatlas.roslin.ed.ac.uk/gwas/allWhites/imputed/data.byGender/imputed.allWhites.", out,".females.chr", chrom, ".csv.gz")
  data <- read_table(url)
  data2<-read.csv(paste0("data/replication/geneatlas/snps.imputed.chr",chrom,"_rsid.csv"))
  data_all<-left_join(data, data2, by = "SNP")
  outcome<-data_all [, c(2,4,5,6,7,8,10,13)]
  colnames(outcome) <- c("A2", "BETA", "SE", "pval", "pos", "A1", "MAF", "rsids")
  outcome[,"chrom"]<-chrom
  outcome <- outcome[!is.na(outcome$rsids), ]
  #Format outcome data
  outcomes_format<-format_data(outcome, 
                               type="outcome", 
                               header=TRUE,  
                               snp= exposures_format$SNP,
                               snp_col = "rsids",
                               beta_col = "BETA",
                               se_col = "SE",
                               eaf_col = "MAF",
                               effect_allele_col = "A2",
                               other_allele_col = "A1",
                               pval_col = "pval",
                               chr_col = "chrom",
                               pos_col = "pos")
  outcomes_format$outcome<- out
  print(nrow(outcomes_format))
  write.csv(outcomes_format,paste0("results/replication/geneatlas/outcomes/outcome_",out,"_",chrom,".csv"), row.names = F)
  outcome_combine<-rbind(outcome_combine,outcomes_format)
}
  write.csv(outcome_combine,paste0("results/replication/geneatlas/outcomes/outcome_",out,"_all.csv"), row.names = F)  
}

#run MR

exposures_format<-read.csv("data/exposures_format_LST.csv")

#Outcome list
outcome_list<-c(
  "c_Block_J40-J47", 
  "c_Block_J30-J39",
  "c_Block_K20-K31",
  "c_Block_K00-K14",
  "c_Block_K55-K64",
  "c_Block_M50-M54",
  "c_Block_M05-M14", 
  "c_Block_N80-N98",
  "c_Block_N70-N77",
  "c_Block_N10-N16")

#Load number of cases and controls
phenotypes<-read.xlsx(xlsxFile ="data/replication/geneatlas/phenotypes.xlsx")
#Individuals of the unaffected sex were excluded from the analysis for these phenotypes 

for (i in c(1:10)){
  print(i)
  out<-outcome_list[i]
  print(out)
  sel<-filter(phenotypes, phenotype==out)
  ncase<-sel$n_cases
  ncase
  ncontrol<-sel$n_controls
  ncontrol
  u<-ncase/(ncase+ncontrol)
  uu<-u*(1-u)
  
  #Transform beta se
  outcomes_format<-read.csv(paste0("results/replication/geneatlas/outcomes/outcome_clinical_",out,"_all.csv"))
  outcomes_format$beta.outcome<-outcomes_format$beta.outcome/uu
  outcomes_format$beta.outcome[1]
  outcomes_format$se.outcome<-outcomes_format$se.outcome/uu
  outcomes_format$se.outcome[1]
  outcomes_format$id.outcome<-outcomes_format$id.outcome[1]
  head(outcomes_format)
  
  #Run MR IVW
  dat <- harmonise_data(exposures_format, outcomes_format, action = 2)
  #Proxies found for 0 out of 17 missing SNPs.
  res <- mr(dat, method_list = c("mr_wald_ratio", "mr_ivw"))
  
  #Transform to or
  res_or<-generate_odds_ratios(res)
  res_or$or_ivwci<-sprintf("%.2f (%.2f, %.2f)", res_or$or, res_or$or_lci95, res_or$or_uci95)
  write.csv(res_or, paste0("results/replication/geneatlas/", out, "_LST.csv"), row.names = FALSE)}

