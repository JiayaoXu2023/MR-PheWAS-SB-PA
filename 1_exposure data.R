#Load packages
library(TwoSampleMR)
library(dplyr)
library(data.table)

#GWAS https://www.ebi.ac.uk/gwas/publications/36071172
data <- read.delim("data/GCST90104339_buildGRCh37.tsv", header = TRUE, sep = "\t")

#Select SNPs with p_value < 5e-8
subdata<-subset(data, p_value<5e-8)

#Rename
names(subdata)[names(subdata)=="rsid"]<-"SNP"
names(subdata)[names(subdata)=="chromosome"]<-"chr_name"
names(subdata)[names(subdata)=="base_pair_location"]<-"chrom_start"
names(subdata)[names(subdata)=="p_value"]<-"pval.exposure"

#Clump
clumped_data<-clump_data(
  subdata,
  clump_kb = 10000,
  clump_r2 = 0.001, 
  clump_p1 = 1,
  bfile = "~/data_maf0.01_rs_ref",
  plink_bin = "~/plink"
)

#Rename
names(clumped_data)[names(clumped_data)=="BETA"]<-"beta.exposure"
names(clumped_data)[names(clumped_data)=="SE"]<-"se.exposure"
names(clumped_data)[names(clumped_data)=="SampleSize"]<-"samplesize.exposure"

#Calculate F statistics and R square
clumped_data$F_stat <- clumped_data$beta.exposure^2 / clumped_data$se.exposure^2
clumped_data$Rsq <- exp_data$F_stat/(clumped_data$samplesize.exposure-2+clumped_data$F_stat)

#Calculate mean F statistics and total R square
F_stat<-mean(clumped_data$F_stat)
Rsq<-sum(clumped_data$Rsq)

#Save 
write.csv(clumped_data,"data/exposures_format.csv", row.names = FALSE)
