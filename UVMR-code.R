####UVMR-code   ####
##install package
install.packages("devtools")
devtools::install_github("MRCIEU/TwoSampleMR")
install.packages("ggplot2") 

#library package
library(TwoSampleMR)
library(ggplot2)

setwd("D:/R")

####Read Exposure Data####
## way 1
exp_dat <- read_exposure_data(
  filename = "All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq", #filename
  sep = "\t",
  snp_col = "SNP",  
  beta_col = "b",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "Freq1.Hapmap",
  pval_col = "p",
  samplesize_col = "N"
)

#Understanding variable names for data
bmi<-read.table(file = "All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq",sep = "\t",header = T)
head(bmi) 

#Screening P-value
exp_dat2 <- exp_dat[exp_dat$pval.exposure < 5e-8,]  

#Variables that exclude chain imbalances
exp_dat3 <- clump_data(exp_dat2)

##way 2：read ID on line
exp_dat4 <- extract_instruments(
  outcomes = 'ieu-b-40'
)


#Save data
write.csv(exp_dat3, file="exp_dat.csv")

##Read outcome data

out_dat <- read_outcome_data(
  snps = exp_dat3$SNP,
  filename = "cad.add.160614.website.txt",
  sep = "\t",
  snp_col = "markername",
  beta_col = "beta",
  se_col = "se_dgc",
  effect_allele_col = "effect_allele",
  other_allele_col = "noneffect_allele",
  eaf_col = "effect_allele_freq",
  pval_col = "p_dgc")

##Read ID on line
out_dat2 <- extract_outcome_data(
  snps = exp_dat3$SNP,
  outcomes = 'ieu-a-7')

#Harmonized effects, combined data
dat <- harmonise_data(
  exposure_dat =  exp_dat3, 
  outcome_dat = out_dat
)

#Save data
write.csv(dat, file="dat.csv")

#UVMR analysis
res <- mr(dat)
res


#heterogeneity test
mr_heterogeneity(dat)

#Finding outliers
run_mr_presso(dat,NbDistribution = 1000)

##pleiotropy test
mr_pleiotropy_test(dat)

res_single <- mr_singlesnp(dat)
res_single

##leave one out
res_loo <- mr_leaveoneout(dat)
res_loo 

##scatterplot 
p1 <- mr_scatter_plot(res, dat)
p1[[1]]
ggsave(p1[[1]], file="res.pdf", width=7, height=7)

#forest plot
p2 <- mr_forest_plot(res_single)
p2[[1]]
ggsave(p2[[1]], file="p2.pdf", width=7, height=7)

##leave-one-out_plot 
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]
ggsave(p3[[1]], file="p3.pdf", width=7, height=7)

##funnel plot 
res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(res_single)
p4[[1]]
ggsave(p4[[1]], file="p4.pdf", width=7, height=7)




