setwd("R")
library(TwoSampleMR)
id_exposure <- c("ebi-a-GCST90000047","ieu-a-1239","ieu-b-4877","ukb-b-19953") 
id_outcome <- "ieu-b-5099"
exposure_dat <- mv_extract_exposures(id_exposure)
outcome_dat <- extract_outcome_data(exposure_dat$SNP, id_outcome) 
mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)
res <- mv_multiple(mvdat)
res_OR<-generate_odds_ratios(res$result)
res_OR
write.table(res_OR, file="MVMR.xls",sep="\t",quote=F)
SummaryStats<-cbind(mvdat[["outcome_beta"]],
                    mvdat[["exposure_beta"]][,1],
                    mvdat[["exposure_beta"]][,2],
                    mvdat[["exposure_beta"]][,3],
                    mvdat[["exposure_beta"]][,4],
                    mvdat[["exposure_se"]][,1],
                    mvdat[["exposure_se"]][,2],
                    mvdat[["exposure_se"]][,3],
                    mvdat[["exposure_se"]][,4],
                    mvdat[["outcome_se"]])
SummaryStats<-data.frame(SummaryStats)
library(MRPRESSO)
mr_presso(BetaOutcome = "X1",
          BetaExposure = c("X2", "X3","X4","X5"), 
          SdOutcome = "X10", 
          SdExposure = c("X6", "X7","X8","X9"),
          OUTLIERtest = TRUE, 
          DISTORTIONtest = TRUE, 
          data = SummaryStats,
          NbDistribution = 1000, 
          SignifThreshold = 0.05)
library(MendelianRandomization)
MRMVInputObject_1<- mr_mvinput(bx = cbind(SummaryStats$X2,SummaryStats$X3,SummaryStats$X4,SummaryStats$X5),
                               bxse = cbind(SummaryStats$X6,SummaryStats$X7,SummaryStats$X8,SummaryStats$X9),
                               by = SummaryStats$X1, 
                               byse = SummaryStats$X10)
MRMVInputObject_1
MRMVObject <- mr_mvivw(MRMVInputObject_1, 
                       model = "default",
                       correl = FALSE,
                       distribution = "normal",
                       alpha = 0.05)
MRMVObject
MRMVObject<-mr_mvegger(
  MRMVInputObject_1,
  orientate = 1,
  correl = FALSE,
  distribution = "normal",
  alpha = 0.05)
MRMVObject
MRMVObject<-mr_mvlasso(
  MRMVInputObject_1,
  orientate = 1,
  distribution = "normal",
  alpha = 0.05,
  lambda = numeric(0)
)
MRMVObject
MRMVObject<-mr_mvmedian(
  MRMVInputObject_1,
  distribution = "normal",
  alpha = 0.05,
  iterations = 10000,
  seed = 314159265
)
MRMVObject
MRMVObject<-mr_mvegger(
  MRMVInputObject_1,
  orientate = 1,
  correl = FALSE,
  distribution = "normal",
  alpha = 0.05)
MRMVObject
write.csv(exposure_dat, file="exposure_dat.csv")
write.csv(outcome_dat, file="outcome_dat.csv")
