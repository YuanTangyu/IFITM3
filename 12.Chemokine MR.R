

library(TwoSampleMR) 
library(ieugwasr)

library(TwoSampleMR)

exposureID="prot-a-389" 
outcomeID="ebi-a-GCST90020071"    
setwd("")      

exposure_dat <- extract_instruments(exposureID, p1=5e-06, clump=TRUE)


outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=outcomeID)

dat <- harmonise_data(exposure_dat, outcome_dat)

outTab=dat[dat$mr_keep=="TRUE",]
write.csv(outTab, file="table.SNP.csv", row.names=F)

mrResult=mr(dat)
mrTab=generate_odds_ratios(mrResult)
write.csv(mrTab, file="table.MRresult.csv", row.names=F)

heterTab=mr_heterogeneity(dat)
write.csv(heterTab, file="table.heterogeneity.csv", row.names=F)

pleioTab=mr_pleiotropy_test(dat)
write.csv(pleioTab, file="table.pleiotropy.csv", row.names=F)

pdf(file="pic.scatter_plot.pdf", width=7.5, height=7)
mr_scatter_plot(mrResult, dat)
dev.off()

res_single=mr_singlesnp(dat)     
pdf(file="pic.forest.pdf", width=7, height=5.5)
mr_forest_plot(res_single)
dev.off()

pdf(file="pic.funnel_plot.pdf", width=7, height=6.5)
mr_funnel_plot(singlesnp_results = res_single)
dev.off()

pdf(file="pic.leaveoneout.pdf", width=7, height=5.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
dev.off()


