
library(VariantAnnotation)
library(gwasglue)
library(TwoSampleMR)

exposureFile="protein.exposure.csv"     
geneFile="modelGene.diff.txt"      
outcomeName="Crohn's disease"     
setwd("")     

rt=read_exposure_data(filename=exposureFile,
                                sep = ",",
                                snp_col = "SNP",
                                beta_col = "beta.exposure",
                                se_col = "se.exposure",
                                pval_col = "pval.exposure",
                                effect_allele_col="effect_allele.exposure",
                                other_allele_col = "other_allele.exposure",
                                eaf_col = "eaf.exposure",
                                phenotype_col = "exposure",
                                id_col = "id.exposure",
                                samplesize_col = "samplesize.exposure",
                                chr_col="chr.exposure", pos_col = "pos.exposure",
                                clump=FALSE)
files=dir()                          
files=grep(".gz$",files,value=T)    
for(outcomeFile in files){
	outcomeID=gsub(".gz", "", outcomeFile)
	outcomeData=read_outcome_data(snps=rt$SNP,
	                 filename=outcomeFile, sep = "\t",
	                 snp_col = "rsids",
	                 beta_col = "beta",
	                 se_col = "sebeta",
	                 effect_allele_col = "alt",
	                 other_allele_col = "ref",
	                 pval_col = "pval",
	                 eaf_col = "af_alt")
	write.csv(outcomeData, file="outcome.csv", row.names=F)
	
	geneRT=read.table(geneFile, header=T, sep="\t", check.names=F, row.names=1)
	sameGene=intersect(row.names(geneRT), unique(rt$exposure))
	outTab=data.frame()
	for(expoName in sameGene){
		i=paste0(expoName, "__", outcomeID)
		singleExposureFile=paste0(i, ".exposure.csv")
		exposure_set=rt[rt$exposure==expoName,]
		if(nrow(exposure_set)>=3){
			write.csv(exposure_set, file=singleExposureFile, row.names=F)
			exposure_dat=read_exposure_data(filename=singleExposureFile,
		                                sep = ",",
		                                snp_col = "SNP",
		                                beta_col = "beta.exposure",
		                                se_col = "se.exposure",
		                                pval_col = "pval.exposure",
		                                effect_allele_col="effect_allele.exposure",
		                                other_allele_col = "other_allele.exposure",
		                                eaf_col = "eaf.exposure",
		                                phenotype_col = "exposure",
		                                samplesize_col = "samplesize.exposure",
		                                chr_col="chr.exposure", pos_col = "pos.exposure",
		                                clump=FALSE)
		    file.remove(singleExposureFile)
			
			outcome_data=read_outcome_data(snps=exposure_dat$SNP,
				             filename="outcome.csv", sep = ",",
				             snp_col = "SNP",
				             beta_col = "beta.outcome",
				             se_col = "se.outcome",
				             effect_allele_col = "effect_allele.outcome",
				             other_allele_col = "other_allele.outcome",
				             pval_col = "pval.outcome",
				             eaf_col = "eaf.outcome")
			
			outcome_data$outcome=outcomeName
			dat=harmonise_data(exposure_dat, outcome_data)
			dat=dat[dat$pval.outcome>5e-08,]
			snpTab=dat[dat$mr_keep=="TRUE",]
				
			if(nrow(snpTab)>=3){
				mrResult=mr(dat)
				mrTab=generate_odds_ratios(mrResult)
				outTab=rbind(outTab, mrTab)
				pleioTab=mr_pleiotropy_test(dat)
			
				if((nrow(mrTab)>=3) && (mrResult$pval[3]<0.05)){
					if(sum(mrTab$or>1)==nrow(mrTab) | sum(mrTab$or<1)==nrow(mrTab)){
						if(as.numeric(pleioTab$pval)>0.05){
							if((geneRT[expoName,"logFC"] * mrResult$b[3])>0){
								write.csv(snpTab, file=paste0(i, ".table.SNP.csv"), row.names=F)
								write.csv(mrTab, file=paste0(i, ".table.MRresult.csv"), row.names=F)

								if(nrow(snpTab)>3){
									presso=run_mr_presso(dat)
									write.csv(presso[[1]]$`MR-PRESSO results`$`Global Test`, file=paste0(i, ".table.MR-PRESSO_Global.csv"))
									write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file=paste0(i, ".table.MR-PRESSO_Outlier.csv"))
								}
								

								heterTab=mr_heterogeneity(dat)
								write.csv(heterTab, file=paste0(i, ".table.heterogeneity.csv"), row.names=F)

								write.csv(pleioTab, file=paste0(i, ".table.pleiotropy.csv"), row.names=F)
							
								pdf(file=paste0(i, ".scatter_plot.pdf"), width=7, height=6.5)
								p1=mr_scatter_plot(mrResult, dat)
								print(p1)
								dev.off()
								
								res_single=mr_singlesnp(dat)    
								pdf(file=paste0(i, ".forest.pdf"), width=6.5, height=5)
								p2=mr_forest_plot(res_single)
								print(p2)
								dev.off()
								
						
								pdf(file=paste0(i, ".funnel_plot.pdf"), width=6.5, height=6)
								p3=mr_funnel_plot(singlesnp_results = res_single)
								print(p3)
								dev.off()
								
				
								pdf(file=paste0(i, ".leaveoneout.pdf"), width=6.5, height=5)
								p4=mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
								print(p4)
								dev.off()
							}
						}
					}
				}
			}
		}
	}

	write.csv(outTab, file=paste0(outcomeID, "__all.MRresult.csv"), row.names=F)
}




