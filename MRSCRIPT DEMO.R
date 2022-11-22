library(TwoSampleMR)
ao <- available_outcomes()
exposure_dat <- extract_instruments(c('ukb-b-13423'))
exposure_dat <- clump_data(exposure_dat)
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('ebi-a-GCST006250'), proxies = 0, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 3)
mr_results <- mr(dat)

#COMPLETE MR RESULTS
mr_singlesnp(dat,single_method = "mr_meta_fixed") #or
mr_res_complete<-mr_singlesnp(dat,parameters = default_parameters(),single_method = "mr_wald_ratio",all_method = c("mr_ivw","mr_ivw_fe","mr_ivw_mre","mr_simple_mode","mr_weighted_mode","mr_simple_median","mr_weighted_median","mr_egger_regression")
                              
                              #Heterogeneity test
                              heterogeneity_test<-mr_heterogeneity(dat, method_list=c("mr_egger_regression","mr_ivw"))
                              pleiotropy_test<-mr_pleiotropy_test(dat)
                              mr_rucker_cooksdistance(dat, parameters = default_parameters())
                              
                              #Produce PLOT
                              mr_scatter_plot(mr_results,dat)
                              
                              res_single <- mr_singlesnp(dat)
                              mr_forest_plot(res_single)
                              
                              res_loo <- mr_leaveoneout(dat)
                              mr_leaveoneout_plot(res_loo)
                              
                              mr_funnel_plot(res_single)
                              
                              #Rename
                              dat$BetaXG<-dat$beta.exposure
                              dat$seBetaXG<-dat$se.exposure
                              dat$BetaYG<-dat$beta.outcome
                              dat$seBetaYG<-dat$se.outcome
                              BetaXG <- dat$BetaXG
                              BetaYG <- dat$BetaYG
                              seBetaXG <- dat$seBetaXG
                              seBetaYG <- dat$seBetaYG
                              
                              BYG <- BetaYG*sign(BetaXG) #absolutenumber
                              BXG <- abs(BetaXG) 
                              
                              BYG2 <- BetaYG*sign(BetaXG) #absolutenumber
                              BXG2 <- abs(BetaXG)
                              
                              #odd ratio
                              odd_ratio<-generate_odds_ratios(mr_res_complete)
                              
                              #r
                              p<-dat$pval.exposure
                              n<-dat$samplesize.exposure
                              get_r_from_pn(p, n)
                              get_r_from_bsen(BXG,seBetaXG,n)
                              
                              #r2
                              r2<-add_rsq(dat)
                              
                              #F-statistics
                              F   = BXG^2/seBetaXG^2
                              mF  = mean(F)
                              
                              #Isq
                              function(y,s){
                                k          = length(y)
                                w          = 1/s^2; sum.w  = sum(w)
                                mu.hat     = sum(y*w)/sum.w
                                Q          = sum(w*(y-mu.hat)^2)
                                Isq        = (Q - (k-1))/Q
                                Isq        = max(0,Isq)
                                return(Isq)
                              }
                              
                              #Isq
                              Isq_unweighted <- Isq(BXG,seBetaXG)
                              Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) 
                              