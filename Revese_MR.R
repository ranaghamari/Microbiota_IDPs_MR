
##### libraries requirement #####

library(readr)
library(dplyr)
library(remotes)
library(TwoSampleMR)
library(devtools)
library(MRInstruments)
library(MRPRESSO)
library(ieugwasr)
library(plyr)
library(genetics.binaRies)
library(tidyverse)

##### Preparing exposure file ####

# Here we use local pannel for clumping:
a <- tophits(id = s_IDPs_id, clump = 0)

# The studies listed below have unknown SNPs which resulted in error during clumping. Thus, they're removed from the study
genetics.binaRies::get_plink_binary()

a01 <- filter(a,a$id!="ubm-b-876")
a01 <- filter(a01,a01$id!="ubm-b-872")
a01 <- filter(a01,a01$id!="ubm-b-1236")
a01 <- filter(a01,a01$id!="ubm-b-518")
a01 <- filter(a01,a01$id!="ubm-b-1286")
a01 <- filter(a01,a01$id!="ubm-b-578")
a01 <- filter(a01,a01$id!="ubm-b-574")
a01 <- filter(a01,a01$id!="ubm-b-608")
a01 <- filter(a01,a01$id!="ubm-b-377")
a01 <- filter(a01,a01$id!="ubm-b-1119")
a01 <- filter(a01,a01$id!="ubm-b-774")
a01 <- filter(a01,a01$id!="ubm-b-2091")
a01 <- filter(a01,a01$id!="ubm-b-1839")
a01 <- filter(a01,a01$id!="ubm-b-1812")

a01 <- a01 %>% arrange((desc(id)))

b <- ld_clump(
  dplyr::tibble(rsid=a01$rsid, pval=a01$pval, id=a01$id),
  plink_bin = "C:/R/Microbiota_IDPs_MR/plink_win64_20231211 (1)/plink.exe",
  bfile = "C:/R/Microbiota_IDPs_MR/1kg.v3/EUR")

b000 <- ld_clump(
  dplyr::tibble(rsid=b$rsid, pval=b$pval, id=b$id),
  plink_bin = "C:/R/Microbiota_IDPs_MR/plink_win64_20231211 (1)/plink.exe",
  bfile = "C:/R/Microbiota_IDPs_MR/1kg.v3/EUR")

idp_exp_dat <- inner_join(b, a01)

write.csv(idp_exp_dat, "C:/R/Microbiota_IDPs_MR/idp_exp_dat.csv")

##### Preparing outcome file ####

idp_snps_clumped <- idp_exp_dat$rsid
microbiota_h2_id <-Microbiota_H2$id

microbiota_out_dat <- extract_outcome_data(
  idp_snps_clumped,
  microbiota_h2_id,
  proxies = TRUE,
  rsq = 0.8,
  align_alleles = 1,
  palindromes = 1,
  maf_threshold = 0.3,
  opengwas_jwt = ieugwasr::get_opengwas_jwt(),
  splitsize = 10000,
  proxy_splitsize = 500
)

write.csv(microbiota_out_dat, "C:/R/Microbiota_IDPs_MR/microbiota_out_dat.csv")

##### Harmonization ####

reverse_harmonized_data <- harmonise_data(
  exposure_dat = idp_exp_dat, 
  outcome_dat = microbiota_out_dat, action = 2)

# Since our "reverse_harmonized_data" variable does not have true samplesize.exposure, we must add it from "IDPs_H2_No_fMRI" variable using left_join:
# First, we select outcome id and samplesize from "IDPs_H2_No_fMRI" file:
IDPs_H2_No_fMRI_only_samplesize <- IDPs_H2_No_fMRI[,c(2,15)]

# Secondly, we change the "id" to "id.exposure" and "sample_size" to "samplesize.exposure":
names(IDPs_H2_No_fMRI_only_samplesize)[1] <- "id.exposure"
names(IDPs_H2_No_fMRI_only_samplesize)[2] <- "samplesize.exposure"

# Thirdly, we remove samplesize.outcome column from "forward_harmonized_data" variable:
reverse_harmonized_data <- reverse_harmonized_data[,-36]

# Next, we left_join "forward_harmonized_data" and "IDPs_H2_No_fMRI_only_samplesize" variables:
reverse_harmonized_data_joined_samplesize <- left_join(reverse_harmonized_data,
                                                       IDPs_H2_No_fMRI_only_samplesize)

reverse_harmonized_data <- reverse_harmonized_data_joined_samplesize

write.csv(reverse_harmonized_data, "C:/R/Microbiota_IDPs_MR/reverse_harmonized_data.csv")

##### Mendelian Randomization ####

reverse_mr <- mr(reverse_harmonized_data,
                 method_list=c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))

write.csv(reverse_mr, "C:/R/Microbiota_IDPs_MR/reverse_mr.csv")

##### Heterogeneity Test ####

reverse_heterogeneity <- mr_heterogeneity(reverse_harmonized_data, 
                                          parameters = default_parameters(),
                                          method_list = subset(mr_method_list(), 
                                                               heterogeneity_test)$obj)

write.csv(reverse_heterogeneity, "C:/R/Microbiota_IDPs_MR/reverse_heterogeneity.csv")

##### Pleiotropy Test ####

reverse_pleiotropy <- mr_pleiotropy_test(reverse_harmonized_data)

write.csv(reverse_pleiotropy, "C:/R/Microbiota_IDPs_MR/reverse_pleiotropy.csv")

##### leave-one-out MR ####

reverse_mr_loo <- mr_leaveoneout(reverse_harmonized_data, 
                                 parameters = default_parameters(),
                                 method = TwoSampleMR::mr_ivw) 

write.csv(reverse_mr_loo, "C:/R/Microbiota_IDPs_MR/reverse_mr_loo.csv")

##### Single SNP MR ####

reverse_mr_single <- mr_singlesnp(reverse_harmonized_data, parameters = default_parameters(),
                          single_method = 'mr_wald_ratio',         
                          all_method = c('mr_ivw', 'mr_egger_regression'))

write.csv(reverse_mr_single, "C:/R/Microbiota_IDPs_MR/reverse_mr_single.csv")

##### sign concordance test ####

reverse_mr_fit_sign <- mr_sign(b_exp = reverse_harmonized_data$beta.exposure,
                               b_out = reverse_harmonized_data$beta.outcome)

##### robust adjusted profile score ####

## This works:
reverse_mr_fit_raps_all <- mr.raps.all(b_exp = reverse_harmonized_data$beta.exposure,
                                       b_out = reverse_harmonized_data$beta.outcome,
                                       se_exp = reverse_harmonized_data$se.exposure,
                                       se_out = reverse_harmonized_data$se.outcome)

##### MR-PRESSO (pleiotropy and outliers removal) ####

reverse_mr_mrpresso <- mr_presso(BetaOutcome = 'beta.outcome',
                                 BetaExposure = 'beta.exposure',
                                 SdOutcome = 'se.outcome',
                                 SdExposure = 'se.exposure',
                                 data = reverse_harmonized_data,              
                                 OUTLIERtest = FALSE,
                                 DISTORTIONtest = FALSE,
                                 SignifThreshold = 0.05,
                                 NbDistribution = 5e6,
                                 seed = NULL) 

##### Directionality Test ####

reverse_harmonized_data$r.outcome <- get_r_from_pn(reverse_harmonized_data$pval.outcome,
                                                   reverse_harmonized_data$samplesize.outcome)

reverse_harmonized_data$r.exposure <- get_r_from_pn(reverse_harmonized_data$pval.exposure,
                                                    reverse_harmonized_data$samplesize.exposure)
                                                  
reverse_mr_steiger <- directionality_test(reverse_harmonized_data)

write.csv(reverse_mr_steiger, "C:/R/Microbiota_IDPs_MR/reverse_mr_steiger.csv")

##### Joining Results ####

# Reverse-MR-Heterogeneity:
reverse_mr_heterogeneity_result <- left_join(reverse_mr,reverse_heterogeneity)

write.csv(reverse_mr_heterogeneity_result, 
          "C:/R/Microbiota_IDPs_MR/reverse_mr_heterogeneity_result.csv")

# Reverse-MR-Pleiotropy:
reverse_pleiotropy_for_join <- reverse_pleiotropy
names(reverse_pleiotropy_for_join)[6] <- "se_pleiotropy"
names(reverse_pleiotropy_for_join)[7] <- "pval_pleiotropy"
reverse_mr_pleiotropy_result <- left_join(reverse_mr,reverse_pleiotropy_for_join)

write.csv(reverse_mr_pleiotropy_result, 
          "C:/R/Microbiota_IDPs_MR/reverse_mr_pleiotropy_result.csv")

# Reverse-MR-LOO:
reverse_mr_loo_for_join <- reverse_mr_loo
names(reverse_mr_loo_for_join)[7] <- "b_LOO"
names(reverse_mr_loo_for_join)[8] <- "se_LOO"
reverse_mr_loo_result <- left_join(reverse_mr,reverse_mr_loo_for_join,
                                   relationship = "many-to-many")

write.csv(reverse_mr_loo_result, 
          "C:/R/Microbiota_IDPs_MR/reverse_mr_loo_result.csv")

# Reverse-MR-SingleSNP:
reverse_mr_single_join <- reverse_mr_single
names(reverse_mr_single_join)[7] <- "b_SingleSNP"
names(reverse_mr_single_join)[8] <- "se_SingleSNP"
reverse_mr_SingleSNP_result <- left_join(reverse_mr,reverse_mr_single_join,
                                         relationship = "many-to-many")

write.csv(reverse_mr_SingleSNP_result, 
          "C:/R/Microbiota_IDPs_MR/reverse_mr_SingleSNP_result.csv")

# Reverse-MR-Directionality:
reverse_mr_directionality_result <- left_join(reverse_mr, reverse_mr_steiger)
write.csv(reverse_mr_directionality_result, 
          "C:/R/Microbiota_IDPs_MR/reverse_mr_directionality_result.csv")

##### Plotting ####

# Scatter plot:
reverse_scatter_plot <- mr_scatter_plot(reverse_mr)

# forest-plot:
reverse_forest_plot <- mr_forest_plot(reverse_mr_single)

# LOO plot: 
reverse_loo_plot <- mr_leaveoneout_plot(reverse_mr_loo)