
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
library(MRMix)
library(mr.raps)

##### Preparing exposure file ####

# Getting all available studies (available outcomes - ao) from IEU GWAS database:
ao <- available_outcomes(access_token=NULL)
ieugwasr::get_opengwas_jwt()

# Select the Kurilshikov A et al. study from available outcomes (ao variable):
microbiota_data <- filter(ao, pmid == 33462485) # it contains all 211 taxa GWAS mentioned in the study
write.csv(microbiota_data, "C:/R/Microbiota_IDPs_MR/microbiota.csv")

# Import Microbiota summary statistics combined with H2 results:
library(readr)
Microbiota_H2 <- read_csv("Microbiota_H2.csv")

# Here we define a variable for filtering Microbiota taxa based on H2 significance. Here we only keep taxa which H2-H2(SE) > 0 
microbiota_H2_index <- c()
for (i in 1:length(Microbiota_H2$id)) {
  microbiota_H2_index[i] <- Microbiota_H2$H2[i] - Microbiota_H2$`H2(SE)`[i]
}

Microbiota_H2 <- cbind(Microbiota_H2, microbiota_H2_index)

Microbiota_H2 <- filter(Microbiota_H2, 
                        Microbiota_H2$microbiota_H2_index > 0 )


# extract the instrument (variants) from Kurilshikov A et al. study:
microbiota_h2_id <-Microbiota_H2$id # Isolate each taxa GWAS id for next step

# Extract significant variants from the Kurilshikov A et al. study at the significant level of 5e-05:
microbiota_exp_dat <- extract_instruments(outcomes = microbiota_h2_id,
                                              p1 = 5e-05,
                                              clump = FALSE,
                                              p2 = 5e-08,
                                              r2 = 0.001,
                                              kb = 10000)

# Clump filtered Kurilshikov A et al. study variants:
microbiota_exp_dat <- clump_data(microbiota_exp_dat, 
                                     clump_kb = 10000,
                                     clump_r2 = 0.001)

# Save final exposure file as .csv file in the project directory:
write.csv(microbiota_exp_dat,
          "C:/R/Microbiota_IDPs_MR/Microbiome_exposure_-5_h2.csv")

##### Preparing outcome file ####

# Isolating IDPs from the Elliott et al. study in 2021
# According to this study, there are 3935 IDPs
idp <- ao[grepl("ubm-b",ao$id),]  

write.csv(idp, "C:/R/Microbiota_IDPs_MR/IDPs.csv")

# Import IDPs joined with H2 and H2(SE):
IDPs_H2 <- read_csv("IDPs_H2.csv")

# Here we define a variable for filtering IDPs based on H2 significance. Here we only keep IDPs which H2-H2(SE) > 0 

H2_index <- c()
for (i in 1:length(IDPs_H2$id)) {
  H2_index[i] <- IDPs_H2$Heritability[i] - IDPs_H2$`Heritability(SE)`[i]
}

IDPs_H2 <- cbind(IDPs_H2, H2_index)

IDPs_H2 <- filter(IDPs_H2, IDPs_H2$H2_index > 0 )

# Remove functional MRI IDPs

IDPs_H2_fMRI <- IDPs_H2[grepl("fMRI",IDPs_H2$trait),]  
IDPs_H2_No_fMRI <- IDPs_H2[!(IDPs_H2$id %in% IDPs_H2_fMRI$id),]

write.csv(IDPs_H2_No_fMRI, "C:/R/Microbiota_IDPs_MR/IDPs_H2_No_fMRI.csv")

# preparing final outcome file by identifying shared SNPs between exposure and outcome study

s_IDPs_id <- IDPs_H2_No_fMRI$id
microbiota_snps <- microbiota_exp_dat$SNP

idp_out_dat <- extract_outcome_data(
  microbiota_snps,
  s_IDPs_id,
  proxies = TRUE,
  rsq = 0.8,
  align_alleles = 1,
  palindromes = 1,
  maf_threshold = 0.3,
  opengwas_jwt = ieugwasr::get_opengwas_jwt(),
  splitsize = 10000,
  proxy_splitsize = 500
)

write.csv(idp_out_dat,"C:/R/Microbiota_IDPs_MR/IDPs_outcome_file.csv")

##### Harmonization ####

forward_harmonized_data <- harmonise_data(
  exposure_dat = microbiota_exp_dat, 
  outcome_dat = idp_out_dat, action = 2)

# Since our "forward_harmonized_data" variable does not have true samplesize.outcome, we must add it from "IDPs_H2_No_fMRI" variable using left_join:
# First, we select outcome id and samplesize from "IDPs_H2_No_fMRI" file:
IDPs_H2_No_fMRI_only_samplesize <- IDPs_H2_No_fMRI[,c(2,15)]

# Secondly, we change the "id" to "id.outcome" and "sample_size" to "samplesize.outcome":
names(IDPs_H2_No_fMRI_only_samplesize)[1] <- "id.outcome"
names(IDPs_H2_No_fMRI_only_samplesize)[2] <- "samplesize.outcome"

# Thirdly, we remove samplesize.outcome column from "forward_harmonized_data" variable:
forward_harmonized_data <- forward_harmonized_data[,-17]

# Next, we left_join "forward_harmonized_data" and "IDPs_H2_No_fMRI_only_samplesize" variables:
forward_harmonized_data_joined_samplesize <- left_join(forward_harmonized_data,
                                                       IDPs_H2_No_fMRI_only_samplesize)

forward_harmonized_data <- forward_harmonized_data_joined_samplesize

# Finally, we save harmonized data:

write.csv(forward_harmonized_data,"C:/R/Microbiota_IDPs_MR/harmonized_data_file.csv")

##### Mendelian Randomization ####

forward_mr <- mr(forward_harmonized_data,
                 method_list=c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))

write.csv(forward_mr,"C:/R/Microbiota_IDPs_MR/forward_mr.csv")

forward_mr <- forward_mr %>% arrange((desc(exposure)))
forward_mr_genus <- forward_mr[grepl("genus",forward_mr$exposure),]
forward_mr_genus <- forward_mr_genus %>% arrange((asc(pval)))

##### Heterogeneity Test ####

forward_heterogeneity <- mr_heterogeneity(forward_harmonized_data, 
                                          parameters = default_parameters(),
                                                       method_list = subset(mr_method_list(), 
                                                                            heterogeneity_test)$obj)

write.csv(forward_heterogeneity,"C:/R/Microbiota_IDPs_MR/forward_heterogeneity.csv")

##### Pleiotropy Test ####

forward_pleiotropy <- mr_pleiotropy_test(forward_harmonized_data)

write.csv(forward_pleiotropy,"C:/R/Microbiota_IDPs_MR/forward_pleiotropy.csv")

##### leave-one-out MR ####

forward_mr_loo <- mr_leaveoneout(forward_harmonized_data, 
                                      parameters = default_parameters(),
                                      method = TwoSampleMR::mr_ivw) 

write.csv(forward_mr_loo,"C:/R/Microbiota_IDPs_MR/forward_mr_loo.csv")

##### Single SNP MR ####

forward_mr_single <- mr_singlesnp(forward_harmonized_data, parameters = default_parameters(),
                                       single_method = 'mr_wald_ratio',         
                                       all_method = c('mr_ivw', 'mr_egger_regression'))

write.csv(forward_mr_single,"C:/R/Microbiota_IDPs_MR/forward_mr_single.csv")

##### sign concordance test ####

forward_mr_fit_sign <- mr_sign(b_exp = forward_harmonized_data$beta.exposure,
                                    b_out = forward_harmonized_data$beta.outcome)

##### robust adjusted profile score ####

## We faced error:
forward_mr_fit_raps <- mr.raps(b_exp = forward_harmonized_data$beta.exposure,
                                b_out = forward_harmonized_data$beta.outcome,
                                se_exp = forward_harmonized_data$se.exposure,
                                se_out = forward_harmonized_data$se.outcome,
                                over.dispersion = TRUE,    
                                loss.function = 'huber',  
                                B = 1000,                 
                                diagnosis = TRUE)

## This works:
forward_mr_fit_raps_all <- mr.raps.all(b_exp = forward_harmonized_data$beta.exposure,
                                        b_out = forward_harmonized_data$beta.outcome,
                                        se_exp = forward_harmonized_data$se.exposure,
                                        se_out = forward_harmonized_data$se.outcome)

##### MR-PRESSO (pleiotropy and outliers removal) ####

forward_mr_mrpresso <- run_mr_presso(forward_harmonized_data, 
                                     NbDistribution = 1000, 
                                     SignifThreshold = 0.05)

##### Directionality Test ####

forward_harmonized_data$r.outcome <- get_r_from_pn(forward_harmonized_data$pval.outcome,
                                                   forward_harmonized_data$samplesize.outcome)

forward_harmonized_data$r.exposure <- get_r_from_pn(forward_harmonized_data$pval.exposure,
                                                    forward_harmonized_data$samplesize.exposure)

forward_mr_steiger <- directionality_test(forward_harmonized_data)

write.csv(forward_mr_steiger,"C:/R/Microbiota_IDPs_MR/forward_mr_steiger.csv")

##### Joining Results ####

# Forward-MR-Heterogeneity:
forward_mr_heterogeneity_result <- left_join(forward_mr,forward_heterogeneity)

write.csv(forward_mr_heterogeneity_result, 
          "C:/R/Microbiota_IDPs_MR/forward_mr_heterogeneity_result.csv")

# Forward-MR-Pleiotropy:
forward_pleiotropy_for_join <- forward_pleiotropy
names(forward_pleiotropy_for_join)[6] <- "se_pleiotropy"
names(forward_pleiotropy_for_join)[7] <- "pval_pleiotropy"
forward_mr_pleiotropy_result <- left_join(forward_mr,forward_pleiotropy_for_join)

write.csv(forward_mr_pleiotropy_result, 
          "C:/R/Microbiota_IDPs_MR/forward_mr_pleiotropy_result.csv")

# Forward-MR-LOO:
forward_mr_loo_for_join <- forward_mr_loo
names(forward_mr_loo_for_join)[7] <- "b_LOO"
names(forward_mr_loo_for_join)[8] <- "se_LOO"
forward_mr_loo_result <- left_join(forward_mr,forward_mr_loo_for_join,
                                   relationship = "many-to-many")

write.csv(forward_mr_loo_result, 
         "C:/R/Microbiota_IDPs_MR/forward_mr_loo_result.csv")

# Forward-MR-SingleSNP:
forward_mr_single_join <- forward_mr_single
names(forward_mr_single_join)[7] <- "b_SingleSNP"
names(forward_mr_single_join)[8] <- "se_SingleSNP"
forward_mr_SingleSNP_result <- left_join(forward_mr,forward_mr_single_join,
                                         relationship = "many-to-many")

write.csv(forward_mr_SingleSNP_result, 
          "C:/R/Microbiota_IDPs_MR/forward_mr_SingleSNP_result.csv")

# Forward-MR-Directionality:
forward_mr_directionality_result <- left_join(forward_mr, forward_mr_steiger)
write.csv(forward_mr_directionality_result, 
          "C:/R/Microbiota_IDPs_MR/forward_mr_directionality_result.csv")

##### Plotting ####

# Scatter plot:
forward_scatter_plot <- mr_scatter_plot(forward_mr)

# forest-plot:
forward_forest_plot <- mr_forest_plot(forward_mr_single)

# LOO plot: 
forward_loo_plot <- mr_leaveoneout_plot(forward_mr_loo)


