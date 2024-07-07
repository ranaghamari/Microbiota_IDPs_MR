
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
library(ggplot2)

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
##### Manhattan plots ####

# We categorize exposures based on IDPs:

reverse_mr_plot <- reverse_mr
idp_category_reverse = c()
idp_category_number_reverse = c()

for (i in 1:nrow(reverse_mr_plot)) {
  if (grepl("QC_", reverse_mr_plot$exposure[i]) &
      grepl("_SWI_", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "QC"
    idp_category_number_reverse[i] <- 1
  } 
  if (grepl("QC_", reverse_mr_plot$exposure[i]) &
      grepl("_FLAIR_", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "QC"
    idp_category_number_reverse[i] <- 1
  } 
  if (grepl("QC_", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "QC"
    idp_category_number_reverse[i] <- 1
  }
  if (grepl("AmygNuclei", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "AmygNuclei"
    idp_category_number_reverse[i] <- 2
  }
  if (grepl("aparc-", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "aparc"
    idp_category_number_reverse[i] <- 3
  }
  if (grepl("aseg_", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "aseg"
    idp_category_number_reverse[i] <- 4
  }
  if (grepl("BA-exvivo", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "BA"
    idp_category_number_reverse[i] <- 5
  }
  if (grepl("Brainstem_", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "Brain Stem"
    idp_category_number_reverse[i] <- 6
  }
  if (grepl("HippSubfield_", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "HippSubfield"
    idp_category_number_reverse[i] <- 7
  }
  if (grepl("_FA_", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "dMRI FA"
    idp_category_number_reverse[i] <- 8
  }
  if (grepl("_L1_", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "dMRI L1"
    idp_category_number_reverse[i] <- 9
  }
  if (grepl("_L2_", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "dMRI L2"
    idp_category_number_reverse[i] <- 10
  }
  if (grepl("_L3_", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "dMRI L3"
    idp_category_number_reverse[i] <- 11
  }
  if (grepl("_MD_", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "dMRI MD"
    idp_category_number_reverse[i] <- 12
  }
  if (grepl("_MO_", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "dMRI MO"
    idp_category_number_reverse[i] <- 13
  }
  if (grepl("_OD_", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "dMRI OD"
    idp_category_number_reverse[i] <- 14
  }
  if (grepl("_ICVF_", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "dMRI ICVF"
    idp_category_number_reverse[i] <- 15
  }
  if (grepl("_ISOVF_", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "dMRI ISOVF"
    idp_category_number_reverse[i] <- 16
  }
  if (grepl("_SWI_", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "SWI"
    idp_category_number_reverse[i] <- 17
  }
  if (grepl("_FAST_", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "cortical, subcortical volume (FAST)"
    idp_category_number_reverse[i] <- 18
  }
  if (grepl("_FIRST_", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "subcortical volume (FIRST)"
    idp_category_number_reverse[i] <- 19
  }
  if (grepl("_SIENAX_", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "SIENAX"
    idp_category_number_reverse[i] <- 20
  }
  if (grepl("ThalamNuclei_", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "ThalamNuclei"
    idp_category_number_reverse[i] <- 21
  }
  if (grepl("wg_", reverse_mr_plot$exposure[i])) {
    idp_category_reverse[i] <- "wg"
    idp_category_number_reverse[i] <- 22
  }
}


exposure_outcome_pair_reverse <- c()
for (i in 1:nrow(reverse_mr_plot)) {
  exposure_outcome_pair_reverse[i] <- paste(reverse_mr_plot$exposure[i],
                                            reverse_mr_plot$outcome[i], sep = " * ")
}

significane_reverse <- c()
for (i in 1:nrow(reverse_mr_plot)) {
  significane_reverse[i] <- -log10(reverse_mr_plot$pval[i])
}

microbiota_category_reverse <- c()
for (i in 1:nrow(reverse_mr_plot)) {
  if (grepl("phylum",reverse_mr_plot$outcome[i])){
    microbiota_category_reverse[i] <- "phylum"
  }
  if (grepl("class",reverse_mr_plot$outcome[i])) {
    microbiota_category_reverse[i] <- "class"
  }
  if (grepl("order",reverse_mr_plot$outcome[i])) {
    microbiota_category_reverse[i] <- "order"
  }
  if (grepl("family",reverse_mr_plot$outcome[i])) {
    microbiota_category_reverse[i] <- "family"
  }
  if (grepl("genus",reverse_mr_plot$outcome[i])) {
    microbiota_category_reverse[i] <- "genus"
  }
}

reverse_mr_plot <- cbind(reverse_mr_plot, idp_category_reverse, 
                         idp_category_number_reverse, exposure_outcome_pair_reverse,
                         significane_reverse,,microbiota_category_reverse)


reverse_plot_qc <- 
  filter(reverse_mr_plot, reverse_mr_plot$idp_category_reverse=="QC") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/1060)) +
  labs(x = "Taxa -QC",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                        legend.title = element_text( size=5), 
                                        legend.text=element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

reverse_plot_AmygNuclei <- 
  filter(reverse_mr_plot, reverse_mr_plot$idp_category_reverse=="AmygNuclei") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/3144)) +
  labs(x = "Taxa - AmygNuclei",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                legend.title = element_text( size=5), 
                                legend.text=element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

reverse_plot_aparc <- 
  filter(reverse_mr_plot, reverse_mr_plot$idp_category_reverse=="aparc") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/97639)) +
  labs(x = "Taxa - aparc",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                legend.title = element_text( size=5), 
                                legend.text=element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

reverse_plot_aseg <- 
  filter(reverse_mr_plot, reverse_mr_plot$idp_category_reverse=="aseg") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/14277)) +
  labs(x = "Taxa - aseg",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                legend.title = element_text( size=5), 
                                legend.text=element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

reverse_plot_BA <- 
  filter(reverse_mr_plot, reverse_mr_plot$idp_category_reverse=="BA") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/10569)) +
  labs(x = "Taxa - BA",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                legend.title = element_text( size=5), 
                                legend.text=element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

reverse_plot_brainstem <- 
  filter(reverse_mr_plot, reverse_mr_plot$idp_category_reverse=="Brain Stem") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/894)) +
  labs(x = "Taxa - Brain Stem",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                legend.title = element_text( size=5), 
                                legend.text=element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

reverse_plot_HippSubfield <- 
  filter(reverse_mr_plot, reverse_mr_plot$idp_category_reverse=="HippSubfield") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/7199)) +
  labs(x = "Taxa - HippSubfield",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                legend.title = element_text( size=5), 
                                legend.text=element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

reverse_plot_FA <- 
  filter(reverse_mr_plot, reverse_mr_plot$idp_category_reverse=="dMRI FA") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/12424)) +
  labs(x = "Taxa - dMRI FA",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                legend.title = element_text( size=5), 
                                legend.text=element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

reverse_plot_l1 <- 
  filter(reverse_mr_plot, reverse_mr_plot$idp_category_reverse=="dMRI L1") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/11884)) +
  labs(x = "Taxa - dMRI L1",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                legend.title = element_text( size=5), 
                                legend.text=element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

reverse_plot_l2 <- 
  filter(reverse_mr_plot, reverse_mr_plot$idp_category_reverse=="dMRI L2") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/12254)) +
  labs(x = "Taxa - dMRI L2",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                legend.title = element_text( size=5), 
                                legend.text=element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

reverse_plot_l3 <- 
  filter(reverse_mr_plot, reverse_mr_plot$idp_category_reverse=="dMRI L3") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/11899)) +
  labs(x = "Taxa - dMRI L3",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                legend.title = element_text( size=5), 
                                legend.text=element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

reverse_plot_md <- 
  filter(reverse_mr_plot, reverse_mr_plot$idp_category_reverse=="dMRI MD") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/11810)) +
  labs(x = "Taxa - dMRI MD",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                legend.title = element_text( size=5), 
                                legend.text=element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

reverse_plot_mo <- 
  filter(reverse_mr_plot, reverse_mr_plot$idp_category_reverse=="dMRI MO") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/8991)) +
  labs(x = "Taxa - dMRI MO",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                legend.title = element_text( size=5), 
                                legend.text=element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

reverse_plot_od <- 
  filter(reverse_mr_plot, reverse_mr_plot$idp_category_reverse=="dMRI OD") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/10590)) +
  labs(x = "Taxa - dMRI OD",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                legend.title = element_text( size=5), 
                                legend.text=element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

reverse_plot_ICVF <- 
  filter(reverse_mr_plot, reverse_mr_plot$idp_category_reverse=="dMRI ICVF") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/13071)) +
  labs(x = "Taxa - dMRI ICVF",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                legend.title = element_text( size=5), 
                                legend.text=element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

reverse_plot_ISOVF <- 
  filter(reverse_mr_plot, reverse_mr_plot$idp_category_reverse=="dMRI ISOVF") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/9462)) +
  labs(x = "Taxa - dMRI ISOV",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                legend.title = element_text( size=5), 
                                legend.text=element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

reverse_plot_SWI <- 
  filter(reverse_mr_plot, reverse_mr_plot$idp_category_reverse=="SWI") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/1529)) +
  labs(x = "Taxa - SWI",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                legend.title = element_text( size=5), 
                                legend.text=element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

reverse_plot_SIENAX <- 
  filter(reverse_mr_plot, reverse_mr_plot$idp_category_reverse=="SIENAX") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/1760)) +
  labs(x = "Taxa - SIENAX",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                legend.title = element_text( size=5), 
                                legend.text=element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

reverse_plot_fast <- 
  filter(reverse_mr_plot, 
         reverse_mr_plot$idp_category_reverse=="cortical, subcortical volume (FAST)") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/17503)) +
  labs(x = "Taxa - cortical, subcortical volume (FAST)",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                legend.title = element_text( size=5), 
                                legend.text=element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

reverse_plot_first <- 
  filter(reverse_mr_plot, 
         reverse_mr_plot$idp_category_reverse=="subcortical volume (FIRST)") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/2345)) +
  labs(x = "Taxa - subcortical volume (FIRST)",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                legend.title = element_text( size=5), 
                                legend.text=element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

reverse_plot_wg <- 
  filter(reverse_mr_plot, 
         reverse_mr_plot$idp_category_reverse=="wg") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/12324)) +
  labs(x = "Taxa - wg",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                legend.title = element_text( size=5), 
                                legend.text=element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

reverse_plot_ThalamNuclei <- 
  filter(reverse_mr_plot, 
         reverse_mr_plot$idp_category_reverse=="ThalamNuclei") %>%
  arrange(microbiota_category_reverse) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = microbiota_category_reverse)) +
  geom_hline(yintercept = -log10(0.05/8303)) +
  labs(x = "Taxa - ThalamNuclei",
       y = "-log(p-value)", 
       colour = "Taxa") + theme(axis.text.x = element_blank(),
                                legend.title = element_text( size=5), 
                                legend.text=element_text(size=5),
                                axis.title.y = element_text(size=5)) + 
  scale_color_manual(values = c("aquamarine2",
                                "blue",
                                "#1b98e0",
                                "blueviolet",
                                "brown1"))

library(gtable)
library(grid)
library(ggpubr)

final_reverse_plot <- ggarrange(reverse_plot_AmygNuclei, reverse_plot_aparc,
                                reverse_plot_aseg, reverse_plot_BA,
                                reverse_plot_brainstem, reverse_plot_FA, 
                                reverse_plot_fast, reverse_plot_first, 
                                reverse_plot_HippSubfield, reverse_plot_ICVF,
                                reverse_plot_ISOVF, reverse_plot_l1,
                                reverse_plot_l2, reverse_plot_l3, 
                                reverse_plot_md, reverse_plot_mo, 
                                reverse_plot_od, reverse_plot_qc, 
                                reverse_plot_SIENAX, reverse_plot_SWI,
                                reverse_plot_ThalamNuclei, reverse_plot_wg,
                                ncol = 3, nrow = 8, 
                                common.legend = TRUE,legend = "bottom")
