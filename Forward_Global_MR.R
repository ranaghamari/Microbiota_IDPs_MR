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
library(gtable)
library(grid)
library(gridExtra)
library(ggpubr)

##### Preparing exposure file ####

Microbiota_H2_LargerThanZero <- read_csv("Microbiota_H2_LargerThanZero.csv")

microbiota_h2_id <- Microbiota_H2_LargerThanZero$id

microbiota_exp_dat <- "Microbiome_exposure_-5_h2.csv"

##### Preparing outcome file ####

global_idps_files = list.files(path = "C:/R/MR_Microbiota_GlobalIDPs/global_IDPs",
                               pattern = ".csv")

idp_outcome_dat <- list()

for (i in 1:length(global_idps_files)) {
  idp_outcome_dat[[i]] <- read_outcome_data(
    snps = microbiota_exp_dat$SNP,
    filename = global_idps_csv[i],
    sep = ",",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "P",
    chr_col = "CHR"
  )
  
  idp_outcome_dat[[i]][,"outcome"] <- toupper(strsplit(global_idps_files[i],"_")[[1]][2])
}

idp_outcome_data <- data.frame()

for (j in 1:length(idp_outcome_dat)) {
  idp_outcome_data <- rbind(idp_outcome_data, idp_outcome_dat[[j]])
}

##### Harmonization ####

forward_harmonized_data <- harmonise_data(
  exposure_dat = microbiota_exp_dat,
  outcome_dat = idp_outcome_data,
  action = 2
)

##### ---- MR analysis ---- ####

## leave-one-out MR

forward_mr_loo <- TwoSampleMR::mr_leaveoneout(forward_harmonized_data, 
                                      parameters = default_parameters(),  
                                      method = TwoSampleMR::mr_ivw)

## single SNP MR

forward_mr_single <- TwoSampleMR::mr_singlesnp(forward_harmonized_data, 
                                               parameters = default_parameters(),
                                               single_method = 'mr_wald_ratio',         
                                               all_method = c('mr_ivw', 'mr_egger_regression'))

## heterogeneity tests

forward_mr_het <- TwoSampleMR::mr_heterogeneity(forward_harmonized_data)

## MR Egger (horizontal pleiotropy test)

forward_mr_egger <- TwoSampleMR::mr_pleiotropy_test(forward_harmonized_data)

## sign concordance test

forward_mr_fit_sign <- TwoSampleMR::mr_sign(b_exp = forward_harmonized_data$beta.exposure,
                                    b_out = forward_harmonized_data$beta.outcome)

## robust adjusted profile score

mr_fit_raps_all <- mr.raps::mr.raps.all(b_exp = forward_harmonized_data$beta.exposure,
                                        b_out = forward_harmonized_data$beta.outcome,
                                        se_exp = forward_harmonized_data$se.exposure,
                                        se_out = forward_harmonized_data$se.outcome)

## two sample MR 

forward_mr <- mr(forward_harmonized_data)

write.csv(forward_mr, "C:/R/MR_Microbiota_GlobalIDPs/forward_mr.csv")

## MR-PRESSO (pleiotropy and outliers removal)

mr_mrpresso <- MRPRESSO:::mr_presso(BetaOutcome = 'beta.outcome',
                                    BetaExposure = 'beta.exposure',
                                    SdOutcome = 'se.outcome',
                                    SdExposure = 'se.exposure',
                                    data = forward_harmonized_data,             
                                    OUTLIERtest = TRUE,
                                    DISTORTIONtest = TRUE,
                                    SignifThreshold = 0.05,
                                    NbDistribution = 5e4,
                                    seed = NULL) 

##### Directionality Test ####

forward_harmonized_data$r.outcome <- get_r_from_pn(forward_harmonized_data$pval.outcome,
                                                   rep(36663,33767))

forward_harmonized_data$r.exposure <- get_r_from_pn(forward_harmonized_data$pval.exposure,
                                                    forward_harmonized_data$samplesize.exposure)

forward_mr_steiger <- directionality_test(forward_harmonized_data)

##### Manhattan plots #### 

forward_mr_phylum <- forward_mr[grepl("phylum",
                                            forward_mr$exposure),]
forward_mr_class <- forward_mr[grepl("class",
                                           forward_mr$exposure),]
forward_mr_order <- forward_mr[grepl("order",
                                           forward_mr$exposure),]
forward_mr_family <- forward_mr[grepl("family",
                                            forward_mr$exposure),]
forward_mr_genus <- forward_mr[grepl("genus",
                                           forward_mr$exposure),]

forward_mr_plot <- list(forward_mr_phylum, forward_mr_class, forward_mr_family,
                        forward_mr_order, forward_mr_genus)

exposure_outcome_pair <- c()
for (m in 1:length(forward_mr_plot)) {
  for (n in 1:nrow(forward_mr_plot[[m]])) {
    exposure_outcome_pair[n] <- paste(forward_mr_plot[[m]]$exposure[n], 
                                      forward_mr_plot[[m]]$outcome[n], sep = " * ")
  }
  forward_mr_plot[[m]] <- cbind(forward_mr_plot[[m]], exposure_outcome_pair)
  exposure_outcome_pair <- c()
}  

significane <- c()
for (k in 1:length(forward_mr_plot)) {
  for (l in 1:nrow(forward_mr_plot[[k]])) {
    significane[l] <- -log10(forward_mr_plot[[k]]$pval[l])
  }
  forward_mr_plot[[k]] <- cbind(forward_mr_plot[[k]], significane)
  significane <- c()
}


forward_plot_phylum <- 
  forward_mr_plot[[1]] %>%
  arrange(outcome) %>%
  mutate(exposure_outcome_pair = as_factor(exposure_outcome_pair))  %>%
  ggplot(aes(x = exposure_outcome_pair, y = significane)) +
  geom_point(aes(colour = outcome)) +
  labs(x = "Phyla",y = "-log (p-value)", colour = "Global IDPs") + 
  theme(axis.text.x = element_blank(), legend.title = element_text(size=15), 
        legend.text=element_text(size=10), axis.title.y = element_text(size = 8)) + 
  scale_color_manual(values = c("aquamarine2",
                                "aquamarine4",
                                "blue",
                                "#1b98e0",
                                "brown3",
                                "coral",
                                "cyan2",
                                "darkgoldenrod1",
                                "darkolivegreen3",
                                "darkred",
                                "gold",
                                "purple",
                                "darkcyan",
                                "maroon"))

forward_plot_class <- 
  forward_mr_plot[[2]] %>%
  arrange(outcome) %>%
  mutate(exposure_outcome_pair = as_factor(exposure_outcome_pair))  %>%
  ggplot(aes(x = exposure_outcome_pair, y = significane)) +
  geom_point(aes(colour = outcome)) +
  labs(x = "Class",y = "-log (p-value)", colour = "Global IDPs") + 
  theme(axis.text.x = element_blank(), legend.title = element_text(size=15), 
        legend.text=element_text(size=10), axis.title.y = element_text(size = 8)) + 
  scale_color_manual(values = c("aquamarine2",
                                "aquamarine4",
                                "blue",
                                "#1b98e0",
                                "brown3",
                                "coral",
                                "cyan2",
                                "darkgoldenrod1",
                                "darkolivegreen3",
                                "darkred",
                                "gold",
                                "purple",
                                "darkcyan",
                                "maroon"))

forward_plot_family <- 
  forward_mr_plot[[3]] %>%
  arrange(outcome) %>%
  mutate(exposure_outcome_pair = as_factor(exposure_outcome_pair))  %>%
  ggplot(aes(x = exposure_outcome_pair, y = significane)) +
  geom_point(aes(colour = outcome)) +
  labs(x = "Family",y = "-log (p-value)", colour = "Global IDPs") + 
  theme(axis.text.x = element_blank(), legend.title = element_text(size=15), 
        legend.text=element_text(size=10), axis.title.y = element_text(size = 8)) + 
  scale_color_manual(values = c("aquamarine2",
                                "aquamarine4",
                                "blue",
                                "#1b98e0",
                                "brown3",
                                "coral",
                                "cyan2",
                                "darkgoldenrod1",
                                "darkolivegreen3",
                                "darkred",
                                "gold",
                                "purple",
                                "darkcyan",
                                "maroon"))

forward_plot_order <- 
  forward_mr_plot[[4]] %>%
  arrange(outcome) %>%
  mutate(exposure_outcome_pair = as_factor(exposure_outcome_pair))  %>%
  ggplot(aes(x = exposure_outcome_pair, y = significane)) +
  geom_point(aes(colour = outcome)) +
  labs(x = "Order",y = "-log (p-value)", colour = "Global IDPs") + 
  theme(axis.text.x = element_blank(), legend.title = element_text(size=15), 
        legend.text=element_text(size=10), axis.title.y = element_text(size = 8)) + 
  scale_color_manual(values = c("aquamarine2",
                                "aquamarine4",
                                "blue",
                                "#1b98e0",
                                "brown3",
                                "coral",
                                "cyan2",
                                "darkgoldenrod1",
                                "darkolivegreen3",
                                "darkred",
                                "gold",
                                "purple",
                                "darkcyan",
                                "maroon"))

forward_plot_genus <- 
  forward_mr_plot[[5]] %>%
  arrange(outcome) %>%
  mutate(exposure_outcome_pair = as_factor(exposure_outcome_pair))  %>%
  ggplot(aes(x = exposure_outcome_pair, y = significane)) +
  geom_point(aes(colour = outcome)) +
  labs(x = "Genra",y = "-log (p-value)", colour = "Global IDPs") + 
  theme(axis.text.x = element_blank(), legend.title = element_text(size=15), 
        legend.text=element_text(size=10), axis.title.y = element_text(size = 8)) + 
  scale_color_manual(values = c("aquamarine2",
                                "aquamarine4",
                                "blue",
                                "#1b98e0",
                                "brown3",
                                "coral",
                                "cyan2",
                                "darkgoldenrod1",
                                "darkolivegreen3",
                                "darkred",
                                "gold",
                                "purple",
                                "darkcyan",
                                "maroon"))


final_forward_plot <- ggarrange(forward_plot_phylum, forward_plot_class, 
                                forward_plot_order, forward_plot_family,
                                forward_plot_genus, ncol = 2, nrow = 3, 
                                common.legend = TRUE, legend = "bottom")

final_forward_plot <- grid.arrange(final_forward_plot,nrow=1,
                                   top=text_grob("Forward MR \n Gut microbiota abundance -> Global Imaging-derived phenotypes"))
