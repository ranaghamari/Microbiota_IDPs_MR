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

idps_exp_files <- list.files("C:/R/MR_Microbiota_GlobalIDPs" ,  pattern = ".csv")

idps_exp_dat <- list()

library(gsubfn)

for (i in 1:length(idps_exp_files)) {
  idps_exp_dat[[i]] <- read_exposure_data(
    idps_exp_files[i],
    sep = ",",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "P",
    chr_col = "CHR"
  )
  idps_exp_dat[[i]] <- filter(idps_exp_dat[[i]],
                              idps_exp_dat[[i]]$pval.exposure<=5*(10^-8))
  
  idps_exp_dat[[i]][,"outcome"] <- toupper(strsplit(idps_exp_dat[i],"_")[[1]][2])
}

idps_exp_data <- data.frame()

for (j in 1:length(idps_exp_dat)) {
  idps_exp_data <- rbind(idps_exp_data, idps_exp_data[[j]])
}

idps_exp_clumped <- ld_clump(
  dplyr::tibble(rsid=idps_exp_data$SNP, 
                pval=idps_exp_data$pval.exposure, 
                id=idps_exp_data$exposure),
  plink_bin = "C:/R/Microbiota_IDPs_MR/plink_win64_20231211 (1)/plink.exe",
  bfile = "C:/R/Microbiota_IDPs_MR/1kg.v3/EUR")

global_idps_exp_dat <- inner_join(idps_exp_clumped, idps_exp_data)

##### Preparing outcome file ####

# The data is available

##### Harmonization ####

reverse_harmonized_data <- harmonise_data(
  exposure_dat = global_idps_exp_dat, 
  outcome_dat = microbiota_out_dat, action = 2)

##### ---- MR analysis ---- ####

## leave-one-out MR

reverse_mr_loo <- TwoSampleMR::mr_leaveoneout(reverse_harmonized_data, 
                                              parameters = default_parameters(),  
                                              method = TwoSampleMR::mr_ivw)

## single SNP MR

reverse_mr_single <- TwoSampleMR::mr_singlesnp(reverse_harmonized_data, 
                                               parameters = default_parameters(),
                                               single_method = 'mr_wald_ratio',         
                                               all_method = c('mr_ivw', 'mr_egger_regression'))

## heterogeneity tests

reverse_mr_het <- TwoSampleMR::mr_heterogeneity(reverse_harmonized_data)

## MR Egger (horizontal pleiotropy test)

reverse_mr_egger <- TwoSampleMR::mr_pleiotropy_test(reverse_harmonized_data)

## sign concordance test

reverse_mr_fit_sign <- TwoSampleMR::mr_sign(b_exp = reverse_harmonized_data$beta.exposure,
                                            b_out = reverse_harmonized_data$beta.outcome)

## robust adjusted profile score

mr_fit_raps_all <- mr.raps::mr.raps.all(b_exp = forward_harmonized_data$beta.exposure,
                                        b_out = forward_harmonized_data$beta.outcome,
                                        se_exp = forward_harmonized_data$se.exposure,
                                        se_out = forward_harmonized_data$se.outcome)

## two sample MR 

reverse_mr <- mr(reverse_harmonized_data)

write.csv(reverse_mr, "C:/R/MR_Microbiota_GlobalIDPs/reverse_mr.csv")

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

reverse_harmonized_data$r.outcome <- get_r_from_pn(reverse_harmonized_data$pval.outcome,
                                                   reverse_harmonized_data$samplesize.outcome)

reverse_harmonized_data$r.exposure <- get_r_from_pn(reverse_harmonized_data$pval.exposure,
                                                    rep(36663,3675))

reverse_mr_steiger <- directionality_test(reverse_harmonized_data)

##### Manhattan plots ####

reverse_mr_phylum <- reverse_mr[grepl("phylum", reverse_mr$outcome),]

reverse_mr_class <- reverse_mr[grepl("class", reverse_mr$outcome),]
                                         
reverse_mr_order <- reverse_mr[grepl("order", reverse_mr$outcome),]
                                          
reverse_mr_family <- reverse_mr[grepl("family", reverse_mr$outcome),]
                                           
reverse_mr_genus <- reverse_mr[grepl("genus",  reverse_mr$outcome),]
                                         
reverse_mr_plot <- list(reverse_mr_phylum, reverse_mr_class, reverse_mr_family,
                        reverse_mr_order, reverse_mr_genus)

exposure_outcome_pair_reverse <- c()
for (m in 1:length(reverse_mr_plot)) {
  for (n in 1:nrow(reverse_mr_plot[[m]])) {
    exposure_outcome_pair_reverse[n] <- paste(reverse_mr_plot[[m]]$exposure[n], 
                                              reverse_mr_plot[[m]]$outcome[n], sep = " * ")
  }
  reverse_mr_plot[[m]] <- cbind(reverse_mr_plot[[m]], exposure_outcome_pair_reverse)
  exposure_outcome_pair_reverse <- c()
}

significane_reverse <- c()
for (k in 1:length(reverse_mr_plot)) {
  for (l in 1:nrow(reverse_mr_plot[[k]])) {
    significane_reverse[l] <- -log10(reverse_mr_plot[[k]]$pval[l])
  }
  reverse_mr_plot[[k]] <- cbind(reverse_mr_plot[[k]], significane_reverse)
  significane_reverse <- c()
}

reverse_plot_phylum <- 
  reverse_mr_plot[[1]] %>%
  arrange(exposure) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = exposure)) +
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
                                "maroon"))

reverse_plot_class <- 
  reverse_mr_plot[[2]] %>%
  arrange(exposure) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = exposure)) +
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
                                "maroon"))

reverse_plot_family <- 
  reverse_mr_plot[[3]] %>%
  arrange(exposure) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = exposure)) +
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
                                "maroon"))

reverse_plot_order <- 
  reverse_mr_plot[[4]] %>%
  arrange(exposure) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = exposure)) +
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
                                "maroon"))

reverse_plot_genus <- 
  reverse_mr_plot[[5]] %>%
  arrange(exposure) %>%
  mutate(exposure_outcome_pair_reverse = as_factor(exposure_outcome_pair_reverse))  %>%
  ggplot(aes(x = exposure_outcome_pair_reverse, y = significane_reverse)) +
  geom_point(aes(colour = exposure)) +
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
                                "maroon"))

library(gtable)
library(grid)
library(ggpubr)

final_reverse_plot <- ggarrange(reverse_plot_phylum, reverse_plot_class, 
                                reverse_plot_order, reverse_plot_family,
                                reverse_plot_genus, ncol = 2, nrow = 3, 
                                common.legend = TRUE, legend = "bottom")

final_reverse_plot <- grid.arrange(final_reverse_plot,nrow=1,
                                   top=text_grob("Reverse MR \n Global Imaging-derived phenotypes -> Gut microbiota abundance"))



