
insinstall.packages("dplyr")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("vegan")
install.packages("devtools")
install.packages("readr")
install.packages("remotes")
install.packages("plyr")

devtools::install_github("gqi/MRMix")
library(remotes)
library(devtools)
remotes::install_github("mrcieu/ieugwasr")
#install_github("MRCIEU/TwoSampleMR")
remotes::install_github("MRCIEU/TwoSampleMR")

install.packages('C:/R/MRCIEU-MRInstruments-0.3.3-16-gefa2ca0.tar.gz', 
                 repos=NULL, type='source')
install.packages("C:/Program Files/R/R-4.3.2/library/MRCIEU-genetics.binaRies-b0324f1/MRCIEU-genetics.binaRies-b0324f1",repos=NULL, type='source')

devtools::install_github("rondolab/MR-PRESSO", force = TRUE)
library(tidyverse)
library(vegan)
library(ggplot2)
