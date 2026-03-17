library(dplyr)

setwd('~/Desktop/chr15_CHRNA5_bladder_project/')

df <- read.csv('CCLE_CHRNA5_Expression_Public_24Q4_All_Cells.csv')

df <- df[,c(1:6,11)]


df <- df %>% mutate(CHRNA5_TPM_max=max(CHRNA5))
df <- df %>% mutate(CHRNA5_TPM_min=min(CHRNA5))
df <- df %>% mutate(CHRNA5_TPM_median=median(CHRNA5))



