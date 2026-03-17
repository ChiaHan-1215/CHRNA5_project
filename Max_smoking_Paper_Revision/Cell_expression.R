## Date: 03172026
## Goal: The Revision of Max C5 smoking paper. Finding possible datase that have lymphocyte cells (EBV +/-) with CHRNA5 expression 



library(dplyr)

setwd('~/Desktop/chr15_CHRNA5_bladder_project/')

df <- read.csv('CCLE_CHRNA5_Expression_Public_24Q4_All_Cells.csv')

df <- df[,c(1:6,11)]


df <- df %>% mutate(CHRNA5_TPM_max=max(CHRNA5))
df <- df %>% mutate(CHRNA5_TPM_min=min(CHRNA5))
df <- df %>% mutate(CHRNA5_TPM_median=median(CHRNA5))



# loading file from CELL ALTAS

cellinfo <- read.delim('Max_and_Auth_C5_smoking_paper/CCLE_and_ALTAS_table_for_CHRNA5_expression_in_Lyc/rna_cell_lines.tsv')
cellexp <- read.delim('Max_and_Auth_C5_smoking_paper/CCLE_and_ALTAS_table_for_CHRNA5_expression_in_Lyc/rna_celline.tsv')

# CHRNA5 extraction 
# Merge data, and select nTPM (Normalized TPM)
cell_c5 <- cellexp %>% filter(Gene.name == "CHRNA5")
cell_c5 <- cell_c5[,c(2,3,6)]
Mer_cell <- inner_join(cell_c5,cellinfo,by='Cell.line')


# Some simple table
# Cell line group
Per_Cancer <- Mer_cell %>% group_by(Primary.disease) %>% summarise(TPM_mean=mean(nTPM))


# Common cell lines 
Comm_cell <- Mer_cell %>% group_by(Cell.line) %>% summarise(TPM_mean=mean(nTPM))



################################################
################################################
################################################

# MAGE dataset from 1000GP
#https://www.internationalgenome.org/data-portal/data-collection/mage_rnaseq
# https://zenodo.org/records/10535719

# Downloading...




