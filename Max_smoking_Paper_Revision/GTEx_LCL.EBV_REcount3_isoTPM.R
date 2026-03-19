## Goal : lm anaylsis of GTEx bladder data with GT 

### Load VcfR 
library(vcfR)
library(biomaRt)
library(dplyr)
library(tidyr)
library(parallel)
library(doParallel)
library(foreach)
library(RNOmni)
library(tibble)

library(recount3)
library(snapcount)
library(megadepth)
library(dplyr)
library(tidyr)

# load cov, for PC studd

cov <- read.delim('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/GTEx_analysis_v10/QTL/eQTL_covariates/Bladder.v10.covariates.txt')
cov <- t(cov) %>% as.data.frame()
names(cov) <- cov[1,]
cov <- cov[-1,]
cov$GTEx_ID <- rownames(cov) %>% gsub('\\.','-',.)


# laod mainfest file 
mani <- read.delim('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/GTEx_analysis_v10/Phenotype_dbGaP/phs000424.v10.pht002743.v10.p2.c1.GTEx_Sample_Attributes.GRU.txt',skip = 10)

mani_bladder <- mani %>% dplyr::filter(SMTSD == 'Cells - EBV-transformed lymphocytes')
mani_bladder <- mani_bladder[,c(1,2)]
names(mani_bladder)[2] <- "Sample_ID"

########################################################################################################################
# laod Recount3 files Junction COUNT extraction
########################################################################################################################

# GTEx 
setwd('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/project_CHRNA5/Recount3_GTEx_junction_count_table/All_CHRAN5_jt/')


lf <- list.files('.',pattern = '.csv')

# As we load only EBV LCL, it's in the blood section 

blood <- read_csv('BLOOD_GTEx_CHRNA5_jc.csv')
names(blood)[1] <- "Sample_ID"
# Merge to see how many sample 
# REmove NA
Mer <- inner_join(mani_bladder, blood, by = "Sample_ID") 
Mer <- Mer[,-1]

# the target isoform in C5, check with SHSY5Y Jct count, only 1bp move so it's ok

#chr15:78589850-78593091  R153  Super-short
#chr15:78589983-78593091  Q197
#chr15:78590099-78593091  W236
#chr15:78590105-78593091	P238
#chr15:78590507-78593091  R372
#chr15:78590637-78593091  E415  Full-length

## Other Exons, so is +1,-1
 
# chr15:78565826-78580810 Ex1-2
# chr15:78580963-78586644 Ex2-3
# chr15:78586690-78588313 Ex3-4
# chr15:78588424-78589804 Ex4-5

df <- Mer
names(df)[1] <- "full_id"


jxn_set <- c(
  "R153"        = "chr15:78589850-78593091",
  "Q197"        = "chr15:78589983-78593091",
  "W236"        = "chr15:78590099-78593091",
  "P238"        = "chr15:78590105-78593091",
  "R372"        = "chr15:78590507-78593091",
  "E415"        = "chr15:78590637-78593091",
  "Ex1-2"       = "chr15:78565826-78580810",
  "Ex2-3"       = "chr15:78580963-78586644",
  "Ex3-4"       = "chr15:78586690-78588313",
  "Ex4-5"       = "chr15:78588424-78589804"
)

df <- df[, c(1, which(names(df) %in% jxn_set)), drop = FALSE]

idx <- match(names(df)[-1], jxn_set)
names(df)[-1] <- names(jxn_set)[idx]

df <- add_column(df, GTEx_ID=gsub("(GTEX-[^-]*)-.*", "\\1", df$full_id), .before = names(df)[1])

df <-df %>%
  mutate(iso_sum = rowSums(across(7:12)))

for (j in names(df)[7:12]){
  df[[paste0('Perc_',j)]] <- round(df[[j]] / df$iso_sum,3)
}

# after sum, have to filter total sum reads > maybe > 10 or > 5?
#df_isoform <- df %>% filter(iso_sum > 5)

########################################################################################################################
# laod TPM files FOR TPM AND GT EXTRACTION
########################################################################################################################

# load the SMOKE and GT and TPM 
# also have EBV cell "Cells_EBV_transformedlymphocytes"
isoTPM_GT <- read.delim('~/Desktop/chr15_CHRNA5_bladder_project/Finalized_GTEx_GT_CHRAN5_isofrom_TPM.tsv')
isoTPM_GT <- isoTPM_GT[,c(1:12,23,24,names(isoTPM_GT) %>% grep("Cells_EBV_",.,value = F))]
# Total CHRNA5 TPM 
totalTPM <-read.csv('~/Desktop/chr15_CHRNA5_bladder_project/GTEx_CHRNA5_TPM_in_alltissue.csv')
totalTPM <- totalTPM[,c(1,35)]
names(totalTPM)[1] <- "GTEx_ID"
# Merged above data

GTEx_meta <- left_join(isoTPM_GT,totalTPM,by='GTEx_ID')


df <- left_join(df,GTEx_meta,by="GTEx_ID")

df <- df[,c(1,2,20:32,3:19,33:ncol(df))]

df[, 16:37] <- lapply(df[, 16:37], function(x) as.numeric(x))

