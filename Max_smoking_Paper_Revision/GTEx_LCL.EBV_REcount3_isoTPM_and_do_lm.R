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

mani_blood <- mani %>% dplyr::filter(SMTSD == 'Cells - EBV-transformed lymphocytes')
mani_blood <- mani_blood[,c(1,2)]
names(mani_blood)[2] <- "Sample_ID"

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
Mer <- inner_join(mani_blood, blood, by = "Sample_ID") 
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
gtex_info <- read.delim('~/Desktop/chr15_CHRNA5_bladder_project/Finalized_GTEx_GT_CHRAN5_isofrom_TPM.tsv')
gtex_info <- gtex_info[,c(1:12,23,24)]


#load CHRNA5 TPMs in GTExv10, correct 
GTExv10.C5.TPM <- read.csv('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/project_CHRNA5/masterFile_GTEx_v10_gene_tpm_CHRNA5_project.csv')
GTExv10.C5.TPM <- GTExv10.C5.TPM %>% select(GTEx_ID,CHRNA5_geneTPM_cells_ebv.transformed_lymphocytes)


######### FIX ##############
######### FIX ##############
######### FIX ##############
# 
# df <- read.delim('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/GTEx_analysis_v10/RNA_Seq/subset_C5.C3.GTExv10.tpm.renamed.txt')
# df <- df %>% add_column(GTEx_ID=paste0(df$gene_id,"_",df$transcript_id,"_TPM"),.before = 1)
# df <- df[,c(1,4:ncol(df))]
# 
# df_trans <- as.data.frame(t(df[, -1]))
# colnames(df_trans) <- df$GTEx_ID
# df_trans$Sample_ID <- rownames(df_trans)
# rownames(df_trans) <- NULL
# df_trans <- df_trans[, c(ncol(df_trans), 1:(ncol(df_trans)-1))]
# df_trans$Sample_ID <- gsub('\\.','-',df_trans$Sample_ID)
# df_trans <- df_trans %>% add_column(GTEx_ID=sub("^([^-]+-[^-]+).*$", "\\1", df_trans$Sample_ID),.after = 1)
# 
# 
# # laod mainfest file 
# mani <- read.delim('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/GTEx_analysis_v10/Phenotype_dbGaP/phs000424.v10.pht002743.v10.p2.c1.GTEx_Sample_Attributes.GRU.txt',skip = 10)
# mani_sub <- mani %>% select(SAMPID,SMTSD)
# mani_sub$SMTSD <- gsub("[[:space:]()/-]+", "_", mani_sub$SMTSD)
# mani_sub$SMTSD <- gsub("_$", "", gsub("[[:space:]()/-]+", "_", mani_sub$SMTSD))
# mani_sub$SMTSD <- tolower(mani_sub$SMTSD)
# names(mani_sub)[1] <- "Sample_ID"
# df_trans <- left_join(df_trans,mani_sub,by="Sample_ID")
# 
# 
# library(dplyr)
# library(tidyr)
# 
# df_wide <- df_trans %>%
#   # 1. Remove Sample_ID because it's too specific for a one-row-per-GTEx_ID summary
#   select(-Sample_ID) %>%
#   # 2. Pivot the data wider
#   pivot_wider(
#     id_cols = GTEx_ID, 
#     names_from = SMTSD, 
#     values_from = starts_with("CHRNA"),
#     names_glue = "{.value}_{SMTSD}"
#   )
# 

# now save 

#write.csv(df_wide,'/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/project_CHRNA5/masterFile_GTEx_v10_isoform_tpm_CHRNA5_project_V2.csv',row.names = F,quote = F)

######### FIXed ##############


## The table needs correct as CHRNA3/5 are assign incorrect 
GTExv10.C5.iso.TPM <- read.csv('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/project_CHRNA5/masterFile_GTEx_v10_isoform_tpm_CHRNA5_project.csv')
GTExv10.C5.iso.TPM <- GTExv10.C5.iso.TPM %>% select(GTEx_ID,grep('CHRNA5_.*_cells_ebv.transformed_lymphocytes',names(GTExv10.C5.iso.TPM),value = T))

# Merge 
gc.Merged <- left_join(GTExv10.C5.TPM,GTExv10.C5.iso.TPM,by="GTEx_ID")

GTEx_meta <- left_join(gtex_info,gc.Merged,by='GTEx_ID')

df <- left_join(df,GTEx_meta,by="GTEx_ID")

df <- df[,c(1,2,20:32,3:19,33:ncol(df))]

df[, 16:37] <- lapply(df[, 16:37], function(x) as.numeric(x))


### Do lm()

# First chkec df factor or not

df$SEX <- factor(df$SEX)
df$RACE <- factor(df$RACE)
df$Smoke_CURRENT_FORMER_NEVER <- factor(df$Smoke_CURRENT_FORMER_NEVER)
str(df)

df.out <- data.frame()
df_count.summary <- data.frame()


inputdf <- df

for (i in grep("_add",names(inputdf),value = T)){
    # i <- "rs16969968_add"
    # 
    snp_id <- gsub("_add|_dom", "", i)   # "rs62286992"
    
    #Strand_Status <- m$Strand_Status
    
    REF <- "G"
    ALT <- "A"
    
    
    for(j in names(inputdf) %>% grep("Perc|CHRNA5",.,value = T)){
      
      # j <- "CHRNA5_geneTPM_cells_ebv.transformed_lymphocytes"
      fmla <- as.formula(paste0(j, "~" , i))
      fmla_adj <- as.formula(paste0(j, "~" , i, " + SEX + AGE + RACE + Smoke_CURRENT_FORMER_NEVER"))
      fmla_int <- as.formula(paste(j, "~", i, " * Smoke_CURRENT_FORMER_NEVER  + SEX + AGE + RACE"))

      
      # Function to extract variable names from a formula
      
      get_vars_from_formula <- function(formula) {
        all.vars(formula)
      }
      
      # List of formulas
      
      formula_list <- list(fmla,fmla_adj,fmla_int)
      
      # Iterate through the list of formulas
      
      for (current_formula in formula_list) {
        # current_formula <- formula_list[[1]]
        # Extract variable names from the formula
        vars <- get_vars_from_formula(current_formula)
        
        # Filter the tar data frame based on the variables in the formula and filter out NA or empty data to get the N_sample
        
        #tar <- inputdf %>% filter(inputdf[[i]] %in% c(0, 1, 2))
        tar <- inputdf %>%
          filter_at(vars, all_vars(!is.na(.) & . != ""))
        
        
        tmp.snp <- tar
        tmp.snp.0a <- tmp.snp[ which(tmp.snp[,i] == 0), ]
        tmp.snp.0b <- as.character(unique(tmp.snp.0a[gsub("_add", "", i)]))
        tmp.snp.1a <- tmp.snp[ which(tmp.snp[,i] == 1), ]
        tmp.snp.1b <- as.character(unique(tmp.snp.1a[gsub("_add", "", i)]))
        tmp.snp.2a <- tmp.snp[ which(tmp.snp[,i] == 2), ]
        tmp.snp.2b <- as.character(unique(tmp.snp.2a[gsub("_add", "", i)]))
        
        out <- tryCatch(lm(current_formula,tar), error=function(e){
          print(paste0('error of tissue ',k))
          return(NA)})
        
        sk_int <- tryCatch(
          grep(":", rownames(coef(summary(out))), value = TRUE),
          error = function(e) {
            print(paste0('Error for tissue ', k))
            return(character(0))
          })
        
        
        # now set up the summary table
        
        tmp.df <- data.frame(
          snp = snp_id,
          #Strand_Status = Strand_Status,
          variable = j,
          model = paste0(formula(current_formula)[2],"~",formula(current_formula)[3]),
          N_sample_not_NA = sum(!is.na(tar[[j]])),
          
          geno_0 = tmp.snp.0b,
          geno_1 = tmp.snp.1b,
          geno_2 = tmp.snp.2b,
          n_0 = nrow(tmp.snp.0a),
          n_1 = nrow(tmp.snp.1a),
          n_2 = nrow(tmp.snp.2a),
          
          ## add TPM ###
          
          exp_mean = round(mean(tar[[j]]),3),
          exp_geno_0 = round(mean(tmp.snp.0a[[j]]),3),
          exp_geno_1 = round(mean(tmp.snp.1a[[j]]),3),
          exp_geno_2 = round(mean(tmp.snp.2a[[j]]),3),
          
          REF = REF,
          ALT = ALT,
          
          
          
          beta = tryCatch(round(coef(summary(out))[2,1],3), error=function(e){
            print(paste0('error of tissue'))
            return(NA)}),
          
          p_val = tryCatch(round(coef(summary(out))[2,4],3), error=function(e){
            print(paste0('error of tissue'))
            return(NA)}),
          
          p_val_interaction = tryCatch(
            ifelse(
              length(coef(summary(out))[sk_int, 4]) == 0, 
              NA, 
              round(coef(summary(out))[sk_int, 4], 3)
            ), 
            error = function(e) {
              print(paste0('Error for tissue'))
              return(NA)})
          
        )
        
        
        
        # bind rows of temporary data frame to the results data frame
        
        df.out <- rbind(df.out, tmp.df)
        
        
      }
    }}


