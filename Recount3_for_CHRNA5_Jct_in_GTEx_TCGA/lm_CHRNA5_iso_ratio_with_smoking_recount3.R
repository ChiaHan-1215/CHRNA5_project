################################################
# 07172024
# the TPM now include CHRAN5 and PSMA4 of total gene TPM and isofrom TPM
################################################
# CHRAN5 and PSMA4 in GTEx  with smoking lm() analysis #
# data: 07/10/2024
# add new variable "Cig_day_coded" as smoking intensity
################################################
# Note: 12172024: Previous is focus on the Brain tissue (BR)
# Now can add other tissue to test 
 


library(tidyverse)
library(dplyr)
library(ggsci)
library(ggrepel)
library(scales)
library(cowplot)
library(ggpubr)
library(rstatix)
library(patchwork)
library(tibble)
library(ggridges)

###############################################

# Function to extract variable names from a formula
get_vars_from_formula <- function(formula) {
  all.vars(formula)
}


##########################################
# 07222024 load the newest dataset
##########################################


# the file is named "GTEx_CHRNA3.CHRAN5.CHRNB4.PSMA4_total.iso.TPM_with_smoke_age_sex_race_SNPs_07232024.tsv"

df <- read_tsv('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/CHRNA5_project/updated_version_07232024/GTEx_CHRNA3.CHRAN5.CHRNB4.PSMA4_total.iso.TPM_with_smoke_age_sex_race_SNPs_07252024.tsv')



### Now we work on GT
# convert SNP to dosage 
# minor freq, effect = 0, need to table or sort the thing 

for (i in names(df)[15:55]){
  # i <- names(df)[19]
  n1 <- paste0(i,"_","add")
  tp1 <- as.data.frame(table(df[,i]) %>% sort())
  tp1$Var1 <- as.character(tp1[[1]])
  tp1 <- tp1 %>% mutate(Var1_copy = Var1) %>% separate(Var1_copy, into = c("n1", "n2"), sep = " ")
  tp1 <- tp1[,c(3,4,5,2)]
  
  SAMEGT <- data.frame()
  for (x in 1:3){
    if(tp1[x,2] == tp1[x,3]){
      SAMEGT <- rbind(tp1[x,],SAMEGT)
      SAMEGT <- arrange(SAMEGT,Freq)
    }}
  tp_left <- anti_join(tp1,SAMEGT)
  df[[n1]] <- car::recode(df[[i]]," SAMEGT$Var1[1]=0 ; SAMEGT$Var1[2]=2 ; tp_left$Var1[1]=1 ")
}

# organized total TPM
#df_f <- df[,c(names(df)[1:14],grep('rs[0-9]+',names(df),value = T) %>% sort(),names(df)[53:160])]
# organized isoform TPM
df_f <- df[,c(names(df)[1:14],grep('rs[0-9]+',names(df),value = T) %>% sort(),names(df)[56:2377])]


df_f[df_f == "No_data"] <- NA
df_f[, 97:ncol(df_f)] <- lapply(df_f[, 97:ncol(df_f)], as.numeric)
 
#names(df_f)[91:ncol(df_f)] <- gsub('\\.\\.\\.|\\.\\.|\\.','_',names(df_f)[91:ncol(df_f)])
#names(df_f)[91:ncol(df_f)] <- gsub('_$','',names(df_f)[91:ncol(df_f)])


# for isoform md 
 names(df_f)[97:ncol(df_f)] <- gsub(' - ','.',names(df_f)[97:ncol(df_f)])
 names(df_f)[97:ncol(df_f)] <- gsub(' ','.',names(df_f)[97:ncol(df_f)])
 names(df_f)[97:ncol(df_f)] <- gsub('\\(|\\)','',names(df_f)[97:ncol(df_f)])
# 


# makig factors of sex and race
df_f$SEX <- as.factor(df_f$SEX)
df_f$RACE <- as.factor(df_f$RACE)
names(df_f)[grep("c-1",names(df_f))] <- gsub('\\.c-1','',names(df_f)[grep("c-1",names(df_f))])
names(df_f)[grep("Cells.EBV-transformed.lymphocytes",names(df_f))] <- gsub('-','.',names(df_f)[grep("Cells.EBV-transformed.lymphocytes",names(df_f))])


# 07252024 make new column of assign new smoke
# Never vs Ever
df_f <- add_column(df_f,smoke_Current.Former_vs_Never=car::recode(df$smoking," 'no'='0';'former'='1';'current'='1'"),.before = 'smoking')
# Current vs no+former
df_f <- add_column(df_f,smoke_Current_vs_No.Former=car::recode(df$smoking," 'no'='0';'former'='0';'current'='1'"),.before = 'smoking')

### load the recount3 Brian
df_f_br <- df_f[,c(1:98,grep("CHRNA5.*geneTPM.*Brain", names(df_f)))]

BR <- read.csv('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/CHRNA5_project/Recount3_GTEx_junction_count_table/CHRNA5_6_isofoms_jc_with_GT_smoking/BRAIN_GTEx_CHRNA5_6_isofrom_jc.csv')
# the srm file from chr15_project_data_clean can trace back which part of brain
BR <- BR[,c(1,2,14:31)]
# change # to NA
BR <- BR %>% mutate(across(everything(), ~ ifelse(. == "#DIV/0!", NA, .)))


################################################
######## Update: 12192024 test muscle ########
################################################

df_f_mus <- df_f[,c(1:98,grep("CHRNA5.*geneTPM.*Muscle", names(df_f)))]
MUS <- read.csv('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/CHRNA5_project/Recount3_GTEx_junction_count_table/CHRNA5_6_isofoms_jc_with_GT_smoking/Muscle_GTEx_CHRNA5_6_isofrom_jc.csv',header = T)
MUS <- MUS[,c(1,2,16:22)]
# No calculating the ratio of isofoms
for (n in grep("V416",names(MUS),value = T)){
  # n <- grep("V416",names(MUS),value = T)[1]
  MUS[[paste0(n,'_rat')]] <- round(MUS[[n]] / MUS$iso_sum,3)
}
# change NaN to NA or other text
MUS <- MUS %>% mutate(across(everything(), ~ifelse(is.nan(.), "sum_is_0", .)))

# load the tissue detail

setwd('~/Desktop/dbGAP_GTEX_and_other_dataset/GTEX_v9/TPM_and_rawcount_stuff/')


# GTEX sample info for skecetal muscle

samid <- read.csv('../All_GTEx_sample_info_sex_age_race_etc/All_sample_info_donor_detailed/GTEx_Analysis_2021-02-11_v9_Annotations_GTEx_Analysis_2021-02-11_v9_Annotations_SampleAttributesDS.csv')
samid <- samid[,c(1,14)]
names(samid) <- c('full_id','tissue')


################################################################
# Updated
################################################################

set1 <- left_join(MUS,samid,by='full_id')
# clean up the name in tissue 
set1$tissue <- gsub(' - ','.',set1$tissue)
set1$tissue <- gsub(' ','.',set1$tissue)
set1$tissue <- gsub('\\(|\\)','',set1$tissue)
set1$tissue[grep("c-1",set1$tissue)] <- gsub('\\.c-1','',set1$tissue[grep("c-1",set1$tissue)])


set1[,c(3:15)] <- lapply(set1[,c(3:15)], function(x) as.numeric(x))
set1[,c(3:15)] <-set1[,c(3:15)] %>% mutate_if(is.numeric, round, 4)

set1 <- set1[,c(1,2,16,3:15)]

# remove tissue is NA. the sample is duplicated
set1 <- set1 %>% filter(!is.na(tissue))
set1 <- left_join(set1,df_f_mus,by='GTEx_ID')
#set1 <- set1[,c(1:21,119:131,22:118)]

################################################################

#### For Brain tissue ####
# set1 <- left_join(BR,samid,by='full_id')
# set1[,c(3:20)] <- lapply(set1[,c(3:20)], function(x) as.numeric(x))
# set1[,c(3:20)] <-set1[,c(3:20)] %>% mutate_if(is.numeric, round, 4)
# 
# 
# set1 <- set1[,c(1,2,21,3:20)]
# # remove tissue is NA. the sample is duplicated
# set1 <- set1 %>% filter(!is.na(tissue))
# set1 <- left_join(set1,df_f_br,by='GTEx_ID')
# set1 <- set1[,c(1:21,119:131,22:118)]

################################################################
################################################################


df_list <- list()
df_list[["all_data"]] <- set1
df_list[["Current_smoker"]] <- set1 %>% filter(set1$smoking == "current")
df_list[["Former_smoker"]] <- set1 %>% filter(set1$smoking == "former")
df_list[["None_smoker"]] <- set1 %>% filter(set1$smoking == "no")


######################## BACKUP plan B ######################## 
######################## ######################## #############
# if we seperate tissue and filter the CHRNA5_sum  > 2
# set1 can also work?
# 07292024 # use df.wide to test first 

# 
# df.wide <- set1 %>% tidyr::pivot_wider(names_from = tissue, values_from = c(names(set1)[4:21]))
# 
# # df.wide[] <- lapply(df.wide, function(x) replace(x, is.na(x), "No_data"))
# # df.wide[] <- lapply(df.wide, function(x) replace(x, is.na(x) | is.null(x), "No_data"))
# 
# # combine the df dataset 
# 
# cb <- df.wide
# cb <- left_join(df.wide,df_f_br,by="GTEx_ID")



# 
# # make list
# df_list <- list()
# df_list[["all_data"]] <- cb 
# df_list[["Current_smoker"]] <- cb %>% filter(cb$smoking == "current")
# df_list[["Former_smoker"]] <- cb %>% filter(cb$smoking == "former")
# df_list[["None_smoker"]] <- cb %>% filter(cb$smoking == "no")

######################## ######################## ######################## 
######################## ######################## ######################## 
######################## ######################## ######################## 

# make ratio of each tissue 
#tissue <- gsub('PSMA4_|CHRNA5_','',names(df_f)[91:ncol(df_f)]) %>% unique() 

# Get tissue unique names
# tissue <- gsub('CHRNA3_ENST.+_|CHRNA3_geneTPM_','',names(df_f)[91:576]) %>% unique()
tissue <- set1$tissue %>% unique()

# check
# names(cb) %>% grep("CHRNA5",.,value = T) %>% gsub('CHRNA5_geneTPM_','',.) %in% tissue


# for isoform 
gender_sp <- c("Ovary" , 'Uterus' , 'Testis' ,  'Vagina',
               "Prostate" , "Cervix.Ectocervix","Fallopian.Tube")

######################################################################################################################################################
######################################################################################################################################################
# Now seperate data by current, former, never and do eqtl
# TPM~SNP + sex + age 
# with interaction at the end
######################################################################################################################################################

#full_set <- df_list$all_data
df.out_all <- data.frame()

for(z in c("all_data","Current_smoker","Former_smoker","None_smoker") ){

  # z <- "all_data"
  # z <- "None_smoker"
  df_sub <- df_list[[z]]
  
  for (i in grep("_add",names(df_sub),value = T)){
    # grep("_add",names(df_sub),value = T)
    
    # i <- "rs8053_add"
    # i <- "rs142774214_add"
    # i <- "rs114205691_add"
    
    #smoking vs ratio with no SNPs
    
    # grep tissue 
    for (k in tissue){
      
      # k <- tissue[1]
      # k <- tissue[5]
      # tissue[-54]
      # k <- names(df_sub)[90:ncol(df_sub)][15]
      # k <- names(df_sub)[90:ncol(df_sub)][69]
      
      # TPM >0 within tissue k, NA keep or not?
      ssd <- df_sub %>%
        filter(tissue == k & .data[[grep(k, names(df_sub), value = TRUE)]] > 0) 
      
      # then filter CHRAN5_sum >=2
      
      ssd <- ssd %>% filter(CH5RNA_sum >= 2)
      
      # now select column for ratio
      
      selected_col <- grep('rat' ,names(ssd),value = T)
      
      for (j in selected_col){
        # j <- c(names(df_sub) %>% grep("rat",.,value = T))[1]
        # j <- selected_col[2]

        if(any(gender_sp %in% k)){
          
          # all types of fmla in the analysis
          fmla <- as.formula(paste(j, "~", i, " + AGE"))
          fmla.Cig <- as.formula(paste(j, "~", i, " + AGE + Cig_day_coded"))
          
          fmla_inter <- as.formula(paste(j, "~", i, " * smoking_coded + AGE"))
          fmla_inter.Cig <- as.formula(paste(j, "~", i, " * smoking_coded + AGE + Cig_day_coded"))
          
          fmla_inter_set2 <- as.formula(paste(j, "~", i, " * smoke_Current.Former_vs_Never + AGE"))
          fmla_inter_set3 <- as.formula(paste(j, "~", i, " * smoke_Current_vs_No.Former + AGE"))
          
          fmla_intens <- as.formula(paste(j, "~", i, " + AGE + smoking_coded"))
          fmla_intens.Cig <- as.formula(paste(j, "~", i, " + AGE + Cig_day_coded + smoking_coded"))
          
          # new model
          fmla_noSNP <- as.formula(paste(j, "~", " AGE + Cig_day_coded + smoking_coded"))
          
          
      
          
        } else { 
          
          fmla <- as.formula(paste(j, "~", i, " + SEX + AGE"))
          fmla.Cig <- as.formula(paste(j, "~", i, " + SEX + AGE + Cig_day_coded"))
          
          fmla_inter <- as.formula(paste(j, "~", i, " * smoking_coded + SEX + AGE"))
          fmla_inter.Cig <- as.formula(paste(j, "~", i, " * smoking_coded + SEX + AGE + Cig_day_coded"))
          
          fmla_inter_set2 <- as.formula(paste(j, "~", i, " * smoke_Current.Former_vs_Never  + SEX + AGE"))
          fmla_inter_set3 <- as.formula(paste(j, "~", i, " * smoke_Current_vs_No.Former + SEX + AGE"))
      
          
          fmla_intens <- as.formula(paste(j, "~", i, " + SEX + AGE + smoking_coded"))
          fmla_intens.Cig <- as.formula(paste(j, "~", i, " + SEX + AGE + Cig_day_coded + smoking_coded"))
          # new md 
          fmla_noSNP <- as.formula(paste(j, "~", " SEX + AGE + Cig_day_coded + smoking_coded"))
          
          }
        
        
          
        
        ############### Here ###
        
        # List of formulas
        formula_list <- list()
        
        if(z=='all_data'){
          formula_list <- list(fmla,fmla.Cig,fmla_inter,fmla_inter.Cig,fmla_intens,fmla_intens.Cig,fmla_inter_set2,fmla_inter_set3,fmla_noSNP)}else{formula_list <- list(fmla,fmla.Cig,fmla_intens,fmla_intens.Cig)}
        
        
        # Iterate through the list of formulas
        for (current_formula in formula_list) {
          # current_formula <- formula_list[[1]]
          # Extract variable names from the formula
          vars <- get_vars_from_formula(current_formula)
          
          # Filter the tar data frame based on the variables in the formula
          
          tar <- ssd[,c(names(ssd)[1:3],names(ssd)[35:49],gsub('_add',"",i),i,j)]
          # 
          tar <- tar %>%
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
          
          smoke_inter_rows <- tryCatch(
            grep(":", rownames(coef(summary(out))), value = TRUE),
            error = function(e) {
              print(paste0('Error for tissue ', k))
              return(character(0))
            }
          )
            
          ###############
          #### TABLE #####
          ###############
          
          tmp.df <- data.frame(
            snp = i,
            tissue=k,
            # add gene id nad isoform id
            variable=gsub(paste0('_',k),'',j),
            
            # add isoform id if ther's isoform, else put "whole gene"
            
            # isoformid=gsub(".*(ENST[0-9]+\\.[0-9]+).*", "\\1", j),
            
            # isoformid=gsub(".*(ENST[0-9]+\\.[0-9]+).*|.*(geneTPM).*", "\\1\\2", j),
            # tissue_gene_expr = j,
            #model = paste0(formula(current_formula)[2],"~",formula(current_formula)[3]),
            #model = paste0(formula(current_formula)[3]),
            model = gsub(i,'add',paste0(formula(current_formula)[3])),
            
            dataset = z,
            N_sample_not_NA = sum(!is.na(tar[[j]])),
            filter=paste0("TPM>0_and_CHRAN5_sum>=2"),
            
            geno_0 = tmp.snp.0b,
            geno_1 = tmp.snp.1b,
            geno_2 = tmp.snp.2b,
            n_0 = nrow(tmp.snp.0a),
            n_1 = nrow(tmp.snp.1a),
            n_2 = nrow(tmp.snp.2a),
            
            ## add TPM ###
            
            variable_mean = round(mean(tar[[j]]),3),
            variable_mean_g0 = round(mean(tmp.snp.0a[[j]]),3),
            variable_mean_g1 = round(mean(tmp.snp.1a[[j]]),3),
            variable_mean_g2 = round(mean(tmp.snp.2a[[j]]),3),
            
            REF= sapply(tmp.snp.2b, function(x) unique(strsplit(x, " ")[[1]])),
            ALT= sapply(tmp.snp.0b, function(x) unique(strsplit(x, " ")[[1]])),
            
            beta = tryCatch(round(coef(summary(out))[2,1],4), error=function(e){
              print(paste0('error of tissue ',k))
              return(NA)}) ,
            p_val = tryCatch({coef(summary(out))[2,4]}, error=function(e){
              print(paste0('error of tissue ',k))
              return(NA)}) ,
            
            p_val_inter_smoking = tryCatch(
              ifelse(
                length(coef(summary(out))[smoke_inter_rows, 4]) == 0, 
                NA, 
                round(coef(summary(out))[smoke_inter_rows, 4], 5)
              ), 
              error = function(e) {
                print(paste0('Error for tissue ', k))
                return(NA)
              }
            )
            
          )
          
          
          rownames(tmp.df) <- NULL
          df.out_all <- rbind(df.out_all, tmp.df)
          
          
        }
        

      }

    }

  }
}



### check p-value formatting 

df.out_all$p_val <- ifelse(df.out_all$p_val < 1e-4, 
                                     sprintf("%.2e", df.out_all$p_val), 
                                     as.numeric(df.out_all$p_val) %>% round(.,4))


#########################################################################################################
#########################################################################################################
#########################################################################################################
###NOW we add ref frequency
#########################################################################################################
#########################################################################################################
#########################################################################################################


# re-format
# or the output df.out_all
#summf <- read.csv('~/Desktop/lm_GTEx.All.tissue_CHRNA5.PSMA4.TPM_vs_SNP.smoke.age.sex_with_maf.TPM_value_07192024.csv')

summf<- df.out_all
names(summf)[1] <- 'rsid'
summf$rsid <- gsub('_add','',summf$rsid)

# freq file
rff <- read.csv('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/CHRNA5_project/files_for_calculating_freq/mapBC_setSNPs_to_test_in_GTEx_07112024.csv')
# select mulians
rff_select <- rff %>% filter(population == 'Multiancestry') %>% select(c("rsid", "RefAllele", "EffectAllele"))


# rename SNPs
rff_select$rsid[which(rff_select$rsid %in% "rs67624739")] <- "rs142774214"
# rs113931022 rename to rs114205691
rff_select$rsid[which(rff_select$rsid %in% "rs113931022")] <- "rs114205691"


unique(rff_select$rsid)

bb <- data.frame(
  rsid=c("rs3813572","rs12916483","rs56117933"),
  RefAllele=c("T","G","T"),
  EffectAllele=c("C","A","C")
)

rff_select <- rbind(rff_select,bb)
# check snps whcih are missing
# rs114205691 is missing

df_multiCancer.out <- left_join(summf, rff_select, by= "rsid")

table(df_multiCancer.out$REF)
table(df_multiCancer.out$ALT)

# df_multiCancer.out$REF <- ifelse(df_multiCancer.out$EffectAllele == "TGGGCGGGGCCAGAGGGAAATAG", "TGGGCGGGGCCAGAGGGAAATAG", df_multiCancer.out$REF)

df_multiCancer.out$beta.update <- round(ifelse(df_multiCancer.out$EffectAllele == df_multiCancer.out$ALT, df_multiCancer.out$beta, df_multiCancer.out$beta*-1),digits = 4)

names(df_multiCancer.out)

df_multiCancer.out$N_geno.0.update <- ifelse(df_multiCancer.out$RefAllele == df_multiCancer.out$REF, df_multiCancer.out$n_0, df_multiCancer.out$n_2)
df_multiCancer.out$N_geno.1.update <- df_multiCancer.out$n_1
df_multiCancer.out$N_geno.2.update <- ifelse(df_multiCancer.out$EffectAllele == df_multiCancer.out$ALT, df_multiCancer.out$n_2, df_multiCancer.out$n_0)
                     
# Gene name need to adjust 
df_multiCancer.out$geno.0.update <- ifelse(df_multiCancer.out$RefAllele == df_multiCancer.out$REF, df_multiCancer.out$geno_0, df_multiCancer.out$geno_2)
df_multiCancer.out$geno.1.update <- df_multiCancer.out$geno_1
df_multiCancer.out$geno.2.update <- ifelse(df_multiCancer.out$EffectAllele == df_multiCancer.out$ALT, df_multiCancer.out$geno_2, df_multiCancer.out$geno_0)
 
# Also TPM value need to adjust accordingly

df_multiCancer.out$variable_mean_g0.update <- ifelse(df_multiCancer.out$RefAllele == df_multiCancer.out$REF, df_multiCancer.out$variable_mean_g0, df_multiCancer.out$variable_mean_g2)
df_multiCancer.out$variable_mean_g1.update <- df_multiCancer.out$variable_mean_g1
df_multiCancer.out$variable_mean_g2.update <- ifelse(df_multiCancer.out$EffectAllele == df_multiCancer.out$ALT, df_multiCancer.out$variable_mean_g2, df_multiCancer.out$variable_mean_g0)


# N is number of sample 
                         
df_multiCancer.out$maf_REF <- round((df_multiCancer.out$N_geno.0.update*2 + df_multiCancer.out$N_geno.1.update) / (df_multiCancer.out$N_sample_not_NA*2),3)
df_multiCancer.out$maf_EFFECT <- round((df_multiCancer.out$N_geno.2.update*2 + df_multiCancer.out$N_geno.1.update) / (df_multiCancer.out$N_sample_not_NA*2),3)

# Now check which SNPs did not have meta analysis 

########################################
# If there are missing REF, use below
########################################

# 
# # Replace NA values in RefAllele and EffectAllele with "Not_included"
# df_multiCancer.out$RefAllele <- ifelse(is.na(df_multiCancer.out$RefAllele), "Not_included", df_multiCancer.out$RefAllele)
# df_multiCancer.out$EffectAllele <- ifelse(is.na(df_multiCancer.out$EffectAllele), "Not_included", df_multiCancer.out$EffectAllele)
# 
# # Update the columns based on the presence of NA in RefAllele and EffectAllele
# df_multiCancer.out <- df_multiCancer.out %>%
#   mutate(
#     beta.update = ifelse(RefAllele == "Not_included" & EffectAllele == "Not_included", beta, beta.update),
#     N_geno.0.update = ifelse(RefAllele == "Not_included" & EffectAllele == "Not_included", n_0, N_geno.0.update),
#     N_geno.1.update = ifelse(RefAllele == "Not_included" & EffectAllele == "Not_included", n_1, N_geno.1.update),
#     N_geno.2.update = ifelse(RefAllele == "Not_included" & EffectAllele == "Not_included", n_2, N_geno.2.update),
#     geno.0.update = ifelse(RefAllele == "Not_included" & EffectAllele == "Not_included", geno_0, geno.0.update),
#     geno.1.update = ifelse(RefAllele == "Not_included" & EffectAllele == "Not_included", geno_1, geno.1.update),
#     geno.2.update = ifelse(RefAllele == "Not_included" & EffectAllele == "Not_included", geno_2, geno.2.update)
#   )

                                        
# final data.frame
names(df_multiCancer.out)


# select the newer version of data 

#df_multiCancer.out_final <- df_multiCancer.out[, c(1:5,17,18,26,27,20:25,19,15:16) ]
# df_multiCancer.out_final <- df_multiCancer.out[, c(1:6,22,23,34,35,25:30,13,31:33,24,20,21) ]
df_multiCancer.out_final <- df_multiCancer.out[, c(1:7,23,24,35,36,26:31,14,32:34,25,21,22) ]


names(df_multiCancer.out_final) <- gsub(".update", "", names(df_multiCancer.out_final))


write.table(df_multiCancer.out_final,'~/Desktop/SET1.csv',col.names = T,row.names = F,sep = ',',quote = F)




                  ############################################################################################################
                  ################################################==== PLOT ===###############################################
                  ############################################################################################################
                  # making Plots for brain tissue for rs1427, 
                  # with CHRNA5 all isoforms 
                  # and with PSMA4 for all tissue
                  # with all variable not NA, not empty. 
                  ############################################################################################################
                  ############################################################################################################
                  ############################################################################################################

# load summary table
st <- read.csv('~/Desktop/SET1.csv')
st$dataset <- car::recode(st$dataset," 'Current_smoker'='current';'Former_smoker'='former';'None_smoker'='no' ")

# for CHRNA5
SummaryTB <-  st %>% filter(rsid=='rs142774214'& grepl("Brain", tissue)) %>% filter(geneid=='CHRNA5')
SummaryTB$model <- gsub("rs[0-9]+_","",SummaryTB$model)
SummaryTB$model <- gsub("Cig_day_coded","Cig_day",SummaryTB$model)
SummaryTB$model <- gsub("smoking_coded","smoking",SummaryTB$model)

####### for PSMA4, now add new 3 snps 
# SummaryTB <- st %>% filter(rsid!='rs142774214') %>% filter(geneid=='PSMA4') %>% filter(grepl("9082.5|44462.11|geneTPM", isoformid))
# SummaryTB$model <- gsub("rs[0-9]+_","",SummaryTB$model)
# SummaryTB$model <- gsub("Cig_day_coded","Cig_day",SummaryTB$model)
# SummaryTB$model <- gsub("smoking_coded","smoking",SummaryTB$model)

###### select model to include in the plot 
s.MD.all <- SummaryTB$model %>% unique() %>% .[c(1,3,7,8)]



# for CHRNA5
target_tissue <- grep("Brain",tissue,value = T)

# for PSMA4 
#target_tissue <- tissue 

df_sub <- df_list[["all_data"]]
  
  for (i in c("rs142774214_add")){
    # i <- "rs142774214_add"
    # i <- "rs8053_add"
    # grep tissue 
    for (k in target_tissue){
      
      # k <- target_tissue[1]
      
      selected_col <- grep(k,names(df_sub), value = TRUE) %>% grep("CHRNA5", ., value = TRUE) %>% setdiff(grep("ENST00000559554.5", ., value = TRUE))
      
      ########### with PSMA4 total and ENST00000559082.5 and ENST00000044462.11 and plot all tissue
      #selected_col <- grep(k ,names(df_sub),value = T) %>% grep("PSMA4",.,value = T) %>% grep("9082.5|44462.11|geneTPM",.,value = T)
      
      combined_plots <- list()# List to store all combined plots for this target_tissue
      
      for (j in selected_col){
        # j <- selected_col[4]
        
        plots <- list()
        
        fmla_intens.Cig <- as.formula(paste(j, "~", i, " + SEX + AGE + Cig_day_coded + smoking_coded"))
        
      
        
        # Extract variable names from the formula
        vars <- get_vars_from_formula(fmla_intens.Cig)
        tar <- df_sub[,c(names(df_sub)[1:16],gsub('_add',"",i),i,selected_col)]
        
        # tar <- tar %>%
        #   filter_at(vars, all_vars(!is.na(.) & . != ""))
        
        tar <- tar %>%
          filter(if_all(all_of(vars), ~ !is.na(.) & . != ""))
        
        # Calculate y-axis limits
        y_max <- max(tar[[j]], na.rm = TRUE) + 0.25
        
        
        # now seperate current former never
        
        for (m in c("all","current","former","no")){
          
          # m <- "Current_smoker"
          ###### NOW plot #######
          
          if (m == "all") {
            ppdata <- tar  # Use the whole dataset
            # for C5
            st.TARGET <- SummaryTB %>% filter(tissue==k & isoformid==gsub('CHRNA5_|_Brain.*','',j) & dataset=='all_data') 
            # for PSMA4
            # st.TARGET <- SummaryTB %>% filter(tissue==k & isoformid==gsub('PSMA4_|','',j) & dataset=='all_data') 
            
            st.TARGET <- st.TARGET[which(st.TARGET$model %in% s.MD.all),]
            
           
            
            annotations <- st.TARGET %>%
              mutate(annotation = paste0(
                "M:", model,
                ", beta=", round(beta, 3),
                ", p-val=", ifelse(p_val < 0.001, format(p_val, scientific = TRUE, digits = 3), round(p_val, 3)),
                ", p.int=", ifelse(p_val_inter_smoking < 0.001, format(p_val_inter_smoking, scientific = TRUE, digits = 3), round(p_val_inter_smoking, 3))
              )) %>%
              select(annotation)
            
            annotation_text <- paste0(annotations$annotation, collapse = "\n")
            
            
          } else {
            ppdata <- tar %>% filter(smoking == m)  # Filter by smoking status
            # for C5
            st.TARGET <- SummaryTB %>% filter(tissue==k & isoformid==gsub('CHRNA5_|_Brain.*','',j) & dataset==m)
            # for PSMA4 
            # st.TARGET <- SummaryTB %>% filter(tissue==k & isoformid==gsub('PSMA4_|','',j) & dataset==m)
            
            
            
            
            annotations <- st.TARGET %>%
              mutate(annotation = paste0(
                "M:", model,
                ", beta=", round(beta, 3),
                ", p-val=", ifelse(p_val < 0.001, format(p_val, scientific = TRUE, digits = 3), round(p_val, 3)))) %>%
              select(annotation)
            
            
            annotation_text <- paste0(annotations$annotation, collapse = "\n")

          }
          

          
          counts <- ppdata %>% group_by(.data[[gsub("_add","",i)]]) %>% summarise(n=n())
          x_labels <- paste(counts[[gsub("_add","",i)]], "\n(n=", counts$n, ")", sep = "")
          
          p <-ppdata %>%
            #ggplot(aes(x = rs8053, y = .data[[j]], fill = rs8053)) +
            ggplot(aes(x = rs142774214, y = .data[[j]], fill = rs142774214)) +
            # geom_violin(width = 1, lwd = 0.2, colour = "black",
            #             show.legend = NA, inherit.aes = TRUE) +
            # geom_violin(trim = FALSE) +
             ## modified
            
            geom_boxplot(
              width=0.5,lwd = 0.5, 
              outlier.shape=NA) + 
            
            scale_fill_manual(values = c("#f4a582", "#92c5de", "#4197EC"),drop=F) +
            
            geom_point(shape = 1, size = 1.5, colour = "black",
                       position = position_jitterdodge(jitter.width = 0.2,dodge.width = 1)) +
            # Add the mean point
            stat_summary(fun = mean, geom = "point",
                         position = position_dodge(width = 1),
                         shape = 19, size = 1.5, colour = "red") + theme_classic() + 
            theme(axis.title.x=element_blank(),
                  axis.text.y = element_text(color = "black",size = 10),
                  axis.text.x = element_text(color = "black",size = 8,angle = 45, hjust = 1), 
                  axis.title.y =element_blank(),legend.position = "none",
                  plot.title = element_text(size = 10, face = "plain", color = "black")) + 
            ggtitle(paste0("Data: ",m,'\nTarget: ',j,"\nsnp: ",gsub('_add','',i),"\n\n",annotation_text)) +
            
            ylim(0, y_max) + 
            scale_x_discrete(labels = x_labels)
          
          
          #print(p)
          plots[[m]] <- p
          
        }
        
        combined_plot <- plots[["all"]] + plots[["current"]] + plots[["former"]] + plots[["no"]] + plot_layout(ncol = 4)
        
        combined_plots[[length(combined_plots) + 1]] <- combined_plot
        
        # Print the combined plot
        #print(combined_plot)
        #ggsave(paste0("~/Desktop/chr15_CHRNA5_bladder_project/plots/","rs142774214_",j,".pdf"), combined_plot, width = 16, height = 8)
        
      }
      mega_plot <- wrap_plots(combined_plots, ncol = 1)
      
      ggsave(paste0("~/Desktop/chr15_CHRNA5_bladder_project/plots_tpm_gt/",gsub('_ENST[0-9].*|_geneTPM.*','',j),"_",gsub("_add","",i),"_",k,".pdf"), mega_plot, width = 26, height = 8 * length(combined_plots))
      
      
      
    }
  }



