################################################
# 07192024,ver3 
# the TPM now include CHRAN5 and PSMA4 of total gene TPM and isofrom TPM
################################################
# CHRAN5 and PSMA4 in GTEx  with smoking lm() analysis #
# data: 07/10/2024
# data: 07/17/2024 add new variable "Cig_day_coded" as smoking intensity
################################################

# load package 

library(tidyverse)
library(dplyr)

###############################################
# # load total gene TPM data with smoking intensity

#df <- read_tsv('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/CHRNA5_project/GTEx_CHRNA5_PSMA4_TPM_with_smoke_age_sex_race_SNPs_Mila_07172024.txt')

# clean data
#df <- df[,-c(18,21,23,26,28,30,32)]
#df <- df[,-c(19,22,24,27,29,31,33)]

# load isoform TPM data
df <- read_tsv('~/Desktop/GTEx_CHRNA5_PSMA4_isoform.TPM_with_smoke_age_sex_race_SNPs_07192024.tsv')


# convert SNP to dosage 
# minor freq, effect = 0, need to table or sort the thing 
# the code below is for SNPs to convert to 0,1,2
# use the freq table to get the GT

# the SNPs are located from column 15 to 52
for (i in names(df)[15:52]){
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

# organized total TPM dataset columns

#df_f <- df[,c(names(df)[1:14],grep('rs[0-9]+',names(df),value = T) %>% sort(),names(df)[53:160])]

# organized isoform TPM dataset columns
df_f <- df[,c(names(df)[1:14],grep('rs[0-9]+',names(df),value = T) %>% sort(),names(df)[53:1402])]


df_f[df_f == "No_data"] <- NA
df_f[, 91:ncol(df_f)] <- lapply(df_f[, 91:ncol(df_f)], as.numeric)
 
#names(df_f)[91:ncol(df_f)] <- gsub('\\.\\.\\.|\\.\\.|\\.','_',names(df_f)[91:ncol(df_f)])
#names(df_f)[91:ncol(df_f)] <- gsub('_$','',names(df_f)[91:ncol(df_f)])


# for isoform md 
 names(df_f)[91:ncol(df_f)] <- gsub(' - ','.',names(df_f)[91:ncol(df_f)])
 names(df_f)[91:ncol(df_f)] <- gsub(' ','.',names(df_f)[91:ncol(df_f)])
 names(df_f)[91:ncol(df_f)] <- gsub('\\(|\\)','',names(df_f)[91:ncol(df_f)])
# 


# makig factors of sex and race
df_f$SEX <- as.factor(df_f$SEX)
df_f$RACE <- as.factor(df_f$RACE)
names(df_f)[grep("c-1",names(df_f))] <- gsub('\\.c-1','',names(df_f)[grep("c-1",names(df_f))])
names(df_f)[grep("Cells.EBV-transformed.lymphocytes",names(df_f))] <- gsub('-','.',names(df_f)[grep("Cells.EBV-transformed.lymphocytes",names(df_f))])

# make list that is eazier to loop
df_list <- list()
df_list[["all_data"]] <- df_f 
df_list[["Current_smoker"]] <- df_f %>% filter(df_f$smoking == "current")
df_list[["Former_smoker"]] <- df_f %>% filter(df_f$smoking == "former")
df_list[["None_smoker"]] <- df_f %>% filter(df_f$smoking == "no")


# make ratio of each tissue 
#tissue <- gsub('PSMA4_|CHRNA5_','',names(df_f)[91:ncol(df_f)]) %>% unique() 

# for isoform md 
tissue <- gsub('PSMA4_ENST.+_|CHRNA5_ENST.+_','',names(df_f)[91:ncol(df_f)]) %>% unique()


####################################################################
# Now seperate data by current, former, never and do eqtl
# TPM~SNP + sex + age 
# with interaction at the end
########################################################################################################################################


for(z in c("all_data","Current_smoker","Former_smoker","None_smoker") ){

  # z <- "all_data"
  # z <- "Current_smoker"
  df_sub <- df_list[[z]]
  
  for (i in "rs8053_add" ){
    # grep("_add",names(df_sub),value = T)
    # i <- "rs8053_add"
    # grep tissue 
    for (k in tissue[-54]){
    
      # k <- tissue[15]
      # Since some tissue only belongs to Female or Male, need to set up the if condition
      # gender_sp <- c("Ovary" , 'Uterus' , 'Testis' ,  'Vagina',"Prostate" , "Cervix_Ectocervix","Fallopian_Tube")
      # 
      
      # for isoform 
      gender_sp <- c("Ovary" , 'Uterus' , 'Testis' ,  'Vagina',
                     "Prostate" , "Cervix.Ectocervix","Fallopian.Tube")
      
      selected_col <- grep(k ,names(df_sub),value = T)
      
      for (j in selected_col){
        # j <- "CHRNA5_Lung"
        # j <- "PSMA4_ENST00000560099.5_Thyroid"
    

        if(any(gender_sp %in% k)){
          
          fmla_all <- as.formula(paste(j, "~", i, " + AGE + smoking_coded + ", i," * smoking_coded"))
          fmla_inter <- as.formula(paste(j, "~", i, " * smoking_coded + AGE"))
          fmla <- as.formula(paste(j, "~", i, " + AGE"))
          
          #### add Cig_per_day as new variable
          
          fmla_all.Cig <- as.formula(paste(j, "~", i, " + AGE + smoking_coded + ", i," * smoking_coded + Cig_day_coded"))
          fmla_inter.Cig <- as.formula(paste(j, "~", i, " * smoking_coded + AGE + Cig_day_coded"))
          fmla.Cig <- as.formula(paste(j, "~", i, " + AGE + Cig_day_coded"))
          
          
        } else { 
          
          fmla_all <- as.formula(paste(j, "~", i, " + SEX + AGE + smoking_coded + ",i," * smoking_coded"))
          fmla_inter <- as.formula(paste(j, "~", i, " * smoking_coded + SEX + AGE"))
          fmla <- as.formula(paste(j, "~", i, " + SEX + AGE"))
          
          #### add Cig_per_day as new variable
          
          fmla_all.Cig <- as.formula(paste(j, "~", i, " + SEX + AGE + smoking_coded + ",i," * smoking_coded  + Cig_day_coded"))
          fmla_inter.Cig <- as.formula(paste(j, "~", i, " * smoking_coded + SEX + AGE + Cig_day_coded"))
          fmla.Cig <- as.formula(paste(j, "~", i, " + SEX + AGE + Cig_day_coded"))
      
          
        }
        
        # Function to extract variable names from a formula
        
        get_vars_from_formula <- function(formula) {
          all.vars(formula)
        }
        
        # List of formulas
        formula_list <- list()
        
        if(z=='all_data'){
          formula_list <- list(fmla_all, fmla_inter, fmla, fmla_all.Cig, fmla_inter.Cig, fmla.Cig)}else{formula_list <- list(fmla,fmla.Cig)}
        
        
        # Iterate through the list of formulas
        for (current_formula in formula_list) {
          # current_formula <- formula_list[[3]]
          # Extract variable names from the formula
          vars <- get_vars_from_formula(current_formula)
          
          # Filter the tar data frame based on the variables in the formula and filter out NA or empty data to get the N_sample
          
          tar <- df_sub[,c(names(df_sub)[1:14],gsub('_add',"",i),i,selected_col)]


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
            
            
          
          # now set up the summary table
          tmp.df <- data.frame(
            snp = i,
            tissue=k,
            tissue_gene_expr = j,
            model = paste0(formula(current_formula)[2],"~",formula(current_formula)[3]),
            
            dataset = z,
            N_sample_not_NA = sum(!is.na(tar[[j]])),
            
            geno_0 = tmp.snp.0b,
            geno_1 = tmp.snp.1b,
            geno_2 = tmp.snp.2b,
            n_0 = nrow(tmp.snp.0a),
            n_1 = nrow(tmp.snp.1a),
            n_2 = nrow(tmp.snp.2a),
            
            ## add TPM ###
            
            TPM_mean = round(mean(tar[[j]]),3),
            TPM_geno_0 = round(mean(tmp.snp.0a[[j]]),3),
            TPM_geno_1 = round(mean(tmp.snp.1a[[j]]),3),
            TPM_geno_2 = round(mean(tmp.snp.2a[[j]]),3),
            
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


###################################################################
###NOW add ref frequency
####################################################################



#load the output df.out_all

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

# check snps whcih are missing
unique(rff_select$rsid)

df_multiCancer.out <- left_join(summf, rff_select, by= "rsid")

table(df_multiCancer.out$REF)
table(df_multiCancer.out$ALT)

# df_multiCancer.out$REF <- ifelse(df_multiCancer.out$EffectAllele == "TGGGCGGGGCCAGAGGGAAATAG", "TGGGCGGGGCCAGAGGGAAATAG", df_multiCancer.out$REF)

# check if the Ref or Effect allele are mathed or not, if not, change to the same pattern
df_multiCancer.out$beta.update <- round(ifelse(df_multiCancer.out$EffectAllele == df_multiCancer.out$ALT, df_multiCancer.out$beta, df_multiCancer.out$beta*-1),digits = 4)

names(df_multiCancer.out)

df_multiCancer.out$N_geno.0.update <- ifelse(df_multiCancer.out$RefAllele == df_multiCancer.out$REF, df_multiCancer.out$n_0, df_multiCancer.out$n_2)
df_multiCancer.out$N_geno.1.update <- df_multiCancer.out$n_1
df_multiCancer.out$N_geno.2.update <- ifelse(df_multiCancer.out$EffectAllele == df_multiCancer.out$ALT, df_multiCancer.out$n_2, df_multiCancer.out$n_0)
                     
# GT adjust 
df_multiCancer.out$geno.0.update <- ifelse(df_multiCancer.out$RefAllele == df_multiCancer.out$REF, df_multiCancer.out$geno_0, df_multiCancer.out$geno_2)
df_multiCancer.out$geno.1.update <- df_multiCancer.out$geno_1
df_multiCancer.out$geno.2.update <- ifelse(df_multiCancer.out$EffectAllele == df_multiCancer.out$ALT, df_multiCancer.out$geno_2, df_multiCancer.out$geno_0)
 
# Also TPM value need to adjust accordingly
df_multiCancer.out$TPM_geno_0.update <- ifelse(df_multiCancer.out$RefAllele == df_multiCancer.out$REF, df_multiCancer.out$TPM_geno_0, df_multiCancer.out$TPM_geno_2)
df_multiCancer.out$TPM_geno_1.update <- df_multiCancer.out$TPM_geno_1
df_multiCancer.out$TPM_geno_2.update <- ifelse(df_multiCancer.out$EffectAllele == df_multiCancer.out$ALT, df_multiCancer.out$TPM_geno_2, df_multiCancer.out$TPM_geno_0)


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
df_multiCancer.out_final <- df_multiCancer.out[, c(1:6,22,23,34,35,25:30,13,31:33,24,20,21) ]

names(df_multiCancer.out_final) <- gsub(".update", "", names(df_multiCancer.out_final))
