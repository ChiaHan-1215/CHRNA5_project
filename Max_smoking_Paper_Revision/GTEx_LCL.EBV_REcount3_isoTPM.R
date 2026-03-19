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
mani_bladder_with_iso <- inner_join(iso,mani_bladder,by = 'Sample_ID')
rownames(mani_bladder_with_iso) <- mani_bladder_with_iso$Sample_ID
mani_bladder_with_iso <- mani_bladder_with_iso[,c(-1,-12)]

########################################################################################################################
# laod Recount3 files
########################################################################################################################

# GTEx
setwd('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/project_CHRNA5/Recount3_GTEx_junction_count_table/All_CHRAN5_jt/')



lf <- list.files('.',pattern = '.csv')




# the target isoform in C5, check with SHSY5Y Jct count, only 1bp move so it's ok

#chr15:78589850-78593091  R153_V416 Super-short
#chr15:78589983-78593091  Q197_V416 
#chr15:78590099-78593091  W236_V416 
#chr15:78590105-78593091	P238_V416 
#chr15:78590507-78593091  R372_V416 
#chr15:78590637-78593091  E415_V416 Full-length

# the PSAM4_C5 isoform 
# chr15.78546699.78580810 PSAM4_C5
# chr15.78565826.78580810 C5_intron1

all_tar <- c( "R153_V416",
              "Q197_V416",
              "W236_V416",
              "P238_V416",
              "R372_V416",
              "E415_V416",
              "PSAM4_C5",
              "C5_intron1")

iso_id <- c(  "R153_V416",
              "Q197_V416",
              "W236_V416",
              "P238_V416",
              "R372_V416",
              "E415_V416")


novel <- c("PSAM4_C5","C5_intron1")

# load the SMOKE and GT and TPM 
# also have EBV cell "Cells_EBV_transformedlymphocytes"
fin <- read.delim('~/Desktop/chr15_CHRNA5_bladder_project/Finalized_GTEx_GT_CHRAN5_isofrom_TPM.tsv')

### Form Here we can get isoform TPM 


# Get GT only, this time we extract rs169xxxx
fin <- fin[,c(1:88)]

df_ls <- list()

# 
# df.out <- data.frame()
# for (i in lf){
#   
#   # i <- lf[16]
#   df <- read.csv(i)
#   names(df)[1] <- "full_id"
#   
#   old_names <- c("chr15.78589850.78593091",
#                  "chr15.78589983.78593091",
#                  "chr15.78590099.78593091",
#                  "chr15.78590105.78593091",
#                  "chr15.78590507.78593091",
#                  "chr15.78590637.78593091")
#   
#   new_names <- c( "R153_V416",
#                   "Q197_V416",
#                   "W236_V416",
#                   "P238_V416",
#                   "R372_V416",
#                   "E415_V416")
#   
#   for (x in 1:length(old_names)) {
#     # x <- 1
#     if (old_names[x] %in% names(df)) {
#       names(df)[grep(old_names[x], names(df))] <- new_names[x]
#     } else {
#       df[[new_names[x]]] <- NA
#     }
#   }
#   
#   
#   #the PSMA4 C5 
#   old_names2 <- c("chr15.78546699.78580810","chr15.78565826.78580810")
#   new_names2 <- c("PSAM4_C5","C5_intron1")
#   
#   
#   
#   for (g in 1:length(old_names2)) {
#     if (old_names2[g] %in% names(df)) {
#       names(df)[grep(old_names2[g], names(df))] <- new_names2[g]
#     } else {
#       df[[new_names2[g]]] <- NA
#     }
#   }
#   
#   
#   
#   
#   df <- df[,c(1,which(names(df) %in% all_tar))]
#   #df <- add_column(df, GTEx_ID=gsub("(GTEX-[^-]*)-.*", "\\1", df$full_id), .before = names(df)[1])
#   df$iso_sum <- rowSums(df[, which(names(df) %in% iso_id)], na.rm = TRUE)
#   df$novel_sum <-  rowSums(df[, which(names(df) %in% novel)], na.rm = TRUE)
#   # after sum, have to filter total sum reads > maybe > 10 or > 5?
#   df_isoform <- df %>% filter(iso_sum > 5)
#   # df_PSAM4_C5_ext <-df %>% filter(novel_sum > 5)
#   
#   
#   
#   
#   for (j in grep('V416',names(df_isoform),value = T)){
#     # j <- iso_id[2]
#     df_isoform[[paste0('Perc_',j)]] <- round(df_isoform[[j]] / df_isoform$iso_sum * 100,2)
#   }
#   
#   # combine somke GT
#   df <- inner_join(fin,df_isoform,by="GTEx_ID")
#   # filter never and current 
#   df <- df %>% filter(Smoke_CURRENT_FORMER_NEVER == "Current" | Smoke_CURRENT_FORMER_NEVER == "Never")
#   
#   
#   
#   # do lm 
#   # grep('_add',names(df_final),value = T)
#   for (y in c("rs142774214_add","rs503464_add","rs55853698_add","rs55781567_add")){
#     
#     # y <- "rs142774214_add"
#     plist <- list()
#     
#     df$Smoke_CURRENT_FORMER_NEVER_bi <- car::recode(df$Smoke_CURRENT_FORMER_NEVER, " 'Current'=2 ; 'Former'=1  ; 'Never'=0 ")
#     
#     df$Smoke_CURRENT_FORMER_NEVER_bi <- as.factor(df$Smoke_CURRENT_FORMER_NEVER_bi)
#     
#     
#     tmp.snp <- df
#     tmp.snp.0a <- tmp.snp[ which(tmp.snp[,y] == 0), ]
#     tmp.snp.0b <- as.character(unique(tmp.snp.0a[gsub("_add", "", y)]))
#     tmp.snp.1a <- tmp.snp[ which(tmp.snp[,y] == 1), ]
#     tmp.snp.1b <- as.character(unique(tmp.snp.1a[gsub("_add", "", y)]))
#     tmp.snp.2a <- tmp.snp[ which(tmp.snp[,y] == 2), ]
#     tmp.snp.2b <- as.character(unique(tmp.snp.2a[gsub("_add", "", y)]))
#     
#     for (k in grep('Perc',names(df),value = T)){
#       # k <- "Perc_R153_V416"
#       
#       # fmla <- as.formula(paste(k, "~", y, " * Smoke_CURRENT_FORMER_NEVER_bi + AGE"))
#       
#       fmla <- as.formula(paste(k, "~", y, " * Smoke_CURRENT_FORMER_NEVER_bi + SEX + AGE"))
#       
#       
#       out <- tryCatch(lm(fmla,df), error=function(e){
#         print(paste0('error of tissue ',k))
#         return(NA)})
#       
#       summary(out)
#       
#       tmp.df <- data.frame(
#         snp = y,
#         variable = k,
#         geno_0 = tmp.snp.0b,
#         geno_1 = tmp.snp.1b,
#         geno_2 = tmp.snp.2b,
#         n_0 = nrow(tmp.snp.0a),
#         n_1 = nrow(tmp.snp.1a),
#         n_2 = nrow(tmp.snp.2a),
#         
#         beta = coef(summary(out))[2,1],
#         p_val = coef(summary(out))[2,4],
#         beta_smoke_inter = NA,  # Initialize to NA
#         p_val_smoke_inter = NA  # Initialize to NA
#       )
#       
#       # Check if there are any values for beta_smoke_inter and p_val_smoke_inter
#       smoke_inter_rows <- grep(":", rownames(coef(summary(out))), value = TRUE)
#       if (length(smoke_inter_rows) > 0) {
#         tmp.df$beta_smoke_inter <- coef(summary(out))[smoke_inter_rows, 1]
#         tmp.df$p_val_smoke_inter <- coef(summary(out))[smoke_inter_rows, 4]
#       }
#       
#       # tmp.df$variable <- paste0(k,"_",tmp.df$variable)
#       
#       
#       
#       
#       rownames(tmp.df) <-NULL
#       
#       df.out <- rbind(df.out, tmp.df)
#       
#       
#       
#       p <- interact_plot(
#         model = lm(fmla,df),
#         pred = !!y,
#         modx = Smoke_CURRENT_FORMER_NEVER_bi,
#         plot.points = TRUE, jitter = 0.125,
#         #x.label = "snp GT",
#         #y.label = short_name,
#         interval = TRUE, int.width = 0.8,
#         color.class = "Qual3",
#         point.size = 2,
#         point.shape = 19,
#         facet.modx =  FALSE)  + theme(legend.position = "top")
#       
#       
#       
#       plist[[paste0(k,"_",'interact_plot')]] <- p
#       
#       
#     }
#     
#     # Assuming plist is your list of plots
#     plots <- lapply(plist, function(x) x)
#     
#     # Now use do.call with plot_grid to combine the plots
#     combined_plot <- do.call(plot_grid, c(plots, list(nrow = 2)))
#     print(combined_plot)
#     
#     
#     
#     
#     
#     
#   }
# }






