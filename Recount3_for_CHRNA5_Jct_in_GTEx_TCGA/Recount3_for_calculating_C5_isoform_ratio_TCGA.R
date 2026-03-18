####
# Date: 12132024
# GOAL: Once get the junction count of all the tissue/cacer, extract the CHRNA5 iso of interest and do the ratio 
###

library(dplyr)
library(tidyr)
library(tibble)

#setwd('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/CHRNA5_project/Recount3_GTEx_junction_count_table/All_CHRAN5_jt/')
# TCGA
setwd('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/TCGA_data/TCGA_CHRNA5_project/TCGA_C5_junctioncount_gencodev26_from_Recount3/')

lf <- list.files('.',pattern = '.csv')

# load C5 TPM of TCGA in Xena data (TCGA_TARGET_GTEx project)

TCGA_TPM <- read.delim('../CHRNA5_geneTPM_from_Xena_TCGA_GTEx_TARGET_project.tsv')

# since it log2(x+0.0001) conver to TPM 

options(scipen = 999)
TCGA_TPM <- TCGA_TPM %>%
  mutate(CHRNA5_TPM = round(2^CHRNA5_TPM - 0.001, 2))


# the target isoform in C5

#chr15:78589850-78593091  R153_V416 Super-short
#chr15:78589983-78593091  Q197_V416 
#chr15:78590099-78593091  W236_V416 
#chr15:78590105-78593091	P238_V416 
#chr15:78590507-78593091  R372_V416 
#chr15:78590637-78593091  E415_V416 Full-length

# the PSAM4_C5 isoform 
# chr15.78546699.78580810 PSAM4_C5
# chr15.78565826.78580810 C5_intron1
					
# all_tar <- c( "R153_V416",
#            "Q197_V416",
#            "W236_V416",
#            "P238_V416",
#            "R372_V416",
#            "E415_V416",
#           "PSAM4_C5",
#           "C5_intron1")

iso_id <- c(  "R153_V416",
              "Q197_V416",
              "W236_V416",
              "P238_V416",
              "R372_V416",
              "E415_V416")


novel <- c("PSAM4_C5","C5_intron1")

# load the SMOKE and GT 

#fin <- read.delim('~/Desktop/chr15_CHRNA5_bladder_project/Finalized_GTEx_GT_CHRAN5_isofrom_TPM.tsv')
#fin <- fin[,c(1:88)]

df_ls <- list()

namelss <- data.frame(
old_names = c("chr15.78589850.78593091",
               "chr15.78589983.78593091",
               "chr15.78590099.78593091",
               "chr15.78590105.78593091",
               "chr15.78590507.78593091",
               "chr15.78590637.78593091"),

new_names = c( "R153_V416",
                "Q197_V416",
                "W236_V416",
                "P238_V416",
                "R372_V416",
                "E415_V416"))

name_map <- setNames(namelss$new_names, namelss$old_names)

df.out <- data.frame()

for (i in lf){
  
  # i <- lf[1]
  df <- read.csv(i)
  # rename the TCGA if "."
  df[,1] <- gsub("\\.","-",df[,1])
  names(df)[1] <- "full_id"
  df <- df %>% mutate(Sample_ID = gsub("^(TCGA-[A-Z0-9]+-[A-Z0-9]+-\\d{2}).*", "\\1", df$full_id),.before = 2)
  
  # Merge TPM result
  df <- left_join(df,TCGA_TPM,by = "Sample_ID")
  ###### make NA as empty value?? ######
  df$CHRNA5_TPM <- as.numeric(df$CHRNA5_TPM)

  # Replace matching names in df
  matching_cols <- intersect(names(df), names(name_map))
  names(df)[match(matching_cols, names(df))] <- name_map[matching_cols]

  

  #the PSMA4 C5 
  #old_names2 <- c("chr15.78546699.78580810","chr15.78565826.78580810")
  #new_names2 <- c("PSAM4_C5","C5_intron1")

# 
#   for (g in 1:length(old_names2)) {
#     if (old_names2[g] %in% names(df)) {
#       names(df)[grep(old_names2[g], names(df))] <- new_names2[g]
#     } else {
#       df[[new_names2[g]]] <- NA
#     }
#   }
  
  
  df <- df[,c(1:2,which(names(df) %in% "CHRNA5_TPM") ,which(names(df) %in% iso_id))]
  
  #df <- add_column(df, GTEx_ID=gsub("(GTEX-[^-]*)-.*", "\\1", df$full_id), .before = names(df)[1])
  
  # ISOFORM SUMMARY is enough
  # No need noval summaty which consider the PSMA4_CHRNA5 fusion
  
  df$iso_sum <- rowSums(df[, which(names(df) %in% iso_id)], na.rm = TRUE)
  
  
  #df$novel_sum <-  rowSums(df[, which(names(df) %in% novel)], na.rm = TRUE)
  # after sum, have to filter total sum reads > maybe > 10 or > 5?
  # df_isoform <- df %>% filter(iso_sum > 5)
  # df_PSAM4_C5_ext <-df %>% filter(novel_sum > 5)
  
  for (j in grep('V416',names(df),value = T)){
    # j <- "R153_V416"
    df[[paste0(j,"_rat")]] <- round(df[[j]] / df$iso_sum * 100,2)
  }
  
  # combine somke GT
  #df <- inner_join(fin,df_isoform,by="GTEx_ID")
  # filter never and current 
  #df <- df %>% filter(Smoke_CURRENT_FORMER_NEVER == "Current" | Smoke_CURRENT_FORMER_NEVER == "Never")
  
  
  write.table(df.out,paste0('~/Desktop/chr15_CHRNA5_bladder_project/Recount3_for_GTEx_and_TCGA/TCGA_JC/',gsub('_GTEx_.+','',i),'_'),quote = F,col.names = T,row.names = F,sep = ',')

}

  
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
#       tmp.snp <- df
#       tmp.snp.0a <- tmp.snp[ which(tmp.snp[,y] == 0), ]
#       tmp.snp.0b <- as.character(unique(tmp.snp.0a[gsub("_add", "", y)]))
#       tmp.snp.1a <- tmp.snp[ which(tmp.snp[,y] == 1), ]
#       tmp.snp.1b <- as.character(unique(tmp.snp.1a[gsub("_add", "", y)]))
#       tmp.snp.2a <- tmp.snp[ which(tmp.snp[,y] == 2), ]
#       tmp.snp.2b <- as.character(unique(tmp.snp.2a[gsub("_add", "", y)]))
#       
#       for (k in grep('Perc',names(df),value = T)){
#         # k <- "Perc_R153_V416"
#                
#           # fmla <- as.formula(paste(k, "~", y, " * Smoke_CURRENT_FORMER_NEVER_bi + AGE"))
#     
#           fmla <- as.formula(paste(k, "~", y, " * Smoke_CURRENT_FORMER_NEVER_bi + SEX + AGE"))
#           
#         
#         out <- tryCatch(lm(fmla,df), error=function(e){
#           print(paste0('error of tissue ',k))
#           return(NA)})
#         
#         summary(out)
#         
#         tmp.df <- data.frame(
#           snp = y,
#           variable = k,
#           geno_0 = tmp.snp.0b,
#           geno_1 = tmp.snp.1b,
#           geno_2 = tmp.snp.2b,
#           n_0 = nrow(tmp.snp.0a),
#           n_1 = nrow(tmp.snp.1a),
#           n_2 = nrow(tmp.snp.2a),
#           
#           beta = coef(summary(out))[2,1],
#           p_val = coef(summary(out))[2,4],
#           beta_smoke_inter = NA,  # Initialize to NA
#           p_val_smoke_inter = NA  # Initialize to NA
#         )
#         
#         # Check if there are any values for beta_smoke_inter and p_val_smoke_inter
#         smoke_inter_rows <- grep(":", rownames(coef(summary(out))), value = TRUE)
#         if (length(smoke_inter_rows) > 0) {
#           tmp.df$beta_smoke_inter <- coef(summary(out))[smoke_inter_rows, 1]
#           tmp.df$p_val_smoke_inter <- coef(summary(out))[smoke_inter_rows, 4]
#         }
#         
#         # tmp.df$variable <- paste0(k,"_",tmp.df$variable)
#         
#         
#         
#         
#         rownames(tmp.df) <-NULL
#         
#         df.out <- rbind(df.out, tmp.df)
#         
#         
#         
#         p <- interact_plot(
#           model = lm(fmla,df),
#           pred = !!y,
#           modx = Smoke_CURRENT_FORMER_NEVER_bi,
#           plot.points = TRUE, jitter = 0.125,
#           #x.label = "snp GT",
#           #y.label = short_name,
#           interval = TRUE, int.width = 0.8,
#           color.class = "Qual3",
#           point.size = 2,
#           point.shape = 19,
#           facet.modx =  FALSE)  + theme(legend.position = "top")
#         
#         
#         
#         plist[[paste0(k,"_",'interact_plot')]] <- p
#         
#         
#       }
#       
#       # Assuming plist is your list of plots
#       plots <- lapply(plist, function(x) x)
#       
#       # Now use do.call with plot_grid to combine the plots
#       combined_plot <- do.call(plot_grid, c(plots, list(nrow = 2)))
#       print(combined_plot)
#       
#       
#       
#       
#     
#     
#   }
#   
#   
#   
#   
#   
#   
# }
# 
# 
# 


