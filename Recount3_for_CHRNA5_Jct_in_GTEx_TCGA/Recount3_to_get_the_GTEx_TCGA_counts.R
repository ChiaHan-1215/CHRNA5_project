library(recount3)
library(snapcount)
library(megadepth)
library(dplyr)
library(tidyr)



# the temp file is located in leec20/.cache/R/recount3
# needd to clean if too many?
# Note the TCGA file need to used in BIouwlf Rstudio, it's too large

human_projects <- available_projects()

# get the TCGA and GTEx 

#gtex <- subset(human_projects, file_source == "gtex" & project_type == "data_sources")
TCGA <- subset(human_projects, file_source == "tcga" & project_type == "data_sources")

list_TCGA <- TCGA$project


for ( i in list_TCGA){
  
  df <- data.frame()
  
  # i <- "BLCA"
  target_tissue  <- TCGA %>% filter(project == i)
  
  rse_jxn_b <- create_rse(
    target_tissue,
    type = "jxn"
  )
  
  # To find the TCGA ID, it's in the package 
  TCGAIDmatch <- rse_jxn_b@colData@listData %>% as.data.frame()
  TCGAIDmatch <- TCGAIDmatch[,c(2,4)]
  
  matrix_data <- as.matrix(assay(rse_jxn_b))
  
  df <- as.data.frame(matrix_data)
  colnames(df) <- TCGAIDmatch$tcga.tcga_barcode[match(colnames(df), TCGAIDmatch$external_id)]
  
  # extract CHRNA5 region , it's hg38 cord
  
  # PSAM4 to CHRAN5 
  # chr15:78538430-78597269
  
  df$POS <- rownames(df)
  df <- df[which(grepl("chr15:",df$POS)),]
  df <- df %>%  separate(POS, into = c("chr", "start", "end"), sep = "[:-]")
  
  
  df$start <- as.numeric(df$start)
  df$end <- as.numeric(df$end)
  df <- df %>% filter(start > 78538430 & end > 78565138 & end < 78597269)
  
  df$location <- paste0('chr15:',df$start,"-",df$end)
  df <- df[,c(grep("location",names(df),value = T),grep("TCGA",names(df),value = T))]
  
  
  df.t <- as.data.frame(t(df))
  names(df.t) <- df.t[1,]
  df.t <- df.t[-1,]
  df.t[] <- lapply(df.t[], as.numeric)
  df.t <- df.t %>% mutate(TCGA_ID=rownames(df.t),.before = 1)
  
  write.table(df.t,paste0("~/Desktop/chr15_CHRNA5_bladder_project/Recount3_for_GTEx_and_TCGA/TCGA_JC/",i,"_C5iso_jc.csv"),col.names = T,row.names = F,sep = ',',quote = F)
  
}

# # Filter 
# target_tissue  <- TCGA %>% filter(project == "ACC")
#target_tissue  <- gtex %>% filter(project == "BLADDER")
# # subset 
# rse_target_tissue <- create_rse(target_tissue)
# 
# # explore
# metadata(rse_target_tissue)
# # raw count # transfrom
# assayNames(rse_target_tissue)
# assay(rse_target_tissue, "counts") <- transform_counts(rse_target_tissue)

# exon

# target_tissue_exon <- create_rse(
#   target_tissue,
#   type = "exon"
# )
# 
# target_tissue_exon
# rowRanges(target_tissue_exon)

 # EXON EXON Junxtion count




# hg38
# end of PSAM4 78548789 , exon8
# CHRNA5 exon1 skipping event from PSMA4 :  chr15 78546699 78580810
# intron1: chr15 78565826 78580810
# intron2: chr15 78580963 78586644


# ####
# 
# sb <- QueryBuilder(compilation="gtex", regions="CHRNA5")
# 
# c5.gene <- query_gene(sb)
# c5.jx <-  query_exon(sb)
# dim(c5.jx)
# t1 <- data.frame(c5.jx@colData)




