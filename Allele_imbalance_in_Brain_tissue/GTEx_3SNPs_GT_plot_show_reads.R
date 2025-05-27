# Note the 3 SNPs are 
#POS_hg38	rsID	REF, R	ALT, A
#78595892	rs10637216	A	ACCCC
#78596058	rs578776	G	A
#78595490	rs660652	G	A

# the rs10637216 is interest, the location need to set to 78595895 for detect INS
# the result are located in Biowulf igvtool_count_GT_rstudio_biowulf_env folder 

library(dplyr)

setwd('/Volumes/data/igvtool_count_GT_rstudio_biowulf_env/gt_res/Green/')

ff <- list.files('.',pattern = '.txt')

library(dplyr)
library(tidyr)
library(ggplot2)

df <- data.frame()

# need to modify column number when use! 

for (x in ff){
  # x <- ff[2]
  tar <- read.delim(x)
  tar$snpid <- gsub('_.*','',x)
  tar$bam_ID <- gsub('.hg38','',tar$bam_ID)
  tar <- tar[,c(1,10,2,3,4,5,6,7,8,9)]
  df <- rbind(df,tar)
}

for (i in 1:nrow(df)){
  df$tot[i] <- sum(df[i,c(4:10)])
}

# remove row of total reads = 0 

for (x in 1:nrow(df)) {
  if (!is.na(df$tot[x]) && df$tot[x] == 0) {
    df <- df[-x, ]
  }
}

# Get all the brain tissue 
tissue_list <- gsub('.+\\.+','',df$bam_ID) %>% sort() %>% unique()

# remove total read count less that 5 read
df <- df %>% filter(tot >= 5)

# laod GT form T-drive 
genotype <- read.delim('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/CHRNA5_project/updated_version_07232024/GTEx_CHRNA3.CHRAN5.CHRNB4.PSMA4_total.iso.TPM_with_smoke_age_sex_race_SNPs_05142025.tsv')
genotype <- genotype[,c(2,15,17,16)]
genotype$rs660652 <- car::recode(genotype$rs660652,"'A G'='G A' ")



df_plot2 <- df %>%
  filter(snpid == "rs578776") %>%
  separate(bam_ID, into = c("id", "tissue"), sep = "\\.", extra = "merge") %>%
  mutate(
    GTEx_ID = sub("(GTEX-[^-]+).*", "\\1", id),
    frac_A  = round(A / tot, 2),
    frac_G  = round(G / tot, 2)
  ) %>%
  pivot_longer(
    cols      = c(frac_A, frac_G),
    names_to  = "allele",
    values_to = "fraction"
  ) %>%
  # add a combined x‐axis label with tissue and total reads
  mutate(
    #label = paste0(tissue, "\nreads: ", tot)
    label = paste0(tissue)
  )

# Table 

df_wide.rs57 <- df_plot2 %>%
  select(GTEx_ID, tissue, snpid, allele, fraction,tot) %>%
  pivot_wider(
    names_from   = c(snpid, allele),
    values_from  = fraction,
    names_sep    = "_"
  )%>%
  rename(
    rs578776_total_read = tot
  )

# Merge all the table 
snps_merge <- df_wide.rs57
snps_merge <- snps_merge[,c(1,2,4,5,3)]

fff <- left_join(snps_merge,genotype,by="GTEx_ID")

############ Newer displaced which makes force show all tissue for alignment all sample#######


# 1. Make tissue a factor with all your desired levels
df_plot_update <- df_plot2 %>%
  mutate(
    tissue = factor(tissue, levels = tissue_list)
  ) %>%
  # 2. Explicitly complete every combo of sample, allele, and tissue
  complete(
    GTEx_ID,
    allele,
    tissue,
    fill = list(fraction = NA)       # or = 0 if you want zero-height bars
  )

# make it NA tissue to display only tissue name instead of NA 
for(x in 1:nrow(df_plot_update)){
  if (is.na(df_plot_update$label[x])){
    df_plot_update$label[x] <-  as.character(df_plot_update$tissue[x])
  }
}


#################################################
#################################################
##### ADD READ COUNT ON the PLOT
##################################################
##################################################

# 1. Build the plot data, adding a 'reads' column
df_plot_reads <- df_plot_update %>%
  mutate(
    reads = case_when(
      allele == "frac_A" ~ A,    # take A‐count when showing frac_A
      allele == "frac_G" ~ G,    # take G‐count when showing frac_G
      TRUE               ~ NA_real_
    )
  )

df_plot_reads$label <- gsub('Brain_','',df_plot_reads$label)


# 2. Base bar plot
green_plot <- ggplot(df_plot_reads, aes(x = label, y = fraction, fill = allele)) +
  geom_col(position = "dodge", color = "black", na.rm = TRUE) +
  scale_fill_manual(values = c("frac_A" = "#FBB4AE", "frac_G" = "#B3CDE3")) +
  facet_wrap(~ GTEx_ID, scales = "free_x", ncol = 3) +
  scale_x_discrete(drop = FALSE) +
  theme_classic() +
  theme(
    panel.border    = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    legend.title    = element_blank(),
    strip.background= element_rect(fill = "#a1a7ae", color = NA),
    strip.text      = element_text(color = "black", face = "bold"),
    plot.margin     = margin(l = 15, t = 5, r = 5, b = 5)
  ) +
  labs(
    title = "Allelic read counts of rs578776 across brain tissues"
  ) 

# 3. Add read‐count labels above each bar
p1 <- green_plot +
  geom_text(
    aes(label = reads, group = allele),
    position = position_dodge(width = 0.9),
    vjust    = -0.35,
    size     = 3
  ) + scale_y_continuous(expand = expansion(mult = 0),limits = c(0, 1.10)) + 
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 50, unit = "pt"))+
  labs(
    x     = NULL,
    y     = NULL) 

### Keep the Green dataset and plot 
green_data <- df_plot_reads
p1

################################################# 
## For yellow colored sample, print every 6 sample in one plot 
#################################################

setwd('/Volumes/data/igvtool_count_GT_rstudio_biowulf_env/gt_res/Yellow/')
ff <- list.files('.',pattern = '.txt')
df <- data.frame()

# need to modify column number when use! 

for (x in ff){
  # x <- ff[2]
  tar <- read.delim(x)
  tar$snpid <- gsub('_.*','',x)
  tar$bam_ID <- gsub('.hg38','',tar$bam_ID)
  tar <- tar[,c(1,10,2,3,4,5,6,7,8,9)]
  df <- rbind(df,tar)
}

for (i in 1:nrow(df)){
  df$tot[i] <- sum(df[i,c(4:10)])
}

# remove row of total reads = 0 

for (x in 1:nrow(df)) {
  if (!is.na(df$tot[x]) && df$tot[x] == 0) {
    df <- df[-x, ]
  }
}
# remove total read count less that 5 read
df <- df %>% filter(tot >= 5)



df_plot2 <- df %>%
  filter(snpid == "rs578776") %>%
  separate(bam_ID, into = c("id", "tissue"), sep = "\\.", extra = "merge") %>%
  mutate(
    GTEx_ID = sub("(GTEX-[^-]+).*", "\\1", id),
    frac_A  = round(A / tot, 2),
    frac_G  = round(G / tot, 2)
  ) %>%
  pivot_longer(
    cols      = c(frac_A, frac_G),
    names_to  = "allele",
    values_to = "fraction"
  ) %>%
  # add a combined x‐axis label with tissue and total reads
  mutate(
    label = paste0(tissue)
    #label = paste0(tissue, "\nreads: ", tot)
  )



# HAVE to select heterozygyes GT just for detection AI 
# The rs578776 AG sample ID:
SoInterest <- c("GTEX-11DXY", "GTEX-11ZUS", "GTEX-12ZZZ", "GTEX-13VXU", "GTEX-15DYW", "GTEX-1EMGI",
                "GTEX-1HSKV", "GTEX-N7MS", "GTEX-R55E", "GTEX-WHSE" )





# 1. Make tissue a factor with all your desired levels
df_plot_update <- df_plot2 %>%
  mutate(
    tissue = factor(tissue, levels = tissue_list)
  ) %>%
  # 2. Explicitly complete every combo of sample, allele, and tissue
  complete(
    GTEx_ID,
    allele,
    tissue,
    fill = list(fraction = NA)       # or = 0 if you want zero-height bars
  ) %>% filter(GTEx_ID %in% SoInterest)

# make it NA tissue to display only tissue name instead of NA 
for(x in 1:nrow(df_plot_update)){
  if (is.na(df_plot_update$label[x])){
    df_plot_update$label[x] <-  as.character(df_plot_update$tissue[x])
  }
}


# 1. Build the plot data, adding a 'reads' column
df_plot_reads_Yellow <- df_plot_update %>%
  mutate(
    reads = case_when(
      allele == "frac_A" ~ A,    # take A‐count when showing frac_A
      allele == "frac_G" ~ G,    # take G‐count when showing frac_G
      TRUE               ~ NA_real_
    )
  )


df_plot_reads_Yellow$label <- gsub('Brain_','',df_plot_reads_Yellow$label)

#

# 2. Base bar plot
Yellow_plot <- ggplot(df_plot_reads_Yellow, aes(x = label, y = fraction, fill = allele)) +
  geom_col(position = "dodge", color = "black", na.rm = TRUE) +
  scale_fill_manual(values = c("frac_A" = "#FBB4AE", "frac_G" = "#B3CDE3")) +
  facet_wrap(~ GTEx_ID, scales = "free_x", ncol = 5) +
  scale_x_discrete(drop = FALSE) +
  theme_classic() +
  theme(
    panel.border    = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    legend.title    = element_blank(),
    strip.background= element_rect(fill = "#a1a7ae", color = NA),
    strip.text      = element_text(color = "black", face = "bold"),
    plot.margin     = margin(l = 15, t = 5, r = 5, b = 5)
  ) +
  labs(
    title = "Allelic read counts of rs578776 across brain tissues"
  ) 

# 3. Add read‐count labels above each bar
p2 <- Yellow_plot +
  geom_text(
    aes(label = reads, group = allele),
    position = position_dodge(width = 0.9),
    vjust    = -0.35,
    size     = 3
  ) + scale_y_continuous(expand = expansion(mult = 0),limits = c(0, 1.10)) + 
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 50, unit = "pt"))+
  labs(
    x     = NULL,
    y     = NULL) 


###
yellow_data <- df_plot_reads_Yellow
# bind the data
data_GY <- rbind(green_data,yellow_data)

### 
# 
# Sample_present <- c("GTEX-13NYB",
#                     "GTEX-1B8L1",
#                     "GTEX-1GN73",
#                     "GTEX-11DXY",
#                     "GTEX-13VXU",
#                     "GTEX-N7MS",
#                     "GTEX-WHSE",
#                     "GTEX-R55E")
# 

# New Set of sample nased on Mila on 05220205
Sample_present <- c("GTEX-13VXU",
                    "GTEX-1B8L1",
                    "GTEX-N7MS",
                    "GTEX-WHSE")

data_GY <- data_GY %>% filter(GTEx_ID %in% Sample_present)

### Next is to keep the tissue want to present

brain_to_keep <- data_GY$tissue %>% unique() %>% as.character()
brain_to_keep <- brain_to_keep[c(2,7,8,3,4)]
brain_to_keep <- gsub('Brain_','',brain_to_keep)

data_GY_final <- data_GY %>% filter(label %in% brain_to_keep)

my_label_order <- brain_to_keep

# Plot Final

ggplot(data_GY_final, aes(x = label, y = fraction, fill = allele)) +
  geom_col(position = "dodge", color = "black", na.rm = TRUE,linewidth = 0) +
  scale_fill_manual(
    values = c(
      "frac_A" = "#f9a9a0",   # pastel red
      "frac_G" = "#7f8fdc"    # pastel blue
    )
  ) +
  facet_wrap(~ GTEx_ID, scales = "free_x", ncol = 2) +
  scale_x_discrete(limits = my_label_order, drop = FALSE) +
  theme_classic() +
  theme(
    panel.border    = element_rect(color = "black", fill = NA, size = 0.25),
    axis.text.x          = element_text(angle = 45, hjust = 1,size = 11),
    legend.title         = element_blank(),
    strip.background     = element_rect(fill = NA, color = NA),
    strip.text           = element_text(color = "black",face = "bold"),
    # make facet labels black
  ) +
  labs(
    x     = NULL,
    y     = NULL,
    title = "Allelic ratios of rs578776 across brain tissues in heterozygous GTEx samples"
  ) + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 50, unit = "pt")) +
  geom_text(
    aes(label = reads, group = allele),
    position = position_dodge(width = 0.9),
    vjust    = -0.35,
    size     = 3
  ) + scale_y_continuous(expand = expansion(mult = 0),limits = c(0, 1.10)) + 
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 50, unit = "pt"))+
  labs(
    x     = NULL,
    y     = NULL) 


