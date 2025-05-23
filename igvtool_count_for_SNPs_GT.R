library(dplyr)

ls.rsnumber <- c('rs660652','rs10637216','rs578776')
chromosme_location <- c('chr15:78595490-78595490','chr15:78595895-78595895','chr15:78596058-78596058')

# Set working directory
setwd('/data/leec20/GTEX_and_gdc_download/GTEx_v8/GTEx_v8_Brain_RNA-seq_BAMs_selected_by_Oscar_for_CHRNA5_project/ID_coded_in_yellow/')

# Get list of BAM files
bamlocation <- list.files(full.names = FALSE, pattern = '\\.bam$')

# Function to count reads at specific positions
autocount <- function(rs_number, chr_location) {
  result <- data.frame()
  
  for (i in bamlocation) {
    tryCatch({
      # Run IGVTools count
      system(paste0('igvtools count -w 1 --bases --query ', chr_location, " ", i, " ",
                    "/data/leec20/igvtool_count_GT_rstudio_biowulf_env/gt_res/", i, ".wig ",
                    "/fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"))
      
      # Read the result wig file
      wig_path <- paste0('/data/leec20/igvtool_count_GT_rstudio_biowulf_env/gt_res/', i, '.wig')
      if (file.size(wig_path) > 0) {
        ori <- read.delim(wig_path, header = FALSE, skip = 3)
        colnames(ori) <- c('pos', 'A', 'C', 'G', 'T', 'N', 'DEL', 'INS')
        ori$bam_ID <- gsub('\\.bam', "", i)
        ori <- ori[, c(9, 1:8)]
        ori <- ori[grep(gsub('chr15:[0-9]+-', '', chr_location), ori$pos), ]
        result <- rbind(result, ori)
      }
    }, error = function(e) {
      cat("ERROR with BAM file", i, ":", conditionMessage(e), "\n")
    })
  }
  
  # Write the result to a file
  output_path <- paste0('/data/leec20/igvtool_count_GT_rstudio_biowulf_env/gt_res/', rs_number, '_Yellow_result.txt')
  write.table(result, output_path, quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)
  
  # Remove temporary wig files
  system('rm -f /data/leec20/igvtool_count_GT_rstudio_biowulf_env/gt_res/*.wig')
}

# Run the function with mapply
mapply(autocount, ls.rsnumber, chromosme_location)
