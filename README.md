# CHRNA5_project

- Goal:
  
  Downloading GTEx RNA_seq BAMs (v8, hg38 alignment) for detecting rs578776 allele imbalance around CHRNA5/CHRNA3 gene.

- Data:
  
  GTEx v8 brain RNA-seq BAM files, located in T-drive: /Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/GTEx_v8_RNA_seq_BAM_hg38/Brain_BAM/

- Code:
  
  `igvtool_count_for_SNPs_GT.R` : Counting BAM files reads in 3 target SNPs: rs578776, rs10637216, rs660652. the output will generate how many read count per each A,T,C,G,N,DEL,INS

  `run_igvtool.sh`: for running the `igvtool_count_for_SNPs_GT.R` rscript in Biowulf 
