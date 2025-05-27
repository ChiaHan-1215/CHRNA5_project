- Goal:

  - GTEx selected brain sample for detecting allele imbalance

  - Scoring read count for each nucleotid using `igvtoolcount`, and do the ratio of each SNPs. DO boxplot and summary table to show the result

  - SNPs:

  ```
  # Notes for rs10637216, the location for detecting insertion is located in `78595865`, due to the short-read cDNA alignment diffcult to align same repeat sequence for   this SNP 

  POS_hg38	rsID	REF, R	ALT, A
  78595892	rs10637216	A	ACCCC
  78596058	rs578776	G	A
  78595490	rs660652	G	A

  ```


- Data: GTEx brain RNA-seq BAM files mentioned at the main README

  - Selected sample ID and tissues

  ```
  GTEX-13VXU
  GTEX-1B8L1
  GTEX-N7MS
  GTEX-WHSE
  selected tissue order:
  All the xxx_ganglia, Cerebellar_Hemixxx, Cerebellum
  No Hypothalamus
  ```

- Code:

  - `GTEx_3SNPs_GT_plot_show_reads.R` : We focused on rs578776, Generating the data summary table of rs578776 read count and showed the result in barplot.
  - script also located in the T-drive `~/Desktop/chr15_CHRNA5_bladder_project/GTEx_v8_Brain_selected_by_Oscar/Allelic_imbalance_analysis`



  
