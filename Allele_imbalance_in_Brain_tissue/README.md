## Goal: For testing imbalance
- File: GTEx brain RNA-seq BAM files
- SNPs:

```
POS_hg38	rsID	REF, R	ALT, A
78595892	rs10637216	A	ACCCC
78596058	rs578776	G	A
78595490	rs660652	G	A

```

- methods:

Scoring read count for each nucleotid using `igvtoolcount`, and do the ratio of each SNPs.
Boxplot and summary table to show the result

- results/script location:
script name `GTEx_3SNPs_GT.R`
`~/Desktop/chr15_CHRNA5_bladder_project/GTEx_v8_Brain_selected_by_Oscar/Allelic_imbalance_analysis`
will be in T-drive

