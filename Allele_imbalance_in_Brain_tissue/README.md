## Goal: For testing imbalance
- File: GTEx brain RNA-seq BAM files
- SNPs:
- **Notes** for rs10637216, the location for detecting insertion is located in `78595865`, due to the short-read cDNA alignment diffcult to align same repeat sequence for this SNP 

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


------------


- GTEx sample select and tissue

```
GTEX-13VXU
GTEX-1B8L1
GTEX-N7MS
GTEX-WHSE

selected tissue order:
All the xxx_ganglia, Cerebellar_Hemixxx, Cerebellum

No Hypothalamus
```
