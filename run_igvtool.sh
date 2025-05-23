#!/bin/bash

# for running igvtoolcount with Rscript in Biowulf

# module load tools
ml R igvtools/2.17.4 

# run the script
Rscript igvtool_count_for_SNPs_GT.R && echo "done"


