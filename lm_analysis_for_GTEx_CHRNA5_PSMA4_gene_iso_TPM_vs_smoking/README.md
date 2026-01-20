
#### Goal: lm analysis


#### Data: Target dataset: GTEx CHRNA5 and PSMA4 whole gene/isoform TPM
- Location:
- The new files are in folder in T-drive `/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/CHRNA5_project/updated_version_07232024`
- Target SNPs: 33 SNPs in the region

### code: 
- Analysis model
    The variable are `TPM`, `SNPs_GT` as coded as numeric dosage SNPs GT as `0,1,2` . 
      Adjust for:
      numeric `AGE` and factor `SEX`
      numeric `smoking_coded` for `0,1,2` as for no, former, current smoker.
      numeric`Cig_day_Coded` for `0,1,2,3,4,5` as smoke intensity.

```
fmla_all <- as.formula(paste(j, "~", i, " + SEX + AGE + smoking_coded + ",i," * smoking_coded"))

fmla_inter <- as.formula(paste(j, "~", i, " * smoking_coded + SEX + AGE"))

fmla <- as.formula(paste(j, "~", i, " + SEX + AGE"))
          
#### add Cig_per_day
          
fmla_all.Cig <- as.formula(paste(j, "~", i, " + SEX + AGE + smoking_coded + ",i," * smoking_coded  + Cig_day_coded"))

fmla_inter.Cig <- as.formula(paste(j, "~", i, " * smoking_coded + SEX + AGE + Cig_day_coded"))

fmla.Cig <- as.formula(paste(j, "~", i, " + SEX + AGE + Cig_day_coded"))
              

```


