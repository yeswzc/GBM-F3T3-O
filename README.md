### R code used to predict UMAP embeding for GBMs with FGFR3-TACC3 fusion.

This code is used to test if tumor samples can be embeded with recently identitified FGFR3-TACC3 fusion (F3T3) positive methylation outlier group (GBM-F3T3-O).

![result](R/result.jpg)


Note: Because minfi generated some NaN beta values due to Methy = Umethy = 0, I have to set the NaN beta values to 0.5 as a neutral methylation value.


Referemce: under submission
