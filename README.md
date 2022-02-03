### R code used to predict UMAP embeding for GBMs with FGFR3-TACC3 fusion.

This code is used to test if tumor samples can be embeded with recently identitified FGFR3-TACC3 fusion (F3T3) positive methylation outlier group (GBM-F3T3-O).

A reference UMAP was trained and can be applied to newly Illumina Infinium MethylationEPIC array profiled GBM tumors. The UMAP was trained from 24 non trivial principal components of 32,000 most variable probes (standard deviation) based on 1000 times permutations. As the sample size for different tumors is very imblanced, using 10,000 or 20,000 probes cannot find the GBM-F3T3-O group well.

![result](R/result.jpg) The figure showed case BL97 (red dot with idat file name indicated) embeded closer to GBM-F3T3-O. 
Note: Because minfi generated some NaN beta values due to Methy = Umethy = 0, I have to set the NaN beta values to 0.5 as a neutral methylation value. Processing multiple samples together will most likely solve the issue.

The reference raw methylation data was process using R minfi package preprocessFunnorm() function. Probes were filtered as recommended by Zhou et al., NAR, 2017.

<br/>
<br/>
<br/>
<br/>


### Referemces:

This study: under submission

Minfi:
Fortin J, Labbe A, Lemire M, Zanke BW, Hudson TJ, Fertig EJ, Greenwood CM, Hansen KD (2014). “Functional normalization of 450k methylation array data improves replication in large cancer studies.” Genome Biology, 15(12), 503. doi: 10.1186/s13059-014-0503-2.
Fortin J, Triche TJ, Hansen KD (2017). “Preprocessing, normalization and integration of the Illumina HumanMethylationEPIC array with minfi.” Bioinformatics, 33(4). doi: 10.1093/bioinformatics/btw691.

Probe filering: 
Zhou W, Laird PW and Shen H, Comprehensive characterization, annotation and innovative use of Infinium DNA Methylation BeadChip probes, Nucleic Acids Research 2017
