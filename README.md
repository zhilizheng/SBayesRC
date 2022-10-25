# SBayesRC
This is the R implementation of SBayesRC. SBayesRC extends SBayesR to incorporate functional annotations and allows joint analysis of all common SNPs in the genome. Similar to SBayesR, SBayesRC only requires summary statistics from GWAS (i.e., marginal SNP effect estimates, standard errors, and GWAS sample size) and LD correlations from a reference sample as input data. In addition to joint SNP effect estimates for deriving PGS, SBayesRC also generates the fine mapping Bayesian statistics of posterior inclusion probabilities (PIP) for SNPs as measures of trait associations, and estimates of functional genetic architecture parameters such as SNP-based heritability and polygenicity associated with the functional annotations. 

# Install
```
## Install by devtools
devtools::install_github("zhilizheng/SBayesRC")
## package under reformating of the output, the functional related outputs will be enhanced soon.
```

## Resources
Download the resources and decompress by unzip:
* [Baseline model 2.2](https://broadinstitute-my.sharepoint.com/:u:/g/personal/zhengzhi_broadinstitute_org/EaFPvEwSuUJIvNRD5qWnnj0BhsM08JWos7DH1Z4aMzZUfg?e=tDCcOb) (decompressed to ANNOT_FILE): functional annotation information for 8M SNPs from baseline model 2.2 ([Márquez-Luna 2021](https://doi.org/10.1038/s41467-021-25171-9))
* LD refernce (decompressed to LD_PATH): UKB EUR, UKB EAS ([link](https://broadinstitute-my.sharepoint.com/:f:/g/personal/zhengzhi_broadinstitute_org/EmN4lg1qO6xEoRQwEh_I148BZFKG3ndcbaDCIWVirahLLw?e=acWEjq)).

We will provided functions to generate LD from customized genotypes soon. 

# How to run
## Tidy the summary
Match the summary data with LD reference and QC.
Tidy functions will do those automatically:
* Check the SNPs with LD reference
* Check the consistency of alleles (A1, A2)
* Filter NA or inf values in the summary data
* Check the difference of allele frequency in summary data to LD reference smaller 0.2.
* Check per-SNP sample size range in mean +- 3SD

```
SBayesRC::tidy(mafile, LD_PATH, output_FILE)
```

Parameters:

mafile: the path to summary statitics followed [COJO format](https://yanglab.westlake.edu.cn/software/gcta/#COJO), example:
```
SNP A1 A2 freq b se p N 
rs1001 A G 0.8493 0.0024 0.0055 0.6653 129850 
rs1002 C G 0.0306 0.0034 0.0115 0.7659 129799 
rs1003 A C 0.5128 0.0045 0.0038 0.2319 129830
...
```
LD_PATH: the path to the downloaded and decompressed LD reference.

output_FILE: the output path.

## Imputation
Impute the missing SNPs from the LD reference. This is an optional step if the summary data covers the LD reference completely. 

```SBayesRC::impute(mafile, LD_PATH, output_FILE)```

Parameters are same as above. mafile is the output of the tidy steps. 


## SBayesRC
Run the method with or without functional annotation. SBayesRC with functional annotations is preferred most of the time. 

```
# With annotation
SBayesRC::sbayesrc(mafile, LD_PATH, output_FILE, fileAnnot=annotation_file)
# Alternative: without annotation
SBayesRC::sbayesrc(mafile, LD_PATH, output_FILE)
```

fileAnnot is the path to annotation file. Other parameters are same above. 

# Example code
This is an example to generate the SNP weights for a summary data (in bash). 

```
ld_folder="YOUR_PATH"    # LD reference path
annot="YOUR_ANNOT_FILE"  # annotation
ma_file="MA_file"        # summary file name
out_prefix="YOUR_PATH"   # output

# Tidy
Rscript -e "SBayesRC::tidy('$ma_file', '$ld_folder', '${out_prefix}_tidy.ma')"

# Impute
Rscript -e "SBayesRC::impute('${out_prefix}_tidy.ma', '$ld_folder', '${out_prefix}_imp.ma')"

# SBayesRC
Rscript -e "SBayesRC::sbayesrc('${out_prefix}_imp.ma', '$ld_folder', '${out_prefix}_sbrc', fileAnnot='$annot')"

```

The outputs are:

* SNP weights (in 0/1/2 genotype scale): ${out_prefix}_sbrc.txt
* Running logs:  ${out_prefix}_sbrc.log
* Parameter estimation: ${out_prefix}_sbrc.rds 
* Other files: Functional heritability enrichments, functional genetic architecture (documents update soon). 
# Citation
Zheng Z, Liu S, Sidorenko, J, Yengo L, Turley P, Ani A, Wang R, Nolt I, Snieder H, Lifelines Cohort Study, Yang J, Wray NR, Goddard ME, Visscher PM, Zeng J. (2022) Leveraging functional genomic annotations and genome coverage to improve polygenic prediction of complex traits within and between ancestries. bioRxiv 2022.10.12.510418; doi: https://doi.org/10.1101/2022.10.12.510418
