# SBayesRC
This is the R implementation of SBayesRC. SBayesRC extends SBayesR to incorporate functional annotations and allows joint analysis of all common SNPs in the genome. Similar to SBayesR, SBayesRC only requires summary statistics from GWAS (i.e., marginal SNP effect estimates, standard errors, and GWAS sample size) and LD correlations from a reference sample as input data. In addition to joint SNP effect estimates for deriving PGS, SBayesRC also generates the fine mapping Bayesian statistics of posterior inclusion probabilities (PIP) for SNPs as measures of trait associations, and estimates of functional genetic architecture parameters such as SNP-based heritability and polygenicity associated with the functional annotations. 

# Install
```
# Install by devtools
devtools::install_github("zhilizheng/SBayesRC")

# If you find difficulties to install from devtools
# Alternative: install in R by downloading the tar.gz from releases
install.packages(c("Rcpp", "data.table", "BH",  "RcppArmadillo", "RcppEigen"))
install.packages("PATH_DOWNLOAD_SBayesRC_version.tar.gz", repos=NULL, type="source")
```

## Resources
Download the resources and decompress by unzip:
* [Baseline model 2.2](https://drive.google.com/drive/folders/1cq364c50vMw1inJBTkeW7ynwyf2W6WIP?usp=sharing) (unzip to ANNOT_FILE): functional annotation information for 8M SNPs from baseline model 2.2 ([Márquez-Luna 2021](https://doi.org/10.1038/s41467-021-25171-9)).  Customized annotation should be provided in the same format with the first two columns as SNP and Intercept (all 1); binary annotation input as 0 and 1 (in the functional category); continous annotation could be input as its raw value. 
* LD refernce (unzip to LD_PATH): [UKB Imputed](https://drive.google.com/drive/folders/1ZTYv_qlbb1EO70VVSSQFaEP9zH7c9KHt?usp=sharing), [UKB HapMap3](https://drive.google.com/drive/folders/1YTnw1cY-TZfAnLjuwF6wsVHdM4DOXA_G?usp=sharing). We suggest to download imputed LD same ancestry as your GWAS summary data. We will integrate functions to generate LD from customized genotypes soon. 

# How to run
The complete code can be copied from "Example code" section below directly, this section could be skipped.
## Tidy the GWAS summary
Match the GWAS summary data with LD reference and perform further QC.
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

mafile: the path to GWAS summary statitics in [COJO format](https://yanglab.westlake.edu.cn/software/gcta/#COJO). The GWAS summary data is the only essential input for SBayesRC in text format, header line is neccessary, example:
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
Run SBayesRC with or without functional annotation. SBayesRC with functional annotations is preferred most of the time. 

```
# With annotation
SBayesRC::sbayesrc(mafile, LD_PATH, output_FILE, fileAnnot=ANNOT_FILE)
# Alternative: without annotation
SBayesRC::sbayesrc(mafile, LD_PATH, output_FILE)
```

fileAnnot is the path to annotation file. Other parameters are same above. 

# Example code
This is an complete example to run SBayesRC for a GWAS summary data (in bash). 

```
ma_file="MA_file"        # GWAS summary data in COJO text format (the only input needed)
out_prefix="YOUR_OUTPUT_PATH"   # output prefix, e.g. "./test"
ld_folder="YOUR_LD_PATH"    # LD reference path (download and unzip from Resources section)
annot="YOUR_ANNOT_FILE"  # Functional annotation (download and unzip from Resources section)

# Tidy
Rscript -e "SBayesRC::tidy('$ma_file', '$ld_folder', '${out_prefix}_tidy.ma')"

# Impute
Rscript -e "SBayesRC::impute('${out_prefix}_tidy.ma', '$ld_folder', '${out_prefix}_imp.ma')"

# SBayesRC
Rscript -e "SBayesRC::sbayesrc('${out_prefix}_imp.ma', '$ld_folder', '${out_prefix}_sbrc', fileAnnot='$annot')"

```

The outputs are:

* SNP weights (${out_prefix}_sbrc.txt).  First 3 columns are for PGS calculation (SNP: SNP id; A1: effect allele; BETA: joint effect on 0/1/2 genotype scale), it can be the input for other tools directly (e.g. PLINK). Other columns:  PIP: posterior inclusion probability of the variant to be causal from MCMC iterations after burn-in; BETAlast: joint effect obtained in last iteration on 0/1/2 scale.
* Running logs(${out_prefix}_sbrc.log): the estimated runtime is included; provide the logs if you would like to report a problem.  
* Parameter estimation (${out_prefix}_sbrc.par): parameter estimations for heritablity (hsq) and number of non-zero effect variants (nnz). Details (${out_prefix}_sbrc.rds) could be loaded by R readRDS, it's a list with names indicates the variables estimated.
* Functional per-SNP heritability enrichments (${out_prefix}.vg.enrich.qt):  Heritability enrichment for each annotation in the MCMC iterations after burn-in. Each row is an output of enrichment (output per 10 iterations); each column is the enrichment indicated in the header line.
* Proportion of variants' effects in a functional annotation belonging to each of the mixture distributions (${out_prefix}.annoJointProb${comp}):  ${comp} 0 zero effec; 1 small effect; 2 median effect; 3 large effect; 4 very large effect. Each row is an output from MCMC (output per 10 iterations); each column is the functional category indicated in the header line.
* Documents for other outputs will be provided in the future. 

# Reference
Zheng Z, Liu S, Sidorenko, J, Yengo L, Turley P, Ani A, Wang R, Nolt I, Snieder H, Lifelines Cohort Study, Yang J, Wray NR, Goddard ME, Visscher PM, Zeng J. (2022) Leveraging functional genomic annotations and genome coverage to improve polygenic prediction of complex traits within and between ancestries. bioRxiv 2022.10.12.510418; doi: https://doi.org/10.1101/2022.10.12.510418

# Bug report
Report an issue by "new issue" tab in github.
