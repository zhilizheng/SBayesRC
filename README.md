## SBayesRC
SBayesRC is a method that incorporates functional genomic annotations with high-density SNPs (> 7 millions) for polygenic prediction. Our manuscript is available [here](https://www.biorxiv.org/content/10.1101/2022.10.12.510418v1). Our method requires summary statistics from GWAS and LD data only. It also estimates functional genetic architecture parameters such as SNP-based heritability and polygenicity associated with the functional annotations. 

This is the R implementation of SBayesRC that produced the results in the manuscript and had a good performance (benchmarked in the manuscript). This version has more experimental flags, but it is user-friendly as demonstrated in the example below.

## Minimal example
A complete script to run SBayesRC without coding in R for Linux or Mac is provided below. Only the top few variables require modification.

```bash
##############################################
# Variables: need to be fixed
ma_file="MA_file"               # GWAS summary in COJO format (the only input)
ld_folder="YOUR_LD_PATH"        # LD reference (download from "Resources")
annot="YOUR_ANNOT_FILE"         # Functional annotation (download from "Resources")
out_prefix="YOUR_OUTPUT_PATH"   # Output prefix, e.g. "./test"
threads=4                       # Number of CPU cores

##############################################
# Code: usually don't need a change in this section
## Note: Flags were documented in the package, use ?function in R to lookup.
## We suggest to run those in multiple jobs (tasks)
export OMP_NUM_THREADS=$threads # Revise the threads

# Tidy: optional step, tidy summary data
## "log2file=TRUE" means the messages will be redirected to a log file 
Rscript -e "SBayesRC::tidy(mafile='$ma_file', LDdir='$ld_folder', \
                  output='${out_prefix}_tidy.ma', log2file=TRUE)"
## Best practice: read the log to check issues in your GWAS summary data.  

# Impute: optional step if your summary data doesn't cover the SNP panel
Rscript -e "SBayesRC::impute(mafile='${out_prefix}_tidy.ma', LDdir='$ld_folder', \
                  output='${out_prefix}_imp.ma', log2file=TRUE)"

# SBayesRC: main function for SBayesRC
Rscript -e "SBayesRC::sbayesrc(mafile='${out_prefix}_imp.ma', LDdir='$ld_folder', \
                  outPrefix='${out_prefix}_sbrc', annot='$annot', log2file=TRUE)"
# Alternative run, SBayesRC without annotation (similar to SBayesR, not recommended)
# Rscript -e "SBayesRC::sbayesrc(mafile='${out_prefix}_imp.ma', LDdir='$ld_folder', \
#                  outPrefix='${out_prefix}_sbrc_noAnnot', log2file=TRUE)"

##############################################
# Polygenic risk score
## Just a toy demo to calculate the polygenic risk score
# genoPrefix="test_chr{CHR}" # {CHR} means multiple genotype file.
## If just one genotype, input the full prefix genoPrefix="test"
# genoCHR="1-22,X" ## means {CHR} expands to 1-22 and X,
## if just one genotype file, input genoCHR=""
# output="test"
# Rscript -e "SBayesRC::prs(weight='${out_prefix}_sbrc.txt', genoPrefix='$genoPrefix', \
#                    out='$output', genoCHR='$genoCHR')"
## test.score.txt is the polygenic risk score

#################################
## SBayesRC multi
## Run each ancestry: summary data and ancestry matched LD,
##    to obtain prs1 and prs2 from the SBayesRC::prs
# prs1="eur.score.txt"
# prs2="eas.score.txt"
# tuneid="tune.id" # two columns FID IID, without header
# pheno="trait.pheno" # three columns FID IID phenotype, without header, only the samples in tuneid are used 
# outPrefix="tuned_eur_eas"
# Rscript -e "SBayesRC::sbrcMulti(prs1='$prs1', prs2='$prs2', \
#             outPrefix='$outPrefix', tuneid='$tuneid', pheno='$pheno')"
## weighted PRS in tuned_eur_eas.score.txt
## Please don't forget to exclude the tuning sample to calculate the prediction accuracy
```
**New**: we also provided docker version for the users who can't install the R package. Note: only x86_64 support for container version currently.
```
# docker image address:  zhiliz/sbayesrc
# We use Apptainer (formerly Singularity) as an example
# The users can also use the container in Docker directly with the mapping volumes
# same as the previsous example
ma_file="MA_file"               # GWAS summary in COJO format (the only input)
ld_folder="YOUR_LD_PATH"        # LD reference (download from "Resources")
annot="YOUR_ANNOT_FILE"         # Functional annotation (download from "Resources")
out_prefix="YOUR_OUTPUT_PATH"   # Output prefix, e.g. "./test"
threads=4                       # Number of CPU cores

# Tidy
apptainer run docker://zhiliz/sbayesrc --ldm-eigen $ld_folder \
    --gwas-summary $ma_file --impute-summary --out ${out_prefix} --threads $threads
# Main: SBayesRC
apptainer run docker://zhiliz/sbayesrc --ldm-eigen $ld_folder \
    --gwas-summary ${out_prefix}.imputed.ma --sbayes RC --annot $annot --out ${out_prefix} \
    --threads $threads

#########################
# PRS
## --geno support both PLINK BED and PGEN format
## test1.txt and test2.txt are weights obtained from SBayesRC
## if only subset SNPs needed in calculation: --extract snp.list
## --keep and --extract are optional flags
# apptainer run docker://zhiliz/sbayesrc --geno test1{CHR} --chr-range 1-22 --score test1.txt --keep keep.id --out prs1
# apptainer run docker://zhiliz/sbayesrc --geno test1{CHR} --chr-range 1-22 --score test2.txt --keep keep.id --out prs2

# SBayesRC-multi
## First: get two PRS from each population with --score
## --keep keep.id is optional flag
# apptainer run docker://zhiliz/sbayesrc  --sbayesrc-multi prs1.score.txt prs2.score.txt --pheno test.pheno --tuneid tune.id --keep keep.id --out weighted_prs

```
### Inputs
Note: We suggested to use the TAB delimited files for GWAS summary statistics and annotation, although from our testing, space and comma also works. However, there were various reports that other delimiter works bad. (For test, use the fread function in data.table R package)

* `ma_file` is the file of GWAS summary statistics with the following COJO format:
```{r, eval=FALSE, indent="   " }
SNP A1 A2 freq b se p N
rs1001 A G 0.8493 0.0024 0.0055 0.6653 129850
rs1002 C G 0.0306 0.0034 0.0115 0.7659 129799
rs1003 A C 0.5128 0.0045 0.0038 0.2319 129830
```
    * SNP: SNP name to match LD reference; A1: effect allele; A2: alternative allele; 
    * freq: frequencies of A1; b: marginal effects (reference to A1); se: sd of marginal effects; 
    * p: p value; N: per-SNP sample size;

* `ld_folder` is a folder contains the eigen-decomposition data for each LD block. We have provided LD information for several ancestries and SNP panel (refer to the "Resources" section for download). If the SNP panel used is significantly different from ours, it is possible to customize the LD (refer to the "Generate LD" section).

* `annot` is the annotation file, with columns being SNP ID, Intercept (a column of one), and annotation values (it's best to use TAB delimited text file). We have provided the baseline model 2.2 (refer to the "Resources section for download). It's easy to customize the annotation in the text format:
```{r, eval=FALSE, indent="   "}
SNP Intercept Coding Conserved CTCF
rs1001 1 0 0 1
rs1002 1 0 1 0
rs1003 1 0 0 0
```

Note: Customized annotation should be provided in the same format with the first two columns as SNP and Intercept (all 1); binary annotation can be input as 0 and 1 (1 means in the functional category); continous annotation can be input as its raw value. 

### Outputs
The outputs share the same name prefix (in the example: ${out\_prefix}\_sbrc):

* SNP weights (PREFIX.txt): First 3 columns are for PGS calculation (SNP: SNP id; A1: effect allele; BETA: joint effect on 0/1/2 genotype scale). Other columns:  PIP: posterior inclusion probability of the variant to be causal; BETAlast: joint effect obtained in last iteration on 0/1/2 scale.

* Running log (PREFIX.log, if set log2file=TRUE): details of running, includes the starting parameters, MCMC status, parameter estimation and estimated runtime.

* Parameter estimation (PREFIX.par): parameter estimations for heritablity (hsq) and number of non-zero effect variants (nnz). The details from MCMC (PREFIX.rds) can be loaded in R (readRDS).

* Functional per-SNP heritability enrichment (PREFIX.AnnoPerSnpHsqEnrichment): per-SNP heritability enrichment for each annotation, averaged value and standard deviation were reported.

* Proportion of causals in a functional category to the mixture distributions (PREFIX.AnnoJointProb, pi):  component 1 zero effec; 2 small effect; 3 medium effect; 4 large effect; 5 very large effect. The pi is the mean value from MCMC iterations.

* MCMC sampling details (PREFIX.mcmcsamples.\*): the value from MCMC iterations after burn-in, at a preset interval (default 10).

### Computational time

Runtime with 4 CPU cores:

* 7 million SNPs with 96 annotation, estimated 8.5 hours, 74GB memory

* 1 million HapMap3 SNPs with 96 annotation, estimated 0.8 hour, 5GB memory


## Install
* Container version. Use container version if your cluster supports Docker or Apptainer (formerly Singularity), image address: `zhiliz/sbayesrc`. It's a good practise to include the version when run, e.g.,  `zhiliz/sbayesrc:0.2.5`, we will maintain the docker version same to R package version; however, the version can be ommited. The image doesn't require any installation, it integrates all depencencies, including R, PLINK and GCTB. If you run multiple jobs with Apptainer, it's good to pre-cache by `apptainer pull docker://zhiliz/sbayesrc`. The default wrapper has the function to SBayesRC, LD generating, PRS, and SBayesRC-Multi, type `apptainer run docker://zhiliz/sbayesrc` for document and examples. The exmaple shown above in this tutorial is also a good start. If you need other parameters that were not provided by our wrapper script, you can pass your own R script directly, by run `apptainer exec docker://zhiliz/sbayesrc Rscript YOUR_SCRIPT`. It's also possible to create your own docker image to include your own script, by *FROM* keyword in Dockerfile and refer to our image. 
  
* Install locally: A valid R is required. 
```r
# Suggest: enable faster backend BLAS for R, e.g. openBlas, MKL
# Run in R to install dependencies
install.packages(c("Rcpp", "data.table", "stringi", "BH",  "RcppEigen"))

# Install SBayesRC package
install.packages("https://github.com/zhilizheng/SBayesRC/releases/download/v0.2.5/SBayesRC_0.2.5.tar.gz",
                 repos=NULL, type="source")

# If R report problem when installing, try alternative version (worse performance and an old version)
## install.packages("https://github.com/zhilizheng/SBayesRC/releases/download/v0.2.0/SBayesRC_0.2.0_comp.tar.gz", repos=NULL, type="source")
```

If you are interested in developing version, try to install by devtools: `devtools::install_github("zhilizheng/SBayesRC")`.

Note: The package was tested under Linux and macOS (x64 and ARM) platform due to availability. The lists of OS we tested: CentOS > 7; Debian > 9; Ubuntu > 20.04; macOS > 11. We find the default R in CentOS7 didn't work good due to an issue in RcppEigen package (gcc 4.8). For users who don't have admin permission and had problem to install the package,  the R from [anaconda](https://www.anaconda.com/products/distribution#Downloads) works great in all available OS (conda install r-base, submit the jobs by full path of /YOUR\_CONDA\_LOCATION/Rscript by `which Rscript`)

## Resources
Download the resources and decompress by "unzip" (.zip) or "tar -xvf" (.tar.xz):

* [Baseline model 2.2](https://sbayes.pctgplots.cloud.edu.au/data/SBayesRC/resources/v2.0/Annotation/annot_baseline2.2.zip): functional annotation information for 8M SNPs from baseline model 2.2 (Reference: Márquez-Luna 2021, DOI: 10.1038/s41467-021-25171-9). 
* LD reference: We provide LD data calculated from different UKB ancestry (EUR, EAS and AFR) in imputed SNPs and HapMap3 SNPs. We suggest to download imputed LD same ancestry as the GWAS summary data. 
    * Imputed SNPs: [EUR](https://sbayes.pctgplots.cloud.edu.au/data/SBayesRC/resources/v2.0/LD/Imputed/ukbEUR_Imputed.zip), [EAS](https://sbayes.pctgplots.cloud.edu.au/data/SBayesRC/resources/v2.0/LD/Imputed/ukbEAS_Imputed.zip) and [AFR](https://sbayes.pctgplots.cloud.edu.au/data/SBayesRC/resources/v2.0/LD/Imputed/ukbAFR_Imputed.zip)
    * HapMap3 SNPs: [EUR](https://sbayes.pctgplots.cloud.edu.au/data/SBayesRC/resources/v2.0/LD/HapMap3/ukbEUR_HM3.zip), [EAS](https://sbayes.pctgplots.cloud.edu.au/data/SBayesRC/resources/v2.0/LD/HapMap3/ukbEAS_HM3.zip) and [AFR](https://sbayes.pctgplots.cloud.edu.au/data/SBayesRC/resources/v2.0/LD/HapMap3/ukbAFR_HM3.zip)
* Example: [A summary data](https://sbayes.pctgplots.cloud.edu.au/data/SBayesRC/resources/v2.0/example/example.tar.xz) as an example to try SBayesRC. You can fill the example code with a HapMap3 LD (ukbEUR\_HM3), Baseline 2.2 annotation, and this summary data. 

If you have downloaded our resources in previous version, don't need to download it again (althrough the format changed). If the links above don't work, a backup here: [Google Drive link](https://drive.google.com/drive/folders/1uxnxDjRJPzo0dTpFnERS5N2NGZX5S-sU?usp=drive_link)

### Summary data and PGS weights
The summary data for the 28 approximately independent UKB traits analysed in the paper can be downloaded here: [summary data](https://sbayes.pctgplots.cloud.edu.au/data/SBayesRC/share/v1.0/summary/), and the corresponding PGS weights can be found here: [PGS weights](https://sbayes.pctgplots.cloud.edu.au/data/SBayesRC/share/v1.0/PGS/). The PGS weights are joint effect estimates derived from ~7 million genome-wide SNPs, therefore it’s important to have matched SNP set between training and validation datasets. This is because if some important SNPs present in the training are missing in the validation, the genetic effects captured by these SNPs will be lost. To maximise the utility of the joint SNP weights, we recommend considering genotype imputation or rerunning SBayesRC with the matched set of SNPs. In addition, note that the PGS weights were estimated using samples of the European ancestry, so they will perform best when applying to individuals of the European ancestry.

## Generate LD
Always try the LD from "Resources" section first. Generate LD If more than 30% SNPs are missing in "tidy" step, or from different SNP panel, different ancestry.

Here is a complete code to generate customized LD:
```bash
##############################################
# Variables: need to be fixed
ma_file="MA_file"                # GWAS summary data in COJO format
genotype="YOUR_GENOTYPE_Prefix"  # genotype prefix as LD reference (PLINK format), with {CHR} to spefify multiple
outDir="YOUR_OUTPUT"             # Output folder that would be created automatically
threads=4                        # Number of CPU cores for eigen decomposition
#---usually don't need change bellow
genoCHR=""                       # If more than 1 genotype file, input range (e.g. "1-22") here.
refblock=""                      # Text file to define LD blocks, by default to use our GRCH37 coordination 
tool="gctb"                      # Command line to run gctb for generating the full LD matrix


##############################################
# Code
# Step1: generate the LD block information and script
# Output $outDir/ldm.info, $outDir/ld.sh, $outDir/snplist/*.snplist
Rscript -e "SBayesRC::LDstep1(mafile='$ma_file', genoPrefix='$genotype', \
            outDir='$outDir', genoCHR='$genoCHR', blockRef='$refblock', log2file=TRUE)"

# Step2: generate each LD matrix for blocks
#  Loop idx from 1 to NUM_BLOCK (591)
#  Submit multiple jobs on your cluster / clouds instead of for loop
#  Input depends on $outDir/ld.sh, $outDir/snplist/$idx.snplist
#  Ouput $outDir/b$idx.ldm.full.info, $outDir/b$idx.ldm.full.bin
for idx in {1..591}; do
    Rscript -e "SBayesRC::LDstep2(outDir='$outDir', blockIndex=$idx, log2file=TRUE)"
done

# Step3: eigen decomposition for each LD block
#  Loop idx from 1 to NUM_BLOCK (591)
#  Submit multiple jobs on your cluster / clouds instead of for loop
#  Input depends on $outDir/ldm.info, $outDir/b$idx.ldm.full.info, $outDir/b$idx.ldm.full.bin
#  Output $outDir/block$block.eigen.bin, $outDir/block$block.eigen.bin.log
export OMP_NUM_THREADS=$threads  # parallel computing supported in this step
for idx in {1..591}; do
    Rscript -e "SBayesRC::LDstep3(outDir='$outDir', blockIndex=$idx, log2file=TRUE)"
done

# Step4: merge LD information
Rscript -e "SBayesRC::LDstep4(outDir='$outDir', log2file=TRUE)"

# Step5: clean if necessary
# Essential for analysis: $outDir/ldm.info, $outDir/snp.info, $outDir/block*.eigen.bin 
# Other files could be removed
# Note: before removing, check all blocks were finished.
```

`refblock` is a text file for block defination with the format:
```{r, eval=FALSE}
Block Chrom StartBP EndBP
1 1 10583 3582736
2 1 3582736 5913893
3 1 5913893 9365199
4 1 9365199 11777841
5 1 11777841 14891511
```
Note: The blocking is based on the coordination, i.e., chromosome (Chrom) and BP [StartBP, EndBP). So the genome build version should be consistent between the block defination and the genotype. 

## News
### v0.2.5
* Enhanced robustness when the meta-ed summary data has error
### v0.2.4
* Fixed bug and add functions to caculate PRS and SBayesRC-multi

### v0.2.3
* Changed the threshold of outlier detection.

### v0.2.2
* Enabled excluding variants by "exclude" flag.

### v0.2.1
* Improved robustness by removing error variants.

### v0.2.0
* Added the tuning step using summary data only
* Performance enhancement
* Added a LD function
* Added a new flag log2file to each function, to output the running log to file
* Updated document with more details (?function in R)

## Citation
Zheng Z, Liu S, Sidorenko J, Wang Y, Lin T, Yengo L, Turley P, Ani A, Wang R, Nolte IM, Snieder H, Lifelines Cohort Study, Yang J, Wray NR, Goddard ME, Visscher PM, Zeng J. (2024) Leveraging functional genomic annotations and genome coverage to improve polygenic prediction of complex traits within and between ancestries. [Nature Genetics, doi:10.1038/s41588-024-01704-y](https://doi.org/10.1038/s41588-024-01704-y)

## Bug report
Report issue to GitHub repository ([https://github.com/zhilizheng/SBayesRC](https://github.com/zhilizheng/SBayesRC)) by using the "Issues" feature.

### Known issues from system:
* Ubuntu 22.04 has a broken openBLAS version 0.3.20.  If you generate the LD and decomposition (or just perform PCA in R for your other analysis), the matrix is incorrect. Our container image is based on Ubuntu 24.04, which has fixed the problem. We tested various LTS version, and found only Ubuntu 22.04 was affected by this openBLAS bug.
* Apptainer will load your local packages in R, if you find the SBayesRC verion is not the one expected, delete the "~/R/x86_64-pc-linux-gnu-library/4.3/SBayesRC" in your local machine (or other libPaths shown in your local R). 
