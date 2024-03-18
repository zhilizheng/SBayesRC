# SBayesRC: PRS function

# License: GPL
# Author: Zhili Zheng <zhilizheng@outlook.com>

#' @title caculate PRS
#' @usage prs(weight, genoPrefix, outprefix)
#' @param weight string, path to the SBayesRC genereated weight (.txt)
#' @param genoPrefix string, path to the genotype, if multiple genotype files, use {CHR} to indicate the chromosome
#' @param outPrefix string, output path
#' @param genoCHR string, if multiple genotype file present, change this to "startCHR-endCHR", e.g., 1-22,X: means from chr1 to chr22, and chrX
#' @param snplist string, path to variant list to be extracted, default "": all variants
#' @param keepid string, path to sample id (FID IID) to be kept, default "": all samples
#' @param scoreFlag string, flag to produce the PRS in PLINK2, default "1 2 3 header center"
#' @param tool string, path to plink2 alpha; default plink2. Please note plink1.9 isn't supported (for pgen format)
#' @param log2file boolean, FALSE: display message on terminal; TRUE: redirect to an output file; default value is FALSE
#' @return none, results in the specified output
#' @export
prs <- function(weight, genoPrefix, outPrefix, genoCHR="", snplist="", keepid="", scoreFlag="1 2 3 header center", tool="plink2", log2file=FALSE){
    message("Calculating PRS from weights and genotype...")

    logfile = paste0(outPrefix, ".log")
    logger.begin(logfile, log2file)
    if(log2file){
        message("Calculating PRS from weights and genotype...")
    } 

    # check if multi CHR
    chrInfo = expandCHR(genoPrefix, genoCHR)
    bMultiCHR = chrInfo$bMultiCHR
    chrs = chrInfo$CHRs
 
    # check the first genotype
    genoFlag = checkGenoFlag(chrInfo)

    if(genoFlag == ""){
        stop("No genotype found")
    }


    if(!file.exists(weight)){
        stop("The SNP weights are not available: ", weight)
    }
    dt = fread(weight)

    snpcol = as.numeric(stringi::stri_split_fixed(scoreFlag, " ", simplify=TRUE)[, 1])
    outFile = paste0(outPrefix, "_tmp")
    snpv = dt[[snpcol]]
    if(snplist != ""){
        if(!file.exists(snplist)){
            stop("snplist doesn't exist: ", snplist)
        }
        snp.extract = fread(snplist, head=FALSE)$V1
        snpv = intersect(snpv, snp.extract)
    }

    if(length(snpv) != 0){
        cat(snpv, file=paste0(outFile, ".snplist"), sep="\n")
    }else{
        stop("No variants exist in weight file or after extracting the snplist")
    }

    keep = ""
    if(keepid != ""){
        if(!file.exists(keepid)){
            stop("keepid not exists: ", keepid)
        }
        keep = paste0(" --keep ", keepid)
    }

    threads = Sys.getenv("OMP_NUM_THREADS")
    if(threads != ""){
        threads =paste0(" --threads ", threads) 
    }

    numMarker = 0
    sumScore = 0
    for(refGeno in chrInfo$genos){
        message("Processing genotype ", refGeno)
        system(paste0(tool, genoFlag, refGeno, " --extract ", outFile, ".snplist",  " --score ", weight, " ", scoreFlag, keep, threads, " --memory 4096 --out ", outFile))
        infile = paste0(outFile, ".sscore")
        dt.in = fread(infile, head=TRUE)
        all_cols = colnames(dt.in)
        colCT = all_cols[grepl("ALLELE_CT", all_cols)]
        numMarker = numMarker + dt.in[[colCT]]
        sumScore = sumScore + dt.in$SCORE1_AVG * dt.in[[colCT]]
    }
    outDT = data.table(FID=dt.in[["#FID"]], IID=dt.in$IID, ALLELE_CT=numMarker, SCORE=sumScore)
    fwrite(outDT, file=paste0(outPrefix, ".score.txt"), quote=F, sep="\t", na="NA")
    message("Done")

    logger.end()
    if(log2file){
        message("Done")
    }

    file.remove(infile)
    file.remove(paste0(outFile, ".log"))
    file.remove(paste0(outFile, ".snplist"))
}
