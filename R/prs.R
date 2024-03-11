# SBayesRC: PRS function

# License: GPL
# Author: Zhili Zheng <zhilizheng@outlook.com>

#' @title caculate PRS
#' @usage prs(weight, genoPrefix, genoCHR, outprefix)
#' @param weight string, path to the SBayesRC genereated weight (.txt)
#' @param genoPrefix string, path to the genotype, if multiple genotype files, use {CHR} to indicate the chromosome
#' @param genoCHR string, if multiple genotype file present, change this to "startCHR-endCHR", e.g., 1-22,X: means from chr1 to chr22, and chrX
#' @param outPrefix string, output path
#' @param snplist string, path to variant list to be extracted, default "": all variants
#' @param keepid string, path to sample id (FID IID) to be kept, default "": all samples
#' @param tool string, path to plink2 alpha; default plink2
#' @param log2file boolean, FALSE: display message on terminal; TRUE: redirect to an output file; default value is FALSE
#' @return none, results in the specified output
#' @export
prs <- function(weight, genoPrefix, genoCHR, outPrefix, snplist="", keepid="", scoreFlag="1 2 3 header center", tool="plink2", log2file=FALSE){
    message("Calculating PRS from weights and genotype...")

    logfile = paste0(outPrefix, ".log")
    logger.begin(logfile, log2file)
    if(log2file){
        message("Calculating PRS from weights and genotype...")
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

    # check if multi CHR
    bMultiCHR = FALSE
    chrs = c()
    if(genoCHR != ""){
        genoCHRnosp = gsub(" ", "", genoCHR, fixed=TRUE)
        genoStrs = c(stringi::stri_split_fixed(genoCHR, ",", simplify=TRUE))
        for(genoStr in genoStrs){
            chrStr = c(stringi::stri_split_fixed(genoCHR, "-", simplify=TRUE))
            if(length(chrStr) == 1){
                chrs = c(chrs, chrStr)
            }else if(length(chrStr) == 2){
                startchr = as.numeric(chrStr[1])
                endchr = as.numeric(chrStr[2])
                chrs = c(chrs, startchr:endchr)
            }else{
                stop("genoCHR doesn't recognize: ", genoCHR)
            }
        }

        if(!grepl("\\{CHR\\}", genoPrefix)){
            stop("genoCHR has the start and end value, however there is no {CHR} in genoPrefix string")
        }
        bMultiCHR = TRUE
    }
    
    if(!bMultiCHR){
        chrs = c("")
    }
 
    # check the first genotype
    geno1 = genoPrefix
    if(bMultiCHR){
        geno1 = gsub("{CHR}", chrs[1], genoPrefix, fixed=TRUE)
    }

    bBfile = TRUE
    exts = c(".bed", ".bim", ".fam")
    for(ext in exts){
        if(file.exists(paste0(geno1, ext))){
            bBfile = bBfile & TRUE 
        }else{
            bBfile = bBfile & FALSE
        }
    }

    bPfile = TRUE
    exts = c(".pgen", ".pvar", ".psam")
    for(ext in exts){
        if(file.exists(paste0(geno1, ext))){
            bPfile = bPfile & TRUE 
        }else{
            bPfile = bPfile & FALSE
        }
    }

    genoFlag = ""
    if(bBfile & bPfile){
        warning("Found two set of genotype, use PLINK BED format")
        genoFlag = " --bfile "
    }else if(bBfile){
        genoFlag = " --bfile "
        message("Found PLINK BED format")
    }else if(bPfile){
        genoFlag = " --pfile "
        message("Found PLINK PGEN format")
    }

    if(genoFlag == ""){
        stop("No genotype found")
    }

    threads = Sys.getenv("OMP_NUM_THREADS")
    if(threads != ""){
        threads =paste0(" --threads ", threads) 
    }

    numMarker = 0
    sumScore = 0
    for(chr in chrs){
        message("Processing genotype ", chr)
        refGeno = gsub("{CHR}", chr, genoPrefix, fixed=TRUE)
        system(paste0(tool, genoFlag, refGeno, " --extract ", outFile, ".snplist",  " --score ", weight, " ", scoreFlag, keep, threads, " --memory 4096 --out ", outFile))
        infile = paste0(outFile, ".sscore")
        dt.in = fread(infile, head=TRUE)
        numMarker = numMarker + dt.in$NMISS_ALLELE_CT
        sumScore = sumScore + dt.in$SCORE1_AVG * dt.in$NMISS_ALLELE_CT
    }
    outDT = data.table(FID=dt.in[["#FID"]], IID=dt.in$IID, NMISS_ALLELE_CT=numMarker, SCORE=sumScore)
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
