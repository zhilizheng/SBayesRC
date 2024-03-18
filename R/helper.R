# Some helper functions
# License: GPL
# Author: Zhili Zheng <zhilizheng@outlook.com>

logger.begin <- function(outPrefix, log2file){
    bLoggerTOKEN <<- log2file
    if(bLoggerTOKEN){
        zzLOGGER <<- file(paste0(outPrefix, ".log"), "wt")
        message("The messages are redirected to ", outPrefix, ".log")
        message("  No output here...")
        sink(zzLOGGER)
        sink(zzLOGGER, type="message")
    }
}

logger.end <- function(){
    if(bLoggerTOKEN){
        sink()
        sink(type="message")
        close(zzLOGGER) 
        message("  Output reactivated.")
    }
}

checkGenoFlag <- function(genoInfo){
    geno1 = genoInfo$genos[1]
    bBfile = TRUE
    bexts = c(".bed", ".bim", ".fam")
    for(ext in bexts){
        if(file.exists(paste0(geno1, ext))){
            bBfile = bBfile & TRUE 
        }else{
            bBfile = bBfile & FALSE
        }
    }

    bPfile = TRUE
    pexts = c(".pgen", ".pvar", ".psam")
    for(ext in pexts){
        if(file.exists(paste0(geno1, ext))){
            bPfile = bPfile & TRUE 
        }else{
            bPfile = bPfile & FALSE
        }
    }

    genoFlag = ""
    checkExts = c()
    if(bBfile & bPfile){
        warning("Found two set of genotype, use PLINK BED format")
        genoFlag = " --bfile "
        checkExts = bexts
    }else if(bBfile){
        genoFlag = " --bfile "
        message("Found PLINK BED format")
        checkExts = bexts
    }else if(bPfile){
        genoFlag = " --pfile "
        message("Found PLINK PGEN format")
        checkExts = pexts
    }else{
        stop("None genotype file found")
    }

    miss = c()
    for(geno1 in genoInfo$genos){
        for(ext in checkExts){
            geno1full = paste0(geno1, ext)
            if(!file.exists(geno1full)){
                miss = c(miss, geno1full)
            }
        }
    }
    if(length(miss) != 0){
        mistr = paste0(miss, collapse="\n") 
        stop("Error: some genotype file specified doesn't exist: \n", mistr)
    }
    
    return(genoFlag)
}

expandCHR <- function(genoPrefix, genoCHR){
     # check if multi CHR
    bMultiCHR = FALSE
    chrs = c()
    if(genoCHR != ""){
        genoCHRnosp = gsub(" ", "", genoCHR, fixed=TRUE)
        genoStrs = c(stringi::stri_split_fixed(genoCHRnosp, ",", simplify=TRUE))
        for(genoStr in genoStrs){
            chrStr = c(stringi::stri_split_fixed(genoStr, "-", simplify=TRUE))
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
    genos = c()
    for(curCHR in chrs){
        refGeno = gsub("{CHR}", curCHR, genoPrefix, fixed=TRUE)
        genos = c(genos, refGeno)
    }


    return(list(bMultiCHR=bMultiCHR, CHRs=chrs, genos=genos))
}
 
