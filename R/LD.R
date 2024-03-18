# SBayesRC calculate the LD
# License: GPL
# Author: Zhili Zheng <zhilizheng@outlook.com>

#' @title LDstep1: generate the LD command line
#' @usage LDstep1(mafile, genoPrefix, outDir)
#' @param mafile string, path to the summary file
#' @param genoPrefix string, path to the genotype, if multiple genotype files, use {CHR} to indicate the chromosome
#' @param outDir string, folder to save the output (the function will create the folder)
#' @param genoCHR string, if multiple genotype file present, change this to "startCHR-endCHR", e.g., 1-22,X: means from chr1 to chr22, and chrX
#' @param blockRef string, reference window file downloaded from our website or customized, the default 4cM window (in GRch37) will be provided if empty
#' @param tool string, the command line to call the underlying tool (currently gctb), make sure the tool can be run in terminal directly
#' @param log2file boolean, FALSE: display message on terminal (default); TRUE: redirect to an output file
#' @return none, results in the specified output
#' @export
LDstep1 <- function(mafile, genoPrefix, outDir, genoCHR="", blockRef="", tool="gctb", log2file=FALSE){
    output = outDir
    message("Step1: prepare the script to generate the LD matrix")
    logger.begin(output, log2file)
    if(log2file){
        message("Step1: prepare the script to generate the LD matrix")
    }

    chrInfo = expandCHR(genoPrefix, genoCHR)
    bMultiCHR = chrInfo$bMultiCHR
    genoFlag = checkGenoFlag(chrInfo)
    if(genoFlag != " --bfile "){
        stop("Only PLINK BED format supporetd")
    }
 
    ma = mafile
    refPos = blockRef
    if(!file.exists(ma)){
        stop("ma summary data is not there: ", ma)
    }

    suma = fread(ma)
    if(!"SNP" %in% colnames(suma)){
        stop("SNP must exist in ma column's name")
    }else{
        message(nrow(suma), " SNPs in summary data")
    }

    if(refPos == ""){
        refPos = system.file("extdata", "ref4cM_v37.pos", package="SBayesRC")
    }

    if(!file.exists(refPos)){
        stop("Invalid refPos: ", refPos)
    }

    valid_poses = fread(refPos)
    if(!all(colnames(valid_poses) %in% c("Block", "Chrom", "StartBP", "EndBP"))){
        stop("refPos is not a valid posotion file with headers: Block  Chrom  StartBP  EndBP")
    }
    setnames(valid_poses, c("blk", "chr", "start", "end"))

    if(dir.exists(output)){
        stop("The output folder: ", output, " already exists, please change it to another one, or remove the old one")
    }

    dir.create(output)

   
    bims = list()
    idx = 1
    for(geno1 in chrInfo$genos){
        bims[[idx]] = fread(paste0(geno1, ".bim"))
        idx = idx + 1
    }

    bims = rbindlist(bims)

    message(nrow(bims), " SNPs in genotype file")

    com_snp = intersect(bims$V2, suma$SNP)

    idx1 = match(com_snp, bims$V2)
    bims_val = bims[idx1]

    idx2 = match(com_snp, suma$SNP)
    ma_val = suma[idx2]
    
    m = nrow(bims_val)
    if(m == 0){
        stop("No SNPs in common between summmary data and genotype")
    }

    if(all(bims_val$V2 == ma_val$SNP)){
        message(m, " SNPs in common")
    }else{
        stop("Strange genotype file and summary data")
    }

    bA1A1 = (bims_val$V5 == ma_val$A1) & (bims_val$V6 == ma_val$A2)
    bA1A2 = (bims_val$V6 == ma_val$A2) & (bims_val$V5 == ma_val$A1)

    bAll = bA1A1 | bA1A2
    
    message(sum(bAll), " SNPs has consistent alleles")

    bims = bims_val[bAll]
    if(nrow(bims) == 0){
        stop(" No SNPs left after matching the alleles")
    }

    for(curBlock in valid_poses$blk){
        message(" matching ", curBlock)
        curPos = valid_poses[blk==curBlock]
        chr = curPos$chr
        start = curPos$start
        end = curPos$end
        bims[V1==chr & V4>=start & V4<end, blk:=curBlock]
    }


    bims3cM = bims
    bims3cM[, newBlk:=blk]
    blks = unique(bims3cM$newBlk)

    outdir = output

    snpdir = file.path(outdir, "snplist")
    dir.create(snpdir)
    for(curBlock in blks){
        curBim = bims3cM[newBlk == curBlock]
        cat(curBim$V2, file=file.path(snpdir, paste0(curBlock, ".snplist")), sep="\n")
    }

    valid_poses = bims3cM[, .N, .(V1, newBlk)]
    setnames(valid_poses, "V1", "chr")

    valid_poses[, out:=file.path(outdir, paste0("b", newBlk))]
    valid_poses[, cmd:=paste0(tool, genoFlag, genoPrefix, " --chr ", chr, " --extract ", file.path(snpdir, paste0(newBlk, ".snplist")), " --make-full-ldm --out ", out, " &> ", out, ".log")]

    if(bMultiCHR){
        valid_poses[, cmd:=stringi::stri_replace_all_fixed(cmd, "{CHR}", chr)]
    }

    cat(valid_poses$cmd, file=file.path(outdir, "ld.sh"), sep="\n")

    bims[, idx:=0:(nrow(bims)-1)]
    setnames(bims, "newBlk", "Block")
    info = bims[, list(Chrom=.SD[1]$V1, StartSnpIdx = head(.SD$idx, 1), StartSnpID=head(.SD$V2, 1),  EndSnpIdx=tail(.SD$idx, 1), EndSnpID=tail(.SD$V2, 1), NumSnps=.N), by="Block"]

    fwrite(info, file=file.path(outdir, "ldm.info"), sep="\t", quote=F, na="NA")

    message("Done.")
    logger.end()
    if(log2file){
        message("Done.")
    }
}

#' @title LDstep2: run LD
#' @usage LDstep2(outDir, blockIndex)
#' @param outDir string, path to LD folder
#' @param blockIndex numeric, block index to run
#' @param log2file boolean, FALSE: display message on terminal; TRUE: redirect to an output file; default FALSE
#' @return none, results in the specified output
#' @export
LDstep2 = function(outDir, blockIndex, log2file=FALSE){
    message("Step2: generate LD matrix")
    output = outDir
    dt.script = fread(file.path(output, "ld.sh"), head=F)
    curLine = dt.script[blockIndex]
    if(grepl("gctb", tolower(curLine$V1))){
        res = paste0(curLine$V10, ".ldm.full.bin")
        if(file.exists(res)){
            message("Block ", blockIndex, " was done before.")
        }else{
            logger.begin(paste0(curLine$V10, "_R"), log2file)
            if(log2file){
                message("Step2: generate LD matrix")
            }
            t_begin = proc.time()
            cmdline = paste(as.matrix(curLine)[1, ], collapse=" ")
            message(" Run: ", cmdline)
            system(cmdline)
            message("Done.")
            print(proc.time() - t_begin)
            logger.end()
            if(log2file){
                message("Done.")
            }
        }
    }else{
        message("can't recognize the tool")
    }
}

#' @title LDstep3: eigen decomposition on LD matrix
#' @usage LDstep3(outDir, blockIndex)
#' @param outDir string, path to the LD folder
#' @param blockIndex numeric, perform analysis on a specified block number
#' @param thresh numeric, variance cutoff threshold, suggested default 0.995
#' @param log2file boolean, FALSE: display message on terminal (default); TRUE: redirect to an output file
#' @return none, results in the specified output
#' @export
LDstep3 = function(outDir, blockIndex, thresh=0.995, log2file=FALSE){
    message("Step3: perform eigen decomposition on LD matrix")
    message("This step must run after step2")
    output = outDir
    dt = fread(file.path(output, "ldm.info"))

    curLine = blockIndex
    act_block = dt[curLine]$Block
    outFile = file.path(output, paste0("block", act_block, ".eigen.bin"))

    if(file.exists(paste0(outFile))){
        message("The eigen decomposition was done on this block before.")
        message("Block index: ", curLine)
        return
    }

    logger.begin(outFile, log2file)
    if(log2file){
        message("Step3: perform eigen decomposition on LD matrix")
        message("This step must run after step2")
    }
 
    t_begin = proc.time()

    filename = file.path(output, paste0("b", act_block, ".ldm.full"))
    message(" Block index: ", curLine)
    message(" Block number: ", act_block)
    message(" Output: ", outFile)

    message(" Read LD information ", filename)

    info = fread(paste0(filename, ".info"))
    nMarker = nrow(info)

    ldfile = file(paste0(filename, ".bin"), "rb")
    R = readBin(ldfile, n = nMarker * nMarker, what=numeric(0), size=4)
    dim(R) = c(nMarker, nMarker)
    close(ldfile)

    message(" Start eigen decomposition: m = ", nrow(R))

    eig = eigen(R, symmetric=TRUE)

    rm(R)
    gc()

    lambda = eig$values
    m = length(lambda)

    selected = (lambda > 0)
    if(length(rle(selected)$values) > 2){
        stop("strange eigen pattern")
    }
    wholeLambda = lambda[selected]

    message("  Cut to variance threshold: ", thresh)
    lambdaSum = sum(wholeLambda)
    k = which(cumsum(wholeLambda/lambdaSum) >= thresh)[1]

    message(" k: ", k, ", m: ", m)
    cFile = file(outFile, "wb")
    writeBin(m, cFile, size=4)
    writeBin(k, cFile, size=4)
    writeBin(lambdaSum, cFile, size=4)
    writeBin(thresh, cFile, size=4)
    message(" Write binary...")
    writeBin(wholeLambda[1:k], cFile, size=4)

    for(idx in 1:k){
        writeBin(eig$vectors[, idx], cFile, size=4)
    }

    close(cFile)

    message("Done.")
    print(proc.time() - t_begin)
    logger.end()
    if(log2file){
        print(proc.time() - t_begin)
        message("Done.")
    }
}

#' @title LDstep4: Merge the LD information
#' @usage LDstep4(outDir)
#' @param outDir string, path to the LD folder
#' @param log2file boolean, FALSE: display message on terminal; TRUE: redirect to an output file; default value is FALSE
#' @return none, results in the specified output
#' @export
LDstep4 <- function(outDir, log2file=FALSE){
    gDir = outDir
    message("Step4: merge SNP information")
    message("If the LD hasn't been generated, please wait.")
    ldm = fread(file.path(gDir, "ldm.info"), head=TRUE)
    ldm[, file:=file.path(gDir, paste0("b", Block, ".ldm.full.info"))]

    infofile = file.path(gDir, "snp.info")
    if(file.exists(infofile)){
        message("snp.info generated before, the LD seems to be completed, nothing need to be done")
        return
    }

    logger.begin(paste0(infofile), log2file)
    if(log2file){
        message("Step4: merge SNP information")
    }

    infos = list()

    for(idx in 1:nrow(ldm)){
        curInfoFile = ldm$file[idx]
        if(file.exists(curInfoFile)){
            infos[[idx]] = fread(curInfoFile)
            infos[[idx]][, Block:=idx]
        }else{
            stop("Info file: ", curInfoFile, " can't be read, please check step2")
        }
    }

    info = rbindlist(infos)
    info[, Index:=0:(nrow(info)-1)]

    #info_nodup = info[!duplicated(ID)]

    fwrite(info[, .(Chrom,ID,Index,GenPos,PhysPos,A1,A2,A1Freq,N, Block)], file=infofile, sep="\t")
    message("snp.info has been generated")
    message("To save the storage, you could remove the ", file.path(gDir, "*.ldm.bin"), " after everything cheked OK")
    message("The LD is completed!")
    logger.end()
    if(log2file){
        message("The LD is completed!")
    }
}
