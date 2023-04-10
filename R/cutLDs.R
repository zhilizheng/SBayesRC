# SBayesRC cut the LD to another cutoff value
# License: GPL
# Author: Zhili Zheng <zhilizheng@outlook.com>


#' @title cut the LD to a different variance threshold to save the storage
#' @usage cutLD(inputDir, outDir)
#' @param inputDir string, path for input LD
#' @param outDir string, path for output LD
#' @param thresh numeric, new threshold to cut, default 0.995
#' @param log2file boolean, FALSE: display message on terminal; TRUE: redirect to an output file; default FALSE
#' @return none, results in the specified output
#' @export
cutLD = function(inputDir, outDir, thresh=1, log2file=FALSE){
    outputDir = outDir
    message("Cut the LD matrix to eigen cutoff ", thresh)
    message(" From: ",inputDir) 
    message(" To: ", outputDir)
    ldmInfoFile = file.path(inputDir, "ldm.info")
    if(file.exists(ldmInfoFile)){
        ldm = fread(ldmInfoFile)
    }else{
        stop("not a valid input folder: ", inputDir)
    }

    if("block" %in% colnames(ldm)){
        setnames(ldm, "block", "Block")
    }
    blocks = ldm$Block
    
    if(file.exists(outputDir)){
        stop("The output folder is there already, can't overwrite")
    }

    logger.begin(outputDir, log2file)
    if(log2file){
        message("Cut the LD matrix to eigen variance cutoff ", thresh)
        message(" From: ",inputDir) 
        message(" To: ", outputDir)
    }

    info = getLDPrefix(inputDir)

    if(info$type < 1){
        stop("Invalid input LD folder")
    }

    if(thresh == 1){
        thresh = info$thresh
    }

    dir.create(outputDir)
    # copy two info file
    #file.copy(ldmInfoFile, file.path(outputDir, "ldm.info"))

    info.snp = fread(file.path(inputDir, "snp.info"))
    infonames = colnames(info.snp)
    if("blk" %in% infonames){
        setnames(info.snp, "blk", "Block")
    }
    if("CHR" %in% infonames){
        setnames(info.snp, "CHR", "Chrom")
    }

    if("SNP" %in% infonames){
        setnames(info.snp, "SNP", "ID")
    }
 
    if("BP" %in% infonames){
        setnames(info.snp, "BP", "PhysPos")
    }
  
    if("POS" %in% infonames){
        setnames(info.snp, "POS", "PhysPos")
    }

    if("GD" %in% infonames){
        setnames(info.snp, "GD", "GenPos")
    }

    if(!"GenPos" %in% infonames){
        info.snp[, GenPos:=0]
    }

    info.snp[, Index:=0:(nrow(info.snp) -1)]
    fwrite(info.snp[, .(Chrom,ID,Index,GenPos,PhysPos,A1,A2,A1Freq,N, Block)], file=file.path(outputDir, "snp.info"), quote=F, sep="\t", na="NA")

    info.ldm = info.snp[, list(Chrom=.SD[1]$Chrom, StartSnpIdx = head(.SD$Index, 1), StartSnpID=head(.SD$ID, 1),  EndSnpIdx=tail(.SD$Index, 1), EndSnpID=tail(.SD$ID, 1), NumSnps=.N), by="Block"]
    fwrite(info.ldm, file=file.path(outputDir, "ldm.info"), sep="\t", quote=F, na="NA")

    cutLDc(info$template, info$type, blocks, outputDir, thresh)

    message("All files (", length(blocks),  ") has been proceeded.")
    logger.end()
    if(log2file){
        message("All files (", length(blocks),  ") has been proceeded.")
    }
}


