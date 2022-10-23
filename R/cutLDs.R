
#' @title cut the LD to a different variance threshold and save the storage
#' @usage cutLDs(inputDir, outputDir, 0.999)
#' @param inputDir string, LD path for input
#' @param outputDir string, LD path for output
#' @param thresh numeric, new threshold to cut
#' @return none, the output is in the output file
#' @export
cutLDs = function(inputDir, outputDir, thresh){
    message("All messages are redirected to ", outputDir, ".log")
    zz = file(paste0(outputDir, ".log"), "wt")
    sink(zz)
    sink(zz, type="message")
 
    ldmInfoFile = paste0(inputDir, "/", "ldm.info")
    if(file.exists(ldmInfoFile)){
        ldm = fread(ldmInfoFile)
    }else{
        stop("not a valid input folder: ", inputDir)
    }

    blocks = ldm$block
    
    guessVars = c(0.9999, 0.9995, 0.999, 0.995)

    curVar = 0
    for(guessVar in guessVars){
        bin1file = paste0(inputDir, "/eig_block", blocks[1], "_var", guessVar, ".bin3")
        if(file.exists(bin1file)){
            curVar = guessVar
            break
        }
    }

    if(curVar == 0){
        stop("Input folder doesn't have a binary format of with a valid threshold, e.g. 0.9999, 0.9995, 0.999, 0.995")
    }

    if(file.exists(outputDir)){
        stop("The output folder is there already, can't overwrite")
    }

    dir.create(outputDir, showWarnings=F)
    # copy two info file
    file.copy(ldmInfoFile, paste0(outputDir, "/ldm.info"))
    file.copy(paste0(inputDir, "/snp.info"),  paste0(outputDir, "/snp.info"))

    # cut blocks
    n = length(blocks)
    for(block in blocks){
        message("Processing block ", block)
        bin1file = paste0(inputDir, "/eig_block", block, "_var", curVar, ".bin3")
        out1file = paste0(outputDir, "/eig_block", block, "_var", thresh, ".bin3")
        cutLD(bin1file, out1file, thresh)
    }

    message("All files (", n,  ") has been proceeded.")

}


