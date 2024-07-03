# SBayesRC main function
# License: GPL
# Author: Zhili Zheng <zhilizheng@outlook.com>

#' @title SBayesRC analysis
#' @usage sbayesrc(mafile, LDdir, outPrefix, annot="AnnotFile_PATH")
#' @param mafile string, path of summary data in COJO format
#' @param LDdir string, path to LD folder
#' @param outPrefix string, output path
#' @param annot string, annotation file path
#' @param log2file boolean, FALSE: display message on terminal; TRUE: redirect to an output file; the default is FALSE
#' @param bTune boolean, perform tuning or not (TRUE, FALSE), default TRUE (yes)
#' @param tuneIter integer, total number of tuning iterations
#' @param tuneBurn integer, number of tuning iterations to burn
#' @param thresh numeric, eigen variance cutoff for LD, default 0.995
#' @param tuneStep numeric vector, eigen variance cutoffs to search
#' @param bTunePrior boolean, use prior estimated from tuning step or not
#' @param niter integer, total number of MCMC iterations
#' @param burn integer, numbmer of MCMC burn-in iterations
#' @param starth2 numeric, heritability value to start iteration
#' @param startPi vector, start proportion of causal in each component 
#' @param gamma vector, gamma of each componnet
#' @param method string, method to use (always sbr_ori)
#' @param sSamVe string, method to re-sample residual (always allMixVe)
#' @param twopq string, method to scale the beta
#' @param bOutDetail boolean, output more details or not
#' @param resam_thresh numeric, threshold to resampling the residual variance
#' @param seed integer, random seed
#' @param outFreq integer, frequencies to calculate the parameters and output
#' @param annoSigmaScale numeric, scale to the annotation Sigma
#' @param bOutBeta boolean, output raw beta or not
#' @param exclude string, path to a file contains the SNPs excluded from modeling, defalut empty;
#' @return none, results in the specified output
#' @examples
#' ## run SBayesRC without annotation
#' #sbayesrc(mafile, LDdir, output)
#' ## run SBayesRC with annotation
#' # sbayesrc(mafile, LDdir, ouput, annot=annot_file_path)
#' @export
sbayesrc = function(mafile, LDdir, outPrefix, annot="", log2file=FALSE,
                    bTune=TRUE, tuneIter=150, tuneBurn=100, thresh=0.995, 
                    tuneStep=c(0.995, 0.99, 0.95, 0.9), bTunePrior=FALSE, 
                    niter=3000, burn=1000, starth2=0.5,
                    startPi=c(0.990, 0.005, 0.003, 0.001, 0.001), 
                    gamma=c(0, 0.001, 0.01, 0.1, 1), 
                    method="sbr_ori", sSamVe="allMixVe", twopq="nbsq",
                    bOutDetail=FALSE, resam_thresh=1.1, seed=22, 
                    outFreq=10, annoSigmaScale=1.0, bOutBeta=TRUE, exclude=''){

    original_options <- options()

    file_summary = mafile
    ld_folder = LDdir
    file_out = outPrefix
    fileAnnot = annot

    cSamVe = "fixVe"
    if(sSamVe == "noReSamVe"){
        cSamVe = "fixVe"
    }else if(sSamVe == "reSamVe"){
        cSamVe = "samVe"
    }else if(sSamVe == "mixSamVe"){
        cSamVe = "mixVe"
    }else if(sSamVe == "allMixVe"){
        cSamVe = "allMixVe"
    }else{
        stop("Unknown resample Ve:", sSamVe)
    }

    bTwopq = TRUE
    if(twopq == "twopq"){
        bTwopq = TRUE
    }else if(twopq == "nbsq"){
        bTwopq = FALSE
    }else{
        stop("twopq is unknown: ", twopq)
    }


    bOri = FALSE
    bAnnot = FALSE
    submethod = ""


    if(fileAnnot != ""){
        if(file.exists(fileAnnot)){
            bAnnot = TRUE
        }else{
            stop("Annotation file: ", fileAnnot, ", is not there")
        }
    }else{
        bAnnot = FALSE
    }


    if(method=="sbr_ori"){
        submethod="sbr"
        bOri = TRUE
    }else if(method=="sbr_robust"){
        submethod="sbr"
        bOri = FALSE
    }else if(method=="sbc"){
        submethod="sbc"
    }else{
        stop("Unknow method: ", method)
    }


    outfile = file_out

    if(file.exists(paste0(outfile, ".txt"))){
        message("The output: ", outfile, " exists")
        message("Don't need to run")
        return
    }

    if(submethod == "sbc"){
        stop("annot can't work in SBayesC yet")
    }

    message("------------------------------------------")
    message("SBayesRC v", packageVersion("SBayesRC"))
    message("Developed by Zhili Zheng and Jian Zeng")
    message("------------------------------------------")

    logger.begin(outfile, log2file)

    if(log2file){
        message("------------------------------------------")
        message("SBayesRC v", packageVersion("SBayesRC"))
        message("Developed by Zhili Zheng and Jian Zeng")
        message("------------------------------------------")
    }

    message(" GWAS summary statistics: ", file_summary)
    message(" Annotation: ", bAnnot)
    if(bAnnot) message(" Annotation file: ", fileAnnot)
    message(" LD input folder: ", ld_folder, ", eigen variance cutoff threshold: ", thresh)
    message(" Tune the best eigen cutoff: ", bTune)
    if(bTune){
        message("   Tuning step: ", paste(tuneStep, collapse=" "))
        message("   Tuning iterations: ", tuneIter)
        message("   Tuning burn-in iterations: ", tuneBurn)
    }
    message(" Start hsq: ", starth2)
    message(" Start Pi: ", paste(startPi, collapse=" "))
    message(" Gamma: ", paste(gamma, collapse=" "))
    nComp = length(startPi)
    if(nComp != length(gamma)){
        stop("The number of component is not consistent between pi and gamma")
    }
    message(" Number of MCMC iteration: ", niter)
    message("    Burn-in iteration: ", burn)
    message(" Use 2pq to scale summary: ", bTwopq)
    message(" Method to resample residual: ", cSamVe)
    message(" Threshold to resample residual: ", resam_thresh)
    message(" Hostname: ", Sys.info()['nodename'])
    message(" Analysis started: ", Sys.time())
    message("------------------------------------------")

    mafile = file_summary
    ma = fread(mafile)

    info = fread(file.path(ld_folder, "snp.info"))
    m = nrow(info)

    if("blk" %in% colnames(info)){
        setnames(info, "blk", "Block")
    }

    # del SNPs 
    rmSNPIndices = as.integer(c())
    if(exclude != ""){
        if(file.exists(exclude)){
            rmSNPs = readLines(exclude)
            message(length(rmSNPs), " SNPs in exclude file [", exclude, "]. ")
            idx1 = match(rmSNPs, info$ID)
            valid_idx = idx1[is.finite(idx1)] - 1
            rmSNPIndices = as.integer(sort(valid_idx))
            message(length(rmSNPIndices), " SNPs will be excluded from the analysis.")
        }else{
            stop("exclude file does not exist: ", exclude)
        }
    }


    # correct order
    idx = match(info$ID, ma$SNP)
    if(any(is.na(idx))){
        stop("Some SNPs are missing in the GWAS summary statistics while existing in LD. Try \"impute\" function to impute.")
    }

    ma_ord = ma[idx]
    if(sum(info$ID != ma_ord$SNP)!=0){
        stop("Some SNPs are missing in the GWAS summary statistics while existing in LD. Try \"impute\" function to impute.")
    }

    # flip allele

    bA1A2 = (info$A1 == ma_ord$A2)
    bA2A1 = (info$A2 == ma_ord$A1)
    bA1A1 = (info$A1 == ma_ord$A1)
    bA2A2 = (info$A2 == ma_ord$A2)

    if(sum(bA1A2 != bA2A1) != 0 | sum(bA1A1 != bA2A2) != 0){
        stop("LD information and GWAS summary statistics have some alleles mismatching, try \"tidy\" function.")
    }
    ma_ord[bA1A2, freq := 1 - freq]
    ma_ord[bA1A2, b := (-b)]
    message(sum(bA1A2), " SNPs flipped allele")

    if("n" %in% colnames(ma_ord)){
        setnames(ma_ord, "n", "N")
    }

    ma_ord[, blk:=info$Block]
    dt.n = ma_ord[, list(N=mean(N)), by=blk][order(blk)]
    n = as.numeric(dt.n$N)

    #n = median(ma_ord$N)
    #n = rep(median(ma_ord$N), nrow(ma_ord))
    ma_ord[, D:=2 * freq * (1 - freq) * N]
    ma_ord[, varps:=D*(N * se^2 + b^2)/N]
    vary = median(ma_ord$varps)
    message("Observed phenotypic variance: ", vary)


    ord_std = c()
    if(bTwopq){
        ord_std = 1.0 / sqrt( 2 * ma_ord$freq * (1 - ma_ord$freq)) 
    }else{
        vary = 1;
        ord_std = sqrt((ma_ord$N) * ma_ord$se^2 + ma_ord$b^2)
    }

    bhat = ma_ord$b / ord_std


    snp = ma_ord$SNP
    a1 = info$A1


    ## create dir
    mcmcDir = paste0(outfile, ".mcmcsamples")
    dir.create(mcmcDir, showWarnings = FALSE)

    ### annot

    annoMat = matrix(0, ncol=0, nrow=0)
    numAnno = 0
    annoStrings = c("")
    tempFileAnnot = file.path(mcmcDir, "SBayesRC.annot.tmp.bin")
    if(bAnnot){
        anno = fread(fileAnnot, head=TRUE)

        colAnno = ncol(anno)
        numAnno = colAnno - 1

        annoMat = matrix(0.0, nrow=nrow(ma_ord), ncol=numAnno)

        annoMat[, 1] = 1

        commSNP = intersect(ma_ord$SNP, anno$SNP)

        idx = match(ma_ord$SNP, anno$SNP)

        conv_anno = anno[idx]


        if(colAnno > 2){
            conv_mat = as.matrix(conv_anno[, 3:colAnno])
            annoMat[, 2:numAnno] = conv_mat
            rm(conv_mat)
        }

        annoStrings = colnames(anno)[2:colAnno]

        annoMat[is.na(annoMat)] = 0

        rm(conv_anno, commSNP, idx, anno)

        tempAnnot = file(tempFileAnnot, "wb")
        for(idx1 in 1:numAnno){
            writeBin(annoMat[, idx1], tempAnnot, size=4) 
        }
        close(tempAnnot)
        rm(annoMat)

    }

    if(bTune){
        message("Generating pseudo-summary...")
        ma_ord[, bhat:=bhat]
        pseudoSuma(ma_ord, ld_folder)
        fwrite(ma_ord[, .(SNP, A1, A2, freq, b, se, p, N, blk, bhat, bhat_t, bhat_v, n_train, n_val, bt_add, factor)], file=file.path(mcmcDir, paste0("SBayesRC_tune_inter.txt")), sep="\t", na="NA", quote=F)

        dt.n = ma_ord[, list(N=mean(n_train)), by=blk][order(blk)]
        n_train = as.numeric(dt.n$N)

        bhat_t = ma_ord$bhat_t
        bhat_v = ma_ord$bhat_v
    }

    rm(ma, ma_ord, info, bA1A1, bA2A2, bA1A2, bA2A1)
    gc()

    if(bTune){
        message("********************************")
        message("Perform tuning...")
        dt.tune = data.table()
        for(tune_thresh in tuneStep){
            curFile = file.path(mcmcDir, paste0("SBayesRC_tune", tune_thresh))
            message("--------------------------------")
            message("Tune with parameter ", tune_thresh)
            res = sbayesr_eigen_joint_annot(tuneIter, tuneBurn, bhat_t, 0, c(""), ld_folder, vary, n_train, gamma, startPi, rmSNPIndices, starth2, tune_thresh, bOri, curFile, cSamVe, resam_thresh, FALSE, outFreq, annoSigmaScale, FALSE)
            saveRDS(res, file=paste0(curFile, ".rds"))

            dt.tune = rbind(dt.tune, data.table(thresh=tune_thresh, r=(t(res$betaMean) %*% bhat_v)[1] / sqrt(sum(res$betaMean^2))))
            rm(res)
            gc()
        }

        dt.tune$rel_r = abs(dt.tune$r / dt.tune[thresh=="0.995"]$r)
        fwrite(dt.tune, file=file.path(mcmcDir, paste0("SBayesRC_tune.txt")), sep="\t", quote=F, na="NA")
        message("Tuning information:")
        options(digits = 3)        
        print(dt.tune)
        options(original_options)
        dt.tune = dt.tune[is.finite(r) & is.finite(rel_r)]

        if(nrow(dt.tune) == 0){
            stop("Invalid tune outputs, all correlations are invalid number. The input data has problem.\nPlease perform QC on summary data.\nYou can still run without tuning by passing bTune=FALSE and set a proper threshold by thresh={value among 0~1}")
        }

        if(all(dt.tune$r < 0)){
            stop("All correlations are negative, this may indicate errors in summary data. \nPlease perform QC on summary data.\nYou can still run without tuning by passing bTune=FALSE and set a proper threshold by thresh={value among 0~1}")
        }

        thresh = max(tuneStep)

        dt.tune.pos = dt.tune[r > 0]
        rel_max = max(dt.tune.pos$rel_r)
        dt.max = dt.tune[rel_r == rel_max]
        max_thresh =  dt.max$thresh
        if(rel_max > 1.25){
            thresh = max_thresh
        }else{
            if(max_thresh != max(tuneStep) & dt.tune[thresh==max(tuneStep)]$r < 0){
                thresh = max_thresh
            }
        }

        if(thresh == min(tuneStep)){
            stop("Warning, the best parameter is the minimumn threshold, we suggest to expand the tuning grid by specify lower tuning value, e.g. tuneStep=c(0.995, 0.9, 0.8, 0.7, 0.6)")
        }

        message("Continue with best eigen variance cutoff: ", thresh)
        if(bTunePrior){
            dt = readRDS(file.path(mcmcDir, paste0("SBayesRC_tune", thresh, ".rds")))
            startPi = dt$pi_hist[nrow(dt$pi_hist), ]
            starth2 = tail(dt$hsq_hist, 1)
            rm(dt)
        }

        runtime = proc.time()
        message("Time elapsed: ", round(runtime["elapsed"] / 3600, 2), " hour(s)")
        message("********************************")
        rm(bhat_t, bhat_v, dt.tune )
        gc()
    }

    message("Start SBayesRC with eigen variance cutoff ", thresh)
    outRes = paste0(outfile, ".rds")
    if(file.exists(outRes)){
        message("The parameter file exists, loading the parameter instead of a re-run: ", outRes)
        res = readRDS(outRes)
    }else{
        res = sbayesr_eigen_joint_annot(niter, burn, bhat, numAnno, annoStrings, ld_folder, vary, n, gamma, startPi, rmSNPIndices, starth2, thresh, bOri, file.path(mcmcDir, "SBayesRC"), cSamVe, resam_thresh, bOutDetail, outFreq, annoSigmaScale, bOutBeta)
        res[["scale"]] = ord_std
        saveRDS(res, file=outRes)
    }

    if(numAnno != 0){
        file.remove(tempFileAnnot)
    }

    message("MCMC cycles completed.")

    out = data.table(SNP=snp, A1=a1, BETA=res$betaMean * ord_std, SE=res$betaSD * ord_std, PIP=(1 - res$pip[,1]),  BETAlast=res$betaLast * ord_std)
    fwrite(out, file=paste0(outfile, ".txt"), sep="\t", na="NA", quote=FALSE)
    message("Use the ", outfile, ".txt, column 1 2 3 to calculate the polygenic risk score.")

    message("------------------------------------------")
    message("Parameter estimation:")

    out.par = data.table(Item=names(res$par), Mean=res$par, SD=res$sd_par)

    mat_nsnps = res$n_comp_hist[(burn+1):niter, -1]
    out.par = rbind(out.par, data.table(Item=paste0("NumSnp", 2:5), Mean=colMeans(mat_nsnps), SD=apply(mat_nsnps, 2, sd)))

    mat_vgs = res$vg_comp_hist[(burn+1):niter, -1]
    out.par = rbind(out.par, data.table(Item=paste0("Vg", 2:5), Mean=colMeans(mat_vgs), SD=apply(mat_vgs, 2, sd)))

    fwrite(out.par, file=paste0(outfile, ".par"), sep="\t", na="NA", quote=FALSE)

    options(digits = 2)
    print(out.par, row.names = FALSE)
    options(original_options)

   if(numAnno != 0){
       # per-SNP heritability
       curDT = fread(file.path(mcmcDir, paste0("SBayesRC.mcmcsamples.AnnoPerSnpHsqEnrichment")), sep=" ")
       curMean = colMeans(curDT)
       curSD = sapply(curDT, sd)
       outDT = data.table(Annotation = names(curMean), Enrich = curMean, SD = curSD)
       fwrite(outDT, file=paste0(outfile, ".AnnoPerSnpHsqEnrichment"), sep="\t")

       # comp pi
       outDTs = data.table()
       for(curComp in 1:nComp){
           curDT = fread(file.path(mcmcDir, paste0("SBayesRC.mcmcsamples.AnnoJointProb_pi", curComp)), sep=" ")
           curMean = colMeans(curDT)
           outDT = data.table(Comp=curComp, Gamma=gamma[curComp])
           outDTs = rbind(outDTs, cbind(outDT, t(curMean)))
       }
       fwrite(outDTs, file=paste0(outfile, ".AnnoJointProb"), sep="\t")
   }

    message("Analysis finished: ", Sys.time())
    runtime = proc.time()
    message("Computational time: ", round(runtime["elapsed"] / 3600, 2), " hour(s)")
    logger.end()
    if(log2file){
        message("Use the ", outfile, ".txt, column 1 2 3 to calculate the polygenic risk score.")
        message("Parameter estimation:")

        options(digits = 2)
        print(out.par, row.names = FALSE)
        options(original_options)

        message("Analysis finished: ", Sys.time())
        message("Computational time: ", round(runtime["elapsed"] / 3600, 2), " hour(s)")
    }
}


pseudoSuma <- function(suma, ldm){
    info = getLDPrefix(ldm)
    suma[, idx:=1:nrow(suma)]
    dt.pos = suma[, list(start=.SD[1]$idx - 1), by=blk];
    suma[, idx:=NULL]

    m = nrow(suma)

    rand_vec = rnorm(m)
    bt_add = getPseudoRand(info$template, info$type, m, dt.pos$blk, dt.pos$start, info$thresh, rand_vec)

    suma[, bt_add:=bt_add]
    suma[, n_train:=round(N*0.9)]
    suma[, n_val:=N - n_train]
    suma[, factor := sqrt(1/n_train - 1/N)]
    suma[, bhat_t := bhat + factor * bt_add]
    suma[, bhat_v := (bhat * N - bhat_t * n_train) / n_val]
}

#' @title extract MCMC beta
#' @usage extractMCMCEff(prefix, idx, wPrefix, scale)
#' @param prefix string, output prefix
#' @param idx, index of MCMC to extract
#' @param wFile, output weight file
#' @param scale, orginal scale (TRUE) or standardized scale (FALSE)
#' @return no return value, will be output to wFile (SNP, A1, BETA)
#' @export
extractMCMCEff <- function(prefix, idx, wFile, scale=TRUE){
    mcmcDir = paste0(prefix, ".mcmcsamples")
    filename = file.path(mcmcDir, "SBayesRC.beta.bin")
    if(!file.exists(filename)){
        stop("cannot find ", filename)
    }

    txtFile = paste0(prefix, ".txt")
    if(!file.exists(txtFile)){
        stop("cannot find ", txtFile)
    }
 
    rmIdx = file.path(mcmcDir, "SBayesRC.rm.snpidx")
    rmIndex = c()
    if(file.exists(rmIdx)){
        rmIndex = as.integer(readLines(rmIdx))
        message(length(rmIndex), " removed in the MCMC iterations")
    }

    fp = file(filename, 'rb')
    nr = readBin(fp, "integer", size=4, 1)
    nc = readBin(fp, "integer", size=4, 1)

    nRow = readBin(fp, "integer", size=4, nr)

    start = seek(fp)

    cRow = c(0, cumsum(nRow)*8)

    if(idx > nr){
        stop(idx, " is larger than number of saved MCMC iterations (", nr, ")")
    }
    seek(fp, start + cRow[idx])
    indices = (readBin(fp, "integer", size=4, nRow[idx])) + 1
    value = readBin(fp, "numeric", size=4, nRow[idx])
    close(fp)

    dt.txt = fread(txtFile)
    dt.use = dt.txt[, .(SNP, A1)]
    dt.use[, BETA:=0]
    dt.use[indices, BETA:=value]
    if(scale){
        rds = readRDS(paste0(prefix, ".rds"))
        dt.use[, BETA:=BETA * rds$scale]
    }
    if(length(rmIndex) != 0){
        dt.use[rmIndex+1, BETA:=0]
    }
    fwrite(dt.use, file=wFile, sep="\t", quote=FALSE, na="NA")
    message("Done")
}
