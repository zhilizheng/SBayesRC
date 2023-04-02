
#' @title SBayesRC with eigen decomposition
#' @usage sbayesrc(ma_path, LD_folder, file_out)
#' @param file_summary string, summary data path, COJO format
#' @param ld_folder string, path to LD folder
#' @param file_out string, output path
#' @param thresh number, eigen cutoff for LD
#' @param fileAnnot sring, annotation file path
#' @param niter number, total number of MCMC iterations
#' @param burn number, numbmer of MCMC burn-in iterations
#' @param starth2 number, start heritability
#' @param startPi vector, start proportion of causal in each component 
#' @param gamma vector, gamma of each componnet
#' @param bTune bool, perform tuning or not
#' @return none, outputs are all redirected the file_out
#' @examples
#' # run SBayesRC without annotation
#' sbayesrc(ma_path, LD_folder, file_out)
#' # run SBayesRC with annotation
#' sbayesrc(ma_path, LD_folder, file_out, fileAnnot=annot_file_path)
#' @export
sbayesrc = function(file_summary, ld_folder, file_out, thresh=0.995, niter=3000, burn=1000, fileAnnot="", 
                    method="sbr_ori", sSamVe="allMixVe", twopq="nbsq",
                    bOutDetail=FALSE, resam_thresh=1.1, 
                    starth2=0.5, startPi=c(0.990, 0.005, 0.003, 0.001, 0.001), gamma=c(0, 0.001, 0.01, 0.1, 1), bTune=TRUE, tuneIter=150, tuneBurn=100, tuneStep=c(0.995, 0.99, 0.95, 0.9), bTunePrior=FALSE, seed=22, outFreq=10, annoSigmaScale=1.0, bOutBeta=FALSE){
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
        message("Don't need to run: ", outfile, " exists")
        return;
    }

    if(submethod == "sbc"){
        stop("annot can't work in SBayesC yet")
    }

    message("------------------------------------------")
    message("SBayesRC v", packageVersion("SBayesRC"))
    message("Developed by Zhili Zheng and Jian Zeng")
    message("------------------------------------------")


    logfile = paste0(outfile, ".log")
    message("Logging to ", logfile)
    message("Please look into log file, the console output will be empty")

    zz = file(logfile, "wt")
    sink(zz)
    sink(zz, type="message")

    message("------------------------------------------")
    message("SBayesRC v", packageVersion("SBayesRC"))
    message("Developed by Zhili Zheng and Jian Zeng")
    message("------------------------------------------")

    message("Node: ", Sys.info()['nodename'])
    #message("Run with method ", submethod, ", ori?", bOri)
    message(" With annotation: ", bAnnot)
    if(bAnnot) message(" Annotation file: ", fileAnnot)
    message(" Summary file: ", file_summary)
    message(" LD input: ", ld_folder, ", var threshold: ", thresh)
    message(" Re-sample Ve: ", cSamVe)
    message(" Start Pi: ", paste(startPi, collapse=" "))
    message(" Gamma for effects: ", paste(gamma, collapse=" "))
    nComp = length(startPi)
    if(nComp != length(gamma)){
        stop("The number of component is not consistent between pi and gamma")
    }
    message(" Number of iteration: ", niter)
    message("   Burn-in iteration: ", burn)
    message(" Start hsq: ", starth2)
    message(" Use 2pq: ", bTwopq)
    #message(" output detail: ", bOutDetail)
    message(" Threshold to re-sample residual: ", resam_thresh)
    message("------------------------------------------")

    mafile = file_summary
    ma = fread(mafile)

    info = fread(file.path(ld_folder, "snp.info"))
    m = nrow(info)

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
        stop("LD information and GWAS summary statistics have some alleles mismatching!")
    }
    ma_ord[bA1A2, freq := 1 - freq]
    ma_ord[bA1A2, b := (-b)]
    message(sum(bA1A2), " SNPs flipped allele")

    if("n" %in% colnames(ma_ord)){
        setnames(ma_ord, "n", "N")
    }

    ma_ord[, blk:=info$blk]
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

    ### annot

    annoMat = matrix(0, ncol=0, nrow=0)
    numAnno = 0
    annoStrings = c("")
    tempFileAnnot = paste0(outfile, ".annot.tmp.bin")
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
        cRand = file(file.path(ld_folder, "rand.bin"), "rb")
        randIdx = sample(1:50, 1)
        seek(cRand, where=(randIdx-1) * m * 4)

        ma_ord$bt_add = readBin(cRand, numeric(), size=4, m)

        ma_ord[, n_train:=round(N*0.9)]
        ma_ord[, n_val:=N - n_train]
        ma_ord[, factor := sqrt(1/n_train - 1/N)]

        ma_ord[, bhat_t := bhat + factor * bt_add]
        ma_ord[, bhat_v := (bhat * N - bhat_t * n_train) / n_val]

        fwrite(ma_ord[, .(SNP, A1, A2, freq, b, se, p, N, blk, bhat_t, bhat_v)], file=paste0(outfile, "_tune_inter.txt"), sep="\t", na="NA", quote=F)

        dt.n = ma_ord[, list(N=mean(N)), by=blk][order(blk)]
        n_train = as.numeric(dt.n$N)

        bhat_t = ma_ord$bhat_t
        bhat_v = ma_ord$bhat_v
        message("Generated pseudo-summary ", randIdx)
    }

    rm(ma, ma_ord, info, bA1A1, bA2A2, bA1A2, bA2A1)
    gc()

    if(bTune){
        message("********************************")
        message("Perform tuning...")
        dt.tune = data.table()
        for(tune_thresh in tuneStep){
            curFile = paste0(outfile, "_tune", tune_thresh)
            message("--------------------------------")
            message("Tune with parameter ", tune_thresh)
            res = sbayesr_eigen_joint_annot(tuneIter, tuneBurn, bhat_t, 0, c(""), ld_folder, vary, n_train, gamma, startPi, starth2, tune_thresh, bOri, curFile, cSamVe, resam_thresh, bOutDetail, outFreq, annoSigmaScale, bOutBeta)
            saveRDS(res, file=paste0(curFile, ".rds"))

            dt.tune = rbind(dt.tune, data.table(thresh=tune_thresh, r=(t(res$betaMean) %*% bhat_v)[1] / sqrt(sum(res$betaMean^2))))
            rm(res)
            gc()
        }

        dt.tune$rel_r = abs(dt.tune$r / dt.tune[thresh=="0.995"]$r)
        fwrite(dt.tune, file=paste0(outfile, "_tune.txt"), sep="\t", quote=F, na="NA")

        message("Tuning Done")
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

        message("Continue with best parameter: ", thresh)
        if(bTunePrior){
            dt = readRDS(paste0(outfile, "_tune", thresh, ".rds"))
            startPi = dt$pi_hist[nrow(dt$pi_hist), ]
            starth2 = tail(dt$hsq_hist, 1)
            rm(dt)
        }

        runtime = proc.time()
        message("Time elapsed:", runtime["elapsed"], " seconds")
        message("********************************")
        rm(bhat_t, bhat_v, dt.tune )
        gc()
    }

    message("SBayesRC with eigen cutoff ", thresh)
    outRes = paste0(outfile, ".rds")
    if(file.exists(outRes)){
        message("The parameter file exists, loading the parameter instead of a re-run: ", outRes)
        res = readRDS(outRes)
    }else{
        res = sbayesr_eigen_joint_annot(niter, burn, bhat, numAnno, annoStrings, ld_folder, vary, n, gamma, startPi, starth2, thresh, bOri, outfile, cSamVe, resam_thresh, bOutDetail, outFreq, annoSigmaScale, bOutBeta)
        saveRDS(res, file=outRes)
    }

    if(numAnno != 0){
        file.remove(tempFileAnnot)
    }

    out = data.table(SNP=snp, A1=a1, BETA=res$betaMean * ord_std, PIP=(1 - res$pip[,1]),  BETAlast=res$betaLast * ord_std)
    fwrite(out, file=paste0(outfile, ".txt"), sep="\t", na="NA", quote=FALSE)

    message("------------------------------------------")
    message("Parameter estimation:")
    print(res$par)

    out.par = data.table(name=names(res$par), value=res$par)
    fwrite(out.par, file=paste0(outfile, ".par"), sep="\t", na="NA", quote=FALSE, col.names=FALSE)

   if(numAnno != 0){
       # per-SNP heritability
       curDT = fread(paste0(outfile, ".vg.enrich.qt"), sep=" ")
       curMean = colMeans(curDT)
       curSD = sapply(curDT, sd)
       outDT = data.table(Annotation = names(curMean), Enrich = curMean, SD = curSD)
       fwrite(outDT, file=paste0(outfile, ".hsq.enrich"), sep="\t")

       # comp pi
       outDTs = data.table()
       for(curComp in 0:(nComp - 1)){
           curDT = fread(paste0(outfile, ".annoJointProb", curComp), sep=" ")
           curMean = colMeans(curDT)
           outDT = data.table(Comp=curComp, Gamma=gamma[curComp+1])
           outDTs = rbind(outDTs, cbind(outDT, t(curMean)))
       }
       fwrite(outDTs, file=paste0(outfile, ".annoJointProb"), sep="\t")

   }

    message("SBayesRC run successfully!")
    runtime = proc.time()
    message("Time elapsed:", runtime["elapsed"], " seconds")
}
