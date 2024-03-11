# SBayesRC: SBayesRC-multi

# License: GPL
# Author: Zhili Zheng <zhilizheng@outlook.com>
#' @title caculate SBayesRC-multi PRS
#' @usage sbrcMulti(prs1, prs2, outPrefix, tuneid, pheno)
#' @param prs1 string, path to the polygenic score set1 from SBayesRC::prs
#' @param prs2 string, path to the polygenic score set2 from SBayesRC::prs
#' @param outPrefix string, output path
#' @param tuneid string, path to sample id (FID IID) for the tunning
#' @param pheno string, path to the phenotype (FID IID pheno; without header)
#' @param keepid string, path to sample id (FID IID) to be kept, default "": all samples
#' @param log2file boolean, FALSE: display message on terminal; TRUE: redirect to an output file; default value is FALSE
#' @return none, results in the specified output
#' @export
sbrcMulti <- function(prs1, prs2, outPrefix, tuneid, pheno, keepid="", log2file=FALSE){
    message("Performing SBayesRC-multi...")
    logfile = paste0(outPrefix, ".log")
    logger.begin(logfile, log2file)
    if(log2file){
        message("Performing SBayesRC-multi...")
    } 

    dt.tune = fread(tuneid, head=FALSE)
    setnames(dt.tune, c("FID", "IID"))
    message(nrow(dt.tune), " tuning samples")
    if(! (file.exists(prs1) & file.exists(prs2))){
        stop("Eighter prs1 or prs2 doesn't exist: ", prs1, ", ", prs2)
    }
    dt.prs1 = fread(prs1)[is.finite(SCORE)]
    dt.prs2 = fread(prs2)[is.finite(SCORE)]

    message(nrow(dt.prs1), " samples in PRS1")
    message(nrow(dt.prs2), " samples in PRS2")

    if(file.exists(keepid)){
        keep.id = fread(keepid, head=FALSE)
        message(nrow(keep.id), " samples in keepid")
        idx1 = match(keep.id[[2]], dt.prs1$IID)
        dt.val1 = dt.prs1[idx1[is.finite(idx1)]]
        message(nrow(dt.val1), " samples remained merging with keepid")

        idx2 = match(keep.id[[2]], dt.prs2$IID)
        dt.val2 = dt.prs2[idx2[is.finite(idx2)]]
        message(nrow(dt.val2), " samples remained merging with keepid")
    }else{
        dt.val1 = dt.prs1
        dt.val2 = dt.prs2
    }

    dt.prs = merge(dt.val1, dt.val2[, .(IID, SCORE)], by="IID")
    message(nrow(dt.prs), " common samples between PRS1 and PRS2")

    dt.pheno.raw = fread(pheno)
    setnames(dt.pheno.raw, c("FID", "IID", "pheno"))
    dt.pheno = dt.pheno.raw[is.finite(pheno)]
    message(nrow(dt.pheno), " samples in phenotype")
    dt.prs.all = merge(dt.prs, dt.pheno[, .(IID, pheno)], by="IID")
    message(nrow(dt.prs.all), " samples after merging all data")

    dt.prs.all.tune = dt.prs.all[IID %in% dt.tune$IID]
    message(nrow(dt.prs.all.tune), " valid tuning samples")

    message("Tuning weights...")
    model = lm(pheno ~ SCORE.x + SCORE.y, data=dt.prs.all.tune)
    dt.prs.all[, SCORE:=predict(model, dt.prs.all[ ,.(SCORE.x, SCORE.y)])]
    fwrite(dt.prs.all[, .(FID, IID, SCORE)], file=paste0(outPrefix, ".score.txt"), sep="\t", na="NA", quote=FALSE)

    coef1 = summary(model)$coefficients[,1]
    winter = paste0("intercept: ", coef1[1])
    w1 = paste0("weight_PRS1: ", coef1[2])
    w2 = paste0("weight_PRS2: ", coef1[3])
    message(winter)
    message(w1)
    message(w2)

    dt.prs.all.tune2 = dt.prs.all[IID %in% dt.tune$IID]
    #print(summary(model))
    r1 = summary(lm(pheno~SCORE.x, data=dt.prs.all.tune2))$adj.r.squared
    r2 = summary(lm(pheno~SCORE.y, data=dt.prs.all.tune2))$adj.r.squared
    r3 = summary(lm(pheno~SCORE, data=dt.prs.all.tune2))$adj.r.squared
    message("Adj R-squared for PRS1 in tuning samples: ", r1)
    message("Adj R-squared for PRS2 in tuning samples: ", r2)
    message("Adj R-squared for weighted PRS in tuning samples: ", r3)
    message("Note: those R-squared was calculated in tuning samples, please don't report as prediction accuracy!")
    message("Use the weighted PRS in ", outPrefix, ".score.txt and exclude the tuning samples to check the prediction accuracy.")
    if(r3 < r1 | r3 < r2){
        warning("Weighted PRS has lower R-squared, the inputs may have problem, please double-check.")
    }
    #print(summary(lm(pheno~SCORE.y, data=dt.prs.all.tune)))

    cat(c(winter, w1, w2), file=paste0(outPrefix, ".weight"), sep="\n")

    message("Done.")
    logger.end()
    if(log2file){
        message("Done.")
    }
}
